!-------------------------------------------------------------------------------------------------
! Authors: Adam Kit and Kostis Papadakis (2023)
! Program: euler_cfd
! Solving the compressible Euler equations in 3 dimensions using the Lax–Friedrichs
! flux method.
! TODO:
! --Add a flux limiter using a high order reconstruction : probably use Lax–Wendroff for
! F_{h}
! --Supported Boundaries: Periodic, Slip Wall, Outflow
!-------------------------------------------------------------------------------------------------

module types_and_kinds
   implicit none
   integer, parameter :: rk = selected_real_kind(8)
   integer, parameter :: ik = selected_int_kind(8)
end module types_and_kinds

module global
   use types_and_kinds
   implicit none
   real(rk), parameter :: gamma = 5.0_rk/3.0_rk, rs = 287.0_rk, cfl = 0.2_rk
end module global

module initialization
   use types_and_kinds
   implicit none
contains

   subroutine init_grid(grid, nx, ny, nz, nGhosts, init_value)
      real(rk), dimension(:, :, :), intent(out) :: grid
      integer(ik), intent(in) :: nx, ny, nz, nGhosts
      real(rk), intent(in) :: init_value
      integer(ik) :: i, j, k
      ! to efficintly traverse the grid in the order of memory layout
      ! loop over 3d array with proper cache efficiency
      do k = nGhosts + 1, nz - nGhosts
         do j = nGhosts + 1, ny - nGhosts
            do i = nGhosts + 1, nx - nGhosts
               grid(i, j, k) = init_value
            end do
         end do
      end do
   end subroutine init_grid

   ! similar to above subroutine, we now intialize the grid with a 3d gaussian centered in middle of grid
   real(rk) function gaussian3d(i, j, k, mu, sigma)
      real(rk), intent(in) :: mu, sigma
      integer(ik), intent(in) :: i, j, k
      real(rk) :: a = 1.0
      !gaussian3d = a*sin(  2*dacos(-1.d0)* ((i-3)/127.0) )
      gaussian3d = a*exp(-((i - mu)**2 + (j - mu)**2 + (k - mu)**2)/(2*sigma**2))
   end function gaussian3d

   ! use the above function to initialize the grid
   subroutine init_grid_gaussian(grid, nx, ny, nz, nGhosts, init_value_mu, init_value_sigma, magnitude, offset)
      real(rk), dimension(:, :, :), intent(out) :: grid
      integer(ik), intent(in) :: nx, ny, nz, nGhosts
      real(rk), intent(in) :: init_value_mu, init_value_sigma, magnitude, offset
      integer(ik) :: i, j, k
      do k = nGhosts + 1, nz - nGhosts
         do j = nGhosts + 1, ny - nGhosts
            do i = nGhosts + 1, nx - nGhosts
               grid(i, j, k) = magnitude*gaussian3d(i, j, k, init_value_mu, init_value_sigma) + offset
            end do
         end do
      end do
   end subroutine init_grid_gaussian
end module initialization

module io
   use types_and_kinds
   implicit none
contains

   !DEPRECATED USE WRITE_STATE
   subroutine write_grid_sad(grid, filename)
      real(rk), dimension(:, :, :), intent(in) :: grid
      character(len=*), intent(in) :: filename
      ! need to use stream as access for the binary (why?)
      ! maybe retarded but name can not be *.bin
      open (1, file=filename, form='unformatted', access='stream', status='replace')
      write (1) grid
      close (1)
   end subroutine write_grid_sad

   !only for double 8 and same 3D shape for all arrays
   !Credits to https://gist.github.com/dmentipl/ed609e7278050bd6de0016c09155e039
   subroutine add_field_h5(id, name, field)
      use hdf5
      use types_and_kinds
      implicit none
      real(rk), dimension(:, :, :), intent(in) :: field
      character(len=*), intent(in) :: name
      integer(hid_t), intent(in) :: id
      integer :: error

      integer, parameter :: ndims = 3
      integer(hsize_t)   :: data_shape(ndims)
      integer(hsize_t)   :: chunk(ndims)
      integer(hid_t)     :: dspace_id
      integer(hid_t)     :: dset_id
      integer(hid_t)     :: dtype_id

      data_shape = shape(field)
      chunk = shape(field)
      dtype_id = H5T_NATIVE_DOUBLE

      call h5screate_f(h5s_scalar_f, dspace_id, error)
      if (error /= 0) then
         write (*, '("cannot create hdf5 dataspace",/)')
         return
      end if

      call h5dcreate_f(id, name, dtype_id, dspace_id, dset_id, error)
      if (error /= 0) then
         write (*, '("cannot create hdf5 dataset",/)')
         return
      end if

      call h5dwrite_f(dset_id, dtype_id, field, data_shape, error)
      if (error /= 0) then
         write (*, '("cannot write to hdf5 file",/)')
         return
      end if

      call h5dclose_f(dset_id, error)
      if (error /= 0) then
         write (*, '("cannot close hdf5 dataset",/)')
         return
      end if

      call h5sclose_f(dspace_id, error)
      if (error /= 0) then
         write (*, '("cannot close hdf5 dataspace",/)')
         return
      end if

   end subroutine add_field_h5

   subroutine create_hdf5file(filename, file_id, error)
      use hdf5
      character(len=*), intent(in)  :: filename
      integer(hid_t), intent(out) :: file_id
      integer, intent(out) :: error
      integer :: filter_info
      integer :: filter_info_both
      logical :: avail

      ! initialise hdf5
      call h5open_f(error)
      if (error /= 0) then
         write (*, '("cannot initialise hdf5",/)')
         return
      end if

      ! check if gzip compression is available.
      call h5zfilter_avail_f(h5z_filter_deflate_f, avail, error)
      if (.not. avail) then
         write (*, '("gzip filter not available.",/)')
         return
      end if
      call h5zget_filter_info_f(h5z_filter_deflate_f, filter_info, error)
      filter_info_both = ior(h5z_filter_encode_enabled_f, h5z_filter_decode_enabled_f)
      if (filter_info /= filter_info_both) then
         write (*, '("gzip filter not available for encoding and decoding.",/)')
         return
      end if

      ! create file
      call h5fcreate_f(filename, h5f_acc_trunc_f, file_id, error)
      if (error /= 0) then
         write (*, '("cannot create hdf5 file",/)')
         return
      end if

   end subroutine create_hdf5file

   !Credits to https://gist.github.com/dmentipl/ed609e7278050bd6de0016c09155e039
   subroutine open_hdf5file(filename, file_id, error)
      use hdf5
      implicit none
      character(len=*), intent(in)  :: filename
      integer(hid_t), intent(out) :: file_id
      integer, intent(out) :: error

      ! initialise hdf5
      call h5open_f(error)
      if (error /= 0) then
         write (*, '("cannot initialise hdf5",/)')
         return
      end if

      ! open file
      call h5fopen_f(filename, h5f_acc_rdwr_f, file_id, error)
      if (error /= 0) then
         write (*, '("cannot open hdf5 file",/)')
         return
      end if

   end subroutine open_hdf5file

   !Credits to https://gist.github.com/dmentipl/ed609e7278050bd6de0016c09155e039
   subroutine close_hdf5file(file_id, error)
      use hdf5
      implicit none
      integer(hid_t), intent(in)  :: file_id
      integer, intent(out) :: error

      ! close file
      call h5fclose_f(file_id, error)
      if (error /= 0) then
         write (*, '("cannot close hdf5 file",/)')
         return
      end if

      ! close hdf5
      call h5close_f(error)
      if (error /= 0) then
         write (*, '("cannot close hdf5",/)')
         return
      end if

   end subroutine close_hdf5file

   subroutine write_state(tstep, rho, vx, vy, vz, p,&
                                 mass,momentum_x,momentum_y,&
                                 momentum_z,energy,temp)
      use hdf5
      use types_and_kinds
      implicit none
      real(rk), dimension(:, :, :), intent(in) :: rho, vx, vy, vz, p,&
                                    mass,momentum_x,momentum_y,&
                                    momentum_z,energy,temp
      integer(ik), intent(in) :: tstep
      character(32) :: fname
      integer(hid_t) :: file_id
      integer:: error
      write (fname, '(a,"000",i7.7,a)') "state.", tstep, ".h5"
      call create_hdf5file(fname, file_id, error)
      call open_hdf5file(fname, file_id, error)
      call add_field_h5(file_id, "rho",rho)
      call add_field_h5(file_id, "vx", vx)
      call add_field_h5(file_id, "vy", vy)
      call add_field_h5(file_id, "vz", vz)
      call add_field_h5(file_id, "p", p)
      call add_field_h5(file_id, "mass",mass)
      call add_field_h5(file_id, "momentum_x", momentum_x)
      call add_field_h5(file_id, "momentum_y", momentum_y)
      call add_field_h5(file_id, "momentum_z", momentum_z)
      call add_field_h5(file_id, "energy", energy)
      call add_field_h5(file_id, "temperature",temp)
      call close_hdf5file(file_id, error)
   end subroutine write_state

   subroutine write_conserved(tstep,mass,momentum_x,momentum_y,momentum_z,energy,temp)
      use hdf5
      use types_and_kinds
      implicit none
      real(rk), dimension(:, :, :), intent(in) :: mass,momentum_x,momentum_y,momentum_z,energy,temp
      integer(ik), intent(in) :: tstep
      character(32) :: fname
      integer(hid_t) :: file_id
      integer:: error
      write (fname, '(a,"000",i7.7,a)') "primitives.", tstep, ".h5"
      call create_hdf5file(fname, file_id, error)
      call open_hdf5file(fname, file_id, error)
      call add_field_h5(file_id, "mass",mass)
      call add_field_h5(file_id, "momentum_x", momentum_x)
      call add_field_h5(file_id, "momentum_y", momentum_y)
      call add_field_h5(file_id, "momentum_z", momentum_z)
      call add_field_h5(file_id, "energy", energy)
      call add_field_h5(file_id, "temperature", temp)
      call close_hdf5file(file_id, error)
   end subroutine write_conserved

end module io

module physics
   use types_and_kinds
   use global
   implicit none
contains

   subroutine conservative(mass, momentum_x, momentum_y, momentum_z, energy, rho, p, vx, vy, vz, temp, ds)
      real(rk), dimension(:, :, :), intent(inout) :: mass, momentum_x, momentum_y, momentum_z, energy, temp
      real(rk), dimension(:, :, :), intent(in) :: rho, vx, vy, vz, p
      real(rk), intent(in) :: ds
      real(rk) :: cell_volume
      cell_volume = ds**3
      mass = rho*cell_volume
      momentum_x = rho*vx*cell_volume
      momentum_y = rho*vy*cell_volume
      momentum_z = rho*vz*cell_volume
      energy = cell_volume*(rho*(vx**2 + vy**2 + vz**2)/2.0_rk + p/(gamma - 1.0_rk))
      temp = p/(rs*rho) - 273.15_rk
   end subroutine conservative

   real(rk) function get_timestep(ds, vx, vy, vz, p, rho)
      real(rk), intent(in) :: ds
      real(rk), intent(in), dimension(:, :, :) :: vx, vy, vz, p, rho
      get_timestep = minval(cfl*ds/sqrt((gamma*p/rho) + (vx**2 + vy**2 + vz**2)))
      ! todo: double check this
      ! get_timestep = minval(cfl*ds/(sqrt(gamma*p/rho) + sqrt(vx**2 + vy**2 + vz**2)) )
   end function get_timestep

   subroutine primitive(mass, momentum_x, momentum_y, momentum_z, energy, rho, p, vx, vy, vz, temp, ds)
      real(rk), dimension(:, :, :), intent(in) :: mass, momentum_x, momentum_y, momentum_z, energy
      real(rk), dimension(:, :, :), intent(inout) :: p, vx, vy, vz, rho, temp
      real(rk), intent(in) :: ds
      real(rk) :: cell_volume
      cell_volume = ds**3
      rho = mass/cell_volume
      vx = momentum_x/(rho*cell_volume)
      vy = momentum_y/(rho*cell_volume)
      vz = momentum_z/(rho*cell_volume)
      p = (energy/cell_volume - 0.5_rk*rho*(vx*vx + vy*vy + vz*vz))*(gamma - 1.0_rk)
      temp = p/(rs*rho) - 273.15_rk
   end subroutine primitive

   subroutine extrapolateprimitves(rho, vx, vy, vz, p, drho_dx, drho_dy, drho_dz, dvx_dx, dvx_dy, &
                                   dvx_dz, dvy_dx, dvy_dy, dvy_dz, dvz_dx, dvz_dy, dvz_dz, &
                                   dp_dx, dp_dy, dp_dz, rho_xtr, vx_xtr, vy_xtr, vz_xtr, p_xtr, dt)
      use types_and_kinds
      implicit none
      real(rk), dimension(:, :, :), intent(in) :: rho, vx, vy, vz, p, drho_dx, drho_dy, drho_dz, &
                                                  dvx_dx, dvx_dy, dvx_dz, dvy_dx, dvy_dy, dvy_dz, dvz_dx, dvz_dy, &
                                                  dvz_dz, dp_dx, dp_dy, dp_dz

      real(rk), dimension(:, :, :), intent(inout) ::rho_xtr, vx_xtr, vy_xtr, vz_xtr, p_xtr
      real(rk), intent(in) :: dt

      rho_xtr = rho - 0.5*dt*(vx*drho_dx + rho*dvx_dx + vy*drho_dy + rho*dvy_dy + vz*drho_dz + rho*dvz_dz); 
      vx_xtr = vx - 0.5*dt*(vx*dvx_dx + vy*dvx_dy + vz*dvx_dz + (1.0/rho)*dp_dx); 
      vy_xtr = vy - 0.5*dt*(vx*dvy_dx + vy*dvy_dy + vz*dvy_dz + (1.0/rho)*dp_dy); 
      vz_xtr = vz - 0.5*dt*(vx*dvz_dx + vy*dvz_dy + vz*dvz_dz + (1.0/rho)*dp_dz); 
      p_xtr = p - 0.5*dt*(gamma*p*(dvx_dx + dvy_dy + dvz_dz) + vx*dp_dx + vy*dp_dy + vz*dp_dz); 
   end subroutine extrapolateprimitves

   subroutine reconstructflux(mass_flux_x, momentum_x_flux_x, momentum_y_flux_x, &
                              momentum_z_flux_x, energy_flux_x, &
                              drho, dvx, dvy, dvz, dp, rho, vx, vy, vz, p, nx, ny, nz, nGhosts, ds, offsets)
      use types_and_kinds
      implicit none
      real(rk), dimension(:, :, :), intent(inout) :: mass_flux_x, momentum_x_flux_x, &
                                                     momentum_y_flux_x, momentum_z_flux_x, energy_flux_x
      real(rk), dimension(:, :, :), intent(in) :: rho, vx, vy, vz, p, drho, dvx, dvy, dvz, dp
      integer(ik), intent(in) :: offsets(3)
      integer(ik), intent(in) :: nx, ny, nz, nGhosts
      real(rk), intent(in) :: ds
      integer(ik) :: i, j, k
      real(rk) :: rho_star, momentum_x_star, momentum_y_star, momentum_z_star, p_star, en_right, en_left, en_star
      real(rk) :: c_l, c_r, c_star
      real(rk) :: rl, rr, vxl, vxr, vyl, vyr, vzl, vzr, pl, pr
      ! start by calculating rho_star, which is average of density
      do k = nGhosts + 1, nz - nGhosts
         do j = nGhosts + 1, ny - nGhosts
            do i = nGhosts + 1, nx - nGhosts

               rl = rho(i + offsets(1), j + offsets(2), k + offsets(3)) - &
                    (drho(i + offsets(1), j + offsets(2), k + offsets(3)))*(ds/2.)

               rr = rho(i, j, k) + (drho(i, j, k))*(ds/2.)

               vxl = vx(i + offsets(1), j + offsets(2), k + offsets(3)) - &
                     (dvx(i + offsets(1), j + offsets(2), k + offsets(3)))*(ds/2.)

               vxr = vx(i, j, k) + (dvx(i, j, k))*(ds/2.)

               vyl = vy(i + offsets(1), j + offsets(2), k + offsets(3)) - &
                     (dvy(i + offsets(1), j + offsets(2), k + offsets(3)))*(ds/2.)

               vyr = vy(i, j, k) + (dvy(i, j, k))*(ds/2.)

               vzl = vz(i + offsets(1), j + offsets(2), k + offsets(3)) - &
                     (dvz(i + offsets(1), j + offsets(2), k + offsets(3)))*(ds/2.)

               vzr = vz(i, j, k) + (dvz(i, j, k))*(ds/2.)

               pl = p(i + offsets(1), j + offsets(2), k + offsets(3)) - &
                    (dp(i + offsets(1), j + offsets(2), k + offsets(3)))*(ds/2.)

               pr = p(i, j, k) + (dp(i, j, k))*(ds/2.)

               en_left = (pl/(gamma - 1.0_rk)) + 0.5_rk*(rl*(vxl**2 + vyl**2 + vzl**2))
               en_right = (pr/(gamma - 1.0_rk)) + 0.5_rk*(rr*(vxr**2 + vyr**2 + vzr**2))

               rho_star = (rl + rr)/2.0_rk
               momentum_x_star = (rl*vxl + rr*vxr)/2.0_rk
               momentum_y_star = (rl*vyl + rr*vyr)/2.0_rk
               momentum_z_star = (rl*vzl + rr*vzr)/2.0_rk
               en_star = 0.5_rk*(en_left + en_right)
               p_star = (gamma - 1.0_rk)*(en_star - 0.5_rk*(momentum_x_star**2 + momentum_y_star**2 + momentum_z_star**2)/rho_star)

               mass_flux_x(i, j, k) = momentum_x_star
               momentum_x_flux_x(i, j, k) = momentum_x_star*momentum_x_star/rho_star + p_star; 
               momentum_y_flux_x(i, j, k) = momentum_x_star*momentum_y_star/rho_star; 
               momentum_z_flux_x(i, j, k) = momentum_x_star*momentum_z_star/rho_star; 
               energy_flux_x(i, j, k) = (en_star + p_star)*(momentum_x_star/rho_star)

               c_l = sqrt(gamma*pl/rl) + abs(vxl)
               c_r = sqrt(gamma*pr/rr) + abs(vxr)
               c_star = max(c_l, c_r)

               mass_flux_x(i, j, k) = momentum_x_star - (0.5_rk*c_star*(rl - rr))
               momentum_x_flux_x(i, j, k) = momentum_x_flux_x(i, j, k) - (c_star*(rl*vxl - rr*vxr))/2.0_rk
               momentum_y_flux_x(i, j, k) = momentum_y_flux_x(i, j, k) - (c_star*(rl*vyl - rr*vyr))/2.0_rk
               momentum_z_flux_x(i, j, k) = momentum_z_flux_x(i, j, k) - (c_star*(rl*vzl - rr*vzr))/2.0_rk
               energy_flux_x(i, j, k) = energy_flux_x(i, j, k) - (c_star*(en_left - en_right))/2.0_rk
            end do
         end do
      end do
   end subroutine reconstructflux

   subroutine addfluxes(mass_flux_x, mass_flux_y, mass_flux_z, &
                        momentum_x_flux_x, momentum_x_flux_y, momentum_x_flux_z, &
                        momentum_y_flux_x, momentum_y_flux_y, momentum_y_flux_z, &
                        momentum_z_flux_x, momentum_z_flux_y, momentum_z_flux_z, &
                        energy_flux_x, energy_flux_y, energy_flux_z, &
                        mass, momentum_x, momentum_y, momentum_z, energy, &
                        nx, ny, nz, nGhosts, dt, ds)
      real(rk), dimension(:, :, :), intent(inout) :: mass_flux_x, mass_flux_y, mass_flux_z, &
                                                     momentum_x_flux_x, momentum_x_flux_y, momentum_x_flux_z, &
                                                     momentum_y_flux_x, momentum_y_flux_y, momentum_y_flux_z, &
                                                     momentum_z_flux_x, momentum_z_flux_y, momentum_z_flux_z, &
                                                     energy_flux_x, energy_flux_y, energy_flux_z
      real(rk), dimension(:, :, :), intent(inout) :: mass, momentum_x, momentum_y, momentum_z, energy
      integer(ik), intent(in) :: nx, ny, nz, nGhosts
      real(rk), intent(in) :: dt, ds
      integer(ik) :: i, j, k

      do k = nGhosts + 1, nz - nGhosts
         do j = nGhosts + 1, ny - nGhosts
            do i = nGhosts + 1, nx - nGhosts
               mass(i, j, k) = mass(i, j, k) - (dt*ds)*(mass_flux_x(i, j, k) - mass_flux_x(i - 1, j, k) + &
                                                        mass_flux_y(i, j, k) - mass_flux_y(i, j - 1, k) + &
                                                        mass_flux_z(i, j, k) - mass_flux_z(i, j, k - 1))
               momentum_x(i, j, k) = momentum_x(i, j, k) - (dt*ds)*(momentum_x_flux_x(i, j, k) - momentum_x_flux_x(i - 1, j, k) + &
                                                                    momentum_x_flux_y(i, j, k) - momentum_x_flux_y(i, j - 1, k) + &
                                                                    momentum_x_flux_z(i, j, k) - momentum_x_flux_z(i, j, k - 1))
               momentum_y(i, j, k) = momentum_y(i, j, k) - (dt*ds)*(momentum_y_flux_x(i, j, k) - momentum_y_flux_x(i - 1, j, k) + &
                                                                    momentum_y_flux_y(i, j, k) - momentum_y_flux_y(i, j - 1, k) + &
                                                                    momentum_y_flux_z(i, j, k) - momentum_y_flux_z(i, j, k - 1))
               momentum_z(i, j, k) = momentum_z(i, j, k) - (dt*ds)*(momentum_z_flux_x(i, j, k) - momentum_z_flux_x(i - 1, j, k) + &
                                                                    momentum_z_flux_y(i, j, k) - momentum_z_flux_y(i, j - 1, k) + &
                                                                    momentum_z_flux_z(i, j, k) - momentum_z_flux_z(i, j, k - 1))
               energy(i, j, k) = energy(i, j, k) - (dt*ds)*(energy_flux_x(i, j, k) - energy_flux_x(i - 1, j, k) + &
                                                            energy_flux_y(i, j, k) - energy_flux_y(i, j - 1, k) + &
                                                            energy_flux_z(i, j, k) - energy_flux_z(i, j, k - 1))
            end do
         end do
      end do
   end subroutine addfluxes
   function maxmod(a, b)
      use types_and_kinds
      implicit none
      real(rk) :: a, b
      real(rk) :: maxmod

      if (abs(a) > abs(b) .and. a*b > 0.d0) then
         maxmod = a
      else if (abs(b) > abs(a) .and. a*b > 0) then
         maxmod = b
      else
         maxmod = 0.d0
      end if

      return
   end function maxmod
   function minmod(a, b)
      use types_and_kinds
      implicit none
      real(rk) :: a, b
      real(rk) :: minmod
      if (abs(a) < abs(b) .and. a*b > 0.d0) then
         minmod = a
      else if (abs(b) < abs(a) .and. a*b > 0) then
         minmod = b
      else
         minmod = 0.d0
      end if
      return
   end function minmod
   real(rk) function vanleerlimiter(a, b)
      real(rk), intent(in) :: a, b
      vanleerlimiter = (sign(1.0d0, a) + sign(1.0d0, b))*min(abs(a), abs(b))/(abs(a) + abs(b) + 1.0d-30)
   end function vanleerlimiter
   subroutine calculate_gradients(grid, dgdx, dgdy, dgdz, nx, ny, nz, nGhosts, ds)
      use types_and_kinds
      implicit none
      real(rk), dimension(:, :, :), intent(in) :: grid
      real(rk), dimension(:, :, :), intent(inout) :: dgdx, dgdy, dgdz
      real(rk), intent(in) :: ds
      integer(ik), intent(in) :: nx, ny, nz, nGhosts
      integer(ik) :: i, j, k
      do k = nGhosts + 1, nz - nGhosts
         do j = nGhosts + 1, ny - nGhosts
            do i = nGhosts + 1, nx - nGhosts
               dgdy(i, j, k) = (grid(i, j + 1, k) - grid(i, j - 1, k))/(2.0_rk*ds)
               dgdx(i, j, k) = (grid(i + 1, j, k) - grid(i - 1, j, k))/(2.0_rk*ds)
               dgdz(i, j, k) = (grid(i, j, k + 1) - grid(i, j, k - 1))/(2.0_rk*ds)
            end do
         end do
      end do
   end subroutine calculate_gradients

   subroutine update_ghosts(grid, nx, ny, nz, nGhosts)
      use types_and_kinds
      implicit none
      real(rk), dimension(:, :, :), intent(inout) :: grid
      integer(ik), intent(in) :: nx, ny, nz, nGhosts
      grid(1, :, :) = grid(3, :, :)
      grid(2, :, :) = grid(3, :, :)
      grid(nx - 1, :, :) = grid(nx - nGhosts, :, :)
      grid(nx, :, :) = grid(nx - nGhosts, :, :)

      grid(:, 1, :) = grid(:, 3, :)
      grid(:, 2, :) = grid(:, 3, :)
      grid(:, nx - 1, :) = grid(:, ny - nGhosts, :)
      grid(:, nx, :) = grid(:, ny - nGhosts, :)

      grid(:, :, 1) = grid(:, :, 3)
      grid(:, :, 2) = grid(:, :, 3)
      grid(:, :, nx - 1) = grid(:, :, nz - nGhosts)
      grid(:, :, nx) = grid(:, :, nz - nGhosts)

      !grid(1,:,:)=grid(nx-nGhosts-1,:,:)
      !grid(2,:,:)=grid(nx-nGhosts,:,:)
      !grid(nx-1,:,:)=grid(3,:,:)
      !grid(nx,:,:)=grid(4,:,:)

      !grid(:,1,:)=grid(:,nx-nGhosts-1,:)
      !grid(:,2,:)=grid(:,nx-nGhosts,:)
      !grid(:,nx-1,:)=grid(:,3,:)
      !grid(:,nx,:)=grid(:,4,:)

      !grid(:,:,1)=grid(:,:,nx-nGhosts-1)
      !grid(:,:,2)=grid(:,:,nx-nGhosts)
      !grid(:,:,nx-1)=grid(:,:,3)
      !grid(:,:,nx)=grid(:,:,4)

   end subroutine update_ghosts

end module physics

program euler_cfd
   use types_and_kinds
   use initialization
   use io
   use physics
   implicit none
   real(rk), dimension(:, :, :), allocatable  :: &
      rho, vx, vy, vz, p, mass, momentum_x, momentum_y, momentum_z, energy, &
      rho_prime, vx_prime, vy_prime, vz_prime, temp, &
      ! fluxes
      mass_flux_x, mass_flux_y, mass_flux_z, &
      momentum_x_flux_x, momentum_x_flux_y, momentum_x_flux_z, &
      momentum_y_flux_x, momentum_y_flux_y, momentum_y_flux_z, &
      momentum_z_flux_x, momentum_z_flux_y, momentum_z_flux_z, &
      energy_flux_x, energy_flux_y, energy_flux_z, &
      drho_dx, drho_dy, drho_dz, &
      dvx_dx, dvx_dy, dvx_dz, &
      dvy_dx, dvy_dy, dvy_dz, &
      dvz_dx, dvz_dy, dvz_dz, &
      dp_dx, dp_dy, dp_dz, rho_xtr, vx_xtr, vy_xtr, vz_xtr, p_xtr

   integer(ik), parameter :: xcells = 128, &
                             ycells = 128, &
                             zcells = 128, &
                             nGhosts = 2
   integer(ik), parameter :: nx = xcells + 2*nGhosts, ny = ycells + 2*nGhosts, nz = zcells + 2*nGhosts
   real(rk), parameter :: ds = 1.0_rk
   real(rk):: dt = 0.0_rk, time = 0.0_rk, time_max = 1.e-2_rk
   integer(ik) :: timestep = 1
   integer(ik) :: shiftx(3), shifty(3), shiftz(3)

   allocate (rho(nx, ny, nz), vx(nx, ny, nz), vy(nx, ny, nz), vz(nx, ny, nz), p(nx, ny, nz), mass(nx, ny, nz), &
             momentum_x(nx, ny, nz), momentum_y(nx, ny, nz), momentum_z(nx, ny, nz), energy(nx, ny, nz), &
             rho_prime(nx, ny, nz), vx_prime(nx, ny, nz), vy_prime(nx, ny, nz), vz_prime(nx, ny, nz), temp(nx, ny, nz), &
             mass_flux_x(nx, ny, nz), mass_flux_y(nx, ny, nz), mass_flux_z(nx, ny, nz), &
             momentum_x_flux_x(nx, ny, nz), momentum_x_flux_y(nx, ny, nz), momentum_x_flux_z(nx, ny, nz), &
             momentum_y_flux_x(nx, ny, nz), momentum_y_flux_y(nx, ny, nz), momentum_y_flux_z(nx, ny, nz), &
             momentum_z_flux_x(nx, ny, nz), momentum_z_flux_y(nx, ny, nz), momentum_z_flux_z(nx, ny, nz), &
             energy_flux_x(nx, ny, nz), energy_flux_y(nx, ny, nz), energy_flux_z(nx, ny, nz), &
             drho_dx(nx, ny, nz), drho_dy(nx, ny, nz), drho_dz(nx, ny, nz), &
             dvx_dx(nx, ny, nz), dvx_dy(nx, ny, nz), dvx_dz(nx, ny, nz), &
             dvy_dx(nx, ny, nz), dvy_dy(nx, ny, nz), dvy_dz(nx, ny, nz), &
             dvz_dx(nx, ny, nz), dvz_dy(nx, ny, nz), dvz_dz(nx, ny, nz), &
             dp_dx(nx, ny, nz), dp_dy(nx, ny, nz), dp_dz(nx, ny, nz), rho_xtr(nx, ny, nz), &
             vx_xtr(nx, ny, nz), vy_xtr(nx, ny, nz), vz_xtr(nx, ny, nz), p_xtr(nx, ny, nz))

   shiftx(:) = 0
   shifty(:) = 0
   shiftz(:) = 0
   shiftx(1) = 1
   shifty(2) = 1
   shiftz(3) = 1
   ! initialize vx, vy, vz to 0
   call init_grid_gaussian(rho, nx, ny, nz, nGhosts, nx/2.0_rk, 8.0_rk, 0.0_rk*1.2_rk, 1.2_rk)
   call init_grid_gaussian(p, nx, ny, nz, nGhosts, nx/2.0_rk, 8.0_rk, 0.0001_rk*101325, 101325.0_rk)
   call init_grid(vx, nx, ny, nz, nGhosts, 0.0_rk)
   call init_grid(vy, nx, ny, nz, nGhosts, 0.0_rk)
   call init_grid(vz, nx, ny, nz, nGhosts, 0.0_rk)
   call init_grid(mass_flux_x, nx, ny, nz, nGhosts, 0.0_rk)
   call init_grid(mass_flux_y, nx, ny, nz, nGhosts, 0.0_rk)
   call init_grid(mass_flux_z, nx, ny, nz, nGhosts, 0.0_rk)
   call init_grid(momentum_x_flux_x, nx, ny, nz, nGhosts, 0.0_rk)
   call init_grid(momentum_x_flux_y, nx, ny, nz, nGhosts, 0.0_rk)
   call init_grid(momentum_x_flux_z, nx, ny, nz, nGhosts, 0.0_rk)
   call init_grid(momentum_y_flux_x, nx, ny, nz, nGhosts, 0.0_rk)
   call init_grid(momentum_y_flux_y, nx, ny, nz, nGhosts, 0.0_rk)
   call init_grid(momentum_y_flux_z, nx, ny, nz, nGhosts, 0.0_rk)
   call init_grid(momentum_z_flux_x, nx, ny, nz, nGhosts, 0.0_rk)
   call init_grid(momentum_z_flux_y, nx, ny, nz, nGhosts, 0.0_rk)
   call init_grid(momentum_z_flux_z, nx, ny, nz, nGhosts, 0.0_rk)
   call init_grid(energy_flux_x, nx, ny, nz, nGhosts, 0.0_rk)
   call init_grid(energy_flux_y, nx, ny, nz, nGhosts, 0.0_rk)
   call init_grid(energy_flux_z, nx, ny, nz, nGhosts, 0.0_rk)

   call update_ghosts(rho, nx, ny, nz, nGhosts)
   call update_ghosts(vx, nx, ny, nz, nGhosts)
   call update_ghosts(vy, nx, ny, nz, nGhosts)
   call update_ghosts(vz, nx, ny, nz, nGhosts)
   call update_ghosts(p, nx, ny, nz, nGhosts)

   call conservative(mass, momentum_x, momentum_y, momentum_z, energy, rho, p, vx, vy, vz, temp, ds)

   !main
   do while (time <= time_max)

      call update_ghosts(rho, nx, ny, nz, nGhosts)
      call update_ghosts(vx, nx, ny, nz, nGhosts)
      call update_ghosts(vy, nx, ny, nz, nGhosts)
      call update_ghosts(vz, nx, ny, nz, nGhosts)
      call update_ghosts(p, nx, ny, nz, nGhosts)
      call primitive(mass, momentum_x, momentum_y, momentum_z, energy, rho, p, vx, vy, vz, temp, ds)
      dt = get_timestep(ds, vx, vy, vz, p, rho)

      call calculate_gradients(rho, drho_dx, drho_dy, drho_dz, nx, ny, nz, nGhosts, ds)
      call calculate_gradients(vx, dvx_dx, dvx_dy, dvx_dz, nx, ny, nz, nGhosts, ds)
      call calculate_gradients(vy, dvy_dx, dvy_dy, dvy_dz, nx, ny, nz, nGhosts, ds)
      call calculate_gradients(vz, dvz_dx, dvz_dy, dvz_dz, nx, ny, nz, nGhosts, ds)
      call calculate_gradients(p, dp_dx, dp_dy, dp_dz, nx, ny, nz, nGhosts, ds)

      call update_ghosts(drho_dx, nx, ny, nz, nGhosts)
      call update_ghosts(drho_dy, nx, ny, nz, nGhosts)
      call update_ghosts(drho_dz, nx, ny, nz, nGhosts)
      call update_ghosts(dvx_dx, nx, ny, nz, nGhosts)
      call update_ghosts(dvx_dy, nx, ny, nz, nGhosts)
      call update_ghosts(dvx_dz, nx, ny, nz, nGhosts)
      call update_ghosts(dvy_dx, nx, ny, nz, nGhosts)
      call update_ghosts(dvy_dy, nx, ny, nz, nGhosts)
      call update_ghosts(dvy_dz, nx, ny, nz, nGhosts)
      call update_ghosts(dvz_dx, nx, ny, nz, nGhosts)
      call update_ghosts(dvz_dy, nx, ny, nz, nGhosts)
      call update_ghosts(dvz_dz, nx, ny, nz, nGhosts)
      call update_ghosts(dp_dx, nx, ny, nz, nGhosts)
      call update_ghosts(dp_dy, nx, ny, nz, nGhosts)
      call update_ghosts(dp_dz, nx, ny, nz, nGhosts)

      call extrapolateprimitves(rho, vx, vy, vz, p, drho_dx, drho_dy, drho_dz, &
                                dvx_dx, dvx_dy, dvx_dz, dvy_dx, dvy_dy, dvy_dz, &
                                dvz_dx, dvz_dy, dvz_dz, dp_dx, dp_dy, dp_dz, &
                                rho_xtr, vx_xtr, vy_xtr, vz_xtr, p_xtr, dt)

      call update_ghosts(rho_xtr, nx, ny, nz, nGhosts)
      call update_ghosts(vx_xtr, nx, ny, nz, nGhosts)
      call update_ghosts(vy_xtr, nx, ny, nz, nGhosts)
      call update_ghosts(vz_xtr, nx, ny, nz, nGhosts)
      call update_ghosts(p_xtr, nx, ny, nz, nGhosts)

      call reconstructflux(mass_flux_x, momentum_x_flux_x, momentum_y_flux_x, momentum_z_flux_x, &
                           energy_flux_x, drho_dx, dvx_dx, dvy_dx, dvz_dx, dp_dx, &
                           rho_xtr, vx_xtr, vy_xtr, vz_xtr, p_xtr, nx, ny, nz, nGhosts, ds, shiftx)

      call reconstructflux(mass_flux_y, momentum_y_flux_y, momentum_x_flux_y, momentum_z_flux_y, energy_flux_y, &
                           drho_dy, dvy_dy, dvx_dy, dvz_dy, dp_dy, &
                           rho_xtr, vy_xtr, vx_xtr, vz_xtr, p_xtr, nx, ny, nz, nGhosts, ds, shifty)

      call reconstructflux(mass_flux_z, momentum_z_flux_z, momentum_y_flux_z, momentum_x_flux_z, energy_flux_z, &
                           drho_dz, dvz_dz, dvy_dz, dvx_dz, dp_dz, &
                           rho_xtr, vz_xtr, vy_xtr, vx_xtr, p_xtr, nx, ny, nz, nGhosts, ds, shiftz)

      call update_ghosts(mass_flux_x, nx, ny, nz, nGhosts)
      call update_ghosts(mass_flux_y, nx, ny, nz, nGhosts)
      call update_ghosts(mass_flux_z, nx, ny, nz, nGhosts)
      call update_ghosts(momentum_x_flux_x, nx, ny, nz, nGhosts)
      call update_ghosts(momentum_y_flux_x, nx, ny, nz, nGhosts)
      call update_ghosts(momentum_z_flux_x, nx, ny, nz, nGhosts)
      call update_ghosts(momentum_x_flux_y, nx, ny, nz, nGhosts)
      call update_ghosts(momentum_y_flux_y, nx, ny, nz, nGhosts)
      call update_ghosts(momentum_z_flux_y, nx, ny, nz, nGhosts)
      call update_ghosts(momentum_x_flux_z, nx, ny, nz, nGhosts)
      call update_ghosts(momentum_y_flux_z, nx, ny, nz, nGhosts)
      call update_ghosts(momentum_z_flux_z, nx, ny, nz, nGhosts)
      call update_ghosts(energy_flux_x, nx, ny, nz, nGhosts)
      call update_ghosts(energy_flux_y, nx, ny, nz, nGhosts)
      call update_ghosts(energy_flux_z, nx, ny, nz, nGhosts)

      call addfluxes(mass_flux_x, mass_flux_y, mass_flux_z, &
                     momentum_x_flux_x, momentum_x_flux_y, momentum_x_flux_z, &
                     momentum_y_flux_x, momentum_y_flux_y, momentum_y_flux_z, &
                     momentum_z_flux_x, momentum_z_flux_y, momentum_z_flux_z, &
                     energy_flux_x, energy_flux_y, energy_flux_z, &
                     mass, momentum_x, momentum_y, momentum_z, energy, &
                     nx, ny, nz, nGhosts, dt, ds)


      call write_state(timestep, rho, vx, vy, vz, p,&
                                 mass,momentum_x,momentum_y,&
                                 momentum_z,energy,temp)
      print *, "Time=", time, "dt=", dt
      timestep = timestep + 1
      time = time + dt
   end do

   deallocate (rho, vx, vy, vz, p, mass, momentum_x, momentum_y, momentum_z, &
               energy, rho_prime, vx_prime, vy_prime, vz_prime, temp, &
               mass_flux_x, mass_flux_y, mass_flux_z, &
               momentum_x_flux_x, momentum_x_flux_y, momentum_x_flux_z, &
               momentum_y_flux_x, momentum_y_flux_y, momentum_y_flux_z, &
               momentum_z_flux_x, momentum_z_flux_y, momentum_z_flux_z, &
               energy_flux_x, energy_flux_y, energy_flux_z, &
               drho_dx, drho_dy, drho_dz, &
               dvx_dx, dvx_dy, dvx_dz, &
               dvy_dx, dvy_dy, dvy_dz, &
               dvz_dx, dvz_dy, dvz_dz, &
               dp_dx, dp_dy, dp_dz, rho_xtr, vx_xtr, vy_xtr, vz_xtr, p_xtr)
end program euler_cfd
