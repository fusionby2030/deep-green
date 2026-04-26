! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
module physics
   use types_and_kinds
   use global
   use omp_lib
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

   real(rk) function compute_timestep(ds, vx, vy, vz, p, rho)
      real(rk), intent(in) :: ds
      real(rk), intent(in), dimension(:, :, :) :: vx, vy, vz, p, rho
      compute_timestep = cfl*minval(ds/(sqrt(gamma*p/rho) + sqrt(vx**2 + vy**2 + vz**2)))
   end function compute_timestep

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

   subroutine extrapolate_primitves(rho, vx, vy, vz, p, drho_dx, drho_dy, drho_dz, dvx_dx, dvx_dy, &
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
   end subroutine extrapolate_primitves

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
      real(rk) :: c_l, c_r, c_star, h
      real(rk) :: rl, rr, vxl, vxr, vyl, vyr, vzl, vzr, pl, pr
      ! start by calculating rho_star, which is average of density
      !$omp parallel do collapse(3)
      do k = nGhosts + 0, nz - nGhosts
         do j = nGhosts + 0, ny - nGhosts
            do i = nGhosts + 0, nx - nGhosts

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

               en_left = 0.5_rk*(rl*(vxl**2 + vyl**2 + vzl**2)) + (pl/(gamma - 1.0_rk))
               en_right = 0.5_rk*(rr*(vxr**2 + vyr**2 + vzr**2)) + (pr/(gamma - 1.0_rk))

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

               if (offsets(3) == 1) then
                  !momentum_x_flux_x(i, j, k) = -1.0*rho_star*g + momentum_x_star*momentum_x_star/rho_star + p_star;
                  !energy_flux_x(i, j, k) = (en_star + p_star)*(momentum_x_star/rho_star)-0.5*rho_star*(vxl+vxr)*g
                  h = (nz - 2.0*nGhosts)*ds - (k - 2.0)*ds
                  momentum_x_flux_x(i, j, k) = momentum_x_flux_x(i, j, k) + 1.0*rho_star*g*h
                  energy_flux_x(i, j, k) = energy_flux_x(i, j, k) + 0.5*rho_star*(vxl + vxr)*g*h
               end if

               c_l = sqrt(gamma*pl/rl) + abs(vxl)
               c_r = sqrt(gamma*pr/rr) + abs(vxr)
               c_star = max(c_l, c_r)

               mass_flux_x(i, j, k) = momentum_x_star - (0.5_rk*c_star*(rl - rr))
               momentum_x_flux_x(i, j, k) = momentum_x_flux_x(i, j, k) - (c_star*(rl*vxl - rr*vxr))/2.0_rk
               momentum_y_flux_x(i, j, k) = momentum_y_flux_x(i, j, k) - (c_star*(rl*vyl - rr*vyr))/2.0_rk
               momentum_z_flux_x(i, j, k) = momentum_z_flux_x(i, j, k) - (c_star*(rl*vzl - rr*vzr))/2.0_rk
               energy_flux_x(i, j, k) = energy_flux_x(i, j, k) - (c_star*(en_left - en_right))/2.0_rk

               mass_flux_x(i, j, k) = ds*mass_flux_x(i, j, k)
               momentum_x_flux_x(i, j, k) = ds*momentum_x_flux_x(i, j, k)
               momentum_y_flux_x(i, j, k) = ds*momentum_y_flux_x(i, j, k)
               momentum_z_flux_x(i, j, k) = ds*momentum_z_flux_x(i, j, k)
               energy_flux_x(i, j, k) = ds*energy_flux_x(i, j, k)
            end do
         end do
      end do
      !$omp end parallel do
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

      !$omp parallel do collapse(3)
      do k = nGhosts, nz - nGhosts
         do j = nGhosts, ny - nGhosts
            do i = nGhosts, nx - nGhosts
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
      !$omp end parallel do
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
      !$omp parallel do collapse(3)
      do k = nGhosts + 0, nz - nGhosts
         do j = nGhosts + 0, ny - nGhosts
            do i = nGhosts + 0, nx - nGhosts
               dgdy(i, j, k) = (grid(i, j + 1, k) - grid(i, j - 1, k))/(2.0_rk*ds)
               dgdx(i, j, k) = (grid(i + 1, j, k) - grid(i - 1, j, k))/(2.0_rk*ds)
               dgdz(i, j, k) = (grid(i, j, k + 1) - grid(i, j, k - 1))/(2.0_rk*ds)
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine calculate_gradients

   subroutine update_ghosts(grid, nx, ny, nz, nGhosts, BCs)
      use types_and_kinds
      implicit none
      real(rk), dimension(:, :, :), intent(inout) :: grid
      integer(4), intent(in) :: BCs(6)
      integer(ik), intent(in) :: nx, ny, nz, nGhosts
      integer(4):: done_x, done_y, done_z

      done_x = 0
      done_y = 0
      done_z = 0

      if (BCs(1) == PERIODIC .and. BCs(2) == PERIODIC) then
         grid(1, :, :) = grid(nx - nGhosts - 1, :, :)
         grid(2, :, :) = grid(nx - nGhosts, :, :)
         grid(nx - 1, :, :) = grid(3, :, :)
         grid(nx, :, :) = grid(4, :, :)
         done_x = 1
      end if

      if (BCs(3) == PERIODIC .and. BCs(4) == PERIODIC) then
         grid(:, 1, :) = grid(:, ny - nGhosts - 1, :)
         grid(:, 2, :) = grid(:, ny - nGhosts, :)
         grid(:, ny - 1, :) = grid(:, 3, :)
         grid(:, ny, :) = grid(:, 4, :)
         done_y = 1
      end if

      if (BCs(5) == PERIODIC .and. BCs(6) == PERIODIC) then
         grid(:, :, 1) = grid(:, :, nz - nGhosts - 1)
         grid(:, :, 2) = grid(:, :, nz - nGhosts)
         grid(:, :, nz - 1) = grid(:, :, 3)
         grid(:, :, nz) = grid(:, :, 4)
         done_z = 1
      end if

      if (done_x == 1 .and. done_y == 1 .and. done_z == 1) then
         return
      end if

      select case (BCs(1))
      case (OUTFLOW)
         grid(1, :, :) = grid(3, :, :)
         grid(2, :, :) = grid(3, :, :)
      end select

      select case (BCs(2))
      case (OUTFLOW)
         grid(nx - 1, :, :) = grid(nx - nGhosts, :, :)
         grid(nx, :, :) = grid(nx - nGhosts, :, :)
      end select

      select case (BCs(3))
      case (OUTFLOW)
         grid(:, 1, :) = grid(:, 3, :)
         grid(:, 2, :) = grid(:, 3, :)
      end select

      select case (BCs(4))
      case (OUTFLOW)
         grid(:, ny - 1, :) = grid(:, ny - nGhosts, :)
         grid(:, ny, :) = grid(:, ny - nGhosts, :)
      end select

      select case (BCs(5))
      case (OUTFLOW)
         grid(:, :, 1) = grid(:, :, 3)
         grid(:, :, 2) = grid(:, :, 3)
      end select

      select case (BCs(6))
      case (OUTFLOW)
         grid(:, :, nz - 1) = grid(:, :, nz - nGhosts)
         grid(:, :, nz) = grid(:, :, nz - nGhosts)
      end select

   end subroutine update_ghosts

   subroutine update_primitive_ghosts(rho, vx, vy, vz, p, nx, ny, nz, nGhosts, BCs)
      use types_and_kinds
      implicit none
      real(rk), dimension(:, :, :), intent(inout) :: rho, vx, vy, vz, p
      integer(ik), intent(in) :: nx, ny, nz, nGhosts
      integer(4), intent(in) :: BCs(6)
      call update_ghosts(rho, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(vx, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(vy, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(vz, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(p, nx, ny, nz, nGhosts, BCs)
      call apply_boundary_conditions(rho, vx, vy, vz, p, nx, ny, nz, nGhosts, BCs)
   end subroutine update_primitive_ghosts

   subroutine update_gradient_ghost(drho_dx, drho_dy, drho_dz, &
                                    dvx_dx, dvx_dy, dvx_dz, dvy_dx, dvy_dy, dvy_dz, dvz_dx, dvz_dy, dvz_dz, &
                                    dp_dx, dp_dy, dp_dz, nx, ny, nz, nGhosts, BCs)
      use types_and_kinds
      implicit none
      real(rk), dimension(:, :, :), intent(inout) :: drho_dx, drho_dy, drho_dz, &
                                                     dvx_dx, dvx_dy, dvx_dz, dvy_dx, dvy_dy, dvy_dz, dvz_dx, dvz_dy, dvz_dz, &
                                                     dp_dx, dp_dy, dp_dz
      integer(ik), intent(in) :: nx, ny, nz, nGhosts
      integer(4), intent(in) :: BCs(6)
      call update_ghosts(drho_dx, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(drho_dy, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(drho_dz, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(dvx_dx, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(dvx_dy, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(dvx_dz, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(dvy_dx, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(dvy_dy, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(dvy_dz, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(dvz_dx, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(dvz_dy, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(dvz_dz, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(dp_dx, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(dp_dy, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(dp_dz, nx, ny, nz, nGhosts, BCs)
   end subroutine update_gradient_ghost

   subroutine update_flux_ghost(mass_flux_x, mass_flux_y, mass_flux_z, &
                                momentum_x_flux_x, momentum_y_flux_x, momentum_z_flux_x, &
                                momentum_x_flux_y, momentum_y_flux_y, momentum_z_flux_y, &
                                momentum_x_flux_z, momentum_y_flux_z, momentum_z_flux_z, &
                                energy_flux_x, energy_flux_y, energy_flux_z, nx, ny, nz, nGhosts, BCs)
      use types_and_kinds
      implicit none
      real(rk), dimension(:, :, :), intent(inout) :: mass_flux_x, mass_flux_y, mass_flux_z, &
                                                     momentum_x_flux_x, momentum_y_flux_x, momentum_z_flux_x, &
                                                     momentum_x_flux_y, momentum_y_flux_y, momentum_z_flux_y, &
                                                     momentum_x_flux_z, momentum_y_flux_z, momentum_z_flux_z, &
                                                     energy_flux_x, energy_flux_y, energy_flux_z
      integer(ik), intent(in) :: nx, ny, nz, nGhosts
      integer(4), intent(in) :: BCs(6)
      call update_ghosts(mass_flux_x, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(mass_flux_y, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(mass_flux_z, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(momentum_x_flux_x, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(momentum_y_flux_x, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(momentum_z_flux_x, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(momentum_x_flux_y, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(momentum_y_flux_y, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(momentum_z_flux_y, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(momentum_x_flux_z, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(momentum_y_flux_z, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(momentum_z_flux_z, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(energy_flux_x, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(energy_flux_y, nx, ny, nz, nGhosts, BCs)
      call update_ghosts(energy_flux_z, nx, ny, nz, nGhosts, BCs)
   end subroutine update_flux_ghost

   !Implements slip wall Boundaries
   subroutine apply_boundary_conditions(rho, vx, vy, vz, p, nx, ny, nz, nGhosts, BCs)
      use types_and_kinds
      implicit none
      real(rk), dimension(:, :, :), intent(inout) :: rho, vx, vy, vz, p
      integer(ik), intent(in) :: nx, ny, nz, nGhosts
      integer(4), intent(in) :: BCs(6)

      select case (BCs(1))
      case (WALL)
         rho(2, :, :) = rho(3, :, :)
         vx(2, :, :) = -vx(3, :, :)
         vy(2, :, :) = vy(3, :, :)
         vz(2, :, :) = vz(3, :, :)
         p(2, :, :) = p(3, :, :)

         rho(1, :, :) = rho(4, :, :)
         vx(1, :, :) = -vx(4, :, :)
         vy(1, :, :) = vy(4, :, :)
         vz(1, :, :) = vz(4, :, :)
         p(1, :, :) = p(4, :, :)
      end select

      select case (BCs(2))
      case (WALL)
         rho(NX - 1, :, :) = rho(NX - 2, :, :)
         vx(NX - 1, :, :) = -vx(NX - 2, :, :)
         vy(NX - 1, :, :) = vy(NX - 2, :, :)
         vz(NX - 1, :, :) = vz(NX - 2, :, :)
         p(NX - 1, :, :) = p(NX - 2, :, :)

         rho(NX, :, :) = rho(NX - 3, :, :)
         vx(NX, :, :) = -vx(NX - 3, :, :)
         vy(NX, :, :) = vy(NX - 3, :, :)
         vz(NX, :, :) = vz(NX - 3, :, :)
         p(NX, :, :) = p(NX - 3, :, :)
      end select

      select case (BCs(3))
      case (WALL)
         rho(:, 2, :) = rho(:, 3, :)
         vx(:, 2, :) = vx(:, 3, :)
         vy(:, 2, :) = -vy(:, 3, :)
         vz(:, 2, :) = vz(:, 3, :)
         p(:, 2, :) = p(:, 3, :)

         rho(:, 1, :) = rho(:, 4, :)
         vx(:, 1, :) = vx(:, 4, :)
         vy(:, 1, :) = -vy(:, 4, :)
         vz(:, 1, :) = vz(:, 4, :)
         p(:, 1, :) = p(:, 4, :)
      end select

      select case (BCs(4))
      case (WALL)
         rho(:, NY - 1, :) = rho(:, NY - 2, :)
         vx(:, NY - 1, :) = vx(:, NY - 2, :)
         vy(:, NY - 1, :) = -vy(:, NY - 2, :)
         vz(:, NY - 1, :) = vz(:, NY - 2, :)
         p(:, NY - 1, :) = p(:, NY - 2, :)

         rho(:, NY, :) = rho(:, NY - 3, :)
         vx(:, NY, :) = vx(:, NY - 3, :)
         vy(:, NY, :) = -vy(:, NY - 3, :)
         vz(:, NY, :) = vz(:, NY - 3, :)
         p(:, NY, :) = p(:, NY - 3, :)
      end select

      select case (BCs(5))
      case (WALL)
         rho(:, :, 2) = rho(:, :, 3)
         vx(:, :, 2) = vx(:, :, 3)
         vy(:, :, 2) = vy(:, :, 3)
         vz(:, :, 2) = -vz(:, :, 3)
         p(:, :, 2) = p(:, :, 3)

         rho(:, :, 1) = rho(:, :, 4)
         vx(:, :, 1) = vx(:, :, 4)
         vy(:, :, 1) = vy(:, :, 4)
         vz(:, :, 1) = -vz(:, :, 4)
         p(:, :, 1) = p(:, :, 4)
      end select

      select case (BCs(6))
      case (WALL)
         rho(:, :, NZ - 1) = rho(:, :, NZ - 2)
         vx(:, :, NZ - 1) = vx(:, :, NZ - 2)
         vy(:, :, NZ - 1) = vy(:, :, NZ - 2)
         vz(:, :, NZ - 1) = -vz(:, :, NZ - 2)
         p(:, :, NZ - 1) = p(:, :, NZ - 2)

         rho(:, :, NZ) = rho(:, :, NZ - 3)
         vx(:, :, NZ) = vx(:, :, NZ - 3)
         vy(:, :, NZ) = vy(:, :, NZ - 3)
         vz(:, :, NZ) = -vz(:, :, NZ - 3)
         p(:, :, NZ) = p(:, :, NZ - 3)
      end select
   end subroutine apply_boundary_conditions
end module physics
