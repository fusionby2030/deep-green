!---------------------------------------------------------------------------------------
! Authors: Adam Kit and Kostis Papadakis (2023)
! Program: euler_cfd
! Solving the compressible Euler equations in 3 dimensions using the Lax–Friedrichs
! flux method. Stable for CFL<0.4.
! TODO:
! --Add a flux limiter using a high order reconstruction :
! probably use Lax–Wendroff for F_{h}
! --Supported Boundaries: Periodic, Slip Wall, Outflow
!
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
!--------------------------------------------------------------------------------------
MODULE types_and_kinds 
    IMPLICIT NONE 
    INTEGER, PARAMETER :: rk = selected_real_kind(4) ! p
    INTEGER, PARAMETER :: ik = selected_int_kind(4)
END MODULE types_and_kinds
MODULE global
    use types_and_kinds
    IMPLICIT NONE 
    REAL(RK), PARAMETER :: c_to_k = 273.15_rk, k_to_c = -273.15_rk
    REAL(rk), PARAMETER :: k_air = 0.024 ! W/(m*K)
    REAL(rk), PARAMETER :: rho_air = 1.2922 ! kg/m^3
    REAL(rk), PARAMETER :: cp_air = 1003.5 ! J/(kg*K)
    REAL(rk), PARAMETER :: alpha_air=1.8408476735367102e-05_rk
    ! REAL(rk), PARAMETER :: alpha_air = k_air / ( rho_air * cp_air )
    REAL(rk), PARAMETER :: thermal_condctivity_plastic = 0.23   ! W/(m*K)
    REAL(RK), PARAMETER :: density_plastic = 1200.0_rk          ! kg/m^3
    REAL(RK), PARAMETER :: specific_heat_plastic = 1200.0_rk    ! J/(kg*K)
    REAL(rk), PARAMETER :: thermal_diffusivity_plastic=1.5972222222222223e-07_rk
    ! REAL(rk), PARAMETER :: thermal_diffusivity_plastic = thermal_condctivity_plastic / (density_plastic * specific_heat_plastic) ! m^3/s
    REAL(RK), PARAMETER :: thermal_conductivity_wood = 0.15_rk ! W/(m*K)
    REAL(RK), PARAMETER :: density_wood = 640.0_rk             ! kg/m^3
    REAL(RK), PARAMETER :: specific_heat_wood = 2805.0_rk      ! J/(kg*K)
    ! REAL(RK), PARAMETER :: thermal_diffusivity_wood = 8.355614973262032e-08_rk
    REAL(RK), PARAMETER :: thermal_diffusivity_wood = thermal_conductivity_wood / (density_wood * specific_heat_wood) ! m^3/s
    REAL(RK), PARAMETER :: wood_thickness = 0.02_rk            ! m
    real(rk), parameter :: pi = 2.d0*dasin(1.d0)
    real(rk), parameter :: gamma = 5.0_rk/3.0_rk, rs = 287.0_rk, cfl = 0.2_rk, g = -9.81_rk
    integer(rk), parameter:: N_VARS = 5, N_FLUX = 15, N_DIMS = 3
    integer(4), parameter:: PERIODIC = 0, OUTFLOW = 1, WALL = 2, INFLOW = 3
END MODULE global
MODULE physics 
    use types_and_kinds
    use global 
    use omp_lib
    IMPLICIT NONE
    CONTAINS
    !-----------------------------------------------------------------------
    ! SUBROUTINE: get_diffusivity_grid
    ! Purpose: Diffusivity grid for the greenhouse, so we do not need to branch in main loop
    ! The diffusivity grid is the same size as the temperature grid
    ! The bottom z-layer is wood 
    ! the rest is plastic 
    !-----------------------------------------------------------------------
    SUBROUTINE get_diffusivity_grid(diffusivity_grid, Nx, Ny, Nz, nGhosts)
        REAL(RK), DIMENSION(:, :, :), INTENT(INOUT) :: diffusivity_grid
        INTEGER(IK), INTENT(IN) :: Nx, Ny, Nz, nGhosts
        diffusivity_grid = alpha_air
        ! plastic border
        diffusivity_grid(nGhosts+1:Nx+1:Nx-1, nGhosts+1:Ny+1, nGhosts+1:Nz+1) = thermal_diffusivity_plastic
        diffusivity_grid(nGhosts+1:Nx+1, nGhosts+1:Ny+1:Ny-1, nGhosts+1:Nz+1) = thermal_diffusivity_plastic
        diffusivity_grid(nGhosts+1:Nx+1, nGhosts+1:Ny+1, nGhosts+1:Nz+1:Nz-1) = thermal_diffusivity_plastic
        ! replace small component where computers are with air 
        diffusivity_grid(1+nGhosts:Nx/2, 1:1+nGhosts, :) = alpha_air
        ! replace bottom z layer as wood 
        diffusivity_grid(nGhosts+1:Nx+1, nGhosts+1:Ny+1, nGhosts+1) = thermal_diffusivity_wood
    END SUBROUTINE get_diffusivity_grid

    !-----------------------------------------------------------------------
    ! SUBROUTINE: get_wall_thickness_grid 
    ! Purpose: Wall thickness grid for the greenhouse, so we do not need to branch in main loop
    ! The wall thickness grid is the same size as the temperature grid
    ! The bottom z-layer is wood, i.e., 0.2 m 
    ! the border is plastic 
    ! the rest is ds, or the grid size
    !-----------------------------------------------------------------------
    SUBROUTINE get_wall_thickness_grid(wall_thickness_grid, Nx, Ny, Nz, nGhosts, ds, wall_thickness)
        REAL(RK), DIMENSION(:, :, :), INTENT(INOUT) :: wall_thickness_grid
        INTEGER(IK), INTENT(IN) :: Nx, Ny, Nz, nGhosts
        REAL(RK), INTENT(IN) :: wall_thickness, ds
        wall_thickness_grid = ds
        wall_thickness_grid(nGhosts+1:Nx+1:Nx-1, nGhosts+1:Ny+1, nGhosts+1:Nz+1) = wall_thickness
        wall_thickness_grid(nGhosts+1:Nx+1, nGhosts+1:Ny+1:Ny-1, nGhosts+1:Nz+1) = wall_thickness
        wall_thickness_grid(nGhosts+1:Nx+1, nGhosts+1:Ny+1, nGhosts+1:Nz+1:Nz-1) = wall_thickness
        ! replace small component where computers are with ds 
        wall_thickness_grid(1+nGhosts:Nx/2, 1:1+nGhosts, :) = ds
        ! replace bottom z layer as wood 
        wall_thickness_grid(nGhosts+1:Nx+1, nGhosts+1:Ny+1, nGhosts+1) = wood_thickness
    END SUBROUTINE get_wall_thickness_grid
    !-----------------------------------------------------------------------
    ! SUBROUTINE: get_source_grid
    ! Purpose: Calculate the source grid for the heat equation
    ! Sources should be in W/m^3
    !-----------------------------------------------------------------------
    SUBROUTINE get_source_grid(sources, Nx, Ny, Nz, ds, nGhosts, computer_power)
        REAL(RK), DIMENSION(:, :, :), INTENT(INOUT) :: sources
        INTEGER(IK), INTENT(IN) :: Nx, Ny, Nz, nGhosts
        REAL(RK), INTENT(IN) :: ds, computer_power
        INTEGER(IK) :: i, j, k
        INTEGER(IK), PARAMETER :: lamp_x = 3, lamp_y = 12, lamp_z = 30
        INTEGER(RK), PARAMETER :: lamp_len_x=8, lamp_len_y=12, lamp_len_z = 5
        REAL(RK) :: lamp_power=500.0_rk, lamp_power_density        
        REAL(RK) :: computer_power_density_per_grid
        lamp_power_density = lamp_power / (real(lamp_len_x*lamp_len_y*lamp_len_z, rk)*ds**3)
        computer_power_density_per_grid = computer_power / ((Ny/2)*ds**3)
        sources = 0.0_rk
        sources(1+nGhosts:Nx/2, 1:1+nGhosts, nGhosts+1:) = computer_power_density_per_grid
        DO k=1, Nz+nGhosts
            DO j=1, Ny+nGhosts
                DO i=1, Nx+nGhosts
                    if (i >= lamp_x .AND. i<lamp_x + lamp_len_x .AND.& 
                        j >= lamp_y .AND. j<lamp_y + lamp_len_y .AND. & 
                        k >= lamp_z .AND. k<lamp_z + lamp_len_z) then
                            sources(i, j, k) =  sources(i, j, k) + lamp_power_density
                    END IF 
                END DO 
            END DO 
        END DO 
    END SUBROUTINE get_source_grid
    !-----------------------------------------------------------------------
    ! SUBROUTINE: conductive_heat_transfer
    ! Purpose: Calculate the new temperature of grid using the heat equation
    ! Direchlet boundary conditions are used for the walls, i.e., ghost points not updated
    !-----------------------------------------------------------------------
    SUBROUTINE conductive_heat_transfer(T, T_new, Nx, Ny, Nz, nGhosts, dt, ds, sources, wall_thickness_grid, diffusivity_grid)
        REAL(RK), DIMENSION(:, :, :), INTENT(IN) :: T, sources, diffusivity_grid, wall_thickness_grid
        REAL(RK), DIMENSION(:, :, :), INTENT(INOUT) :: T_new
        INTEGER(IK), INTENT(IN) :: Nx, Ny, Nz, nGhosts
        REAL(RK), INTENT(IN) :: dt, ds
        INTEGER :: i, j, k
        DO k = 1+Nghosts, Nz+nGhosts
            DO j = 1+Nghosts, Ny+nGhosts
                DO i = 1+Nghosts, Nx+nGhosts
                    T_new(i, j, k) = T(i, j, k) + (diffusivity_grid(i, j, k)*dt)*( & 
                                                    (T(i+1, j, k) - 2*T(i, j, k) + T(i-1, j, k) + &
                                                   T(i, j+1, k) - 2*T(i, j, k) + T(i, j-1, k) + &
                                                   T(i, j, k+1) - 2*T(i, j, k) + T(i, j, k-1))/(wall_thickness_grid(i,j,k)**2) + &
                                                    sources(i, j, k)* (ds**2/ (k_air)))
                END DO
            END DO
        END DO
    END SUBROUTINE conductive_heat_transfer
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
      do k = nGhosts + 0_ik, nz - nGhosts
         do j = nGhosts + 0_ik, ny - nGhosts
            do i = nGhosts + 0_ik, nx - nGhosts

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
      vanleerlimiter = (sign(1.0_rk, a) + sign(1.0_rk, b))*min(abs(a), abs(b))/(abs(a) + abs(b) + 1.0d-30)
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
SUBROUTINE run_conductive(T, Nx, Ny, Nz, nGhosts, dt, ds, sources, wall_thickness, t_max, converged)
    use types_and_kinds
    use global 
    ! use initialization 
    use physics 
    IMPLICIT NONE 
    REAL(rk), INTENT(IN) :: sources(:, :, :)
    INTEGER(ik), INTENT(IN) :: Nx, Ny, Nz, nGhosts
    REAL(rk), INTENT(INOUT), DIMENSION(:, :, :) :: T
    REAL(rk), ALLOCATABLE, DIMENSION(:, :, :) :: T_new, ds_grid, diffusivity_grid
    REAL(rk), INTENT(IN) :: dt, ds, wall_thickness, t_max
    INTEGER(ik), INTENT(INOUT) :: converged
    INTEGER :: timestep=0
    REAL :: epsilon = 0.000001
    allocate(T_new(Nx+2*nGhosts, Ny+2*nGhosts, Nz+2*nGhosts),  & 
            diffusivity_grid(Nx+2*nGhosts, Ny+2*nGhosts, Nz+2*nGhosts), & 
            ds_grid(Nx+2*nGhosts, Ny+2*nGhosts, Nz+2*nGhosts))
    
    ! call init_temperature_grid(T, Nx, Ny, Nz, nGhosts, outside_temperature)
    ! CALL get_source_grid(sources, Nx, Ny, Nz, dx, nGhosts, computer_power)
    CALL get_diffusivity_grid(diffusivity_grid, Nx, Ny, Nz, nGhosts)
    call get_wall_thickness_grid(ds_grid, Nx, Ny, Nz, nGhosts, ds, wall_thickness)
    print *, size(diffusivity_grid), size(ds_grid)
    T = T + c_to_k
    T_new = T 
    do while (timestep*dt < t_max) 
        timestep = timestep + 1
        CALL conductive_heat_transfer(T, T_new, Nx, Ny, Nz, nGhosts, dt, ds, sources, &
                             ds_grid, diffusivity_grid)
        if (maxval(abs(T_new-T)) < epsilon) then 
            converged = 1_ik 
            print *, "converged, exiting"
            exit 
        end if
        T = T_new
        timestep = timestep + 1 
        converged = 0_ik
    end do 
    T = T + k_to_c
    print *, "Finised Conductive", maxval(T)
    print *, maxval(sources)
    deallocate(T_new, ds_grid, diffusivity_grid)
END SUBROUTINE run_conductive


subroutine euler_driver(rho,vx,vy,vz,p,xcells,ycells,zcells,nGhosts,ds,time,time_max)
    use types_and_kinds
    use global
    use physics
    implicit none
    integer(ik),intent(in) :: xcells , ycells , zcells, nGhosts
    real(rk), dimension(:, :, :), intent(inout) :: rho,vx,vy,vz,p
    real(rk),intent(in) :: ds 
    real(rk),intent(inout) :: time 
    real(rk),intent(in)::time_max 
    integer(4):: BCs(6)
    ! integer(4),intent(in):: BCs(6)

    real(rk), dimension(:, :, :), allocatable  :: &
          mass, momentum_x, momentum_y, momentum_z, energy, &
          rho_prime, vx_prime, vy_prime, vz_prime, temp, &
          mass_flux_x, mass_flux_y, mass_flux_z, &
          momentum_x_flux_x, momentum_x_flux_y, momentum_x_flux_z, &
          momentum_y_flux_x, momentum_y_flux_y, momentum_y_flux_z, &
          momentum_z_flux_x, momentum_z_flux_y, momentum_z_flux_z, &
          energy_flux_x, energy_flux_y, energy_flux_z, &
          drho_dx, drho_dy, drho_dz,dvx_dx, dvx_dy, dvx_dz, &
          dvy_dx, dvy_dy, dvy_dz,dvz_dx, dvz_dy, dvz_dz, &
          dp_dx, dp_dy, dp_dz, rho_xtr,vx_xtr, vy_xtr, vz_xtr, p_xtr
 
    real(rk)::dt
    integer(ik):: nx ,ny , nz
    integer(ik) :: timestep = 0
    integer(ik) :: shiftx(3), shifty(3), shiftz(3)
    
    BCs(1) = WALL !OUTFLOW !x-
    BCs(2) = WALL !x+
    BCs(3) = WALL !y-
    BCs(4) = WALL !y-
    BCs(5) = WALL !z-
    BCs(6) = WALL !z+
 
    nx = xcells + 2_ik**nGhosts
    ny = ycells + 2_ik*nGhosts
    nz = zcells + 2_ik*nGhosts
 
 
    allocate ( mass(nx, ny, nz), &
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
 
 
    call update_primitive_ghosts(rho, vx, vy, vz, p, nx, ny, nz, nGhosts, BCs)
    call conservative(mass, momentum_x, momentum_y, momentum_z, energy, rho, p, vx, vy, vz, temp, ds)
    do while (time <= time_max)
 
       call apply_boundary_conditions(rho, vx, vy, vz, p, nx, ny, nz, nGhosts, BCs)
 
       call update_primitive_ghosts(rho, vx, vy, vz, p, nx, ny, nz, nGhosts, BCs)
 
       call primitive(mass, momentum_x, momentum_y, momentum_z, energy, rho, p, vx, vy, vz, temp, ds)
 
       call update_primitive_ghosts(rho, vx, vy, vz, p, nx, ny, nz, nGhosts, BCs)
 
       dt = compute_timestep(ds, vx, vy, vz, p, rho)
 
       call calculate_gradients(rho, drho_dx, drho_dy, drho_dz, nx, ny, nz, nGhosts, ds)
       call calculate_gradients(vx, dvx_dx, dvx_dy, dvx_dz, nx, ny, nz, nGhosts, ds)
       call calculate_gradients(vy, dvy_dx, dvy_dy, dvy_dz, nx, ny, nz, nGhosts, ds)
       call calculate_gradients(vz, dvz_dx, dvz_dy, dvz_dz, nx, ny, nz, nGhosts, ds)
       call calculate_gradients(p, dp_dx, dp_dy, dp_dz, nx, ny, nz, nGhosts, ds)
 
       call update_gradient_ghost(drho_dx, drho_dy, drho_dz, dvx_dx, dvx_dy, dvx_dz, dvy_dx, &
                                  dvy_dy, dvy_dz, dvz_dx, dvz_dy, dvz_dz, dp_dx, dp_dy, dp_dz, nx, ny, nz, nGhosts, BCs)
 
       call extrapolate_primitves(rho, vx, vy, vz, p, drho_dx, drho_dy, drho_dz, &
                                  dvx_dx, dvx_dy, dvx_dz, dvy_dx, dvy_dy, dvy_dz, &
                                  dvz_dx, dvz_dy, dvz_dz, dp_dx, dp_dy, dp_dz, &
                                  rho_xtr, vx_xtr, vy_xtr, vz_xtr, p_xtr, dt)
 
       !Apply BCs if we have anything other than Periodic or Outflow!
       call apply_boundary_conditions(rho_xtr, vx_xtr, vy_xtr, vz_xtr, p_xtr, nx, ny, nz, nGhosts, BCs)
 
       call reconstructflux(mass_flux_x, momentum_x_flux_x, momentum_y_flux_x, momentum_z_flux_x, &
                            energy_flux_x, drho_dx, dvx_dx, dvy_dx, dvz_dx, dp_dx, &
                            rho_xtr, vx_xtr, vy_xtr, vz_xtr, p_xtr, nx, ny, nz, nGhosts, ds, shiftx)
 
       call reconstructflux(mass_flux_y, momentum_y_flux_y, momentum_x_flux_y, momentum_z_flux_y, energy_flux_y, &
                            drho_dy, dvy_dy, dvx_dy, dvz_dy, dp_dy, &
                            rho_xtr, vy_xtr, vx_xtr, vz_xtr, p_xtr, nx, ny, nz, nGhosts, ds, shifty)
 
       call reconstructflux(mass_flux_z, momentum_z_flux_z, momentum_x_flux_z, momentum_y_flux_z, energy_flux_z, &
                            drho_dz, dvz_dz, dvx_dz, dvy_dz, dp_dz, &
                            rho_xtr, vz_xtr, vx_xtr, vy_xtr, p_xtr, nx, ny, nz, nGhosts, ds, shiftz)
 
       call update_flux_ghost(mass_flux_x, mass_flux_y, mass_flux_z, momentum_x_flux_x, momentum_y_flux_x, momentum_z_flux_x, &
                              momentum_x_flux_y, momentum_y_flux_y, momentum_z_flux_y, momentum_x_flux_z, momentum_y_flux_z, &
                              momentum_z_flux_z, energy_flux_x, energy_flux_y, energy_flux_z, nx, ny, nz, nGhosts, BCs)
 
       call addfluxes(mass_flux_x, mass_flux_y, mass_flux_z, &
                      momentum_x_flux_x, momentum_x_flux_y, momentum_x_flux_z, &
                      momentum_y_flux_x, momentum_y_flux_y, momentum_y_flux_z, &
                      momentum_z_flux_x, momentum_z_flux_y, momentum_z_flux_z, &
                      energy_flux_x, energy_flux_y, energy_flux_z, &
                      mass, momentum_x, momentum_y, momentum_z, energy, &
                      nx, ny, nz, nGhosts, dt, ds)
       time = time + dt
       timestep = timestep + 1_ik
       if (mod(timestep, 1000) == 0) then
          print *, time, dt, timestep
       end if
       ! print *, time, dt, timestep
    end do
 
    deallocate (mass, momentum_x, momentum_y, momentum_z, &
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
 
 end subroutine euler_driver