!-----------------------------------------------------------------------
! Program: euler_cfd
! Purpose: Main program for the Euler CFD solver
! 4 steps:
! 1. Initialize grids
! 2. Compute mass, momentum, and energy
! 3. Time step from CFL condition
! 4. Update primitive variables
! 5. Compute fluxes
! 6. Update conserved variables with fluxes
!-----------------------------------------------------------------------

MODULE types_and_kinds
   IMPLICIT NONE
   INTEGER, PARAMETER :: rk = selected_real_kind(8) ! p
   INTEGER, PARAMETER :: ik = selected_int_kind(8)
END MODULE types_and_kinds

MODULE global
   use types_and_kinds
   IMPLICIT NONE
   REAL(RK), PARAMETER :: gamma = 5./3., Rs = 287.0_rk, cfl = 0.1_rk ! rydeberg cons

END MODULE global

MODULE initialization
   use types_and_kinds
   IMPLICIT NONE
CONTAINS
   SUBROUTINE init_grid(grid, Nx, Ny, Nz, nGhosts,init_value)
      REAL(RK), dimension(:, :, :), intent(out) :: grid
      INTEGER(ik), intent(in) :: Nx, Ny, Nz,nGhosts
      REAL(RK), intent(in) :: init_value
      INTEGER(ik) :: i, j, k
      ! to efficintly traverse the grid in the order of memory layout
      ! loop over 3D array with proper cache efficiency
      DO k = nGhosts+1, Nz-nGhosts
         DO j = nGhosts+1, Ny-nGhosts
            DO i = nGhosts+1, Nx-nGhosts
               grid(i, j, k) = init_value
            END DO
         END DO
      END DO
   END SUBROUTINE init_grid
   ! similar to above subroutine, we now intialize the grid with a 3D gaussian centered in middle of grid
   REAL(RK) FUNCTION GAUSSIAN3D(i, j, k, mu, sigma)
      REAL(RK), intent(in) :: mu, sigma
      INTEGER(ik), intent(in) :: i, j, k
      REAL(RK) :: A = 1.0
      !GAUSSIAN3D = A*sin(  2*DACOS(-1.D0)* ((i-3)/127.0) )
      GAUSSIAN3D = A*exp(-((i - mu)**2 + (j - mu)**2 + (k - mu)**2)/(2*sigma**2))
   END FUNCTION GAUSSIAN3D
   ! use the above function to initialize the grid
   SUBROUTINE init_grid_gaussian(grid, Nx, Ny, Nz,nGhosts, init_value_mu, init_value_sigma, magnitude, offset)
      REAL(RK), dimension(:, :, :), intent(out) :: grid
      INTEGER(ik), intent(in) :: Nx, Ny, Nz,nGhosts
      REAL(RK), intent(in) :: init_value_mu, init_value_sigma, magnitude, offset
      INTEGER(ik) :: i, j, k
      DO k = nGhosts+1, Nz-nGhosts
         DO j = nGhosts+1, Ny-nGhosts
            DO i = nGhosts+1, Nx-nGhosts
               grid(i, j, k) = magnitude*GAUSSIAN3D(i, j, k, init_value_mu, init_value_sigma) + offset
            END DO
         END DO
      END DO
   END SUBROUTINE init_grid_gaussian
END MODULE initialization

MODULE io
   use types_and_kinds
   IMPLICIT NONE
CONTAINS
   SUBROUTINE write_grid(grid, filename)
      REAL(RK), dimension(:, :, :), intent(in) :: grid
      CHARACTER(len=*), intent(in) :: filename
      ! need to use stream as access for the binary (WHY?)
      ! maybe retarded but name can not be *.bin
      open (1, file=filename, form='unformatted', access='stream', status='replace')
      write (1) grid
      close (1)
   END SUBROUTINE write_grid
END MODULE io

MODULE physics
   use types_and_kinds
   use global
   IMPLICIT NONE
CONTAINS
   SUBROUTINE conservative(mass, momentum_x, momentum_y, momentum_z, energy, rho, p, vx, vy, vz, temp, ds)
      REAL(RK), dimension(:, :, :), intent(inout) :: mass, momentum_x, momentum_y, momentum_z, energy, temp
      REAL(RK), dimension(:, :, :), intent(in) :: rho, vx, vy, vz, p
      REAL(RK), intent(in) :: ds
      REAL(RK) :: cell_volume
      cell_volume = ds**3
      mass = rho*cell_volume
      momentum_x = rho*vx*cell_volume
      momentum_y = rho*vy*cell_volume
      momentum_z = rho*vz*cell_volume
      energy = cell_volume*(rho*(vx**2 + vy**2 + vz**2)/2.0_rk + p/(gamma - 1.0_rk))
      temp = p/(Rs*rho) - 273.15_rk
   END SUBROUTINE conservative
   REAL(RK) FUNCTION get_timestep(ds, vx, vy, vz, p, rho)
      REAL(RK), intent(in) :: ds
      REAL(RK), intent(in), dimension(:, :, :) :: vx, vy, vz, p, rho
      get_timestep = minval(cfl*ds/sqrt((gamma*p/rho) + (vx**2 + vy**2 + vz**2)))
      ! TODO: DOUBLE CHECK THIS
      ! get_timestep = minval(cfl*ds/(sqrt(gamma*p/rho) + sqrt(vx**2 + vy**2 + vz**2)) )
   END FUNCTION get_timestep
   SUBROUTINE primitive(mass, momentum_x, momentum_y, momentum_z, energy, rho, p, vx, vy, vz, temp, ds)
      REAL(RK), dimension(:, :, :), intent(in) :: mass, momentum_x, momentum_y, momentum_z, energy
      REAL(RK), dimension(:, :, :), intent(inout) :: p, vx, vy, vz, rho, temp
      REAL(RK), intent(in) :: ds
      REAL(RK) :: cell_volume
      cell_volume = ds**3
      rho = mass/cell_volume
      vx = momentum_x/(rho*cell_volume)
      vy = momentum_y/(rho*cell_volume)
      vz = momentum_z/(rho*cell_volume)
      p = (energy/cell_volume - 0.5_rk*rho*(vx*vx + vy*vy + vz*vz))*(gamma - 1.0_rk)
      temp = p/(Rs*rho) - 273.15_rk
   END SUBROUTINE primitive

   SUBROUTINE reconstructFlux(mass_flux_x, momentum_x_flux_x, momentum_y_flux_x, &
                                momentum_z_flux_x, energy_flux_x, &
                                drho, dvx, dvy, dvz, dp, rho, vx, vy, vz, p, Nx, Ny, Nz,nGhosts, ds,offsets)
  REAL(RK), dimension(:, :, :), intent(inout) :: mass_flux_x, momentum_x_flux_x, momentum_y_flux_x, momentum_z_flux_x, energy_flux_x
      REAL(RK), dimension(:, :, :), intent(in) :: rho, vx, vy, vz, p, drho, dvx, dvy, dvz, dp
      INTEGER(IK), intent(in) :: offsets(3)
      INTEGER(IK), intent(in) :: Nx, Ny, Nz,nGhosts
      REAL(RK), intent(in) :: ds
      INTEGER(IK) :: i, j, k
      REAL(RK) :: rho_star, momentum_x_star, momentum_y_star, momentum_z_star, p_star, en_right, en_left, en_star
      REAL(RK) :: C_L, C_R, C_star
      REAL(RK) :: RL, RR, VXL, VXR, VYL, VYR, VZL, VZR, PL, PR
      ! start by calculating rho_star, which is average of density
      DO k = nGhosts+1, Nz-nGhosts
         DO j = nGhosts+1, Ny-nGhosts
            DO i = nGhosts+1, Nx-nGhosts

               RL = rho(i + offsets(1), j+ offsets(2), k+ offsets(3)) - (drho(i + offsets(1), j+ offsets(2), k+ offsets(3)))*(ds/2.)
               RR = rho(i, j, k) + (drho(i, j, k))*(ds/2.)

               VXL = vx(i + offsets(1), j+ offsets(2), k+ offsets(3)) - (dvx(i + offsets(1), j+ offsets(2), k+ offsets(3)))*(ds/2.)
               VXR = vx(i, j, k) + (dvx(i, j, k))*(ds/2.)

               VYL = vy(i + offsets(1), j+ offsets(2), k+ offsets(3)) - (dvy(i + offsets(1), j+ offsets(2), k+ offsets(3)))*(ds/2.)
               VYR = vy(i, j, k) + (dvy(i, j, k))*(ds/2.)

               VZL = vz(i + offsets(1), j+ offsets(2), k+ offsets(3)) - (dvz(i + offsets(1), j+ offsets(2), k+ offsets(3)))*(ds/2.)
               VZR = vz(i, j, k) + (dvz(i, j, k))*(ds/2.)

               PL = p(i + offsets(1), j+ offsets(2), k+ offsets(3)) - (dp(i + offsets(1), j+ offsets(2), k+ offsets(3)))*(ds/2.)
               PR = p(i, j, k) + (dp(i, j, k))*(ds/2.)

               en_left = (PL/(gamma - 1.0_rk)) + 0.5_rk*(RL*(VXL**2 + VYL**2 + VZL**2))
               en_right = (PR/(gamma - 1.0_rk)) + 0.5_rk*(RR*(VXR**2 + VYR**2 + VZR**2))

               rho_star = (RL + RR)/2.0_rk
               momentum_x_star = (RL*VXL + RR*VXR)/2.0_rk
               momentum_y_star = (RL*VYL + RR*VYR)/2.0_rk
               momentum_z_star = (RL*VZL + RR*VZR)/2.0_rk
               en_star = 0.5_rk*(en_left + en_right)
               p_star = (gamma - 1.0_rk)*(en_star - 0.5_rk*(momentum_x_star**2 + momentum_y_star**2 + momentum_z_star**2)/rho_star)

               mass_flux_x(i, j, k) = momentum_x_star
               momentum_x_flux_x(i, j, k) = momentum_x_star*momentum_x_star/rho_star + p_star; 
               momentum_y_flux_x(i, j, k) = momentum_x_star*momentum_y_star/rho_star; 
               momentum_z_flux_x(i, j, k) = momentum_x_star*momentum_z_star/rho_star; 
               energy_flux_x(i, j, k) = (en_star + p_star)*(momentum_x_star/rho_star)

               C_L = sqrt(gamma*PL/RL) + abs(VXL)
               C_R = sqrt(gamma*PR/RR) + abs(VXR)
               C_star = max(C_L, C_R)

               mass_flux_x(i, j, k) =  momentum_x_star- (0.5_rk*C_star*(RL - RR))
               momentum_x_flux_x(i, j, k) = momentum_x_flux_x(i, j, k) - (C_star*(RL*VXL - RR*VXR))/2.0_rk
               momentum_y_flux_x(i, j, k) = momentum_y_flux_x(i, j, k) - (C_star*(RL*VYL - RR*VYR))/2.0_rk
               momentum_z_flux_x(i, j, k) = momentum_z_flux_x(i, j, k) - (C_star*(RL*VZL - RR*VZR))/2.0_rk
               energy_flux_x(i, j, k) = energy_flux_x(i, j, k) - (C_star*(en_left - en_right))/2.0_rk
            END DO
         END DO
      END DO
   END SUBROUTINE reconstructFlux

   SUBROUTINE addFluxes(mass_flux_x, mass_flux_y, mass_flux_z, &
                        momentum_x_flux_x, momentum_x_flux_y, momentum_x_flux_z, &
                        momentum_y_flux_x, momentum_y_flux_y, momentum_y_flux_z, &
                        momentum_z_flux_x, momentum_z_flux_y, momentum_z_flux_z, &
                        energy_flux_x, energy_flux_y, energy_flux_z, &
                        mass, momentum_x, momentum_y, momentum_z, energy, &
                        Nx, Ny, Nz,nGhosts, dt, ds)
      REAL(RK), dimension(:, :, :), intent(inout) :: mass_flux_x, mass_flux_y, mass_flux_z, &
                                                     momentum_x_flux_x, momentum_x_flux_y, momentum_x_flux_z, &
                                                     momentum_y_flux_x, momentum_y_flux_y, momentum_y_flux_z, &
                                                     momentum_z_flux_x, momentum_z_flux_y, momentum_z_flux_z, &
                                                     energy_flux_x, energy_flux_y, energy_flux_z
      REAL(RK), dimension(:, :, :), intent(inout) :: mass, momentum_x, momentum_y, momentum_z, energy
      INTEGER(IK), intent(in) :: Nx, Ny, Nz,nGhosts
      REAL(RK), intent(in) :: dt, ds
      INTEGER(IK) :: i, j, k

      DO k = nGhosts+1, Nz-nGhosts
         DO j = nGhosts+1, Ny-nGhosts
            DO i = nGhosts+1, Nx-nGhosts
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
   END SUBROUTINE addFluxes
   function maxmod(a,b)
     use types_and_kinds
     implicit none
     REAL(RK) :: a, b
     REAL(RK) :: maxmod

     if (abs(a) > abs(b) .and. a*b > 0.d0) then
        maxmod = a
     else if (abs(b) > abs(a) .and. a*b > 0) then
        maxmod = b
     else
        maxmod = 0.d0
     endif

     return 
   end function maxmod
   function minmod(a,b)
     use types_and_kinds
     implicit none
     REAL(RK) :: a, b
     REAL(RK) :: minmod
       if (abs(a) < abs(b) .and. a*b > 0.d0) then
     minmod = a
     else if (abs(b) < abs(a) .and. a*b > 0) then
        minmod = b
     else
        minmod = 0.d0
     endif
     return 
   end function minmod
   REAL (RK) function vanleerlimiter(a, b)
    REAL (RK), intent(in) :: a, b
    vanleerlimiter = (sign(1.0d0, a) + sign(1.0d0, b)) * min(abs(a), abs(b)) / (abs(a) + abs(b) + 1.0d-30)
  end function vanleerlimiter
   SUBROUTINE calculate_gradients(Grid, dGdx, dGdy, dGdz, Nx, Ny, Nz,nGhosts,ds)
      use types_and_kinds
      implicit none
      REAL(RK), dimension(:, :, :), intent(in) :: Grid
      REAL(RK), dimension(:, :, :), intent(inout) :: dGdx, dGdy, dGdz
      REAL(RK), intent(in) :: ds
      REAL (RK):: s1,s2
      INTEGER(ik), intent(in) :: Nx, Ny, Nz,nGhosts
      INTEGER(ik) :: i, j, k
      DO k = nGhosts+1, Nz-nGhosts
         DO j = nGhosts+1, Ny-nGhosts
            DO i = nGhosts+1, Nx-nGhosts

            
               !Naive
               dGdy(i, j, k) = (Grid(i, j+1, k) - Grid(i, j-1, k))/(2.0_rk*ds)
               dGdx(i, j, k) = (Grid(i+1, j, k) - Grid(i-1, j, k))/(2.0_rk*ds)
               dGdz(i, j, k) = (Grid(i, j, k+1) - Grid(i, j, k-1))/(2.0_rk*ds)


               !!Van Leer 
               !dGdx(i,j,k)=vanleerlimiter(Grid(i+1,j,k) - Grid(i,j,k), Grid(i,j,k) - Grid(i-1,j,k))/(2.0_rk*ds)
               !dGdy(i,j,k)=vanleerlimiter(Grid(i,j+1,k) - Grid(i,j,k), Grid(i,j,k) - Grid(i,j-1,k))/(2.0_rk*ds)
               !dGdz(i,j,k)=vanleerlimiter(Grid(i,j,k+1) - Grid(i,j,k), Grid(i,j,k) - Grid(i,j,k-1))/(2.0_rk*ds)

               !Mimond only
               !dGdx(i,j,k)=minmod(Grid(i+1,j,k) - Grid(i,j,k), Grid(i,j,k) - Grid(i-1,j,k))/(2.0*ds)
               !dGdy(i,j,k)=minmod(Grid(i,j+1,k) - Grid(i,j,k), Grid(i,j,k) - Grid(i,j-1,k))/(2.0*ds)
               !dGdz(i,j,k)=minmod(Grid(i,j,k+1) - Grid(i,j,k), Grid(i,j,k) - Grid(i,j,k+1))/(2.0*ds)
               !dGdx(i, j, k)=minmod(minmod(2.0_rk*(Grid(i,j,k)-Grid(i-1,j,k))/ds,2.d0*(Grid(i+1,j,k)-Grid(i,j,k))/ds),0.5d0*(Grid(i+1,j,k) - Grid(i-1,j,k))/ds)
               !dGdy(i, j, k)=minmod(minmod(2.0_rk*(Grid(i,j,k)-Grid(i,j-1,k))/ds,2.d0*(Grid(i,j+1,k)-Grid(i,j,k))/ds),0.5d0*(Grid(i,j+1,k) - Grid(i,j-1,k))/ds)
               !dGdz(i, j, k)=minmod(minmod(2.0_rk*(Grid(i,j,k)-Grid(i,j,k-1))/ds,2.d0*(Grid(i,j,k+1)-Grid(i,j,k))/ds),0.5d0*(Grid(i,j,k+1) - Grid(i,j,k-1))/ds)

               !!SuperBee
               !s1=  minmod((Grid(i+1,j,k) - Grid(i,j,k))/ds,2.d0*(Grid(i,j,k) - Grid(i-1,j,k))/ds)
               !s2 = minmod(2.d0*(Grid(i+1,j,k) - Grid(i,j,k))/ds,(Grid(i,j,k) - Grid(i-1,j,k))/ds)
               !dGdx(i,j,k) = maxmod(s1, s2)

               !s1=  minmod((Grid(i,j+1,k) - Grid(i,j,k))/ds,2.d0*(Grid(i,j,k) - Grid(i,j-1,k))/ds)
               !s2 = minmod(2.d0*(Grid(i,j+1,k) - Grid(i,j,k))/ds,(Grid(i,j,k) - Grid(i,j-1,k))/ds)
               !dGdy(i,j,k) = maxmod(s1, s2)

               !s1=  minmod((Grid(i,j,k+1) - Grid(i,j,k))/ds,2.d0*(Grid(i,j,k) - Grid(i,j,k-1))/ds)
               !s2 = minmod(2.d0*(Grid(i,j,k+1) - Grid(i,j,k))/ds,(Grid(i,j,k) - Grid(i,j,k-1))/ds)
               !dGdz(i,j,k) = maxmod(s1, s2)



            END DO
         END DO
      END DO
   END SUBROUTINE calculate_gradients

   SUBROUTINE update_ghosts(Grid, Nx, Ny, Nz,nGhosts)
      use types_and_kinds
      implicit none
      REAL(RK), dimension(:, :, :), intent(inout) :: Grid
      INTEGER(ik), intent(in) :: Nx, Ny, Nz,nGhosts
      Grid(1,:,:)=Grid(Nx-nGhosts-2,:,:)
      Grid(2,:,:)=Grid(Nx-nGhosts-1,:,:)
      Grid(Nx-1,:,:)=Grid(nGhosts+2,:,:)
      Grid(Nx,:,:)=Grid(nGhosts+3,:,:)

      Grid(:,1,:)=Grid(:,Nx-nGhosts-2,:)
      Grid(:,2,:)=Grid(:,Nx-nGhosts-1,:)
      Grid(:,Nx-1,:)=Grid(:,nGhosts+2,:)
      Grid(:,Nx,:)=Grid(:,nGhosts+3,:)

      Grid(:,:,1)=Grid(:,:,Nx-nGhosts-2)
      Grid(:,:,2)=Grid(:,:,Nx-nGhosts-1)
      Grid(:,:,Nx-1)=Grid(:,:,nGhosts+2)
      Grid(:,:,Nx)=Grid(:,:,nGhosts+3)

   END SUBROUTINE update_ghosts

   SUBROUTINE applyPeriodicBCs(Grid, Nx, Ny, Nz,nGhosts)
   use types_and_kinds
   implicit none
   REAL(RK), dimension(:, :, :), intent(inout) :: Grid
   INTEGER(ik), intent(in) :: Nx, Ny, Nz,nGhosts

   ! Periodic boundary conditions in the X direction
   Grid(2, :, :) = Grid(Nx-1, :, :)
   Grid(Nx-1, :, :) = Grid(2, :, :)

   ! Periodic boundary conditions in the Y direction
   Grid(:, 2, :) = Grid(:, Ny-1, :)
   Grid(:, Ny-1, :) = Grid(:, 2, :)

   ! Periodic boundary conditions in the Z direction
   Grid(:, :, 2) = Grid(:, :, Nz-1)
   Grid(:, :, Nz-1) = Grid(:, :, 2)
   END SUBROUTINE applyPeriodicBCs

   SUBROUTINE applySlipBoundary(Vx,Vy,Vz, Nx, Ny, Nz,nGhosts,Normal)
   use types_and_kinds
   implicit none
   REAL(RK), dimension(:, :, :), intent(inout) :: Vx,Vy,Vz
   Real(RK) ::Normal(3)
   INTEGER(ik), intent(in) :: Nx, Ny, Nz,nGhosts
   END SUBROUTINE applySlipBoundary

END MODULE physics

PROGRAM EULER_CFD
   USE types_and_kinds
   use initialization
   use io
   use physics
   IMPLICIT NONE
   REAL(RK), dimension(:, :, :), allocatable  :: rho, vx, vy, vz, p, mass, momentum_x, momentum_y, momentum_z, energy, &
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
                                                 dp_dx, dp_dy, dp_dz

   INTEGER(ik), PARAMETER :: Xcells = 128, Ycells = 128, Zcells = 128, nGhosts = 2
   INTEGER(ik), PARAMETER :: Nx = Xcells+2*nGhosts, Ny = Ycells+2*nGhosts, Nz = Zcells+2*nGhosts
   REAL(rk), PARAMETER :: ds = 1.0_rk/NX
   REAL(rk):: dt = -1.0_rk, time, time_max
   INTEGER(ik) :: timestep = 1, timestep_max
   INTEGER(ik) :: shiftX(3) , shiftY(3), shiftZ(3)
   CHARACTER(72)::filename
   ! Initialize grids
   allocate (rho(Nx, Ny, Nz), vx(Nx, Ny, Nz), vy(Nx, Ny, Nz), vz(Nx, Ny, Nz), p(Nx, Ny, Nz), mass(Nx, Ny, Nz), &
             momentum_x(Nx, Ny, Nz), momentum_y(Nx, Ny, Nz), momentum_z(Nx, Ny, Nz), energy(Nx, Ny, Nz), &
             rho_prime(Nx, Ny, Nz), vx_prime(Nx, Ny, Nz), vy_prime(Nx, Ny, Nz), vz_prime(Nx, Ny, Nz), temp(Nx, Ny, Nz), &
             mass_flux_x(Nx, Ny, Nz), mass_flux_y(Nx, Ny, Nz), mass_flux_z(Nx, Ny, Nz), &
             momentum_x_flux_x(Nx, Ny, Nz), momentum_x_flux_y(Nx, Ny, Nz), momentum_x_flux_z(Nx, Ny, Nz), &
             momentum_y_flux_x(Nx, Ny, Nz), momentum_y_flux_y(Nx, Ny, Nz), momentum_y_flux_z(Nx, Ny, Nz), &
             momentum_z_flux_x(Nx, Ny, Nz), momentum_z_flux_y(Nx, Ny, Nz), momentum_z_flux_z(Nx, Ny, Nz), &
             energy_flux_x(Nx, Ny, Nz), energy_flux_y(Nx, Ny, Nz), energy_flux_z(Nx, Ny, Nz), &
             drho_dx(Nx, Ny, Nz), drho_dy(Nx, Ny, Nz), drho_dz(Nx, Ny, Nz), &
             dvx_dx(Nx, Ny, Nz), dvx_dy(Nx, Ny, Nz), dvx_dz(Nx, Ny, Nz), &
             dvy_dx(Nx, Ny, Nz), dvy_dy(Nx, Ny, Nz), dvy_dz(Nx, Ny, Nz), &
             dvz_dx(Nx, Ny, Nz), dvz_dy(Nx, Ny, Nz), dvz_dz(Nx, Ny, Nz), &
             dp_dx(Nx, Ny, Nz), dp_dy(Nx, Ny, Nz), dp_dz(Nx, Ny, Nz))

    shiftX(:)=0
    shiftY(:)=0
    shiftZ(:)=0
    shiftX(1)=1
    shiftY(2)=1
    shiftZ(3)=1
   ! initialize vx, vy, vz to 0
   !call init_grid_gaussian(rho, Nx, Ny, Nz, NX/2.0_rk, 8.0_rk, 0.01_rk, 1.2_rk)
   !call init_grid(rho, Nx, Ny, Nz, 1.2_rk) ! pascal
   call init_grid(p, Nx, Ny, Nz,nGhosts, 101325.0_rk) ! pascal
   call init_grid_gaussian(rho, Nx, Ny, Nz,nGhosts, NX/2.0_rk, 8.0_rk, 0.000_rk*1.2_rk, 1.2_rk)
   call init_grid_gaussian(p, Nx, Ny, Nz,nGhosts, NX/2.0_rk, 8.0_rk, 0.01_rk*101325, 101325.0_rk)
   call init_grid(vx, Nx, Ny, Nz,nGhosts, 0.0_rk)
   call init_grid(vy, Nx, Ny, Nz,nGhosts, 0.0_rk)
   call init_grid(vz, Nx, Ny, Nz,nGhosts, 0.0_rk)
   call init_grid(mass_flux_x, Nx, Ny, Nz,nGhosts, 0.0_rk)
   call init_grid(mass_flux_y, Nx, Ny, Nz,nGhosts, 0.0_rk)
   call init_grid(mass_flux_z, Nx, Ny, Nz,nGhosts, 0.0_rk)
   call init_grid(momentum_x_flux_x, Nx, Ny, Nz,nGhosts, 0.0_rk)
   call init_grid(momentum_x_flux_y, Nx, Ny, Nz,nGhosts, 0.0_rk)
   call init_grid(momentum_x_flux_z, Nx, Ny, Nz,nGhosts, 0.0_rk)
   call init_grid(momentum_y_flux_x, Nx, Ny, Nz,nGhosts, 0.0_rk)
   call init_grid(momentum_y_flux_y, Nx, Ny, Nz,nGhosts, 0.0_rk)
   call init_grid(momentum_y_flux_z, Nx, Ny, Nz,nGhosts, 0.0_rk)
   call init_grid(momentum_z_flux_x, Nx, Ny, Nz,nGhosts, 0.0_rk)
   call init_grid(momentum_z_flux_y, Nx, Ny, Nz,nGhosts, 0.0_rk)
   call init_grid(momentum_z_flux_z, Nx, Ny, Nz,nGhosts, 0.0_rk)
   call init_grid(energy_flux_x, Nx, Ny, Nz,nGhosts, 0.0_rk)
   call init_grid(energy_flux_y, Nx, Ny, Nz,nGhosts, 0.0_rk)
   call init_grid(energy_flux_z, Nx, Ny, Nz,nGhosts, 0.0_rk)


   call update_ghosts(rho,Nx,Ny,Nz,nGhosts)
   call update_ghosts(vx,Nx,Ny,Nz,nGhosts)
   call update_ghosts(vy,Nx,Ny,Nz,nGhosts)
   call update_ghosts(vz,Nx,Ny,Nz,nGhosts)
   call update_ghosts(p,Nx,Ny,Nz,nGhosts)

   call conservative(mass, momentum_x, momentum_y, momentum_z, energy, rho, p, vx, vy, vz, temp, ds)
   
   time = 0_rk
   !MAIN
   do while (timestep < 1005)


      call primitive(mass, momentum_x, momentum_y, momentum_z, energy, rho, p, vx, vy, vz, temp, ds)

   call update_ghosts(rho,Nx,Ny,Nz,nGhosts)
   call update_ghosts(vx,Nx,Ny,Nz,nGhosts)
   call update_ghosts(vy,Nx,Ny,Nz,nGhosts)
   call update_ghosts(vz,Nx,Ny,Nz,nGhosts)
   call update_ghosts(p,Nx,Ny,Nz,nGhosts)
      rho(3,:,:)=rho(Nx-2,:,:)
      vx(3,:,:)=vx(Nx-2,:,:)
      vy(3,:,:)=vy(Nx-2,:,:)
      vz(3,:,:)=vz(Nx-2,:,:)
      p(3,:,:)=p(Nx-2,:,:)

      rho(:,3,:)=rho(:,Nx-2,:)
      vx(:,3,:)=vx(:,Nx-2,:)
      vy(:,3,:)=vy(:,Nx-2,:)
      vz(:,3,:)=vz(:,Nx-2,:)
      p(:,3,:)=p(:,Nx-2,:)

      rho(:,:,3)=rho(:,:,Nx-2)
      vx(:,:,3)=vx(:,:,Nx-2)
      vy(:,:,3)=vy(:,:,Nx-2)
      vz(:,:,3)=vz(:,:,Nx-2)
      p(:,:,3)=p(:,:, Nx-2)

      dt =2e-7;!get_timestep(ds, vx, vy, vz, p, rho)

      call calculate_gradients(rho, drho_dx, drho_dy, drho_dz, Nx, Ny, Nz,nGhosts, ds)
      call calculate_gradients(vx, dvx_dx, dvx_dy, dvx_dz, Nx, Ny, Nz,nGhosts, ds)
      call calculate_gradients(vy, dvy_dx, dvy_dy, dvy_dz, Nx, Ny, Nz,nGhosts, ds)
      call calculate_gradients(vz, dvz_dx, dvz_dy, dvz_dz, Nx, Ny, Nz,nGhosts, ds)
      call calculate_gradients(p, dp_dx, dp_dy, dp_dz, Nx, Ny, Nz, nGhosts,ds)

      call update_ghosts(drho_dx,Nx,Ny,Nz,nGhosts)
      call update_ghosts(drho_dy,Nx,Ny,Nz,nGhosts)
      call update_ghosts(drho_dz,Nx,Ny,Nz,nGhosts)
      call update_ghosts(dvx_dx,Nx,Ny,Nz,nGhosts)
      call update_ghosts(dvx_dy,Nx,Ny,Nz,nGhosts)
      call update_ghosts(dvx_dz,Nx,Ny,Nz,nGhosts)

      call update_ghosts(dvy_dx,Nx,Ny,Nz,nGhosts)
      call update_ghosts(dvy_dy,Nx,Ny,Nz,nGhosts)
      call update_ghosts(dvy_dz,Nx,Ny,Nz,nGhosts)

      call update_ghosts(dvz_dx,Nx,Ny,Nz,nGhosts)
      call update_ghosts(dvz_dy,Nx,Ny,Nz,nGhosts)
      call update_ghosts(dvz_dz,Nx,Ny,Nz,nGhosts)

      call update_ghosts(dp_dx,Nx,Ny,Nz,nGhosts)
      call update_ghosts(dp_dy,Nx,Ny,Nz,nGhosts)
      call update_ghosts(dp_dz,Nx,Ny,Nz,nGhosts)


      call reconstructFlux(mass_flux_x, momentum_x_flux_x, momentum_y_flux_x, momentum_z_flux_x, &
                             energy_flux_x, drho_dx, dvx_dx, dvy_dx, dvz_dx, dp_dx, &
                             rho, vx, vy, vz, p, Nx, Ny, Nz,nGhosts ,ds,shiftX)

      call reconstructFlux(mass_flux_y, momentum_y_flux_y, momentum_x_flux_y, momentum_z_flux_y, energy_flux_y, &
                             drho_dy, dvy_dy, dvx_dy, dvz_dy, dp_dy, &
                             rho, vy, vx, vz, p, Nx, Ny, Nz,nGhosts ,ds,shiftY)

      call reconstructFlux(mass_flux_z, momentum_z_flux_z, momentum_y_flux_z, momentum_x_flux_z, energy_flux_z, &
                             drho_dz, dvz_dz, dvy_dz, dvx_dz, dp_dz, &
                             rho, vz, vy, vx, p, Nx, Ny, Nz,nGhosts ,ds,shiftZ)


      call update_ghosts(mass_flux_x,Nx,Ny,Nz,nGhosts)
      call update_ghosts(mass_flux_y,Nx,Ny,Nz,nGhosts)
      call update_ghosts(mass_flux_z,Nx,Ny,Nz,nGhosts)

      call update_ghosts(momentum_x_flux_x,Nx,Ny,Nz,nGhosts)
      call update_ghosts(momentum_y_flux_x,Nx,Ny,Nz,nGhosts)
      call update_ghosts(momentum_z_flux_x,Nx,Ny,Nz,nGhosts)
      call update_ghosts(momentum_x_flux_y,Nx,Ny,Nz,nGhosts)
      call update_ghosts(momentum_y_flux_y,Nx,Ny,Nz,nGhosts)
      call update_ghosts(momentum_z_flux_y,Nx,Ny,Nz,nGhosts)
      call update_ghosts(momentum_x_flux_z,Nx,Ny,Nz,nGhosts)
      call update_ghosts(momentum_y_flux_z,Nx,Ny,Nz,nGhosts)
      call update_ghosts(momentum_z_flux_z,Nx,Ny,Nz,nGhosts)


      call update_ghosts(energy_flux_x,Nx,Ny,Nz,nGhosts)
      call update_ghosts(energy_flux_y,Nx,Ny,Nz,nGhosts)
      call update_ghosts(energy_flux_z,Nx,Ny,Nz,nGhosts)

      call addFluxes(mass_flux_x, mass_flux_y, mass_flux_z, &
                     momentum_x_flux_x, momentum_x_flux_y, momentum_x_flux_z, &
                     momentum_y_flux_x, momentum_y_flux_y, momentum_y_flux_z, &
                     momentum_z_flux_x, momentum_z_flux_y, momentum_z_flux_z, &
                     energy_flux_x, energy_flux_y, energy_flux_z, &
                     mass, momentum_x, momentum_y, momentum_z, energy, &
                     Nx, Ny, Nz,nGhosts, dt, ds)

      call update_ghosts(mass,Nx,Ny,Nz,nGhosts)
      call update_ghosts(momentum_x,Nx,Ny,Nz,nGhosts)
      call update_ghosts(momentum_y,Nx,Ny,Nz,nGhosts)
      call update_ghosts(momentum_z,Nx,Ny,Nz,nGhosts)
      call update_ghosts(energy,Nx,Ny,Nz,nGhosts)




      write (filename, '(A,"000",i7.7,A)') "mass.", timestep, ".dat"
      call write_grid(mass, filename)
      write (filename, '(A,"000",i7.7,A)') "rho.", timestep, ".dat"
      call write_grid(rho, filename)
      write (filename, '(A,"000",i7.7,A)') "vx.", timestep, ".dat"
      call write_grid(vx, filename)
      write (filename, '(A,"000",i7.7,A)') "vy.", timestep, ".dat"
      call write_grid(vy, filename)
      write (filename, '(A,"000",i7.7,A)') "vz.", timestep, ".dat"
      call write_grid(vz, filename)
      write (filename, '(A,"000",i7.7,A)') "pressure.", timestep, ".dat"
      call write_grid(p, filename)
      write (filename, '(A,"000",i7.7,A)') "mass_flux_x.", timestep, ".dat"
      call write_grid(mass_flux_x, filename)
      write (filename, '(A,"000",i7.7,A)') "mass_flux_y.", timestep, ".dat"
      call write_grid(mass_flux_y, filename)
      write (filename, '(A,"000",i7.7,A)') "mass_flux_z.", timestep, ".dat"
      call write_grid(mass_flux_z, filename)
      write (filename, '(A,"000",i7.7,A)') "momentum_x_flux_x.", timestep, ".dat"
      call write_grid(momentum_x_flux_x, filename)
      write (filename, '(A,"000",i7.7,A)') "momentum_x_flux_y.", timestep, ".dat"
      call write_grid(momentum_x_flux_y, filename)
      write (filename, '(A,"000",i7.7,A)') "momentum_x_flux_z.", timestep, ".dat"
      call write_grid(momentum_x_flux_z, filename)
      write (filename, '(A,"000",i7.7,A)') "dvx_dx.", timestep, ".dat"
      call write_grid(dvx_dx, filename)
      write (filename, '(A,"000",i7.7,A)') "drho_dx.", timestep, ".dat"
      call write_grid(drho_dx, filename)
      write (filename, '(A,"000",i7.7,A)') "momentum_x.", timestep, ".dat"
      call write_grid(momentum_x, filename)
      write (filename, '(A,"000",i7.7,A)') "momentum_y.", timestep, ".dat"
      call write_grid(momentum_y, filename)
      write (filename, '(A,"000",i7.7,A)') "momentum_z.", timestep, ".dat"
      call write_grid(momentum_z, filename)

      print *, "tstep=", timestep, "dt=", dt
      timestep = timestep + 1
      time = time + dt
   end do

  deallocate (rho, vx, vy, vz, p, mass, momentum_x, momentum_y, momentum_z, energy, rho_prime, vx_prime, vy_prime, vz_prime, temp, &
               mass_flux_x, mass_flux_y, mass_flux_z, &
               momentum_x_flux_x, momentum_x_flux_y, momentum_x_flux_z, momentum_y_flux_x, momentum_y_flux_y, momentum_y_flux_z, &
               momentum_z_flux_x, momentum_z_flux_y, momentum_z_flux_z, energy_flux_x, energy_flux_y, energy_flux_z, &
               drho_dx, drho_dy, drho_dz, &
               dvx_dx, dvx_dy, dvx_dz, &
               dvy_dx, dvy_dy, dvy_dz, &
               dvz_dx, dvz_dy, dvz_dz, &
               dp_dx, dp_dy, dp_dz)
END PROGRAM EULER_CFD
