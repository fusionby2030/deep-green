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
    REAL(RK), PARAMETER :: gamma=1.4_rk, Rs=287.0_rk, cfl=0.4_rk ! rydeberg cons

END MODULE global

MODULE initialization
    use types_and_kinds
    IMPLICIT NONE 
    CONTAINS
    SUBROUTINE init_grid(grid, Nx, Ny, Nz, init_value)
        REAL(RK), dimension(:, :, :), intent(out) :: grid
        INTEGER(ik), intent(in) :: Nx, Ny, Nz
        REAL(RK), intent(in) :: init_value
        INTEGER(ik) :: i, j, k
        ! to efficintly traverse the grid in the order of memory layout
        ! loop over 3D array with proper cache efficiency 
        DO k = 1, Nz
            DO j = 1, Ny
                DO i = 1, Nx
                    grid(i, j, k) = init_value
                END DO
            END DO
        END DO
    END SUBROUTINE init_grid
    ! similar to above subroutine, we now intialize the grid with a 3D gaussian centered in middle of grid
    REAL(RK) FUNCTION GAUSSIAN3D(i, j, k, mu, sigma)
        REAL(RK), intent(in) :: mu, sigma
        INTEGER(ik), intent(in) :: i, j, k
        REAL(RK) :: A=1.0
        GAUSSIAN3D = A*exp(-((i-mu)**2 + (j-mu)**2 + (k-mu)**2)/(2*sigma**2))
    END FUNCTION GAUSSIAN3D
    ! use the above function to initialize the grid
    SUBROUTINE init_grid_gaussian(grid, Nx, Ny, Nz, init_value_mu, init_value_sigma, magnitude, offset)
        REAL(RK), dimension(:, :, :), intent(out) :: grid
        INTEGER(ik), intent(in) :: Nx, Ny, Nz
        REAL(RK), intent(in) :: init_value_mu, init_value_sigma, magnitude, offset
        INTEGER(ik) :: i, j, k
        DO k = 1, Nz
            DO j = 1, Ny
                DO i = 1, Nx
                    ! realistic values
                    ! 0.01, 1.2_rk
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
        open(1, file=filename, form='unformatted', access='stream', status='replace')
        write(1) grid 
        close(1) 
    END SUBROUTINE write_grid
END MODULE io

MODULE physics
    use types_and_kinds
    use global
    IMPLICIT NONE
    CONTAINS
    SUBROUTINE getconserved(mass, mom_x, mom_y, mom_z, energy, rho, p, vx, vy, vz, temp, ds)
        REAL(RK), dimension(:, :, :), intent(inout) :: mass, mom_x, mom_y, mom_z, energy, temp
        REAL(RK), dimension(:, :, :), intent(in) :: rho, vx, vy, vz, p
        REAL(RK), intent(in) :: ds
        REAL(RK) :: cell_volume
        cell_volume=ds**3
        mass = rho*cell_volume
        mom_x = rho*vx*cell_volume
        mom_y = rho*vy*cell_volume
        mom_z = rho*vz*cell_volume
        energy = rho*(vx**2 + vy**2 + vz**2)/2.0_rk + p/(gamma-1.0_rk)
        temp = p/(Rs*rho) - 273.15_rk
    END SUBROUTINE getconserved
    REAL(RK) FUNCTION get_timestep(ds, vx, vy, vz, p, rho)
        REAL(RK), intent(in) :: ds
        REAL(RK), intent(in), dimension(:, :, :) :: vx, vy, vz, p, rho
        get_timestep = minval(cfl*ds/sqrt((gamma*p/rho) + (vx**2 + vy**2 + vz**2)) )
        ! TODO: DOUBLE CHECK THIS
        ! get_timestep = minval(cfl*ds/(sqrt(gamma*p/rho) + sqrt(vx**2 + vy**2 + vz**2)) )
    END FUNCTION get_timestep
    SUBROUTINE getprimitives(mom_x, mom_y, mom_z, energy, rho, p, vx, vy, vz, temp, ds)
        REAL(RK), dimension(:, :, :), intent(in) :: mom_x, mom_y, mom_z, energy
        REAL(RK), dimension(:, :, :), intent(inout) :: p, vx, vy, vz, rho, temp
        REAL(RK), intent(in) :: ds
        REAL(RK) :: cell_volume
        cell_volume=ds**3
        vx = mom_x/(rho*cell_volume)
        vy = mom_y/(rho*cell_volume)
        vz = mom_z/(rho*cell_volume)
        p = (gamma-1.0_rk)*(energy/cell_volume - rho*(vx**2 + vy**2 + vz**2)/2.0_rk)
        temp = p/(Rs*rho) - 273.15_rk
    END SUBROUTINE getprimitives
    ! calculate rusanov flux
    SUBROUTINE getflux_x(mass_flux_x, mom_x_flux_x, mom_y_flux_x, mom_z_flux_x, energy_flux_x, rho, vx, vy, vz, p, Nx, Ny, Nz)
        REAL(RK), dimension(:, :, :), intent(inout) :: mass_flux_x, mom_x_flux_x, mom_y_flux_x, mom_z_flux_x, energy_flux_x
        REAL(RK), dimension(:, :, :), intent(in) :: rho, vx, vy, vz, p
        INTEGER(IK), intent(in) :: Nx, Ny, Nz
        ! REAL(RK), intent(in) :: ds
        INTEGER(IK) :: i, j, k
        REAL(RK) :: rho_star, mom_x_star, mom_y_star, mom_z_star, p_star, en_right, en_left, en_star
        REAL(RK) :: C_L, C_R, C_star
        ! start by calculating rho_star, which is average of density 
        DO k = 1, Nz
            DO j = 1, Ny
                DO i = 1, Nx-1
                    rho_star = (rho(i, j, k) + rho(i+1, j, k))/2.0_rk
                    mom_x_star = (rho(i, j, k)*vx(i,j,k) + rho(i+1, j, k)*vx(i+1,j,k))/2.0_rk
                    mom_y_star = (rho(i, j, k)*vy(i,j,k) + rho(i+1, j, k)*vy(i+1,j,k))/2.0_rk
                    mom_z_star = (rho(i, j, k)*vz(i,j,k) + rho(i+1, j, k)*vz(i+1,j,k))/2.0_rk
                    en_left = (p(i, j, k)/(gamma-1.0_rk)) + & 
                               (rho(i, j, k)*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))/2.0_rk
                    en_right = (p(i+1, j, k)/(gamma-1.0_rk)) + & 
                               (rho(i+1, j, k)*(vx(i+1,j,k)**2 + vy(i+1,j,k)**2 + vz(i+1,j,k)**2))/2.0_rk

                    mom_x_flux_x(i, j, k) = (rho(i, j, k)*vx(i,j,k)**2 + p(i,j,k) + rho(i+1,j,k)*vx(i+1,j,k)**2 +p(i+1,j,k) )/2.0_rk
                    mom_y_flux_x(i, j, k) = (rho(i, j, k)*vy(i,j,k)**2 + p(i,j,k) + rho(i+1,j,k)*vy(i+1,j,k)**2 +p(i+1,j,k) )/2.0_rk
                    mom_z_flux_x(i, j, k) = (rho(i, j, k)*vz(i,j,k)**2 + p(i,j,k) + rho(i+1,j,k)*vz(i+1,j,k)**2 +p(i+1,j,k) )/2.0_rk
                    en_star = (en_left + en_right)
                    p_star = (gamma-1.0_rk)*(en_star - ((mom_x_star**2 + mom_y_star**2 + mom_z_star**2))/(rho_star*2.0_rk))
                    
                    mass_flux_x(i, j, k) = mom_x_star
                    
                    energy_flux_x(i, j, k) = (en_star + p_star )*(mom_x_star/rho_star)

                    C_L = sqrt(gamma*p(i,j,k)/rho(i,j,k)) + abs(vx(i,j,k))
                    C_R = sqrt(gamma*p(i+1,j,k)/rho(i+1,j,k)) + abs(vx(i+1,j,k))
                    C_star = max(C_L, C_R)
                    ! diffusive term
                    mass_flux_x(i, j, k) = mass_flux_x(i,j,k) -(C_star*(rho(i,j,k) - rho(i+1,j,k)))/2.0_rk
                    mom_x_flux_x(i, j, k) = mom_x_flux_x(i,j,k) -(C_star*(rho(i,j,k)*vx(i,j,k) - rho(i+1,j,k)*vx(i+1,j,k)))/2.0_rk
                    mom_y_flux_x(i, j, k) = mom_y_flux_x(i,j,k) -(C_star*(rho(i,j,k)*vy(i,j,k) - rho(i+1,j,k)*vy(i+1,j,k)))/2.0_rk
                    mom_z_flux_x(i, j, k) = mom_z_flux_x(i,j,k) -(C_star*(rho(i,j,k)*vz(i,j,k) - rho(i+1,j,k)*vz(i+1,j,k)))/2.0_rk
                    energy_flux_x(i,j,k) = energy_flux_x(i,j,k) + (C_star*(en_left - en_right))/2.0_rk
                    
                END DO
            END DO
        END DO
    END SUBROUTINE getflux_x
    SUBROUTINE getflux_y(mass_flux_y, mom_x_flux_y, mom_y_flux_y, mom_z_flux_y, energy_flux_y, rho, vx, vy, vz, p, Nx, Ny, Nz)
        REAL(RK), dimension(:, :, :), intent(inout) :: mass_flux_y, mom_x_flux_y, mom_y_flux_y, mom_z_flux_y, energy_flux_y
        REAL(RK), dimension(:, :, :), intent(in) :: rho, vx, vy, vz, p
        INTEGER(IK), intent(in) :: Nx, Ny, Nz
        ! REAL(RK), intent(in) :: ds
        INTEGER(IK) :: i, j, k
        REAL(RK) :: rho_star, mom_x_star, mom_y_star, mom_z_star, p_star, en_right, en_left, en_star
        REAL(RK) :: C_L, C_R, C_star
        ! start by calculating rho_star, which is average of density 
        DO k = 1, Nz
            DO j = 1, Ny-1
                DO i = 1, Nx
                    rho_star = (rho(i, j, k) + rho(i, j+1, k))/2.0_rk
                    mom_x_star = (rho(i, j, k)*vx(i,j,k) + rho(i, j+1, k)*vx(i,j+1,k))/2.0_rk
                    mom_y_star = (rho(i, j, k)*vy(i,j,k) + rho(i, j+1, k)*vy(i,j+1,k))/2.0_rk
                    mom_z_star = (rho(i, j, k)*vz(i,j,k) + rho(i, j+1, k)*vz(i,j+1,k))/2.0_rk
                    en_left = (p(i, j, k)/(gamma-1.0_rk)) + & 
                               (rho(i, j, k)*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))/2.0_rk
                    en_right = (p(i, j+1, k)/(gamma-1.0_rk)) + & 
                               (rho(i, j+1, k)*(vx(i,j+1,k)**2 + vy(i,j+1,k)**2 + vz(i,j+1,k)**2))/2.0_rk

                    mom_x_flux_y(i, j, k) = (rho(i, j, k)*vx(i,j,k)**2 + p(i,j,k) + rho(i,j+1,k)*vx(i,j+1,k)**2 +p(i,j+1,k) )/2.0_rk
                    mom_y_flux_y(i, j, k) = (rho(i, j, k)*vy(i,j,k)**2 + p(i,j,k) + rho(i,j+1,k)*vy(i,j+1,k)**2 +p(i,j+1,k) )/2.0_rk
                    mom_z_flux_y(i, j, k) = (rho(i, j, k)*vz(i,j,k)**2 + p(i,j,k) + rho(i,j+1,k)*vz(i,j+1,k)**2 +p(i,j+1,k) )/2.0_rk
                    en_star = (en_left + en_right)
                    p_star = (gamma-1.0_rk)*(en_star - ((mom_x_star**2 + mom_y_star**2 + mom_z_star**2))/(rho_star*2.0_rk))
                    
                    mass_flux_y(i, j, k) = mom_x_star
                    
                    energy_flux_y(i, j, k) = (en_star + p_star )*(mom_x_star/rho_star) 
                    
                    C_L = sqrt(gamma*p(i,j,k)/rho(i,j,k)) + abs(vx(i,j,k))
                    C_R = sqrt(gamma*p(i,j+1,k)/rho(i,j+1,k)) + abs(vx(i,j+1,k))
                    C_star = max(C_L, C_R)
                    ! diffusive term
                    mass_flux_y(i, j, k) = mass_flux_y(i,j,k) -(C_star*(rho(i,j,k) - rho(i,j+1,k)))/2.0_rk
                    mom_x_flux_y(i, j, k) = mom_x_flux_y(i,j,k) -(C_star*(rho(i,j,k)*vx(i,j,k) - rho(i,j+1,k)*vx(i,j+1,k)))/2.0_rk
                    mom_y_flux_y(i, j, k) = mom_y_flux_y(i,j,k) -(C_star*(rho(i,j,k)*vy(i,j,k) - rho(i,j+1,k)*vy(i,j+1,k)))/2.0_rk
                    mom_z_flux_y(i, j, k) = mom_z_flux_y(i,j,k) -(C_star*(rho(i,j,k)*vz(i,j,k) - rho(i,j+1,k)*vz(i,j+1,k)))/2.0_rk
                    energy_flux_y(i,j,k) = energy_flux_y(i,j,k) + (C_star*(en_left - en_right))/2.0_rk
                    
                END DO
            END DO
        END DO
    END SUBROUTINE getflux_y
    SUBROUTINE getflux_z(mass_flux_z, mom_x_flux_z, mom_y_flux_z, mom_z_flux_z, energy_flux_z, rho, vx, vy, vz, p, Nx, Ny, Nz)
        use omp_lib
        REAL(RK), dimension(:, :, :), intent(inout) :: mass_flux_z, mom_x_flux_z, mom_y_flux_z, mom_z_flux_z, energy_flux_z
        REAL(RK), dimension(:, :, :), intent(in) :: rho, vx, vy, vz, p
        INTEGER(IK), intent(in) :: Nx, Ny, Nz
        ! REAL(RK), intent(in) :: ds
        INTEGER(IK) :: i, j, k
        REAL(RK) :: rho_star, mom_x_star, mom_y_star, mom_z_star, p_star, en_right, en_left, en_star
        REAL(RK) :: C_L, C_R, C_star
        ! start by calculating rho_star, which is average of density 
        !!$OMP PARALLEL  PRIVATE(i, j, k, rho_star, mom_x_star, mom_y_star, mom_z_star, p_star, en_right, en_left, en_star, C_L, C_R, C_star)
        !!$OMP DO COLLAPSE(3)
        DO k = 1, Nz-1
            DO j = 1, Ny
                DO i = 1, Nx
                    rho_star = (rho(i, j, k) + rho(i, j, k+1))/2.0_rk
                    mom_x_star = (rho(i, j, k)*vx(i,j,k) + rho(i, j, k+1)*vx(i,j,k+1))/2.0_rk
                    mom_y_star = (rho(i, j, k)*vy(i,j,k) + rho(i, j, k+1)*vy(i,j,k+1))/2.0_rk
                    mom_z_star = (rho(i, j, k)*vz(i,j,k) + rho(i, j, k+1)*vz(i,j,k+1))/2.0_rk
                    en_left = (p(i, j, k)/(gamma-1.0_rk)) + & 
                               (rho(i, j, k)*(vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2))/2.0_rk
                    en_right = (p(i, j, k+1)/(gamma-1.0_rk)) + & 
                               (rho(i, j, k+1)*(vx(i,j,k+1)**2 + vy(i,j,k+1)**2 + vz(i,j,k+1)**2))/2.0_rk

                    mom_x_flux_z(i, j, k) = (rho(i, j, k)*vx(i,j,k)**2 + p(i,j,k) + rho(i,j,k+1)*vx(i,j,k+1)**2 +p(i,j,k+1) )/2.0_rk
                    mom_y_flux_z(i, j, k) = (rho(i, j, k)*vy(i,j,k)**2 + p(i,j,k) + rho(i,j,k+1)*vy(i,j,k+1)**2 +p(i,j,k+1) )/2.0_rk
                    mom_z_flux_z(i, j, k) = (rho(i, j, k)*vz(i,j,k)**2 + p(i,j,k) + rho(i,j,k+1)*vz(i,j,k+1)**2 +p(i,j,k+1) )/2.0_rk
                    en_star = (en_left + en_right)
                    p_star = (gamma-1.0_rk)*(en_star - ((mom_x_star**2 + mom_y_star**2 + mom_z_star**2))/(rho_star*2.0_rk))
                    
                    mass_flux_z(i, j, k) = mom_x_star
                    
                    energy_flux_z(i, j, k) = (en_star + p_star )*(mom_x_star/rho_star) 
                    
                    ! rusanov fluxes
                    ! signal speeds
                    C_L = sqrt(gamma*p(i,j,k)/rho(i,j,k)) + abs(vx(i,j,k))
                    C_R = sqrt(gamma*p(i,j,k+1)/rho(i,j,k+1)) + abs(vx(i,j,k+1))
                    C_star = max(C_L, C_R)
                    ! diffusive term
                    mass_flux_z(i, j, k) = mass_flux_z(i,j,k) -(C_star*(rho(i,j,k) - rho(i,j,k+1)))/2.0_rk
                    mom_x_flux_z(i, j, k) = mom_x_flux_z(i,j,k) -(C_star*(rho(i,j,k)*vx(i,j,k) - rho(i,j,k+1)*vx(i,j,k+1)))/2.0_rk
                    mom_y_flux_z(i, j, k) = mom_y_flux_z(i,j,k) -(C_star*(rho(i,j,k)*vy(i,j,k) - rho(i,j,k+1)*vy(i,j,k+1)))/2.0_rk
                    mom_z_flux_z(i, j, k) = mom_z_flux_z(i,j,k) -(C_star*(rho(i,j,k)*vz(i,j,k) - rho(i,j,k+1)*vz(i,j,k+1)))/2.0_rk
                    energy_flux_z(i,j,k) = energy_flux_z(i,j,k) + (C_star*(en_left - en_right))/2.0_rk
                END DO
            END DO
        END DO
        !!$OMP END PARALLEL DO
    END SUBROUTINE getflux_z
    SUBROUTINE applyfluxes(mass_flux_x, mass_flux_y, mass_flux_z, &
        mom_x_flux_x, mom_x_flux_y, mom_x_flux_z, &
        mom_y_flux_x, mom_y_flux_y, mom_y_flux_z, &
        mom_z_flux_x, mom_z_flux_y, mom_z_flux_z, &
        energy_flux_x, energy_flux_y, energy_flux_z, & 
        mass, mom_x, mom_y, mom_z, energy, & 
        Nx, Ny, Nz, dt, ds)
        REAL(RK), dimension(:, :, :), intent(inout) :: mass_flux_x, mass_flux_y, mass_flux_z, &
                                                       mom_x_flux_x, mom_x_flux_y, mom_x_flux_z, &
                                                       mom_y_flux_x, mom_y_flux_y, mom_y_flux_z, &
                                                       mom_z_flux_x, mom_z_flux_y, mom_z_flux_z, &
                                                       energy_flux_x, energy_flux_y, energy_flux_z, & 
                                                       mass, mom_x, mom_y, mom_z, energy
        INTEGER(IK), intent(in) :: Nx, Ny, Nz
        REAL(RK), intent(in) :: dt, ds
        INTEGER(IK) :: i, j, k

        do k=2, Nz-1 
            do j=2, Ny-1
                do i=2, Nx-1
                    ! update mass
                    mass(i,j,k) = mass(i,j,k) - (dt/ds)*(mass_flux_x(i,j,k) - mass_flux_x(i-1,j,k) + &
                                                          mass_flux_y(i,j,k) - mass_flux_y(i,j-1,k) + &
                                                          mass_flux_z(i,j,k) - mass_flux_z(i,j,k-1))
                    ! update momentum
                    mom_x(i,j,k) = mom_x(i,j,k) - (dt/ds)*(mom_x_flux_x(i,j,k) - mom_x_flux_x(i-1,j,k) + &
                                                            mom_x_flux_y(i,j,k) - mom_x_flux_y(i,j-1,k) + &
                                                            mom_x_flux_z(i,j,k) - mom_x_flux_z(i,j,k-1))
                    mom_y(i,j,k) = mom_y(i,j,k) - (dt/ds)*(mom_y_flux_x(i,j,k) - mom_y_flux_x(i-1,j,k) + &
                                                            mom_y_flux_y(i,j,k) - mom_y_flux_y(i,j-1,k) + &
                                                            mom_y_flux_z(i,j,k) - mom_y_flux_z(i,j,k-1))
                    mom_z(i,j,k) = mom_z(i,j,k) - (dt/ds)*(mom_z_flux_x(i,j,k) - mom_z_flux_x(i-1,j,k) + &
                                                            mom_z_flux_y(i,j,k) - mom_z_flux_y(i,j-1,k) + &
                                                            mom_z_flux_z(i,j,k) - mom_z_flux_z(i,j,k-1))
                    ! update energy
                    energy(i,j,k) = energy(i,j,k) - (dt/ds)*(energy_flux_x(i,j,k) - energy_flux_x(i-1,j,k) + &
                                                              energy_flux_y(i,j,k) - energy_flux_y(i,j-1,k) + &
                                                              energy_flux_z(i,j,k) - energy_flux_z(i,j,k-1))  
                end do 
            end do 
        end do
    END SUBROUTINE applyfluxes
END MODULE physics

PROGRAM EULER_CFD
    USE types_and_kinds
    use initialization
    use io 
    use physics
    IMPLICIT NONE  
    REAL(RK), dimension(:, :, :), allocatable  :: rho, vx, vy, vz, p, mass, mom_x, mom_y, mom_z, energy, & 
                                                    rho_prime, vx_prime, vy_prime, vz_prime, temp, & 
                                                    ! fluxes
                                                    mass_flux_x, mass_flux_y, mass_flux_z, &
                                                    mom_x_flux_x, mom_x_flux_y, mom_x_flux_z, &
                                                    mom_y_flux_x, mom_y_flux_y, mom_y_flux_z, &
                                                    mom_z_flux_x, mom_z_flux_y, mom_z_flux_z, &
                                                    energy_flux_x, energy_flux_y, energy_flux_z

    
    INTEGER(ik), PARAMETER :: Nx=128, Ny=128, Nz=128
    REAL(rk), PARAMETER :: ds=0.05_rk / NX
    REAL(rk):: dt=-1.0_rk
    INTEGER(ik) :: timestep=1
    ! Initialize grids
    allocate(rho(Nx, Ny, Nz), vx(Nx, Ny, Nz), vy(Nx, Ny, Nz), vz(Nx, Ny, Nz), p(Nx, Ny, Nz), mass(Nx, Ny, Nz), &
             mom_x(Nx, Ny, Nz), mom_y(Nx, Ny, Nz), mom_z(Nx, Ny, Nz), energy(Nx, Ny, Nz), &
             rho_prime(Nx, Ny, Nz), vx_prime(Nx, Ny, Nz), vy_prime(Nx, Ny, Nz), vz_prime(Nx, Ny, Nz), temp(Nx, Ny, Nz),& 
             mass_flux_x(Nx, Ny, Nz), mass_flux_y(Nx, Ny, Nz), mass_flux_z(Nx, Ny, Nz), &
             mom_x_flux_x(Nx, Ny, Nz), mom_x_flux_y(Nx, Ny, Nz), mom_x_flux_z(Nx, Ny, Nz), &
             mom_y_flux_x(Nx, Ny, Nz), mom_y_flux_y(Nx, Ny, Nz), mom_y_flux_z(Nx, Ny, Nz), &
             mom_z_flux_x(Nx, Ny, Nz), mom_z_flux_y(Nx, Ny, Nz), mom_z_flux_z(Nx, Ny, Nz), &
             energy_flux_x(Nx, Ny, Nz), energy_flux_y(Nx, Ny, Nz), energy_flux_z(Nx, Ny, Nz))
    
    
    ! initialize vx, vy, vz to 0 
    call init_grid(vx, Nx, Ny, Nz, 200.0_rk)
    ! call init_grid_gaussian(vx, Nx, Ny, Nz, NX/2.0_rk, 8.0_rk, 1.0_rk, 0.0_rk) 
    call init_grid(vy, Nx, Ny, Nz, 0.0_rk)
    call init_grid(vz, Nx, Ny, Nz, 0.0_rk)
    ! initialize p to 1
    call init_grid(p, Nx, Ny, Nz, 101325.0_rk) ! pascal

    ! initialize rho with gaussian 
    call init_grid_gaussian(rho, Nx, Ny, Nz, NX/2.0_rk, 8.0_rk, 0.1_rk, 1.2_rk)
    ! initialize all fluxes to 0 
    call init_grid(mass_flux_x, Nx, Ny, Nz, 0.0_rk)
    call init_grid(mass_flux_y, Nx, Ny, Nz, 0.0_rk)
    call init_grid(mass_flux_z, Nx, Ny, Nz, 0.0_rk)
    call init_grid(mom_x_flux_x, Nx, Ny, Nz, 0.0_rk)
    call init_grid(mom_x_flux_y, Nx, Ny, Nz, 0.0_rk)
    call init_grid(mom_x_flux_z, Nx, Ny, Nz, 0.0_rk)
    call init_grid(mom_y_flux_x, Nx, Ny, Nz, 0.0_rk)
    call init_grid(mom_y_flux_y, Nx, Ny, Nz, 0.0_rk)
    call init_grid(mom_y_flux_z, Nx, Ny, Nz, 0.0_rk)
    call init_grid(mom_z_flux_x, Nx, Ny, Nz, 0.0_rk)
    call init_grid(mom_z_flux_y, Nx, Ny, Nz, 0.0_rk)
    call init_grid(mom_z_flux_z, Nx, Ny, Nz, 0.0_rk)
    call init_grid(energy_flux_x, Nx, Ny, Nz, 0.0_rk)
    call init_grid(energy_flux_y, Nx, Ny, Nz, 0.0_rk)
    call init_grid(energy_flux_z, Nx, Ny, Nz, 0.0_rk)

    call getconserved(mass, mom_x, mom_y, mom_z, energy, rho, p, vx, vy, vz, temp, ds)
    call write_grid(mass, 'mass0.dat')
    call write_grid(temp, 'temp0.dat')
    call write_grid(rho, 'rho0.dat')
    call write_grid(p, 'p0.dat')
    dt =  get_timestep(ds, vx, vy, vz, p, rho)
    print *, timestep, dt, size(rho)

    do  while (timestep < 1000) 
        dt =  1E-4!get_timestep(ds, vx, vy, vz, p, rho)
        call getprimitives(mom_x, mom_y, mom_z, energy, rho, p, vx, vy, vz, temp, ds)
        call getflux_x(mass_flux_x, mom_x_flux_x, mom_y_flux_x, mom_z_flux_x, energy_flux_x, rho, vx, vy, vz, p, Nx, Ny, Nz)
        call getflux_y(mass_flux_y, mom_x_flux_y, mom_y_flux_y, mom_z_flux_y, energy_flux_y, rho, vx, vy, vz, p, Nx, Ny, Nz)
        call getflux_z(mass_flux_z, mom_x_flux_z, mom_y_flux_z, mom_z_flux_z, energy_flux_z, rho, vx, vy, vz, p, Nx, Ny, Nz)
        call applyfluxes(mass_flux_x, mass_flux_y, mass_flux_z, &
                         mom_x_flux_x, mom_x_flux_y, mom_x_flux_z, &
                         mom_y_flux_x, mom_y_flux_y, mom_y_flux_z, &
                         mom_z_flux_x, mom_z_flux_y, mom_z_flux_z, &
                         energy_flux_x, energy_flux_y, energy_flux_z, & 
                         mass, mom_x, mom_y, mom_z, energy, & 
                         Nx, Ny, Nz, dt, ds)
        call write_grid(mass_flux_x, 'mass_flux_x.dat')
        call write_grid(mass_flux_y, 'mass_flux_y.dat')
        call write_grid(mass_flux_z, 'mass_flux_z.dat')
        call write_grid(rho, 'rho.dat')
        call write_grid(temp, 'temp.dat')
        call write_grid(p, 'pressure.dat')
        timestep = timestep + 1
        print *, timestep, dt, size(rho)
    end do 

    deallocate(rho, vx, vy, vz, p, mass, mom_x, mom_y, mom_z, energy, rho_prime, vx_prime, vy_prime, vz_prime, temp, &
               mass_flux_x, mass_flux_y, mass_flux_z, &
               mom_x_flux_x, mom_x_flux_y, mom_x_flux_z, mom_y_flux_x, mom_y_flux_y, mom_y_flux_z, &
               mom_z_flux_x, mom_z_flux_y, mom_z_flux_z, energy_flux_x, energy_flux_y, energy_flux_z)
    

END PROGRAM EULER_CFD