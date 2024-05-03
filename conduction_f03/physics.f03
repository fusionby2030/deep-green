MODULE physics
    use types_and_kinds
    use globals 
    use simgrid 
    IMPLICIT NONE
    public :: step
CONTAINS
    SUBROUTINE step(sim, grid, dt)
        TYPE(Simulation), INTENT(IN) :: sim 
        TYPE(SimulationGrid), pointer :: grid
        REAL(RK), INTENT(IN) :: dt
        INTEGER :: i, j, k
        integer :: Nx, Ny, Nz
        nx = sim%Nx - 2*sim%Nghosts
        ny = sim%Ny - 2*sim%Nghosts
        nz = sim%Nz - 2*sim%Nghosts
        DO k=sim%Nghosts+1,nz+sim%Nghosts
            DO j=sim%Nghosts+1,ny+sim%Nghosts
                DO i=sim%Nghosts+1,nx+sim%Nghosts-1
                    grid%primitives%TNEW(i,j,k) = grid%primitives%T(i,j,k) + & 
                                    (dt*grid%coefficients%alpha(i,j,k) /(sim%ds**2))*  & 
                                    (grid%primitives%T(i+1,j,k) - 2.0*grid%primitives%T(i,j,k) + grid%primitives%T(i-1,j,k)  + & 
                                    grid%primitives%T(i,j+1,k) - 2.0*grid%primitives%T(i,j,k) + grid%primitives%T(i,j-1,k) + &
                                    grid%primitives%T(i,j,k+1) - 2.0*grid%primitives%T(i,j,k) + grid%primitives%T(i,j,k-1)) + &
                                    grid%sources%S1(i,j,k)*(sim%ds**2 / thermal_condctivity_air)
                END DO
            END DO
        END DO
    END SUBROUTINE step

    REAL(RK) FUNCTION CFL_CONDITION(sim, grid, dt)
        TYPE(Simulation), INTENT(IN) :: sim 
        TYPE(SimulationGrid), pointer :: grid
        REAL(RK), INTENT(IN) :: dt
        REAL(RK) :: max_diffusion 
        max_diffusion = maxval(grid%coefficients%alpha)
        CFL_CONDITION = max_diffusion * dt / sim%ds / sim%ds
    END FUNCTION CFL_CONDITION
END MODULE physics