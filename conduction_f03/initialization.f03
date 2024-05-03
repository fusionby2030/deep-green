MODULE initialization
    use types_and_kinds
    use simgrid 
    use globals 
    IMPLICIT NONE 
CONTAINS
    SUBROUTINE setup_pipe(sim, grd) ! grid
        INTEGER(IK), PARAMETER :: nx=200, ny=2, nz=20
        INTEGER(IK), PARAMETER :: nghosts=1
        REAL(RK) :: ds
        type(Simulation), intent(inout) :: sim 
        type(SimulationGrid), intent(inout):: grd
        ds = 1.0 / REAL(nx, RK)
        sim = Simulation(nx, ny, nz, nghosts, ds)
        grd = SimulationGrid(sim)
        call newsim(sim, grd)
        grd%primitives%T = 2.0  + c_to_k
        grd%primitives%T(1:nx+2*nghosts, 1:ny+2*nghosts, 1:nz+2*nghosts) = -5.0 + c_to_k
        grd%primitives%TNEW = 0.0
        grd%sources%S1 = 20.0 
        grd%coefficients%alpha = thermal_diffusivity_air
    END SUBROUTINE  setup_pipe

    SUBROUTINE setup_analytical_1(sim, grid) 
        INTEGER(IK), PARAMETER :: nx=100, ny=1, nz=1
        INTEGER(IK), PARAMETER :: nghosts=1
        REAL(RK) :: ds
        INTEGER(IK) :: i 
        type(Simulation), intent(inout) :: sim
        type(SimulationGrid), intent(inout) :: grid
        ds = 2.0 / REAL(nx, RK)
        sim = Simulation(nx+2*nghosts, ny+2*nghosts, nz+2*nghosts, nghosts, ds)
        grid = SimulationGrid(sim)
        call newsim(sim, grid)
        ! domain , 0 < x < 2
        ! u_t = alpha * u_xx, alpha = 3.0
        ! u(0, t) = u(2, 0) = 0
        ! u(x, 0) = 50 
        ! solution: u(x, t) = (200/ pi)sum_k=0^inf 
        ! (1/(2k+1)) * sin((2k+1)pi x / 2) * exp(-3*(2*k+1)^2 \pi^2 t / 4)
        
        grid%sources%S1 = 0.0
        grid%coefficients%alpha = 3.0_rk
        ! set the middle layer to 50, and the rest to 0
        grid%primitives%T = 0.0
        do i=0, nx
            grid%primitives%T(i+2, 2, 2) = 50.0_rk
        end do
        grid%primitives%TNEW = 0.0_rk 
    END SUBROUTINE setup_analytical_1
    SUBROUTINE setup_analytical_2(sim, grid)
        INTEGER(IK), PARAMETER :: nx=100, ny=1, nz=1
        INTEGER(IK), PARAMETER :: nghosts=1
        INTEGER(IK) :: i 
        REAL(RK) :: ds
        type(Simulation), intent(inout) :: sim
        type(SimulationGrid), intent(inout) :: grid
        ds = 80 / REAL(nx, RK)
        sim = Simulation(nx+2*nghosts, ny+2*nghosts, nz+2*nghosts, nghosts, ds)
        grid = SimulationGrid(sim)
        call newsim(sim, grid)
        ! domain , 0 < x < 80 cm
        ! u_t = c^2 * u_xx, alpha = c^2 = 1.158 cm^2 / s 
        ! u(0, t) = u(nx, 0) = 0
        ! u(x, 0) = 100 sin((pi / 80) x) 
        ! solution: u(x, t) = 100 exp(-0.001768 t) sin((pi / 80) x)
        ! grid%primitives%T =
        grid%primitives%T = 0.0         
        do i=0, nx
            grid%primitives%T(i+2, 2, 2) = 100.0 * sin((pi / 80.0) * (i) * ds)
        end do
        grid%sources%S1 = 0.0
        grid%coefficients%alpha = 1.158
        grid%primitives%TNEW = grid%primitives%T
    END SUBROUTINE setup_analytical_2
END MODULE initialization 