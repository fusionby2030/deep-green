PROGRAM conductive_heat 
    use types_and_kinds
    use globals
    use initialization
    use simgrid
    use io 
    use physics
IMPLICIT NONE 
    type(Simulation) :: sim 
    type(SimulationGrid), pointer :: grd_ptr
    type(SimulationGrid), TARGET :: grd
    CHARACTER(LEN=50) :: arg
    REAL(RK) :: total_time = 1.0_rk ! 60.0_rk * 60.0_rk * 24.0_rk
    INTEGER(IK) :: max_steps = 10000000
    REAL(RK) :: t_out = 0.01_rk
    INTEGER(IK) :: i
    call get_command_argument(1, demoname)
    SELECT case (demoname)
        case ('PIPE')
            call setup_pipe(sim, grd)
        case ('ANALYTICAL_1')
            call setup_analytical_1(sim, grd)
        case ('ANALYTICAL_2')
            call setup_analytical_2(sim, grd)
        case ('GREENHOUSE')
            print *, "Not implemented yet!!!"
            stop
        case default 
            print *, "Unknown demo"
            stop
    END SELECT 
    grd_ptr => grd
    write ( *, '(a)' ) ' '
    WRITE (*, '(A, A)') "Simulation: ", demoname
    WRITE (*, '(A, I4, I4, I4, A, F8.4)') "Grid size: (dx dy dz)", sim%nx, sim%ny, sim%nz, " with ds ", sim%ds
    call compute(sim, grd_ptr, total_time, max_steps, t_out)
CONTAINS 
    SUBROUTINE compute(info, grid, total_time, max_steps, t_out)
        TYPE(Simulation), INTENT(IN) :: info 
        TYPE(SimulationGrid),  pointer :: grid
        REAL(RK), INTENT(IN) :: total_time, t_out
        INTEGER(IK), INTENT(IN) :: max_steps
        INTEGER(IK) :: step_num=0, wstep=0
        REAL(RK) :: step_time, dt=0.00001_rk, wtime=0.0_rk, cfl, converged
        call statewrite(step_time, step_num, grid)

        cfl = cfl_condition(sim, grid, dt) 
        
        write ( *, '(a,g14.6)' ) '  CFL stability criterion value = ', cfl
        if (cfl > 0.5_rk) then
            print *, "CFL condition not met (> 0.5)"
            stop
        end if
        
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Entering Loop...'
        do step_num = 1, max_steps
            step_time = step_num * dt
            write (*, "(A, I8, A, F12.4, A, F8.6)") "Computing step: ", step_num, " Time: ", step_time, " with dt ", dt
            call step(info, grid, dt)
            
            if (wtime >= t_out  .or. wstep == 0) then
                wstep = wstep + 1
                call statewrite(step_time, wstep, grid)
                wtime = 0.0_rk
            end if
            wtime = wtime + dt 
            
            if (step_time >= total_time) exit
            
            converged = sum(abs(grid%primitives%T - grid%primitives%TNEW))
            if (converged < 1.0e-6) then
                print *, "Converged at step: ", step_num
                exit
            end if
            grid%primitives%T = grid%primitives%TNEW
        end do
    call cleansim(sim, grid)
    END SUBROUTINE compute
END PROGRAM conductive_heat