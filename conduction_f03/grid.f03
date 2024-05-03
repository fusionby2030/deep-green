MODULE SIMGRID 
    use types_and_kinds
    IMPLICIT NONE 

    TYPE, PUBLIC :: Simulation 
        INTEGER(IK) :: nx, ny, nz, nghosts
        REAL(RK) :: ds
        ! RUN TYPE 
        ! Boundaries 
        CONTAINS 
    
    END TYPE Simulation

    TYPE :: primitives
        REAL(RK), DIMENSION(:, :, :), ALLOCATABLE:: T, TNEW
    END TYPE primitives

    TYPE :: sources 
        REAL(RK), DIMENSION(:, :, :), ALLOCATABLE :: S1
    END TYPE sources

    TYPE :: coefficients
        REAL(RK), DIMENSION(:, :, :), ALLOCATABLE :: alpha
    END TYPE coefficients

    TYPE, PUBLIC :: SimulationGrid
        type(Simulation) :: info
        type(primitives) :: primitives
        type(sources) :: sources
        type(coefficients) :: coefficients
        ! apply boundary condition
    END TYPE SimulationGrid

CONTAINS 
    SUBROUTINE allocate_all(grid, prim, sourc, coef)
        type(SimulationGrid), intent(in) :: grid
        type(primitives) :: prim
        type(sources) :: sourc
        type(coefficients) :: coef
        INTEGER(IK) :: nx, ny, nz
        nx = grid%info%nx
        ny = grid%info%ny
        nz = grid%info%nz
        ALLOCATE(prim%T(nx, ny, nz), prim%TNEW(nx, ny, nz), sourc%S1(nx, ny, nz), coef%alpha(nx, ny, nz))
    END SUBROUTINE allocate_all
    
    SUBROUTINE deallocate_all(prim, sourc, coef)
        type(primitives), INTENT(INOUT) :: prim
        type(sources), INTENT(INOUT) :: sourc
        type(coefficients), INTENT(INOUT) :: coef
        DEALLOCATE(prim%T, sourc%S1, coef%alpha)
    END SUBROUTINE deallocate_all

    SUBROUTINE cleansim(info, grid)
        type(Simulation), intent(in) :: info
        type(SimulationGrid), intent(INOUT) :: grid
        print *, "Deallocating Simulation Grid..."
        call deallocate_all(grid%primitives, grid%sources, grid%coefficients)
    END SUBROUTINE cleansim

    SUBROUTINE newsim(info, grid)
        type(Simulation), intent(inout) :: info
        type(SimulationGrid), intent(inout) :: grid
        ! info%nx = info%nx + 2*info%nghosts
        ! info%ny = info%ny + 2*info%nghosts
        ! info%nz = info%nz + 2*info%nghosts
        call allocate_all(grid, grid%primitives, grid%sources, grid%coefficients)
        grid%info = info
    END SUBROUTINE newsim

END MODULE SIMGRID