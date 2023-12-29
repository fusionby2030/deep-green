!-----------------------------------------------------------------------
! Program: heat transport with conductive heat transfer
! Purpose: Solve heat equation in 3D to simulate air flow in greenhouse
! Execution: gfortran heat_equation.f90 -fsanitize=address -Werror -Wpedantic -Wno-tabs -g3 -ggdb -Wall -O3 && ./a.out 
! Authors: Adam Kit and Kostis Papadakis (2023)
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
    INTEGER, PARAMETER :: rk = selected_real_kind(8) ! p
    INTEGER, PARAMETER :: ik = selected_int_kind(8)
END MODULE types_and_kinds

MODULE global
    use types_and_kinds
    IMPLICIT NONE 
    REAL(RK), PARAMETER :: c_to_k = 273.15_rk, k_to_c = -273.15_rk
    REAL(rk), PARAMETER :: thermal_condctivity_air = 0.024 ! W/(m*K)
    ! http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/thrcn.html
    REAL(rk), PARAMETER :: density_air = 1.2922 ! kg/m^3
    ! https://www.earthdata.nasa.gov/topics/atmosphere/atmospheric-pressure/air-mass-density
    REAL(rk), PARAMETER :: specific_heat_air = 1003.5 ! J/(kg*K)
    ! https://en.wikipedia.org/wiki/Table_of_specific_heat_capacities
    REAL(rk), PARAMETER :: thermal_diffusivity_air = thermal_condctivity_air / (density_air * specific_heat_air) ! m^3/s

    REAL(rk), PARAMETER :: thermal_condctivity_plastic = 0.23   ! W/(m*K)
    REAL(RK), PARAMETER :: density_plastic = 1200.0_rk          ! kg/m^3
    REAL(RK), PARAMETER :: specific_heat_plastic = 1200.0_rk    ! J/(kg*K)
    REAL(rk), PARAMETER :: thermal_diffusivity_plastic = thermal_condctivity_plastic / (density_plastic * specific_heat_plastic) ! m^3/s
    ! https://en.wikipedia.org/wiki/Polycarbonate

    ! yellow pine
    REAL(RK), PARAMETER :: thermal_conductivity_wood = 0.15_rk ! W/(m*K)
    REAL(RK), PARAMETER :: density_wood = 640.0_rk             ! kg/m^3
    REAL(RK), PARAMETER :: specific_heat_wood = 2805.0_rk      ! J/(kg*K)
    REAL(RK), PARAMETER :: thermal_diffusivity_wood = thermal_conductivity_wood / (density_wood * specific_heat_wood) ! m^3/s
    REAL(RK), PARAMETER :: wood_thickness = 0.02_rk            ! m
END MODULE global

MODULE io
    use types_and_kinds
    IMPLICIT NONE 
    CONTAINS
    SUBROUTINE write_grid(grid, filename)
        REAL(RK), DIMENSION(:, :, :), intent(in) :: grid
        CHARACTER(len=*), INTENT(in) :: filename
        ! need to use stream as access for the binary (WHY?)
        ! maybe retarded but name can not be *.bin
        open(1, file=filename, form='unformatted', access='stream', status='replace')
        write(1) grid 
        close(1) 
    END SUBROUTINE write_grid

    SUBROUTINE open_inputfile(file_path, file_unit, iostat)
        !! Check whether file exists, with consitent error message
        !! return the file unit
        CHARACTER(len=*),  INTENT(in)  :: file_path
        INTEGER,  INTENT(out) :: file_unit, iostat

        inquire (file=file_path, iostat=iostat)
        if (iostat /= 0) then
            write (*, '(3a)') 'Error: file "', trim(file_path), '" not found!'
        end if
        open (action='read', file=file_path, iostat=iostat, newunit=file_unit)
    END SUBROUTINE open_inputfile
    SUBROUTINE close_inputfile(file_path, file_unit, iostat)
        !! Check the reading was OK
        !! return error line IF not
        !! close the unit
        character(len=*),  intent(in)  :: file_path
        character(len=1000) :: line
        integer,  intent(in) :: file_unit, iostat

        if (iostat /= 0) then
            write (*, '(2a)') 'Error reading file :"', trim(file_path)
            write (*, '(a, i0)') 'iostat was:"', iostat
            backspace(file_unit)
            read(file_unit,fmt='(A)') line
            write(*,'(A)') &
                'Invalid line : '//trim(line)
        end if
        close (file_unit)   
    END SUBROUTINE close_inputfile
    
    SUBROUTINE read_inputs(filepath, & 
                            Nx, Ny, Nz, nGhosts, wall_thickness, computer_power,outside_temperature, & 
                            dt, ds, simrealtime, save_interval) 
    
        CHARACTER(len=*), INTENT(IN) :: filepath 
        INTEGER(ik), INTENT(INOUT) :: Nx, Ny, Nz,nGhosts, save_interval
        REAL(RK), INTENT(INOUT) :: wall_thickness, ds, dt, simrealtime, outside_temperature
        REAL(RK) :: computer_power
        integer :: file_unit, iostat
        namelist /GREENHOUSE/ Nx, Ny, Nz,nGhosts, wall_thickness, computer_power, outside_temperature
        namelist /SIMULATION/ dt, ds, simrealtime, save_interval 
        
        call open_inputfile(filepath, file_unit, iostat)
        if (iostat /= 0) then
            stop
        end if
        read (nml=GREENHOUSE, iostat=iostat, unit=file_unit)
        read (nml=SIMULATION, iostat=iostat, unit=file_unit)
        call close_inputfile(filepath, file_unit, iostat)
        if (iostat /= 0) then
            stop
        end if
    END SUBROUTINE read_inputs

END MODULE io

MODULE initialization
    use types_and_kinds
    IMPLICIT NONE
    CONTAINS 
    !-----------------------------------------------------------------------
    ! FUNCTION: GAUSSIAN3D
    ! Purpose: Return value of 3D gaussian at point (i, j, k)
    !        centered at mu, with standard deviation sigma
    !-----------------------------------------------------------------------
    REAL(RK) FUNCTION GAUSSIAN3D(i, j, k, mu, sigma)
        REAL(RK), intent(in) :: mu, sigma
        INTEGER(ik), intent(in) :: i, j, k
        GAUSSIAN3D = exp(-((i-mu)**2 + (j-mu)**2 + (k-mu)**2)/(2*sigma**2))
    END FUNCTION GAUSSIAN3D
    !-----------------------------------------------------------------------
    ! SUBROUTINE: init_grid_gaussian
    ! Purpose: Initialize grid with a gaussian of some magnitude and offset
    !         centered at init_value_mu, with standard deviation init_value_sigma
    !         magnitude is the amplitude of gaussian, offset is offset from 0 for rest of points
    !         this is done using index coordinates, not physical coordinates
    !-----------------------------------------------------------------------
    SUBROUTINE init_grid_gaussian(grid, Nx, Ny, Nz, init_value_mu, init_value_sigma, magnitude, offset)
        REAL(RK), dimension(:, :, :), intent(out) :: grid
        INTEGER(ik), intent(in) :: Nx, Ny, Nz
        REAL(RK), intent(in) :: init_value_mu, init_value_sigma, magnitude, offset
        INTEGER(ik) :: i, j, k
        DO k = 1, Nz
            DO j = 1, Ny
                DO i = 1, Nx
                    grid(i, j, k) = magnitude*GAUSSIAN3D(i, j, k, init_value_mu, init_value_sigma) + offset
                END DO
            END DO
        END DO
    END SUBROUTINE init_grid_gaussian
END MODULE initialization

MODULE physics 
    use types_and_kinds
    use global 
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
        diffusivity_grid = thermal_diffusivity_air
        ! plastic border
        diffusivity_grid(nGhosts+1:Nx+1:Nx-1, nGhosts+1:Ny+1, nGhosts+1:Nz+1) = thermal_diffusivity_plastic
        diffusivity_grid(nGhosts+1:Nx+1, nGhosts+1:Ny+1:Ny-1, nGhosts+1:Nz+1) = thermal_diffusivity_plastic
        diffusivity_grid(nGhosts+1:Nx+1, nGhosts+1:Ny+1, nGhosts+1:Nz+1:Nz-1) = thermal_diffusivity_plastic
        ! replace small component where computers are with air 
        diffusivity_grid(1+nGhosts:Nx/2, 1:1+nGhosts, :) = thermal_diffusivity_air

        ! replace bottom z layer as wood 
        diffusivity_grid(nGhosts+1:Nx+1, nGhosts+1:Ny+1, nGhosts+1) = thermal_diffusivity_wood
        print *, thermal_diffusivity_wood, thermal_diffusivity_plastic, thermal_diffusivity_air
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
        ! sources(lamp_x:lamp_x+lamp_len_x, lamp_y:lamp_y+lamp_len_y, lamp_z:lamp_z+lamp_len_z) = lamp_power_density
        
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
        INTEGER(IK) :: i, j, k
        DO k = 1+Nghosts, Nz+nGhosts
            DO j = 1+Nghosts, Ny+nGhosts
                DO i = 1+Nghosts, Nx+nGhosts
                    T_new(i, j, k) = T(i, j, k) + (diffusivity_grid(i, j, k)*dt)*( & 
                                                    (T(i+1, j, k) - 2*T(i, j, k) + T(i-1, j, k) + &
                                                   T(i, j+1, k) - 2*T(i, j, k) + T(i, j-1, k) + &
                                                   T(i, j, k+1) - 2*T(i, j, k) + T(i, j, k-1))/(wall_thickness_grid(i,j,k)**2) + &
                                                    sources(i, j, k)* (ds**2/ (thermal_condctivity_air)))
                    !if (i == 1+nGhosts .OR. i==Nx+nGhosts .OR.& 
                    !    j==1+nGhosts .OR. j==Ny+nGhosts .OR.&
                    !    k==1+nGhosts .OR. k==Nz+nGhosts) then
                    !        T_new(i, j, k) = T(i, j, k) + (thermal_diffusivity_plastic*dt)*( & 
                    !                            (T(i+1, j, k) - 2*T(i, j, k) + T(i-1, j, k) + &
                    !                               T(i, j+1, k) - 2*T(i, j, k) + T(i, j-1, k) + &
                    !                            T(i, j, k+1) - 2*T(i, j, k) + T(i, j, k-1))/(wall_thickness**2) + & 
                    !                            sources(i, j, k)* (ds**2/ (thermal_condctivity_air)))! + &
                    !else 
                    !T_new(i, j, k) = T(i, j, k) + (thermal_diffusivity_air*dt)*( & 
                    !                                (T(i+1, j, k) - 2*T(i, j, k) + T(i-1, j, k) + &
                    !                               T(i, j+1, k) - 2*T(i, j, k) + T(i, j-1, k) + &
                    !                               T(i, j, k+1) - 2*T(i, j, k) + T(i, j, k-1))/(ds**2) + &
                    !                                sources(i, j, k)* (ds**2/ (thermal_condctivity_air)))
                    !end if                            
                END DO
            END DO
        END DO
    END SUBROUTINE conductive_heat_transfer
END MODULE physics
PROGRAM conductive_heat
    use types_and_kinds
    use global 
    use io 
    use initialization
    use physics 
    IMPLICIT NONE 
    REAL(RK), DIMENSION(:, :, :), ALLOCATABLE :: T, T_new, sources, ds_grid, diffusivity_grid 
    INTEGER(IK) :: Nx, Ny, Nz, nGhosts
    INTEGER(ik) :: timestep=1, saveinterval
    REAL(RK) :: dx, dt, wall_thickness, simrealtime, computer_power, outside_temperature
    REAL(RK), PARAMETER :: epsilon=1.5e-4 ! convergence criteria
    character(len=1000) :: filepath
    filepath = '/home/kitadam/Uni/2023/Sci_Computing_II/final_project/deep-green/conduction/run/conductive_input.nml'
    call read_inputs(filepath, Nx, Ny, Nz,nGhosts, wall_thickness, computer_power, outside_temperature,& 
                               dt, dx, simrealtime, saveinterval)
    ALLOCATE(T(Nx+2*nGhosts, Ny+2*nGhosts, Nz+2*nGhosts), T_new(Nx+2*nGhosts, Ny+2*nGhosts, Nz+2*nGhosts), & 
                sources(Nx+2*nGhosts, Ny+2*nGhosts, Nz+2*nGhosts), ds_grid(Nx+2*nGhosts, Ny+2*nGhosts, Nz+2*nGhosts), &
                diffusivity_grid(Nx+2*nGhosts, Ny+2*nGhosts, Nz+2*nGhosts))
    
    ! Initial temperature around outside temperature
    T = outside_temperature + 2.0_rk + c_to_k
    ! T at nghost points is outside temperature 
    T(1:Nx+2*nGhosts:Nx+nGhosts, :, :) = outside_temperature + c_to_k
    T(:, 1:Ny+2*nGhosts:Ny+nGhosts, :) = outside_temperature + c_to_k
    T(:, :, 1:Nz+2*nGhosts:Nz+nGhosts) = outside_temperature + c_to_k

    ! Temperature of plastic 
    T(2:Nx+1:Nx-1, 2:Ny+1, 2:Nz+1) = outside_temperature + 1.0_rk + c_to_k
    T(2:Nx+1, 2:Ny+1:Ny-1, 2:Nz+1) = outside_temperature + 1.0_rk + c_to_k
    T(2:Nx+1, 2:Ny+1, 2:Nz+1:Nz-1) = outside_temperature + 1.0_rk + c_to_k
    
    T_new = T
    CALL get_source_grid(sources, Nx, Ny, Nz, dx, nGhosts, computer_power)
    CALL get_diffusivity_grid(diffusivity_grid, Nx, Ny, Nz, nGhosts)
    call get_wall_thickness_grid(ds_grid, Nx, Ny, Nz, nGhosts, dx, wall_thickness)
    CALL write_grid(T + k_to_c, 'T0.dat')
    CALL write_grid(sources, 'sources.dat')
    CALL write_grid(ds_grid, 'thickness.dat' )
    CALL write_grid(diffusivity_grid, 'diffusivity.dat')
    do while (timestep*dt < simrealtime) 
        timestep = timestep + 1
        CALL conductive_heat_transfer(T, T_new, Nx, Ny, Nz, nGhosts, dt, dx, sources, ds_grid, diffusivity_grid)
        if (maxval(abs(T_new-T)) < epsilon) then 
            print *, "converged, exiting"
            exit 
        end if 
        
        if (mod(timestep, saveinterval) == 0) then 
            print *, timestep, timestep*dt / (60_rk*60_rk), maxval(abs(T_new-T))
        end if 
        T = T_new
        
    end do 
    CALL write_grid(T + k_to_c, 'T.dat')
    DEALLOCATE(T, T_new, sources, ds_grid, diffusivity_grid)

END PROGRAM conductive_heat
