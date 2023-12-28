SUBROUTINE run_conductive(T, Nx, Ny, Nz, nGhosts, dt, ds, sources, wall_thickness, t_max)
    REAL(4), PARAMETER :: c_to_k = 273.15, k_to_c = -273.15
    REAL(4), PARAMETER :: thermal_condctivity_air = 0.024 ! W/(m*K) ! http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/thrcn.html
    REAL(4), PARAMETER :: density_air = 1.2922 ! kg/m^3 ! https://www.earthdata.nasa.gov/topics/atmosphere/atmospheric-pressure/air-mass-density
    REAL(4), PARAMETER :: specific_heat_air = 1003.5 ! J/(kg*K) ! https://en.wikipedia.org/wiki/Table_of_specific_heat_capacities
    REAL(4), PARAMETER :: thermal_diffusivity_air = thermal_condctivity_air / (density_air * specific_heat_air) ! m^3/s
    REAL(4), PARAMETER :: thermal_condctivity_plastic = 0.23   ! W/(m*K)
    REAL(4), PARAMETER :: density_plastic = 1200.0          ! kg/m^3
    REAL(4), PARAMETER :: specific_heat_plastic = 1200.0    ! J/(kg*K)
    REAL(4), PARAMETER :: thermal_diffusivity_plastic = thermal_condctivity_plastic / (density_plastic * specific_heat_plastic)
    REAL(4), INTENT(IN) :: sources(:, :, :)
    INTEGER(4), INTENT(IN) :: Nx, Ny, Nz, nGhosts
    REAL(4), INTENT(INOUT), DIMENSION(:, :, :) :: T
    REAL(4), ALLOCATABLE, DIMENSION(:, :, :) :: T_new
    REAL(4), INTENT(IN) :: dt, ds, wall_thickness, t_max
    INTEGER :: i, j, k, timestep=0
    REAL(4) :: epsilon = 0.000001
    allocate(T_new(Nx+2*nGhosts, Ny+2*nGhosts, Nz+2*nGhosts))

    T = T + c_to_k
    T_new = T 
    DO WHILE(timestep*dt < t_max)
        
        DO k = 1+Nghosts, Nz+nGhosts
            DO j = 1+Nghosts, Ny+nGhosts
                DO i = 1+Nghosts, Nx+nGhosts
                    ! T(i, j, k) = k*Nz*Ny + j*Nz + i
                    ! if (sources(i,j,k) > 0.0) then 
                    !     print *, "source", i, j, k, sources(i,j,k)
                    ! end if 
                    if (i == 1+nGhosts .OR. i==Nx+nGhosts .OR.& 
                        j==1+nGhosts .OR. j==Ny+nGhosts .OR.&
                        k==1+nGhosts .OR. k==Nz+nGhosts) then
                            T_new(i, j, k) = T(i, j, k) + (thermal_diffusivity_plastic*dt)*( & 
                                                (T(i+1, j, k) - 2*T(i, j, k) + T(i-1, j, k) + &
                                                    T(i, j+1, k) - 2*T(i, j, k) + T(i, j-1, k) + &
                                                T(i, j, k+1) - 2*T(i, j, k) + T(i, j, k-1))/(wall_thickness**2) + & 
                                                sources(i, j, k)* (ds**2/ (thermal_condctivity_air)))! + &
                        
                    else 
                    T_new(i, j, k) = T(i, j, k) +(thermal_diffusivity_air*dt)*( & 
                                                    (T(i+1, j, k) - 2*T(i, j, k) + T(i-1, j, k) + &
                                                    T(i, j+1, k) - 2*T(i, j, k) + T(i, j-1, k) + &
                                                    T(i, j, k+1) - 2*T(i, j, k) + T(i, j, k-1))/(ds**2) + &
                                                    sources(i, j, k)* (ds**2/ (thermal_condctivity_air)))

                    end if                            
                END DO
            END DO
        END DO
        if (maxval(abs(T_new-T)) < epsilon) then 
            print *, "converged, exiting"
            exit 
        end if 
        T = T_new
        timestep = timestep + 1
    END DO    
    T = T + k_to_c
    print *, "Finished Conductive", maxval(T)
    print *, maxval(sources)
    deallocate(T_new)
END SUBROUTINE run_conductive