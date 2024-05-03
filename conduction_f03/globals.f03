module globals 
    use types_and_kinds 
    IMPLICIT NONE 
    ! GLOBAL VARIABLES
    REAL(RK), PARAMETER :: c_to_k = 273.15_rk, k_to_c = -273.15_rk 
    REAL(RK), PARAMETER :: rho_air = 1.225_RK ! kg/m^3
    REAL(rk), PARAMETER :: thermal_condctivity_air = 0.024 ! W/(m*K)
    ! http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/thrcn.html
    REAL(rk), PARAMETER :: density_air = 1.2922 ! kg/m^3
    ! https://www.earthdata.nasa.gov/topics/atmosphere/atmospheric-pressure/air-mass-density
    REAL(rk), PARAMETER :: specific_heat_air = 1003.5 ! J/(kg*K)
    ! https://en.wikipedia.org/wiki/Table_of_specific_heat_capacities
    REAL(rk), PARAMETER :: thermal_diffusivity_air = thermal_condctivity_air / (density_air * specific_heat_air) ! m^3/s
    REAL(RK), PARAMETER :: pi = 3.14159265358979323846_rk
end module globals 