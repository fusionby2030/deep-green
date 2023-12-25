module global
   use types_and_kinds
   implicit none
   real(rk), parameter :: gamma = 5.0_rk/3.0_rk, rs = 287.0_rk, cfl = 0.2_rk
   integer(rk), parameter:: N_VARS = 5, N_FLUX = 15, N_DIMS = 3
   integer(4), parameter:: PERIODIC = 0, OUTFLOW = 1, WALL = 2, INFLOW = 3
end module global
