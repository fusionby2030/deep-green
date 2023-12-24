module initialization
   use types_and_kinds
   implicit none
contains

   subroutine init_grid(grid, nx, ny, nz, nGhosts, init_value)
      real(rk), dimension(:, :, :), intent(out) :: grid
      integer(ik), intent(in) :: nx, ny, nz, nGhosts
      real(rk), intent(in) :: init_value
      integer(ik) :: i, j, k
      ! to efficintly traverse the grid in the order of memory layout
      ! loop over 3d array with proper cache efficiency
      do k = nGhosts + 1, nz - nGhosts
         do j = nGhosts + 1, ny - nGhosts
            do i = nGhosts + 1, nx - nGhosts
               grid(i, j, k) = init_value
            end do
         end do
      end do
   end subroutine init_grid

   ! similar to above subroutine, we now intialize the grid with a 3d gaussian centered in middle of grid
   real(rk) function gaussian3d(i, j, k, mu, sigma)
      real(rk), intent(in) :: mu, sigma
      integer(ik), intent(in) :: i, j, k
      real(rk) :: a = 1.0
      !gaussian3d = a*sin(  2*dacos(-1.d0)* ((i-3)/127.0) )
      gaussian3d = a*exp(-((i - mu)**2 + (j - mu)**2 + (k - mu)**2)/(2*sigma**2))
   end function gaussian3d

   ! use the above function to initialize the grid
   subroutine init_grid_gaussian(grid, nx, ny, nz, nGhosts, init_value_mu, init_value_sigma, magnitude, offset)
      real(rk), dimension(:, :, :), intent(out) :: grid
      integer(ik), intent(in) :: nx, ny, nz, nGhosts
      real(rk), intent(in) :: init_value_mu, init_value_sigma, magnitude, offset
      integer(ik) :: i, j, k
      do k = nGhosts + 1, nz - nGhosts
         do j = nGhosts + 1, ny - nGhosts
            do i = nGhosts + 1, nx - nGhosts
               grid(i, j, k) = magnitude*gaussian3d(i, j, k, init_value_mu, init_value_sigma) + offset
            end do
         end do
      end do
   end subroutine init_grid_gaussian
end module initialization

