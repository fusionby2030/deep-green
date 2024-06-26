!---------------------------------------------------------------------------------------
! Authors: Adam Kit and Kostis Papadakis (2023)
! Program: euler_cfd
! Solving the compressible Euler equations in 3 dimensions using the Lax–Friedrichs
! flux method. Stable for CFL<0.4.
! TODO:
! --Add a flux limiter using a high order reconstruction :
! probably use Lax–Wendroff for F_{h}
! --Supported Boundaries: Periodic, Slip Wall, Outflow
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

program euler_cfd
   use types_and_kinds
   use initialization
   use physics
   use io
   implicit none
   real(rk), dimension(:, :, :), allocatable  :: &
      rho, vx, vy, vz, p, mass, momentum_x, momentum_y, momentum_z, energy, &
      rho_prime, vx_prime, vy_prime, vz_prime, temp, &
      ! fluxes
      mass_flux_x, mass_flux_y, mass_flux_z, &
      momentum_x_flux_x, momentum_x_flux_y, momentum_x_flux_z, &
      momentum_y_flux_x, momentum_y_flux_y, momentum_y_flux_z, &
      momentum_z_flux_x, momentum_z_flux_y, momentum_z_flux_z, &
      energy_flux_x, energy_flux_y, energy_flux_z, &
      !gradients
      drho_dx, drho_dy, drho_dz, &
      dvx_dx, dvx_dy, dvx_dz, &
      dvy_dx, dvy_dy, dvy_dz, &
      dvz_dx, dvz_dy, dvz_dz, &
      dp_dx, dp_dy, dp_dz, rho_xtr, &
      ! extrapolated primitives
      vx_xtr, vy_xtr, vz_xtr, p_xtr

   integer(ik), parameter :: xcells = 128, &
                             ycells = 6, &
                             zcells = 256, &
                             nGhosts = 2
   integer(ik), parameter :: nx = xcells + 2*nGhosts, ny = ycells + 2*nGhosts, nz = zcells + 2*nGhosts
   real(rk), parameter :: ds = 1.0
   real(rk), parameter :: tout = 0.1_rk
   real(rk):: dt = 0.0_rk, time = 0.0_rk, write_time = 0.0_rk, time_max = 20.5_rk
   real(rk):: mag, cx = xcells/2, cy = ycells/2, cz = zcells/4! zcells/4 for thermal bubble
   integer(ik) :: timestep = 0, nWrites = 0, i, j, k
   integer(ik) :: shiftx(3), shifty(3), shiftz(3)
   integer(4):: BCs(6)

   !Set  boundary conditions
   BCs(1) = WALL !OUTFLOW !x-
   BCs(2) = WALL !x+
   BCs(3) = PERIODIC !y-
   BCs(4) = PERIODIC !y-
   BCs(5) = WALL !z-
   BCs(6) = WALL !z+

   allocate (rho(nx, ny, nz), vx(nx, ny, nz), vy(nx, ny, nz), vz(nx, ny, nz), p(nx, ny, nz), mass(nx, ny, nz), &
             momentum_x(nx, ny, nz), momentum_y(nx, ny, nz), momentum_z(nx, ny, nz), energy(nx, ny, nz), &
             rho_prime(nx, ny, nz), vx_prime(nx, ny, nz), vy_prime(nx, ny, nz), vz_prime(nx, ny, nz), temp(nx, ny, nz), &
             mass_flux_x(nx, ny, nz), mass_flux_y(nx, ny, nz), mass_flux_z(nx, ny, nz), &
             momentum_x_flux_x(nx, ny, nz), momentum_x_flux_y(nx, ny, nz), momentum_x_flux_z(nx, ny, nz), &
             momentum_y_flux_x(nx, ny, nz), momentum_y_flux_y(nx, ny, nz), momentum_y_flux_z(nx, ny, nz), &
             momentum_z_flux_x(nx, ny, nz), momentum_z_flux_y(nx, ny, nz), momentum_z_flux_z(nx, ny, nz), &
             energy_flux_x(nx, ny, nz), energy_flux_y(nx, ny, nz), energy_flux_z(nx, ny, nz), &
             drho_dx(nx, ny, nz), drho_dy(nx, ny, nz), drho_dz(nx, ny, nz), &
             dvx_dx(nx, ny, nz), dvx_dy(nx, ny, nz), dvx_dz(nx, ny, nz), &
             dvy_dx(nx, ny, nz), dvy_dy(nx, ny, nz), dvy_dz(nx, ny, nz), &
             dvz_dx(nx, ny, nz), dvz_dy(nx, ny, nz), dvz_dz(nx, ny, nz), &
             dp_dx(nx, ny, nz), dp_dy(nx, ny, nz), dp_dz(nx, ny, nz), rho_xtr(nx, ny, nz), &
             vx_xtr(nx, ny, nz), vy_xtr(nx, ny, nz), vz_xtr(nx, ny, nz), p_xtr(nx, ny, nz))

   shiftx(:) = 0
   shifty(:) = 0
   shiftz(:) = 0
   shiftx(1) = 1
   shifty(2) = 1
   shiftz(3) = 1
   !Set an initial state
   !call init_Uniform(rho, vx, vy, vz, p,101000.0_rk, 293.0_rk)
   !call init_Equilibrium(rho,vx,vy,vz,p, nx, ny, nz, nGhosts,ds,101000.0_rk,293.0_rk)
   ! call init_Kelvin_Helmholtz(rho,vx,vy,vz,p, nx, ny, nz, nGhosts,ds)
   ! call init_Thermal_Rising_Bubble(rho, vx, vy, vz, p, nx, ny, nz, nGhosts, ds, cx, cy, cz, 32._rk, 28._rk, 101000.0_rk, 293.0_rk)
   call init_adam_bubble(rho, vx, vy, vz, p, nx, ny, nz, nGhosts, ds, cx, cy, cz, 16._rk, 14._rk, 101000.0_rk, 293.0_rk)

   call update_primitive_ghosts(rho, vx, vy, vz, p, nx, ny, nz, nGhosts, BCs)
   call conservative(mass, momentum_x, momentum_y, momentum_z, energy, rho, p, vx, vy, vz, temp, ds)
   call write_state(time, nWrites, rho, vx, vy, vz, p, mass, momentum_x, momentum_y, momentum_z, energy, temp)
   nWrites = 1; 
   !main
   do while (time <= time_max)

      call apply_boundary_conditions(rho, vx, vy, vz, p, nx, ny, nz, nGhosts, BCs)

      call update_primitive_ghosts(rho, vx, vy, vz, p, nx, ny, nz, nGhosts, BCs)

      call primitive(mass, momentum_x, momentum_y, momentum_z, energy, rho, p, vx, vy, vz, temp, ds)

      call update_primitive_ghosts(rho, vx, vy, vz, p, nx, ny, nz, nGhosts, BCs)

      dt = compute_timestep(ds, vx, vy, vz, p, rho)

      call calculate_gradients(rho, drho_dx, drho_dy, drho_dz, nx, ny, nz, nGhosts, ds)
      call calculate_gradients(vx, dvx_dx, dvx_dy, dvx_dz, nx, ny, nz, nGhosts, ds)
      call calculate_gradients(vy, dvy_dx, dvy_dy, dvy_dz, nx, ny, nz, nGhosts, ds)
      call calculate_gradients(vz, dvz_dx, dvz_dy, dvz_dz, nx, ny, nz, nGhosts, ds)
      call calculate_gradients(p, dp_dx, dp_dy, dp_dz, nx, ny, nz, nGhosts, ds)

      call update_gradient_ghost(drho_dx, drho_dy, drho_dz, dvx_dx, dvx_dy, dvx_dz, dvy_dx, &
                                 dvy_dy, dvy_dz, dvz_dx, dvz_dy, dvz_dz, dp_dx, dp_dy, dp_dz, nx, ny, nz, nGhosts, BCs)

      call extrapolate_primitves(rho, vx, vy, vz, p, drho_dx, drho_dy, drho_dz, &
                                 dvx_dx, dvx_dy, dvx_dz, dvy_dx, dvy_dy, dvy_dz, &
                                 dvz_dx, dvz_dy, dvz_dz, dp_dx, dp_dy, dp_dz, &
                                 rho_xtr, vx_xtr, vy_xtr, vz_xtr, p_xtr, dt)

      !Apply BCs if we have anything other than Periodic or Outflow!
      call apply_boundary_conditions(rho_xtr, vx_xtr, vy_xtr, vz_xtr, p_xtr, nx, ny, nz, nGhosts, BCs)

      call reconstructflux(mass_flux_x, momentum_x_flux_x, momentum_y_flux_x, momentum_z_flux_x, &
                           energy_flux_x, drho_dx, dvx_dx, dvy_dx, dvz_dx, dp_dx, &
                           rho_xtr, vx_xtr, vy_xtr, vz_xtr, p_xtr, nx, ny, nz, nGhosts, ds, shiftx)

      call reconstructflux(mass_flux_y, momentum_y_flux_y, momentum_x_flux_y, momentum_z_flux_y, energy_flux_y, &
                           drho_dy, dvy_dy, dvx_dy, dvz_dy, dp_dy, &
                           rho_xtr, vy_xtr, vx_xtr, vz_xtr, p_xtr, nx, ny, nz, nGhosts, ds, shifty)

      call reconstructflux(mass_flux_z, momentum_z_flux_z, momentum_x_flux_z, momentum_y_flux_z, energy_flux_z, &
                           drho_dz, dvz_dz, dvx_dz, dvy_dz, dp_dz, &
                           rho_xtr, vz_xtr, vx_xtr, vy_xtr, p_xtr, nx, ny, nz, nGhosts, ds, shiftz)

      call update_flux_ghost(mass_flux_x, mass_flux_y, mass_flux_z, momentum_x_flux_x, momentum_y_flux_x, momentum_z_flux_x, &
                             momentum_x_flux_y, momentum_y_flux_y, momentum_z_flux_y, momentum_x_flux_z, momentum_y_flux_z, &
                             momentum_z_flux_z, energy_flux_x, energy_flux_y, energy_flux_z, nx, ny, nz, nGhosts, BCs)

      call addfluxes(mass_flux_x, mass_flux_y, mass_flux_z, &
                     momentum_x_flux_x, momentum_x_flux_y, momentum_x_flux_z, &
                     momentum_y_flux_x, momentum_y_flux_y, momentum_y_flux_z, &
                     momentum_z_flux_x, momentum_z_flux_y, momentum_z_flux_z, &
                     energy_flux_x, energy_flux_y, energy_flux_z, &
                     mass, momentum_x, momentum_y, momentum_z, energy, &
                     nx, ny, nz, nGhosts, dt, ds)

      timestep = timestep + 1
      if (write_time > tout) then
         print *, "Performing IO..."
         call write_state(time, nWrites, rho, vx, vy, vz, p, mass, momentum_x, momentum_y, momentum_z, energy, temp)
         print *, "IO done!"
         nWrites = nWrites + 1
         write_time = 0.0_rk
      end if

      print *, "Time=", time, "s | Timestep=", timestep, "| dt=", dt, "s"
      time = time + dt
      write_time = write_time + dt
   end do
   print *, "Run done go have fun now!"

   deallocate (rho, vx, vy, vz, p, mass, momentum_x, momentum_y, momentum_z, &
               energy, rho_prime, vx_prime, vy_prime, vz_prime, temp, &
               mass_flux_x, mass_flux_y, mass_flux_z, &
               momentum_x_flux_x, momentum_x_flux_y, momentum_x_flux_z, &
               momentum_y_flux_x, momentum_y_flux_y, momentum_y_flux_z, &
               momentum_z_flux_x, momentum_z_flux_y, momentum_z_flux_z, &
               energy_flux_x, energy_flux_y, energy_flux_z, &
               drho_dx, drho_dy, drho_dz, &
               dvx_dx, dvx_dy, dvx_dz, &
               dvy_dx, dvy_dy, dvy_dz, &
               dvz_dx, dvz_dy, dvz_dz, &
               dp_dx, dp_dy, dp_dz, rho_xtr, vx_xtr, vy_xtr, vz_xtr, p_xtr)
end program euler_cfd
