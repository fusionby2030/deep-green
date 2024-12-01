""" 
Warp implementation of the CFD code. 
"""

import warp as wp 
import numpy as np 
import h5py 



# @wp.kernel 
# def calc_cfl_condition(rho: wp.array3d(dtype=wp.float32), 
#                        ux: wp.array3d(dtype=wp.float32),
#                        uy: wp.array3d(dtype=wp.float32),
#                        uz: wp.array3d(dtype=wp.float32),
#                        p: wp.array3d(dtype=wp.float32),
#                        ds: wp.float32, 
#                        dt: wp.array(dtype=wp.float32)):
#     i, j, k = wp.tid()
#     # Does'nt matter because we will loop and only have on thread id 
#     inner = wp.float32(1e10)
#     for _i in range(NGHOST, NX - NGHOST):
#         for _j in range(NGHOST, NY - NGHOST):
#             for _k in range(NGHOST, NZ - NGHOST):
#                 temp = ds / (wp.pow( (GAMMA*p[_i, _j, _k] / rho[_i, _j, _k]), 0.5) + wp.pow(ux[_i, _j, _k]*ux[_i, _j, _k] + uy[_i, _j, _k]*uy[_i, _j, _k] + uz[_i, _j, _k]*uz[_i, _j, _k], 0.5))
#                 inner = wp.min(temp, inner)
#     dt[0] = inner*CFL



@wp.func
def get_jump_index(loc: wp.int32, n: wp.int32): 
    if loc == 0: 
        return n - NGHOST - 2
    if loc == 1: 
        return n - NGHOST - 1
    if loc == n - 1: 
        return NGHOST
    if loc == n - 2:
        return NGHOST + 1
    
@wp.kernel 
def apply_primitive_bdnry_conds_x(rho: wp.array3d(dtype=wp.float32), 
                                ux: wp.array3d(dtype=wp.float32),
                                uy: wp.array3d(dtype=wp.float32),
                                uz: wp.array3d(dtype=wp.float32),
                                p: wp.array3d(dtype=wp.float32),
                                ): 
    
    """ Just periodic boundary conditions """
    i, j, k = wp.tid()

    """ Skipping interior cells"""
    if i >= NGHOST and i < NX - NGHOST:
        return 
    # Define a jump index
    jump_index = get_jump_index(i, NX)

    """ Applying boundary conditions """
    rho[i, j, k] = rho[jump_index, j, k]
    ux[i, j, k] = ux[jump_index, j, k]
    uy[i, j, k] = uy[jump_index, j, k]
    uz[i, j, k] = uz[jump_index, j, k]
    p[i, j, k] = p[jump_index, j, k]

@wp.kernel 
def apply_primitive_grads_bdnry_conds_x(drho: wp.array4d(dtype=wp.float32), 
                                dux: wp.array4d(dtype=wp.float32),
                                duy: wp.array4d(dtype=wp.float32),
                                duz: wp.array4d(dtype=wp.float32),
                                dp: wp.array4d(dtype=wp.float32),
                                ): 
    
    """ Just periodic boundary conditions """
    i, j, k = wp.tid()

    """ Skipping interior cells"""
    if i >= NGHOST and i < NX - NGHOST:
        return 
    # Define a jump index
    jump_index = get_jump_index(i, NX)

    """ Applying boundary conditions """
    for l in range(3): 
        drho[i, j, k, l] = drho[jump_index, j, k, l]
        dux[i, j, k, l] = dux[jump_index, j, k, l]
        duy[i, j, k, l] = duy[jump_index, j, k, l]
        duz[i, j, k, l] = duz[jump_index, j, k, l]
        dp[i, j, k, l] = dp[jump_index, j, k, l]

@wp.kernel 
def apply_primitive_grads_bdnry_conds_y(drho: wp.array4d(dtype=wp.float32), 
                                dux: wp.array4d(dtype=wp.float32),
                                duy: wp.array4d(dtype=wp.float32),
                                duz: wp.array4d(dtype=wp.float32),
                                dp: wp.array4d(dtype=wp.float32),
                                ): 
    
    """ Just periodic boundary conditions """
    i, j, k = wp.tid()

    """ Skipping interior cells"""
    if j >= NGHOST and j < NY - NGHOST:
        return

    # Define a jump index
    jump_index = get_jump_index(j, NY)

    """ Applying boundary conditions """
    for l in range(3): 
        drho[i, j, k, l] = drho[i, jump_index, k, l]
        dux[i, j, k, l] = dux[i, jump_index, k, l]
        duy[i, j, k, l] = duy[i, jump_index, k, l]
        duz[i, j, k, l] = duz[i, jump_index, k, l]
        dp[i, j, k, l] = dp[i, jump_index, k, l]
    
@wp.kernel 
def apply_primitive_grads_bdnry_conds_z(drho: wp.array4d(dtype=wp.float32),
                                dux: wp.array4d(dtype=wp.float32),
                                duy: wp.array4d(dtype=wp.float32),
                                duz: wp.array4d(dtype=wp.float32),
                                dp: wp.array4d(dtype=wp.float32),
                                ): 
    
    """ Just periodic boundary conditions """
    i, j, k = wp.tid()

    """ Skipping interior cells"""
    if k >= NGHOST and k < NZ - NGHOST:
        return

    # Define a jump index
    jump_index = get_jump_index(k, NZ)

    """ Applying boundary conditions """
    for l in range(3): 
        drho[i, j, k, l] = drho[i, j, jump_index, l]
        dux[i, j, k, l] = dux[i, j, jump_index, l]
        duy[i, j, k, l] = duy[i, j, jump_index, l]
        duz[i, j, k, l] = duz[i, j, jump_index, l]
        dp[i, j, k, l] = dp[i, j, jump_index, l]    
    

@wp.kernel
def apply_primitive_bdnry_conds_y(rho: wp.array3d(dtype=wp.float32), 
                                ux: wp.array3d(dtype=wp.float32),
                                uy: wp.array3d(dtype=wp.float32),
                                uz: wp.array3d(dtype=wp.float32),
                                p: wp.array3d(dtype=wp.float32),
                                ): 
    i, j, k = wp.tid()

    if j >= NGHOST and j < NY - NGHOST:
        return
    
    jump_index = get_jump_index(j, NY)


    rho[i, j, k] = rho[i, jump_index, k]
    ux[i, j, k] = ux[i, jump_index, k]
    uy[i, j, k] = uy[i, jump_index, k]
    uz[i, j, k] = uz[i, jump_index, k]
    p[i, j, k] = p[i, jump_index, k]

@wp.kernel
def apply_primitive_bdnry_conds_z(rho: wp.array3d(dtype=wp.float32), 
                                ux: wp.array3d(dtype=wp.float32),
                                uy: wp.array3d(dtype=wp.float32),
                                uz: wp.array3d(dtype=wp.float32),
                                p: wp.array3d(dtype=wp.float32),
                                ): 
    i, j, k = wp.tid()

    if k >= NGHOST and k < NZ - NGHOST:
        return
    
    jump_index = get_jump_index(k, NZ)

    rho[i, j, k] = rho[i, j, jump_index]
    ux[i, j, k] = ux[i, j, jump_index]
    uy[i, j, k] = uy[i, j, jump_index]
    uz[i, j, k] = uz[i, j, jump_index]
    p[i, j, k] = p[i, j, jump_index]

@wp.kernel 
def extrapolate_primitives_time(rho_extrp: wp.array3d(dtype=wp.float32),
                                ux_extrp: wp.array3d(dtype=wp.float32),
                                uy_extrp: wp.array3d(dtype=wp.float32),
                                uz_extrp: wp.array3d(dtype=wp.float32),
                                p_extrp: wp.array3d(dtype=wp.float32),
                                rho: wp.array3d(dtype=wp.float32),
                                ux: wp.array3d(dtype=wp.float32),
                                uy: wp.array3d(dtype=wp.float32),
                                uz: wp.array3d(dtype=wp.float32),
                                p: wp.array3d(dtype=wp.float32),
                                drho: wp.array4d(dtype=wp.float32),
                                dux: wp.array4d(dtype=wp.float32),
                                duy: wp.array4d(dtype=wp.float32),
                                duz: wp.array4d(dtype=wp.float32),
                                dp: wp.array4d(dtype=wp.float32),
                                dt: wp.float32,): 
    i, j, k = wp.tid()

    drhodx = drho[i, j, k, 0]
    drhody = drho[i, j, k, 1]
    drhodz = drho[i, j, k, 2]

    duxdx = dux[i, j, k, 0]
    duxdy = dux[i, j, k, 1]
    duxdz = dux[i, j, k, 2]

    duydx = duy[i, j, k, 0]
    duydy = duy[i, j, k, 1]
    duydz = duy[i, j, k, 2]

    duzdx = duz[i, j, k, 0]
    duzdy = duz[i, j, k, 1]
    duzdz = duz[i, j, k, 2]

    dpdx = dp[i, j, k, 0]
    dpdy = dp[i, j, k, 1]
    dpdz = dp[i, j, k, 2]

    rho_extrp[i,j,k] = rho[i,j,k] - 0.5 * dt * (ux[i,j,k]*drhodx + rho[i,j,k]*duxdx + uy[i,j,k]*drhody + rho[i,j,k]*duydy + uz[i,j,k]*drhodz + rho[i,j,k]*duzdz)
    ux_extrp[i,j,k]  = ux[i,j,k]  - 0.5 * dt * (ux[i,j,k]*duxdx + uy[i,j,k]*duxdy + uz[i,j,k]*duxdz + (1.0 / rho[i,j,k]) * dpdx)
    uy_extrp[i,j,k]  = uy[i,j,k]  - 0.5 * dt * (ux[i,j,k]*duydx + uy[i,j,k]*duydy + uz[i,j,k]*duydz + (1.0 / rho[i,j,k]) * dpdy)
    uz_extrp[i,j,k]  = uz[i,j,k]  - 0.5 * dt * (ux[i,j,k]*duzdx + uy[i,j,k]*duzdy + uz[i,j,k]*duzdz + (1.0 / rho[i,j,k]) * dpdz)
    p_extrp[i,j,k]   = p[i,j,k]   - 0.5 * dt * (GAMMA * p[i,j,k]* (duxdx + duydy + duzdz) + ux[i,j,k]*dpdx + uy[i,j,k]*dpdy + uz[i,j,k]*dpdz)

@wp.kernel
def calculate_primitives(
    rho: wp.array3d(dtype=wp.float32), 
    ux: wp.array3d(dtype=wp.float32),
    uy: wp.array3d(dtype=wp.float32),
    uz: wp.array3d(dtype=wp.float32),
    p: wp.array3d(dtype=wp.float32),
    mass: wp.array3d(dtype=wp.float32),
    etot: wp.array3d(dtype=wp.float32),
    momentum_x: wp.array3d(dtype=wp.float32),
    momentum_y: wp.array3d(dtype=wp.float32),
    momentum_z: wp.array3d(dtype=wp.float32),
    ds: wp.float32
): 
    cell_volume = ds*ds*ds
    i, j, k = wp.tid() 
    # for i in range(1, NX-1): 
    #     for j in range(1, NY-1): 
    #         for k in range(1, NZ-1): 
    rho[i, j, k] = mass[i, j, k] / cell_volume
    ux[i, j, k] = momentum_x[i, j, k] / (rho[i, j, k] * cell_volume)
    uy[i, j, k] = momentum_y[i, j, k] / (rho[i, j, k] * cell_volume)
    uz[i, j, k] = momentum_z[i, j, k] / (rho[i, j, k] * cell_volume)
    p[i, j, k]  = (etot[i,j,k] / cell_volume - 0.5 * rho[i,j,k] * (ux[i,j,k]*ux[i,j,k] + uy[i,j,k]*uy[i,j,k] + uz[i,j,k]*uz[i,j,k])) * (GAMMA - 1.0)

@wp.kernel 
def calc_conserved(
    rho: wp.array3d(dtype=wp.float32), 
    ux: wp.array3d(dtype=wp.float32),
    uy: wp.array3d(dtype=wp.float32),
    uz: wp.array3d(dtype=wp.float32),
    p: wp.array3d(dtype=wp.float32),
    mass: wp.array3d(dtype=wp.float32),
    etot: wp.array3d(dtype=wp.float32),
    momentum_x: wp.array3d(dtype=wp.float32),
    momentum_y: wp.array3d(dtype=wp.float32),
    momentum_z: wp.array3d(dtype=wp.float32),
    ds: wp.float32
): 
    i, j, k = wp.tid() 
    cell_volume = ds*ds*ds
    # for i in range(NX): 
    #     for j in range(NY): 
    #         for k in range(NZ): 
    mass[i, j, k] = rho[i, j, k] * cell_volume
    momentum_x[i, j, k] = rho[i, j, k] * ux[i, j, k] * cell_volume
    momentum_y[i, j, k] = rho[i, j, k] * uy[i, j, k] * cell_volume
    momentum_z[i, j, k] = rho[i, j, k] * uz[i, j, k] * cell_volume
    etot[i, j, k] = 0.5 * rho[i, j, k] * (ux[i, j, k]*ux[i, j, k] + uy[i, j, k]*uy[i, j, k] + uz[i, j, k]*uz[i, j, k]) * cell_volume + p[i, j, k] / (GAMMA - 1.0) * cell_volume 

@wp.kernel 
def calc_gradients(grid_mat: wp.array3d(dtype=wp.float32), 
                   grid_mat_grad: wp.array4d(dtype=wp.float32),
                   ds: wp.float32):
    i, j, k = wp.tid()
    # Check inside the domain 
    if i >= NGHOST and i < NX - NGHOST: 
        if j >= NGHOST and j < NY - NGHOST: 
            if k >= NGHOST and k < NZ - NGHOST: 
                grid_mat_grad[i, j, k, 0] = (grid_mat[i+1, j, k] - grid_mat[i-1, j, k]) / (2.0*ds)
                grid_mat_grad[i, j, k, 1] = (grid_mat[i, j+1, k] - grid_mat[i, j-1, k]) / (2.0*ds)
                grid_mat_grad[i, j, k, 2] = (grid_mat[i, j, k+1] - grid_mat[i, j, k-1]) / (2.0*ds)
    


@wp.kernel 
def calc_flux_reconstruction(mass_flux: wp.array4d(dtype=wp.float32),
                             momentum_x_flux: wp.array4d(dtype=wp.float32),
                            momentum_y_flux: wp.array4d(dtype=wp.float32),
                            momentum_z_flux: wp.array4d(dtype=wp.float32),
                            etot_flux: wp.array4d(dtype=wp.float32),
                            drho: wp.array4d(dtype=wp.float32),
                            dux: wp.array4d(dtype=wp.float32),
                            duy: wp.array4d(dtype=wp.float32),
                            duz: wp.array4d(dtype=wp.float32),
                            dp: wp.array4d(dtype=wp.float32),
                            rho: wp.array3d(dtype=wp.float32),
                            ux: wp.array3d(dtype=wp.float32),
                            uy: wp.array3d(dtype=wp.float32),
                            uz: wp.array3d(dtype=wp.float32),
                            p: wp.array3d(dtype=wp.float32),
                            offsets: wp.array1d(dtype=wp.int32),
                            dim: wp.int32,
                            ds: wp.float32):
    i, j, k = wp.tid()

    # mass_flux_view = mass_flux[i, j, k, dim]
    # momentum_x_flux_view = momentum_x_flux[i, j, k, dim]
    # momentum_y_flux_view = momentum_y_flux[i, j, k, dim]
    # momentum_z_flux_view = momentum_z_flux[i, j, k, dim]
    # etot_flux_view = etot_flux[i, j, k, dim]

    if not (i >= NGHOST and i < NX - NGHOST):
        return
    if not (j >= NGHOST and j < NY - NGHOST):
        return
    if not (k >= NGHOST and k < NZ - NGHOST):
        return

    rho_l = rho[i + offsets[0], j + offsets[1], k + offsets[2]] - (drho[i + offsets[0], j + offsets[1], k + offsets[2], dim] * (ds / 2.0))
    rho_r = rho[i, j, k] + (drho[i, j, k, dim] * (ds / 2.0))

    ux_l = ux[i + offsets[0], j + offsets[1], k + offsets[2]] - (dux[i + offsets[0], j + offsets[1], k + offsets[2], dim] * (ds / 2.0))
    ux_r = ux[i, j, k] + (dux[i, j, k, dim] * (ds / 2.0))

    uy_l = uy[i + offsets[0], j + offsets[1], k + offsets[2]] - (duy[i + offsets[0], j + offsets[1], k + offsets[2], dim] * (ds / 2.0))
    uy_r = uy[i, j, k] + (duy[i, j, k, dim] * (ds / 2.0))

    uz_l = uz[i + offsets[0], j + offsets[1], k + offsets[2]] - (duz[i + offsets[0], j + offsets[1], k + offsets[2], dim] * (ds / 2.0))
    uz_r = uz[i, j, k] + (duz[i, j, k, dim] * (ds / 2.0))

    p_l = p[i + offsets[0], j + offsets[1], k + offsets[2]] - (dp[i + offsets[0], j + offsets[1], k + offsets[2], dim] * (ds / 2.0))
    p_r = p[i, j, k] + (dp[i, j, k, dim] * (ds / 2.0))

    en_l = 0.5 * (rho_l * (ux_l*ux_l + uy_l*uy_l + uz_l*uz_l)) + (p_l / (GAMMA - 1.0))
    en_r = 0.5 * (rho_r * (ux_r*ux_r + uy_r*uy_r + uz_r*uz_r)) + (p_r / (GAMMA - 1.0))

    # TODO: can use wp.lerp here, but this would require finding real space
    rho_star = 0.5 * (rho_l + rho_r)
    momentum_x_star = 0.5 * (rho_l * ux_l + rho_r * ux_r)
    momentum_y_star = 0.5 * (rho_l * uy_l + rho_r * uy_r)
    momentum_z_star = 0.5 * (rho_l * uz_l + rho_r * uz_r)
    etot_star = 0.5 * (en_l + en_r)
    p_star    = (GAMMA - 1.0) * (etot_star - 0.5 * (momentum_x_star*momentum_x_star + momentum_y_star*momentum_y_star + momentum_z_star*momentum_z_star) / rho_star)


    
    mass_flux[i, j, k, dim] = momentum_x_star # - 0.5 * c_star * (rho_l - rho_r)
    momentum_x_flux[i, j, k, dim] = (momentum_x_star*momentum_x_star / rho_star) + p_star # - 0.5 * c_star * (rho_l * ux_l - rho_r * ux_r)
    momentum_y_flux[i, j, k, dim] = (momentum_x_star*momentum_y_star / rho_star) # - 0.5 * c_star * (rho_l * uy_l - rho_r * uy_r)
    momentum_z_flux[i, j, k, dim] = (momentum_x_star*momentum_z_star / rho_star) # - 0.5 * c_star * (rho_l * uz_l - rho_r * uz_r)
    etot_flux[i, j, k, dim] = (etot_star + p_star) * momentum_x_star / rho_star # - 0.5 * c_star * (en_l - en_r)


    # TODO GRAVITY 
    if dim == 2 and GRAVITY: 
        pass 

    c_l = wp.pow(GAMMA * p_l / rho_l, 0.5) + wp.abs(ux_l)
    c_r = wp.pow(GAMMA * p_r / rho_r, 0.5) + wp.abs(ux_r)
    c_star = wp.max(c_l, c_r)


    mass_flux[i, j, k, dim] = momentum_x_star - 0.5 * c_star*(rho_l - rho_r) 
    momentum_x_flux[i,j,k,dim] -= c_star * (rho_l * ux_l - rho_r * ux_r) / 2.0
    momentum_y_flux[i,j,k,dim] -= c_star * (rho_l * uy_l - rho_r * uy_r) / 2.0
    momentum_z_flux[i,j,k,dim] -= c_star * (rho_l * uz_l - rho_r * uz_r) / 2.0
    etot_flux[i,j,k,dim] -= c_star * (en_l - en_r) / 2.0

    mass_flux[i,j,k, dim] *= ds
    momentum_x_flux[i, j, k, dim] *= ds
    momentum_y_flux[i, j, k, dim] *= ds
    momentum_z_flux[i, j, k, dim] *= ds
    etot_flux[i, j, k, dim] *= ds

@wp.kernel 
def add_fluxes(mass: wp.array3d(dtype=wp.float32),
               momentum_x: wp.array3d(dtype=wp.float32),
                momentum_y: wp.array3d(dtype=wp.float32),
                momentum_z: wp.array3d(dtype=wp.float32),
                etot: wp.array3d(dtype=wp.float32),
                mass_flux: wp.array4d(dtype=wp.float32),
                momentum_x_flux: wp.array4d(dtype=wp.float32),
                momentum_y_flux: wp.array4d(dtype=wp.float32),
                momentum_z_flux: wp.array4d(dtype=wp.float32),
                etot_flux: wp.array4d(dtype=wp.float32),
                dt: wp.float32,
                ds: wp.float32):
    i, j, k = wp.tid()

    if not (i >= NGHOST and i < NX - NGHOST):
        return
    if not (j >= NGHOST and j < NY - NGHOST):
        return
    if not (k >= NGHOST and k < NZ - NGHOST):
        return

    mass[i, j, k] -= dt * ds * (mass_flux[i, j, k, 0] - mass_flux[i-1, j, k, 0] + mass_flux[i, j, k, 1] - mass_flux[i, j-1, k, 1] + mass_flux[i, j, k, 2] - mass_flux[i, j, k-1, 2]) 
    momentum_x[i, j, k] -= dt * ds * (momentum_x_flux[i, j, k, 0] - momentum_x_flux[i-1, j, k, 0] + momentum_x_flux[i, j, k, 1] - momentum_x_flux[i, j-1, k, 1] + momentum_x_flux[i, j, k, 2] - momentum_x_flux[i, j, k-1, 2])
    momentum_y[i, j, k] -= dt * ds * (momentum_y_flux[i, j, k, 0] - momentum_y_flux[i-1, j, k, 0] + momentum_y_flux[i, j, k, 1] - momentum_y_flux[i, j-1, k, 1] + momentum_y_flux[i, j, k, 2] - momentum_y_flux[i, j, k-1, 2])
    momentum_z[i, j, k] -= dt * ds * (momentum_z_flux[i, j, k, 0] - momentum_z_flux[i-1, j, k, 0] + momentum_z_flux[i, j, k, 1] - momentum_z_flux[i, j-1, k, 1] + momentum_z_flux[i, j, k, 2] - momentum_z_flux[i, j, k-1, 2])
    etot[i, j, k] -= dt * ds * (etot_flux[i, j, k, 0] - etot_flux[i-1, j, k, 0] + etot_flux[i, j, k, 1] - etot_flux[i, j-1, k, 1] + etot_flux[i, j, k, 2] - etot_flux[i, j, k-1, 2])


kappa = wp.constant(0.1)
s = wp.constant(0.05 / wp.sqrt(2.0))

@wp.kernel 
def initialize_khi(rho: wp.array3d(dtype=wp.float32), 
                   ux: wp.array3d(dtype=wp.float32),
                   uy: wp.array3d(dtype=wp.float32),
                   uz: wp.array3d(dtype=wp.float32),
                   p: wp.array3d(dtype=wp.float32), 
                   ds: wp.float32):
    
    i, j, k = wp.tid() # THe thread id on the kernel. 



    p[i,j,k] = 2.5 
    uy[i,j,k] = 0.0

    xs = wp.float32(i) * ds
    ys = 0.0 
    zs = wp.float32(k) * ds
    
    rho[i,j,k] = 1.0 
    ux[i,j,k] = -0.5
    if wp.abs(zs - ds*wp.float32(NZ) / 2.0) <= ds*wp.float32(NZ) / 4.0: 
        rho[i,j,k] += 1.0 # wp.float32(wp.abs(zs - ds*NZ / 2.0) < ds*NZ / 4.0)
        ux[i,j,k]  += 1.0
    
    uz[i,j,k]  = kappa * wp.sin(4.0 * wp.pi * xs) * (wp.exp(-wp.pow(zs - ds*wp.float32(NZ) / 4.0, 2.0) / (2.0*wp.pow(s, 2.0))) + wp.exp(-wp.pow(zs - 3.0*ds*wp.float32(NZ) / 4.0, 2.0) / (2.0*wp.pow(s, 2.0))))


def save_data(filename, 
              rho, ux, uy, uz, p, 
              mass, momentum_x, momentum_y, momentum_z, etot, 
              mass_flux, momentum_x_flux, momentum_y_flux, momentum_z_flux, etot_flux):
    with h5py.File(filename, 'w') as f: 
        f.create_dataset('rho', data=rho)
        f.create_dataset('ux', data=ux)
        f.create_dataset('uy', data=uy)
        f.create_dataset('uz', data=uz)
        f.create_dataset('p', data=p)
        f.create_dataset('mass', data=mass)
        f.create_dataset('momentum_x', data=momentum_x)
        f.create_dataset('momentum_y', data=momentum_y)
        f.create_dataset('momentum_z', data=momentum_z)
        f.create_dataset('etot', data=etot)
        for i, name in zip(range(3), ['_x', '_y', '_z']):
            f.create_dataset(f'mass_flux{name}', data=mass_flux[:,:,:,i])
            f.create_dataset(f'momentum_x_flux{name}', data=momentum_x_flux[:,:,:,i])
            f.create_dataset(f'momentum_y_flux{name}', data=momentum_y_flux[:,:,:,i])
            f.create_dataset(f'momentum_z_flux{name}', data=momentum_z_flux[:,:,:,i])
            f.create_dataset(f'etot_flux{name}', data=etot_flux[:,:,:,i]) 
            

if __name__ == "__main__":

    CFL   = wp.constant(0.2)
    GAMMA = wp.constant(5.0/3.0)
    NX, NY, NZ = 512, 4, 512
    NGHOST = 2
    GRAVITY = False
    WRITING = True

    NX, NY, NZ = wp.constant(NX + 2*NGHOST), wp.constant(NY + 2*NGHOST), wp.constant(NZ + 2*NGHOST)
    LX, LY, LZ = wp.constant(1.0), wp.constant(1.0), wp.constant(1.0)
    ds = wp.float32(LX / NX)

    CURRSIMTIME, TOTALSIMTIME = 0.0,  1.0
    WSTEP = 0.0025 
    WOUT  = 0.0
    # dt = wp.float32(0.5)
    device = wp.get_device()
    dt = wp.constant(1E-4)

    with wp.ScopedDevice(device): 


        """ Primitives"""
        rho = wp.array3d(np.ones((NX, NY, NZ)), dtype=wp.float32)
        ux = wp.array3d(np.zeros((NX, NY, NZ)), dtype=wp.float32)
        uy = wp.array3d(np.zeros((NX, NY, NZ)), dtype=wp.float32)
        uz = wp.array3d(np.zeros((NX, NY, NZ)), dtype=wp.float32)
        p = wp.array3d(np.ones((NX, NY, NZ)), dtype=wp.float32)

        """ Gradients """
        drho = wp.array4d(np.zeros((NX, NY, NZ, 3)), dtype=wp.float32)
        dux = wp.array4d(np.zeros((NX, NY, NZ, 3)), dtype=wp.float32)
        duy = wp.array4d(np.zeros((NX, NY, NZ, 3)), dtype=wp.float32)
        duz = wp.array4d(np.zeros((NX, NY, NZ, 3)), dtype=wp.float32)
        dp = wp.array4d(np.zeros((NX, NY, NZ, 3)), dtype=wp.float32)

        """ Extrapolation in time primtives """
        rho_extrp = wp.array3d(np.ones((NX, NY, NZ)), dtype=wp.float32)
        ux_extrp = wp.array3d(np.zeros((NX, NY, NZ)), dtype=wp.float32)
        uy_extrp = wp.array3d(np.zeros((NX, NY, NZ)), dtype=wp.float32)
        uz_extrp = wp.array3d(np.zeros((NX, NY, NZ)), dtype=wp.float32)
        p_extrp = wp.array3d(np.ones((NX, NY, NZ)), dtype=wp.float32)

        """ Conserved variables"""
        mass = wp.array3d(np.ones((NX, NY, NZ)), dtype=wp.float32)
        momentum_x = wp.array3d(np.zeros((NX, NY, NZ)), dtype=wp.float32)
        momentum_y = wp.array3d(np.zeros((NX, NY, NZ)), dtype=wp.float32)
        momentum_z = wp.array3d(np.zeros((NX, NY, NZ)), dtype=wp.float32)
        etot = wp.array3d(np.ones((NX, NY, NZ)), dtype=wp.float32)

        """ Fluxes """


        mass_flux = wp.array4d(np.zeros((NX, NY, NZ, 3)), dtype=wp.float32)
        momentum_x_flux = wp.array4d(np.zeros((NX, NY, NZ, 3)), dtype=wp.float32)
        momentum_y_flux = wp.array4d(np.zeros((NX, NY, NZ, 3)), dtype=wp.float32)
        momentum_z_flux = wp.array4d(np.zeros((NX, NY, NZ, 3)), dtype=wp.float32)
        etot_flux = wp.array4d(np.zeros((NX, NY, NZ, 3)), dtype=wp.float32)

        wp.launch(kernel=initialize_khi,
                    dim=(NX, NY, NZ),
                    inputs=(rho, ux, uy, uz, p, ds))

        """ COMPUTE """


        # Calculate the conseved quantities
        wp.launch(kernel=calc_conserved,
                    dim=(NX, NY, NZ),
                    inputs=(rho, ux, uy, uz, p, mass, etot, momentum_x, momentum_y, momentum_z, ds))

        index = 0
        fname = f"state_{index:07d}.h5"
        # save_data(fname, rho, ux, uy, uz, p, mass, momentum_x, momentum_y, momentum_z, etot, mass_flux, momentum_x_flux, momentum_y_flux, momentum_z_flux, etot_flux)

        index += 1
        # Loop: 
        while CURRSIMTIME < TOTALSIMTIME:     
            # if index >= 5: 
            #     break 
            # -- Calculate primitives 

            wp.launch(kernel=calculate_primitives, 
                    dim=(NX, NY, NZ), 
                    inputs=(rho, ux, uy, uz, p, mass, etot, momentum_x, momentum_y, momentum_z, ds))
            # -- Update primtiive ghost ceslls 
            # TODO: This should just be using the indicies of the ghost cells
            wp.launch(kernel=apply_primitive_bdnry_conds_x, 
                    dim=(NX, NY, NZ), 
                    inputs=(rho, ux, uy, uz, p))
            wp.launch(kernel=apply_primitive_bdnry_conds_y,
                    dim=(NX, NY, NZ),
                    inputs=(rho, ux, uy, uz, p))
            wp.launch(kernel=apply_primitive_bdnry_conds_z,
                    dim=(NX, NY, NZ),
                    inputs=(rho, ux, uy, uz, p))
            # -- Calculate CFL dt TODO this is retarded!

            # wp.launch(kernel = calc_cfl_condition, 
            #         dim=1, 
            #         inputs=(rho, ux, uy, uz, p, ds,dt),
            #         )
            # -- Calculate gradients of (rho, vx, vy, vz, p)
            wp.launch(kernel=calc_gradients, 
                    dim=(NX, NY, NZ), 
                    inputs=(rho, drho, ds))
            wp.launch(kernel=calc_gradients,
                        dim=(NX, NY, NZ),
                        inputs=(ux, dux, ds))
            wp.launch(kernel=calc_gradients,
                        dim=(NX, NY, NZ),
                        inputs=(uy, duy, ds))
            wp.launch(kernel=calc_gradients,
                        dim=(NX, NY, NZ),
                        inputs=(uz, duz, ds))
            wp.launch(kernel=calc_gradients,
                        dim=(NX, NY, NZ),
                        inputs=(p, dp, ds))
            
            # -- Update gradients ghost cells 
            wp.launch(kernel=apply_primitive_grads_bdnry_conds_x,
                        dim=(NX, NY, NZ),
                        inputs=(drho, dux, duy, duz, dp))
            wp.launch(kernel=apply_primitive_grads_bdnry_conds_y,
                        dim=(NX, NY, NZ),
                        inputs=(drho, dux, duy, duz, dp))
            wp.launch(kernel=apply_primitive_grads_bdnry_conds_z,
                        dim=(NX, NY, NZ),
                        inputs=(drho, dux, duy, duz, dp))

            # -- Extrapolate primitives
            wp.launch(kernel=extrapolate_primitives_time, 
                        dim=(NX, NY, NZ), 
                        inputs=(rho_extrp, ux_extrp, uy_extrp, uz_extrp, p_extrp, rho, ux, uy, uz, p, drho, dux, duy, duz, dp, dt))
            
            # -- Update primitive ghost cells
            wp.launch(kernel=apply_primitive_bdnry_conds_x,
                        dim=(NX, NY, NZ),
                        inputs=(rho_extrp, ux_extrp, uy_extrp, uz_extrp, p_extrp))
            wp.launch(kernel=apply_primitive_bdnry_conds_y,
                        dim=(NX, NY, NZ),
                        inputs=(rho_extrp, ux_extrp, uy_extrp, uz_extrp, p_extrp))
            wp.launch(kernel=apply_primitive_bdnry_conds_z,
                        dim=(NX, NY, NZ),
                        inputs=(rho_extrp, ux_extrp, uy_extrp, uz_extrp, p_extrp))
            
            # -- Calculate fluxes (in all directions)

            for offset, dim in zip([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 1, 2]): 
                offset = wp.array(offset, dtype=wp.int32)
                dim = int(dim)
                if dim == 0: 
                    inputs = (mass_flux, momentum_x_flux, momentum_y_flux, momentum_z_flux, etot_flux, drho, dux, duy, duz, dp, rho, ux, uy, uz, p, offset, dim, ds)
                elif dim == 1:
                    inputs = (mass_flux, momentum_y_flux, momentum_x_flux, momentum_z_flux, etot_flux, drho, duy, dux, duz, dp, rho, uy, ux, uz, p, offset, dim, ds)
                elif dim == 2:
                    inputs = (mass_flux, momentum_z_flux, momentum_x_flux, momentum_y_flux, etot_flux, drho, duz, dux, duy, dp, rho, uz, ux, uy, p, offset, dim, ds)
                wp.launch(kernel=calc_flux_reconstruction, 
                        dim=(NX, NY, NZ),
                        inputs=inputs)
                

            # -- update fluxes ghost cells
            wp.launch(kernel=apply_primitive_grads_bdnry_conds_x, 
                        dim=(NX, NY, NZ),
                        inputs=(mass_flux, momentum_x_flux, momentum_y_flux, momentum_z_flux, etot_flux))

            wp.launch(kernel=apply_primitive_grads_bdnry_conds_y,
                        dim=(NX, NY, NZ),
                        inputs=(mass_flux, momentum_x_flux, momentum_y_flux, momentum_z_flux, etot_flux))
            wp.launch(kernel=apply_primitive_grads_bdnry_conds_z,
                        dim=(NX, NY, NZ),
                        inputs=(mass_flux, momentum_x_flux, momentum_y_flux, momentum_z_flux, etot_flux))
            
            # -- Add fluxes to the conserved quantities
            wp.launch(kernel=add_fluxes, 
                        dim=(NX, NY, NZ),
                        inputs=(mass, momentum_x, momentum_y, momentum_z, etot, mass_flux, momentum_x_flux, momentum_y_flux, momentum_z_flux, etot_flux, dt, ds))
            


            CURRSIMTIME += dt
            WOUT += dt 
            

            if WOUT > WSTEP and WRITING: 
                # Fname should be state_XXXXXXX.h5 where XXXXXX is the current time index
                fname = f"state_{index:07d}.h5"
                save_data(fname, rho, ux, uy, uz, p, mass, momentum_x, momentum_y, momentum_z, etot, mass_flux, momentum_x_flux, momentum_y_flux, momentum_z_flux, etot_flux)
                WOUT = 0.0
                index += 1 
            print(f"{CURRSIMTIME:.10f} | {dt:.10f}")