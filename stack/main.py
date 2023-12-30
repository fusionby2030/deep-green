import toolbox as tb 
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import scienceplots 
import os 
plt.style.use(['science'])
def plot_grid(grid): 
    stylegrid = grid.transpose()
    vmin = min(stylegrid.min(), stylegrid.min())
    vmax = max(stylegrid.max(), stylegrid.max())
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    fig, axs = plt.subplots(1, 4, figsize=(8, 8), dpi=200, sharey='row', sharex='row')
    # axs = axs.ravel()
    for n, Z_index in enumerate([int(i*Nz//4) for i in range(4)]): 
        cax = axs[n].imshow(stylegrid[Z_index+1, :, :], extent=[0, Nx*ds, Ny*ds, 0], norm=norm)
        axs[n].set_title(f'Z={Z_index*ds:.2}m')#  @ {simrealtime/(60.0*60.0):.3}h')
    fig.colorbar(cax, ax=axs.ravel().tolist(), location='right')
    plt.show() 

def initialize_temperature(T: np.ndarray, outside_temperature: float, Nx: int, Ny: int, Nz: int, nGhosts: int, from_save=False) -> np.ndarray: 
    if not from_save: 
        T = np.full((Nx+2*nGhosts, Ny+2*nGhosts, Nz+2*nGhosts), outside_temperature, dtype='f4') # initialize temperature grid 
        T[nGhosts:Nx+nGhosts, nGhosts:Ny+nGhosts, nGhosts:Nz+nGhosts] = outside_temperature + 1.0 # Set the plastic shell to 1.0
        T[nGhosts+1:Nx+nGhosts-1, nGhosts+1:Ny+nGhosts-1, nGhosts+1:Nz+nGhosts-1] = outside_temperature + 2.0 # Set internal 
    else: 
        T[:nGhosts, :, :] = outside_temperature
        T[-nGhosts:, :, :] = outside_temperature
        T[:, :nGhosts, :] = outside_temperature
        T[:, -nGhosts:, :] = outside_temperature
        T[:, :, :nGhosts] = outside_temperature
        T[:, :, -nGhosts:] = outside_temperature
    return T

def initialize_sources(sources: np.ndarray, computer_power: float, lamp_power: float, 
                       lamp_len_x: int, lamp_len_y: int, lamp_len_z: 
                       int, lamp_x: int, lamp_y: int, lamp_z: int, 
                       Nx: int, Ny: int, Nz: int, nGhosts: int) -> np.ndarray: 
    lamp_power_density = lamp_power / ((lamp_len_x*lamp_len_y*lamp_len_z)*ds**3)
    computer_power_density_per_grid = float(computer_power / ((Ny/2)*ds**3))

    sources[nGhosts:Nx//2, :nGhosts+1,  nGhosts:] = computer_power_density_per_grid

    x_range = np.arange(Nx + 2*nGhosts)
    y_range = np.arange(Ny + 2*nGhosts)
    z_range = np.arange(Nz + 2*nGhosts)

    # Using broadcasting to create a 3D boolean mask
    lamp_mask = ((x_range[:, None, None] >= lamp_x) & (x_range[:, None, None] < lamp_x + lamp_len_x) &
                (y_range[None, :, None] >= lamp_y) & (y_range[None, :, None] < lamp_y + lamp_len_y) &
                (z_range[None, None, :] >= lamp_z) & (z_range[None, None, :] < lamp_z + lamp_len_z))

    # Apply the lamp power density to the lamp area
    sources[lamp_mask] += lamp_power_density
    return sources

def set_all_arrays_as_fortran(T, sources, rho, p, vx, vy, vz) -> list: 
    vx = np.asfortranarray(vx, dtype='f4')
    vy = np.asfortranarray(vy, dtype='f4')
    vz = np.asfortranarray(vz, dtype='f4')
    rho = np.asfortranarray(rho, dtype='f4')
    p = np.asfortranarray(p, dtype='f4')
    T = np.asfortranarray(T, dtype='f4')
    sources = np.asfortranarray(sources, dtype='f4')
    return T, sources, rho, p, vx, vy, vz

from pyevtk.hl import gridToVTK
def save_npy_and_vtk(data, name, i, results_dir='./results/'): 
    fname = results_dir + f'{name}_{i}'
    np.save(fname + '.npy', data)
    x = np.arange(data.shape[0]+1)
    y = np.arange(data.shape[1]+1)
    z = np.arange(data.shape[2]+1)

    gridToVTK(fname, x, y, z, cellData={'data':data.copy()})

results_dir = './results_solverho/'
Nx, Ny, Nz, nGhosts = int(20), int(36), int(40), int(2)
wall_thickness = float(0.002) # meters

computer_power = float(400.0) # Watts
outside_temperature = float(-5.0)

lamp_power = float(500.0) # Watts
lamp_len_x, lamp_len_y, lamp_len_z = 8, 12, 5 # grid points 
lamp_x, lamp_y, lamp_z = 3, 12, 30 # grid points

dt = float(0.1)
ds = float(0.05)
t_max = float(10*60*60.0)
c_to_k = float(273.15)
rs = 287.0

T = np.zeros((Nx+2*nGhosts, Ny+2*nGhosts, Nz+2*nGhosts), dtype='f4', order='F') # initialize temperature grid  
sources = np.zeros((Nx+2*nGhosts, Ny+2*nGhosts, Nz+2*nGhosts), dtype='f4', order='F')
T = initialize_temperature(T, outside_temperature, Nx, Ny, Nz, nGhosts)
sources = initialize_sources(sources, computer_power, lamp_power, lamp_len_x, lamp_len_y, lamp_len_z, lamp_x, lamp_y, lamp_z, Nx, Ny, Nz, nGhosts)

rho = 1.225 # kg/m^3
rho = np.full((Nx+2*nGhosts, Ny+2*nGhosts, Nz+2*nGhosts), rho, dtype='f4', order='F')
p = (T+c_to_k)*rs*rho # Pa

time = 0.0
time_max = 1.0

vx = np.zeros((Nx+2*nGhosts, Ny+2*nGhosts, Nz+2*nGhosts), dtype='f4', order='F')
vy = np.zeros((Nx+2*nGhosts, Ny+2*nGhosts, Nz+2*nGhosts), dtype='f4', order='F')
vz = np.zeros((Nx+2*nGhosts, Ny+2*nGhosts, Nz+2*nGhosts), dtype='f4', order='F')


for i in range(10): 
    print(f'Iteration {i}')
    all_passed = True
    for arr, name in [(rho, 'rho'), (vx, 'vx'), (vy, 'vy'), (vz, 'vz'), (p, 'p'), (T, 'T')]:
        fname = results_dir + f'{name}_{i}.npy'
        if os.path.exists(fname):
            arr = np.load(fname)
        else: 
            all_passed = False 
    if all_passed: 
        print('All arrays loaded, skipping iteration')
        continue
    T, sources, rho, p, vx, vy, vz = set_all_arrays_as_fortran(T, sources, rho, p, vx, vy, vz)  
    
    print('Condctive step')
    tb.run_conductive(T, Nx, Ny, Nz, nGhosts, dt, ds, sources, wall_thickness, 60*60*3.0)
    print('Convective step')
    # p = (T+c_to_k)*rs*rho # Pa
    rho = p / ((T+c_to_k)*rs)

    tb.euler_driver(rho,vx,vy,vz,p,Nx,Ny,Nz,nGhosts,ds,time,time_max)
    for arr, name in [(rho, 'rho'), (vx, 'vx'), (vy, 'vy'), (vz, 'vz'), (p, 'p'), (T, 'T')]:
        save_npy_and_vtk(arr, name, i, results_dir=results_dir)

    T = p/(rs*rho) - c_to_k
    # re initialize temperature
    # T = initialize_temperature(T, outside_temperature, Nx, Ny, Nz, nGhosts, from_save=True)
    # p = (T+c_to_k)*rs*rho # Pa
    T, sources, rho, p, vx, vy, vz = set_all_arrays_as_fortran(T, sources, rho, p, vx, vy, vz)  
    
plot_grid(T)   
    
    
