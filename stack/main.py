import toolbox as tb 
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import scienceplots 
plt.style.use(['science'])

def plot_grid(grid):
    vmin = min(grid.min(), grid.min())
    vmax = max(grid.max(), grid.max())
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    fig, axs = plt.subplots(2, 2, figsize=(8, 8), dpi=200)
    axs = axs.ravel()
    for n, Z_index in enumerate([int(i*Nz//4) for i in range(4)]): 
        cax = axs[n].imshow(grid[:, :, Z_index+1], extent=[0, Nx*ds, Ny*ds, 0], norm=norm, aspect='equal')
        axs[n].set_title(f'Z={Z_index*ds:.2}m')#  @ {simrealtime/(60.0*60.0):.3}h')
    fig.colorbar(cax, ax=axs.ravel().tolist(), location='right')

    
    plt.show()

# python3 -m numpy.f2py  -c physics_toolbox.f90 -m toolbox --fcompiler=gnu95 --f90flags=-O3
""" TODO: READ FROM INPUT FILE """
Nx = int(20)
Ny = int(36)
Nz = int(40)
nGhosts = int(2)
wall_thickness = float(0.005) # meters
computer_power = float(10000.0) # Watts
lamp_power = float(200.0) # Watts
lamp_len_x, lamp_len_y, lamp_len_z = 3, 12, 30 # grid points 
lamp_x, lamp_y, lamp_z = 3, 12, 5 # grid points
outside_temperature = float(0.0)

dt = float(0.1)
ds = float(0.05)
t_max = float(10*60*60.0)

T = np.full((Nx+2*nGhosts, Ny+2*nGhosts, Nz+2*nGhosts), outside_temperature, dtype='f4') # initialize temperature grid 
T[nGhosts:Nx+nGhosts, nGhosts:Ny+nGhosts, nGhosts:Nz+nGhosts] = outside_temperature + 1.0 # Set the plastic shell to 1.0
T[nGhosts+1:Nx+nGhosts-1, nGhosts+1:Ny+nGhosts-1, nGhosts+1:Nz+nGhosts-1] = outside_temperature + 2.0 # Set internal 

lamp_power_density = lamp_power / ((lamp_len_x*lamp_len_y*lamp_len_z)*ds**3)
computer_power_density_per_grid = float(computer_power / ((Ny/2)*ds**3))

sources = np.zeros((Nx+2*nGhosts, Ny+2*nGhosts, Nz+2*nGhosts), dtype='f4')
sources[nGhosts:nGhosts+1, nGhosts:Ny//2, :] = computer_power_density_per_grid

T = np.asfortranarray(T, dtype='f4')

for k in range(Nz+2*nGhosts):
    for j in range(Ny+2*nGhosts):
        for i in range(Nx+2*nGhosts):
            if (i >= lamp_x & i < lamp_x+lamp_len_x & 
                j >= lamp_y & j < lamp_y+lamp_len_y & 
                k >= lamp_z & k < lamp_z+lamp_len_z):
                sources[i, j, k] = sources[i, j, k] + lamp_power_density    
# plot_grid(sources)
sources = np.asfortranarray(sources, dtype='f4')
tb.run_conductive(T, Nx, Ny, Nz, nGhosts, dt, ds, sources, wall_thickness, t_max)

plot_grid(T)

