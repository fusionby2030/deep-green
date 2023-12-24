import numpy as np 
import sys 
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import scienceplots 
plt.style.use(['science'])
print(sys.argv, len(sys.argv))

import f90nml 

with open('./run/conductive_input.nml', 'r') as file: 
    nml = f90nml.read(file)
Nx, Ny, Nz, nGhosts = nml['greenhouse']['Nx'], nml['greenhouse']['Ny'], nml['greenhouse']['Nz'], nml['greenhouse']['nGhosts']
ds, dt,simrealtime = nml['simulation']['ds'], nml['simulation']['dt'], nml['simulation']['simrealtime']

if len(sys.argv) == 2:
    grid = np.fromfile(sys.argv[1], dtype=np.float64)
    grid = grid.reshape((Nz+2*nGhosts, Ny+2*nGhosts, Nx+2*nGhosts))
    vmin = min(grid.min(), grid.min())
    vmax = max(grid.max(), grid.max())
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    fig, axs = plt.subplots(2, 2, figsize=(8, 8), dpi=200)
    axs = axs.ravel()
    for n, Z_index in enumerate([int(i*Nz//4) for i in range(4)]): 
        cax = axs[n].imshow(grid[Z_index+1, :, :], extent=[0, Nx*ds, Ny*ds, 0], norm=norm)
        axs[n].set_title(f'Z={Z_index*ds:.2}m')#  @ {simrealtime/(60.0*60.0):.3}h')
    fig.colorbar(cax, ax=axs.ravel().tolist(), location='right')
    plt.show()

elif len(sys.argv) == 3:
    grid0 = np.fromfile(sys.argv[1], dtype=np.float64)
    grid0 = grid0.reshape((128, 128, 128)).T

    gridt = np.fromfile(sys.argv[2], dtype=np.float64)
    gridt = gridt.reshape((128, 128, 128)).T

    vmin = min(grid0.min(), gridt.min())
    vmax = max(grid0.max(), gridt.max())
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    fig, axs= plt.subplots(1, 2, figsize=(10, 5))
    cax = axs[0].imshow(grid0[:, :, 64], norm=norm)
    cax2 = axs[1].imshow(gridt[:, :, 64], norm=norm)
    fig.colorbar(cax2, ax=axs.ravel().tolist(), location='right')
    plt.show()

elif len(sys.argv) == 4:
    """ plot a 1d slice of the grid for 6 different z coordinates"""

    grid0 = np.fromfile(sys.argv[1], dtype=np.float64)
    grid0 = grid0.reshape((128, 128, 128)).T

    gridt = np.fromfile(sys.argv[2], dtype=np.float64)
    gridt = gridt.reshape((128, 128, 128)).T

    vmin = min(grid0.min(), gridt.min())
    vmax = max(grid0.max(), gridt.max())
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    fig, axs= plt.subplots(2, 3, figsize=(10, 5))
    axs = axs.ravel()
    for n, i in enumerate([i*(128//5) for i in range(6)]):

        cax = axs[n].imshow(gridt[:, :, i], norm=norm)
        axs[n].set_title(f"z = {i}")

    fig.colorbar(cax, ax=axs.ravel().tolist(), location='right')
    plt.show() 
