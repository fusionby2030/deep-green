import numpy as np 
import sys 
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import scienceplots 
plt.style.use(['science'])



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

grid = np.load(sys.argv[1])
print(grid.shape)
grid = grid.T

vmin = min(grid.min(), grid.min())
vmax = max(grid.max(), grid.max())
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
fig, axs = plt.subplots(1, 4, figsize=(8, 8), dpi=200, sharey='row', sharex='row')
# axs = axs.ravel()
for n, Z_index in enumerate([int(i*Nz//4) for i in range(4)]): 
    cax = axs[n].imshow(grid[ Z_index+1, :,:], extent=[0, Nx*ds, Ny*ds, 0], norm=norm)
    axs[n].set_title(f'Z={Z_index*ds:.2}m')#  @ {simrealtime/(60.0*60.0):.3}h')
# axs[0, n].axvline(ds*(Ny//7))
    # axs[1, n].hist(grid[Z_index+1, :, :].flatten(), bins=100, density=True)
    # axs[1, n].axvline(np.mean(grid[Z_index+1, :, :]), color='k', linestyle='dashed', linewidth=1)
    # axs[1, n].axvline(np.median(grid[Z_index+1, Ny//4, Nz//4]), color='r', linestyle='dashed', linewidth=1)
fig.colorbar(cax, ax=axs.ravel().tolist(), location='right')
plt.show()