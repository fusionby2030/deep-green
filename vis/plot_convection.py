import numpy as np 
import sys ,os
import matplotlib.pyplot as plt
from multiprocessing import Pool
import h5py as h5
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.style.use('dark_background')


available=[ "rho", "vx", "vy", "vz", "p", "mass", "momentum_x", "momentum_y", "momentum_x", "energy", "temperature"]
units={
        "rho": "kg m^{-3}",
        "vx" : "m s^{-1}",
        "vy" : "m s^{-1}",
        "vz" : "m s^{-1}",
        "p" :  "kg m s^{-1}", 
        "mass" :"kg",
        "momentum_x" :"kg m s^{-1}",
        "momentum_y": "kg m s^{-1}",
        "momentum_x": "kg m s^{-1}",
        "energy" : "J m^{-3}",
        "temperature":"C"
        }

def plotFile(input):
    var,file,cnt=input
    data=h5.File(file)[var][:]
    data=data[2:-2,2:-2,2:-2]
    nx,ny,nz=np.shape(data)
    print(nx,ny,nz)
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.suptitle(f"Y and Z slices of {var}")
    im1=ax1.imshow(data[:,:,nz//2],cmap="gist_heat")
    im2=ax2.imshow(data[:,ny//2,:],cmap="gist_heat")
    # ax1.set_title(f"{var}, Z = 0 slice")
    # ax2.set_title(f"{var}, Y = 0 slice")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax2.set_xlabel("x")
    ax2.set_ylabel("z")
    divider1 = make_axes_locatable(ax1)
    divider2 = make_axes_locatable(ax2)
    cax1 = divider1.append_axes('right', size='5%', pad=0.05)
    cax2 = divider2.append_axes('right', size='5%', pad=0.05)
    cb1=fig.colorbar(im1, cax=cax1)
    cb2=fig.colorbar(im2, cax=cax2)
    cb1.ax.set_title(f"${units[var]}$ ")
    cb2.ax.set_title(f"${units[var]}$ ")
    plt.tight_layout()
    plt.savefig(var+"_"+str(cnt).zfill(7)+".png")
    plt.close()


if (len(sys.argv)<3):
    print(f"Usage: python3 {sys.argv[0]} <var> <file sequence>")
    sys.exit(1)
re=6378137.0
var=sys.argv[1]
if ( not var in available):
    print(f"Invalid var {var} requested!")
    print("\tvar options : rho, vx, vy, vz, p, mass, momentum_x, momentum_y, momentum_x, energy, temperature")
    sys.exit(1)

files=sys.argv[2::]
pool=Pool(8)
index=np.arange(0,len(files))
pool.map(plotFile,zip(np.repeat(var,len(index)),files,index))
