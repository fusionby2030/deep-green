import numpy as np 
import sys ,os
import matplotlib.pyplot as plt
from multiprocessing import Pool
import h5py as h5
plt.style.use('dark_background')


def plotFile(input):
    var,file,cnt=input
    data=h5.File(file)[var][:]
    data=data[2:-2,2:-2,2:-2]
    nx,ny,nz=np.shape(data)
    print(nx,ny,nz)
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.suptitle('Horizontally stacked subplots')
    ax1.imshow(data[:,:,nz//2],cmap="gist_heat")
    ax2.imshow(data[:,ny//2,:],cmap="gist_heat")
    ax1.set_title(f"{var}, Z = 0 slice")
    ax2.set_title(f"{var}, Y = 0 slice")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax2.set_xlabel("x")
    ax2.set_ylabel("z")
    plt.tight_layout()
    plt.savefig(var+"_"+str(cnt).zfill(7)+".png")
    plt.close()


if (len(sys.argv)<3):
    print(f"Usage: python3 {sys.argv[0]} <var> <file sequence>")
    print("\tvar options : rho, vx, vy, vz, p, mass, momentum_x, momentum_y, momentum_x, energy, temperature")
    sys.exit(1)
re=6378137.0
var=sys.argv[1]
files=sys.argv[2::]
pool=Pool(8)
index=np.arange(0,len(files))
pool.map(plotFile,zip(np.repeat(var,len(index)),files,index))
