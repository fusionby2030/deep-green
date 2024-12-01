import h5py
import sys 
import matplotlib.pyplot as plt
from multiprocessing import Pool
import numpy as np

def read_cfd_data(filename, query: str=None) -> dict: 
    with h5py.File(filename, 'r') as f:
        data = {}
        for key in f.keys():
            data[key] = f[key][()]

    if query: 
        return data[query]
    return data


def find_limits(input): 
    var, file, cnt = input 
    data=h5py.File(file)[var][:]
    data = data[2:-2,2:-2,2:-2]
    vmin = np.min(data)
    vmax = np.max(data)
    return (vmin,vmax)

def plot_file(inp): 
    var, fname, cnt = inp 
    data = read_cfd_data(fname, var)
    print(data.shape, cnt)

    nx,ny,nz=np.shape(data)
    data=data# [2:-2,2:-2,2:-2]
    # nx-=4
    # ny-=4
    # nz-=4
    data= np.flip(data,0)
    
    im=plt.imshow(np.flip(data[:,ny//2,:].T),cmap='gnuplot',interpolation='bilinear') # , vmin=vmin, vmax=vmax) # ,vmin=0.0,vmax=2.0)

    plt.axis('off')
    plt.colorbar()
    # plt.colorbar(im,fraction=0.046, pad=0.04)
    plt.xlabel("x")
    plt.ylabel("z")
    #plt.title(f"{var} [${units[var]}$] ")
    plt.tight_layout()
    var_short=var.split("/")[-1]
    plt.title(str(cnt))
    plt.savefig(var_short+"_"+str(cnt).zfill(7)+".png",dpi=200)
    plt.close()
    
        
# fname = "/home/akadam/dev/deep-green/src/py/state_0000000.h5"


available = ["rho", "ux", "uy", "uz", "p", "mass", "etot", "momentum_x", "momentum_y", "momentum_z", "mass_flux_x", "mass_flux_z", "mass_flux_y"]


var=sys.argv[1]
if ( not var in available):
    print(f"Invalid var {var} requested!")
    print("\tvar options : rho, vx, vy, vz, p, mass, momentum_x, momentum_y, momentum_x, etot")
    sys.exit(1)


files=sys.argv[2::]
pool=Pool(8)
index=np.arange(0,len(files))

vmins, vmaxes = zip(*pool.map(find_limits,zip(np.repeat(var,len(index)),files,index)))
vmin = np.min(vmins)
vmax = np.max(vmaxes)
print(vmin, vmax)

pool.close() 
pool.join()


pool=Pool(8)
print(var)
pool.map(plot_file,zip(np.repeat(var,len(index)),files,index))

# data = read_cfd_data(fname)
# print(data.keys())