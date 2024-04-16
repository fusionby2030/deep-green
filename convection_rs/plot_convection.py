import numpy as np 
import sys ,os
import matplotlib.pyplot as plt
from multiprocessing import Pool
import h5py as h5
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors
plt.style.use('dark_background')


available=[ "primitives/rho", "primitives/vx", "primitives/vy", "primitives/vz", "primitives/pressure", "conserved/mass","conserved/energy"]
units={
        "primitives/rho": "\\frac{kg}{m^{3}}",
        "primitives/vx" : "m s^{-1}",
        "primitives/vy" : "m s^{-1}",
        "primitives/vz" : "m s^{-1}",
        "primitives/pressure" :  "kg m s^{-1}", 
        "conserved/mass" :"kg",
        "conserved/energy" : "J m^{-3}",
        # "temperature":"C"
        }

limits={
        "primitives/rho": [0.7,1.6],
        "primitives/vx" : [None,None],
        "primitives/vy" : [None,None],
        "primitives/vz" : [None,None],
        "primitives/pressure" :[None,None], 
        "conserved/mass" :[None,None],
        "conserved/energy" :[ None,None,]
        }        
def plotFile(input):
    var,file,cnt=input
    data=h5.File(file)[var][:]
    print(data.shape)
    nz,ny,nx=np.shape(data)
    data=data[2:-2,2:-2,2:-2]
    nx-=4
    ny-=4
    nz-=4
    print(np.shape(data))
    data= np.flip(data,0)
    im=plt.imshow(np.flip(data[:,ny//2,:].T),cmap='bone_r',interpolation='bilinear',vmin=0.3,vmax=1.6)
    plt.axis('off')
    # plt.colorbar()
    # plt.colorbar(im,fraction=0.046, pad=0.04)
    plt.xlabel("x")
    plt.ylabel("z")
    #plt.title(f"{var} [${units[var]}$] ")
    plt.tight_layout()
    var_short=var.split("/")[-1]
    plt.savefig(var_short+"_"+str(cnt).zfill(7)+".png",dpi=200)
    plt.close()


if (len(sys.argv)<3):
    print(f"Usage: python3 {sys.argv[0]} <var> <file sequence>")
    sys.exit(1)

var=sys.argv[1]
if ( not var in available):
    print(f"Invalid var {var} requested!")
    print("\tvar options : rho, vx, vy, vz, p, mass, momentum_x, momentum_y, momentum_x, energy, temperature")
    sys.exit(1)

files=sys.argv[2::]
pool=Pool(8)
index=np.arange(0,len(files))
pool.map(plotFile,zip(np.repeat(var,len(index)),files,index))
