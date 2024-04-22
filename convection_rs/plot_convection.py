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
        "primitives/rho": [None, None], # [0.0,2.0],
        "primitives/vx" : [None,None],
        "primitives/vy" : [None,None],
        "primitives/vz" : [None,None],
        "primitives/pressure" :[84000,89000], 
        "conserved/mass" :[None,None],
        "conserved/energy" :[ None,None,]
        }        
def plotFile(input):
    var,file,cnt=input
    data=h5.File(file)[var][:]
    print(data.shape)
    nx,ny,nz=np.shape(data)
    data=data[2:-2,2:-2,2:-2]
    nx-=4
    ny-=4
    nz-=4
    data= np.flip(data,0)
    print(np.shape(data))
    if demo == 'RB': 
        # ny // 2
        # im = plt.imshow(np.flip(data[:, :, nz // 4].T),cmap='gnuplot',interpolation='bilinear', vmin=vmin, vmax=vmax)
        im = plt.contourf(data[:, :, nz - nz // 6].T, cmap='gnuplot', levels=100)
    elif demo == 'SHOCK':
        im=plt.imshow(np.flip(data[:,0,:].T),cmap='gnuplot',interpolation='bilinear', vmin=vmin, vmax=vmax) # ,vmin=0.0,vmax=2.0)
    else: 
        im=plt.imshow(np.flip(data[:,ny//2,:].T),cmap='gnuplot',interpolation='bilinear', vmin=vmin, vmax=vmax) # ,vmin=0.0,vmax=2.0)
        
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
    # plt.show()

def find_limits(input): 
    var, file, cnt = input 
    data=h5.File(file)[var][:]
    data = data[2:-2,2:-2,2:-2]
    vmin = np.min(data)
    vmax = np.max(data)
    return (vmin,vmax)

if (len(sys.argv)<3):
    print(f"Usage: python3 {sys.argv[0]} <var> <file sequence>")
    sys.exit(1)

var=sys.argv[1]
if ( not var in available):
    print(f"Invalid var {var} requested!")
    print("\tvar options : rho, vx, vy, vz, p, mass, momentum_x, momentum_y, momentum_x, energy, temperature")
    sys.exit(1)

demo=sys.argv[2]

# vmin, vmax = limits[var]
files=sys.argv[3::]
pool=Pool(8)
index=np.arange(0,len(files))

# find limits of the data
vmins, vmaxes = zip(*pool.map(find_limits,zip(np.repeat(var,len(index)),files,index)))
vmin = np.min(vmins)
vmax = np.max(vmaxes)
print(vmin, vmax)

pool.close() 
pool.join() 

pool=Pool(8)
print(var)
pool.map(plotFile,zip(np.repeat(var,len(index)),files,index))
