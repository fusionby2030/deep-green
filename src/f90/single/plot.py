import numpy as np
import sys
import matplotlib.pyplot as plt
from multiprocessing import Pool
import h5py as h5

def plotfile(input):
    var, file, cnt = input
    data = h5.File(file)[var][:]
    name = var + "_" + str(cnt).zfill(7) + ".png"
    nx, ny, nz = data.shape
    print(data.shape)

    im = plt.imshow(data[:, ny//2, :].T, cmap="gnuplot", interpolation=None)
    plt.colorbar()
    plt.tight_layout()
    plt.gca().set_aspect("equal")
    plt.savefig(name, dpi=200)
    plt.close()


var = sys.argv[1]
files = sys.argv[2::]

# pool = Pool(4)
index = np.arange(0, len(files))
print(files)
for idx, file in zip(index, files):
    plotfile((var, file, idx))
# pool.map(plotfile, zip(np.repeat(var, len(index)), files, index))
