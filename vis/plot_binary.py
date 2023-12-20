import numpy as np 
import sys 

#  skip first byte in np.fromfile 
#  https://stackoverflow.com/questions/7396849/reading-binary-file-with-numpy-fromfile  
grid = np.fromfile(sys.argv[1], dtype=np.float64)
grid = grid.reshape((128, 128, 128))
import matplotlib.pyplot as plt
fig = plt.figure()
plt.imshow(grid[:, :, 64])
plt.colorbar()
plt.show()
