import numpy as np 
import sys ,os

#  skip first byte in np.fromfile 
#  https://stackoverflow.com/questions/7396849/reading-binary-file-with-numpy-fromfile  
for file in sys.argv[1::]:
    grid = np.fromfile(file, dtype=np.float64)
    grid = grid.reshape((128, 128, 128),order="F")
    import matplotlib.pyplot as plt
    fig = plt.figure()
    # plt.imshow(grid[:,:,63],cmap='jet')
    # plt.colorbar()
    print(grid[:,63,63])
    plt.plot(grid[:,63,63])
    # plt.savefig(os.path.basename(file)+'.png')
    plt.savefig('tmp.png')
