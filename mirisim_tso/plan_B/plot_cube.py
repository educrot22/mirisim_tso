# PLOT CUBE
#
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
#################
def plot_cube(cube, my_title):
    nz, ny, nx = cube.shape
    print(cube.shape, cube.min(), cube.max())
    zm, ym, xm = np.unravel_index(np.argmax(cube), cube.shape)
    print("maximum ", zm, ym, xm, cube[zm, ym, xm] )
    tt = np.arange(nz)*47.7*u.second
    yy = np.array([1,1,1,0])
    xx = np.array([8,7,2,8])
    n = len(xx)
    plt.figure()
    for i in np.arange(n):
        x = xx[i] -8 + xm
        y = yy[i] -1 + ym
        etiquette = "x={},y={},flux={:4.2f}".format(x, y, cube[-1,y, x])
        plt.plot(cube[:,y,x], label=etiquette)
    plt.legend()
    plt.xlabel('time second')
    plt.ylabel('cube '+str(cube.unit))
    plt.title(my_title)
    plt.savefig('plot_'+my_title+'.png')
    return
