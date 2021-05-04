# CHECK POISSON NOISE IN MIRISIM TSO

import numpy as np
import matplotlib.pyplot as plt
#
#from mirisim_tso import effects
###from mirisim_tso.effects import poisson_noise
from add_poisson_noise import add_poisson_noise as poisson_noise
###

cube = np.arange(5*4*3).reshape([1,3,4,5]) +1.
cube = cube*100
mask=np.repeat(False, 4*5).reshape([4,5])
print('cube', cube.min(), cube.max(), cube.shape)
# cube 100.0 6000.0 (1, 3, 4, 5)

diff = cube[0,1:,:,:] - cube[0,0:-1,:,:]
print('diff', diff.min(), diff.max(), diff.shape)
# diff 2000.0 2000.0 (2, 4, 5)


result = np.zeros([10000, 3, 4, 5])
for i in np.arange(10000):
    result[i,:,:,:] = poisson_noise(cube, mask)

# standard deviation = np.sqrt(2000*5.5)/5.5  = 19.069251784911845
ss = np.sqrt(diff.mean()*5.5)/5.5


result_m = result.mean(axis=0)
print('result_m', result_m.min(), result_m.max(), result_m.shape)
#print('statistics', result.mean(), result.std())

coeff = np.polyfit(cube.flatten(), result_m.flatten(),1)
print("linear fit", coeff)
#  array([1.00000425, 0.0397222 ])
x = np.array([0, cube.max()])
y = coeff[0]*x + coeff[1]

plt.figure()
plt.plot(cube.flatten(), result_m.flatten(), '.', label="points")
plt.plot(x, y, label="linear fit slope={0:5.3f}  intercept={1:5.3f}".format(coeff[0], coeff[1]) )
plt.xlabel('input cube')
plt.ylabel('output cube')
plt.legend()
plt.title('test of mirsim tso poisson_noise rene')
plt.savefig('plot_mirsim tso poisson_noise_rene.png')

result_s = result.std(axis=0)


