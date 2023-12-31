# IPython log file

import matplotlib.pyplot as plt
import numpy as np
import os
import time
from scipy.ndimage import gaussian_filter1d
from astropy import units as u
#
from cascade.utilities.readwrite_spectra import plot_spectrum_2d
from cascade.utilities.readwrite_spectra import write_spectra_dir
from cascade.utilities.readwrite_spectra import read_spectra_dir

################
input_dir='/Users/gastaud/test_dap/cascade/simulation_LRS/test_generic3/WASP-43b_delivery0/SPECTRA'
# '/Users/gastaud/test_dap/cascade/simulation_LRS/test_generic3/wasp43b_bckg/SPECTRA/'

output_dir = 'data_10mu_norm/'
output_pattern = 'wasp43b'
################

mx, ferror, wave1d, times = read_spectra_dir(input_dir, pattern='*.fits', verbose=True)
print(mx.shape, ferror.shape, wave1d.shape, times.shape)
flux2d = mx.data
print(type(flux2d))
nw, nt = flux2d.shape
i6 = np.argmin(np.abs(wave1d.value-6))
i7 = np.argmin(np.abs(wave1d.value-7))
i10 = np.argmin(np.abs(wave1d.value-10))
i11 = np.argmin(np.abs(wave1d.value-11))
times.min(), times.max()
ttn= times-times.min()
ttn = ttn.to(u.second)

plt.figure()
plt.title('spectrum')
plt.plot(wave1d, flux2d[:,0])
plt.xlabel('wavelength micron')
plt.ylabel('flux '+str(flux2d.unit))

iw, it = np.unravel_index(flux2d.argmax(), flux2d.shape)
plt.figure()
plt.title('light curve')
plt.plot(ttn, flux2d[0,:], label="iw={}, w={:.2f}".format(0,wave1d[0]))
plt.plot(ttn, flux2d[iw,:], label="iw={}, w={:.2f}".format(iw,wave1d[iw]))
plt.xlabel('time in second')
plt.ylabel('flux '+str(flux2d.unit))
plt.legend()

plt.figure()
ligne6 = gaussian_filter1d(flux2d[i6,:].value, 3)
ligne7 = gaussian_filter1d(flux2d[i7,:].value, 3)
plt.plot(ttn, ligne6/ligne6.max(), label='wavelength 6. micron')
plt.plot(ttn, ligne7/ligne7.max(), label='wavelength 7. micron')
#plt.plot(ttn, flux2d[i10,:]/flux2d[i10,:].max(), label='wavelength 10. micron')
#plt.plot(ttn, flux2d[i11,:]/flux2d[i11,:].max(), label='wavelength 11. micron')
plt.legend()
plt.xlabel('time in second')
plt.ylabel('smoothed normalised flux')
plt.title('wasp43-b smoothed normalised light curve')
plt.savefig('plot_normalised_lightcurve67_bck.png')

lignes = gaussian_filter1d(flux2d.value, 3)
depth = (lignes.max(1) - lignes.min(1))/lignes.max(1)

plt.figure()
plt.title('wasp43-b smoothed normalised light curve  corrected bck')
plt.plot(wave1d, depth)
plt.xlabel('wavelength '+str(wave1d.unit))
plt.ylabel('flux method min/max')
plt.savefig('plot_depth_minmax_bck.png')

plt.figure()
plt.plot(wave1d)

plt.figure()
for i in np.arange(0, 61, 5):
    plt.plot(ttn, flux2d[i,:], label=str(i))
    #plt.plot(ttn, gaussian_filter1d(flux2d[i,:], 2), label=str(i))
plt.legend()

print('wave1d[60]', wave1d[60])
# wave1d[60] 10.772284507751465 micron

plt.figure()
for i in np.arange(60, nw, 10):
    plt.plot(ttn, gaussian_filter1d(flux2d[i,:], 5), label=str(i))
plt.legend()

mxc = mx[i10:,:]
ferrorc = ferror[i10:,:]
wave1dc = wave1d[i10:]
print(mxc.shape, ferrorc.shape, wave1dc.shape, times.shape)
# (142, 150) (142, 150) (142,) (150,)

#mask = np.full( small_flux.shape, False)

print(output_dir, output_pattern, wave1dc.shape, mxc.data.shape, ferrorc.shape, mxc.mask.shape, times.shape)

write_spectra_dir(output_dir, output_pattern, wave1dc, mxc.data, ferrorc, mxc.mask, times, verbose=False, overwrite=False)

