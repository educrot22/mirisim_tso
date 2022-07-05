
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
import os
import sys
sys.path.append('/home/gastaud/exoplanets_simulation/outils_bis/')
#
from read_uncal_dir import read_uncal_dir
from readwrite_generic import write_generic_spectra
from extract_spectra_with_mask import extract_spectra_with_mask
#
input_dir='/home/gastaud/exoplanets_simulation/outils_bis/'
filename='jw01033005001_04103_00001-rgimage.fits'
#########
# create the slopes from the raw data , here the uncal, can be the segmented files 
input_dir = '/feynman/work/projets/webb-commissioning/Uncal220608/'
pattern = 'jw01033005001_04103_00001-seg00*_mirimage_uncal.fits'
cube = read_uncal_dir(input_dir, pattern=pattern,out_file=filename verbose=True)
#####
###  read the inputs
cube = fits.getdata(filename)
header = fits.getheader(filename)
print(cube.shape)
mask = fits.getdata(os.path.join(input_dir, 'mask.fits'))
wavelength = fits.getdata(os.path.join(input_dir, 'wavelength_v3.fits.gz'))
tint = header['tint']
# for the output of the official pipeline
# tint = header['EFFINTTM']
print(tint)
ii = np.where(mask[:,36] == 1)
ii = ii[0]
plt.figure()
plt.plot(wavelength)
plt.plot(ii, wavelength[ii], 'o')
print(wavelength[ii].min(), wavelength[ii].max())
print('5.022635 12.008529')
print(ii.min(), ii.max())
print('148 382')

### EXTRACT SPAECT
flux2d_bkg = extract_spectra_with_mask(cube, mask)
flux2d_bkg.shape
ymax = np.argmax(flux2d_bkg[-1,:])
print(ymax)
nt, ny, nx = cube.shape
tt = np.arange(nt)*tint
plt.figure()
plt.plot(tt, flux2d_bkg[:, ymax])
plt.title('light curve for wavelength'+str(wavelength[ymax]))
plt.xlabel('time s')
plt.ylabel('DN/s')
plt.savefig('plot_light_curve_wave_min.png')
pattern='spectres_LD168-9_maison_'
flux2d_bkg.min(), flux2d_bkg.max()
ferror = np.sqrt(flux2d_bkg*tint*5.5)/tint/5.5
## beware ! different convention for the axes
flux = flux2d_bkg*(u.adu/u.second)
flux = flux.T
ferror = ferror*(u.adu/u.second)
ferror = ferror.T

mask2d = np.zeros(flux.shape)
for i in np.arange(nt):mask2d[:,i]=mask[:,36]

write_generic_spectra('./', pattern, wavelength*u.micron, flux, ferror, mask2d, tt, verbose=True)

## verification
from readwrite_generic import read_generic_spectra

toto = read_generic_spectra('./', pattern=pattern+'*.fits', verbose=True)
waves, flux_out, ferror_out, mask_out, times_out = toto
