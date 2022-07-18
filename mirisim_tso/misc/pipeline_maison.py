# From raw data, ramps in DN, create spectra one dimension array, flux versus wavelength
# the raw data can be :
#         simulated data,  format segmented files as produced by merge_sim_files.py
#         real data, from uncal.fits, as produced by the official pipeline
#
# Author R Gastaud, CEA Saclay
# date : creation July, 5th 2022
#        July, 6th  2022  beware Jeroen convention for mask zero is good
#
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
import os
import sys
sys.path.append('/home/gastaud/exoplanets_simulation/outils_bis/')
#
from read_uncal_dir import read_uncal_dir
from extract_spectra_with_mask import extract_spectra_with_mask
from readwrite_generic import write_generic_spectra
from readwrite_generic import read_generic_spectra

#
cal_dir = '/home/gastaud/exoplanets_simulation/outils_bis/'
input_dir = '/feynman/work/projets/webb-commissioning/Uncal220608/'
in_pattern = 'jw01033005001_04103_00001-seg00*_mirimage_uncal.fits'
slopes_file='jw01033005001_04103_00001-rgimage.fits'
plot = True
out_pattern='spectres_LD168-9_maison_'
out_dir = './'
#
###########  STAGE 1 : SLOPES  ###########
# create the slopes from the raw data , here the uncal, can be the segmented files
cube = read_uncal_dir(input_dir, pattern=in_pattern,out_file=slopes_file, verbose=True)
##
###########  STAGE 2 : SPECTRA  ###########
cube = fits.getdata(slopes_file)
header = fits.getheader(slopes_file)
print(cube.shape)
mask = fits.getdata(os.path.join(cal_dir, 'mask.fits'))
wavelength = fits.getdata(os.path.join(cal_dir, 'wavelength_v3.fits.gz'))
tint = header['tint']
# for the output of the official pipeline
# tint = header['EFFINTTM']
## this is just for verification
if (plot):
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

### EXTRACT SPECTRA
flux2d_bkg = extract_spectra_with_mask(cube, mask)

ymax = np.argmax(flux2d_bkg[-1,:])
print(ymax, flux2d_bkg.shape)
nt, ny, nx = cube.shape
tt = np.arange(nt)*tint
if(plot):
    plt.figure()
    plt.plot(tt, flux2d_bkg[:, ymax])
    plt.title('light curve for wavelength'+str(wavelength[ymax]))
    plt.xlabel('time s')
    plt.ylabel('DN/s')
    plt.savefig('plot_light_curve_wave_min.png')

flux2d_bkg.min(), flux2d_bkg.max()
ferror = np.sqrt(flux2d_bkg*tint*5.5)/tint/5.5
## beware ! different convention for the axes for the python function write
flux = flux2d_bkg*(u.adu/u.second)
flux = flux.T
ferror = ferror*(u.adu/u.second)
ferror = ferror.T

# and different convention for mask : zero is good ! (unix convention)
mask2d = np.zeros(flux.shape)
for i in np.arange(nt):mask2d[:,i] = 1-mask[:,36]

write_generic_spectra(out_dir, out_pattern, wavelength*u.micron, flux, ferror, mask2d, tt, verbose=True)

## verification
waves, flux_out, ferror_out, mask_out, times_out = read_generic_spectra(out_dir, pattern=out_pattern+'*.fits', verbose=True)
