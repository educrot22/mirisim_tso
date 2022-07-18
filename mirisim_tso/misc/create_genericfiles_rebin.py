# Create generic files for CASCADE
# each file contains one calibrated spectrum, with Flux, Wavelength, Flux_error and Mask
#  each variable is a numpy array 1D
# in the FITS header you have the time

import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.io import ascii

# sys.path.append('/Users/gastaud/achren/simulation/outils/')
from readwrite_x1dints import read_x2dints
# from CML
from readwrite_generic import write_generic_spectra
##############
# To be customised
time_file = 'times.dat'
input_file = '../WASP_80b_transit_x2dints.fits'
out_dir='rebin_generic/'
pattern='wasp80_transit'
period=3.06785234*u.day
step = 10
###############################
out_flux, out_wavelength, out_error = read_x2dints(input_file)
print(out_flux.shape, out_wavelength.shape, out_error.shape)
nw, nt = out_flux.shape
#print(nw, n_file, out_flux.shape)

mask = np.isfinite(out_flux)
mask = np.logical_not(mask)
print('mask', mask.shape, mask.dtype, mask[200,1000])

table_time = ascii.read(time_file)
phase = table_time['phase']

tt = phase*period.value
print('time', tt.min(), tt.max(), tt.shape)

##########  rebin ##############
rebin_nt = int(np.floor(nt/step))


rebin_out_flux = np.zeros([nw, rebin_nt] )*out_flux.unit
rebin_out_error = np.zeros([nw, rebin_nt] )*out_error.unit
rebin_tt = np.zeros([rebin_nt] )
for i in np.arange(rebin_nt):
    rebin_out_flux[:, i] = out_flux[:,i*step:i*step+step].mean(axis=1)
    rebin_tt[i] = tt[i*step:i*step+step].mean()
    rebin_out_error[:, i] = out_error[:,i*step:i*step+step].mean(axis=1)/np.sqrt(step)
index = step*np.arange(rebin_nt)
rebin_mask = mask[:, index]

os.mkdir(out_dir)
write_generic_spectra(out_dir, pattern, out_wavelength, rebin_out_flux, rebin_out_error, rebin_mask, rebin_tt)
