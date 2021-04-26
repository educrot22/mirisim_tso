#!/usr/bin/env python
#  -*- coding: utf-8 -*-

"""
Run custom response_drift2d to add drift to the illumination models.
The input of this step is the output of MIRISIM LRS
(not the output of kiss_lrs)

Beware mirisim:
(nb_integrations, nb_frames, nb_y, nb_x) = original_ramp.shape
(nb_frames, nb_y, nb_x)

see mirisim_tso.effects

"""

import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
import glob
import os
from astropy.io import fits
#
from mirisim_tso.misc.response_drift2d import response_drift2d
from mirisim_tso.misc.read_illum_models_dir import read_illum_models_dir

######### CUSTOMISATION
verbose = True
in_dir='data_output_mirisim/'
pattern = 'SED*.fits'
out_dir = 'data_output_mirisim_tso/'
#t0 = 7848.  ## *u.second   ## about 2 hours
t0 = 0  ## *u.second   ## about 2 hours
sampling_time = 10.335*10 ## *u.second
## time and orbital phase should be read from the file time.dat or the header of each file
######### end of customisation

#### read the data in a cube ##############
flux_elec, wavelength = read_illum_models_dir(in_dir, pattern=pattern)
print(flux_elec.max(), flux_elec.min(), flux_elec.shape, flux_elec.unit)

#### convert from elecon/second to DN/second
gain = 5.5*u.electron/u.adu
flux_dn = flux_elec/gain
print(flux_dn.unit, flux_dn.max(), flux_dn.min(), flux_dn.shape)
ii = np.where(flux_dn > 5000*u.adu/u.second)
print("percentage of high flux not changed", len(ii[0])/flux_dn.size*100)

plt.figure()
plt.title('first spectral image')
plt.imshow(flux_dn[0,:,:].value, origin='lower')
#plt.imshow(flux_dn[0,:,:].value, origin='lower')

###  compute the mean image
image = flux_dn.value.mean(axis=0)
nt, ny, nx = flux_dn.shape
t = t0 +np.arange(nt)*sampling_time
## time and orbital phase should be read from the file time.dat or the header of each file

### call response_drift2d
flux_difference = response_drift2d(image, t)

flux_dn2 = flux_dn + flux_difference*flux_dn.unit

data_files = sorted(glob.glob(os.path.join(in_dir, pattern)))
#data_files = data_files[0:2]
n_file = len(data_files)
new_index_int = np.arange(n_file)
for i in np.arange(n_file):
    hdul = fits.open(data_files[i])
    #hdul[0].header['phase'] = phase[i]
    #hdul[0].header['oldfile'] = data_files[i]
    hdul['intensity'].header['BUNIT']   =  str(flux_dn2.unit) # 'DN/s' #  old electron/s'
    hdul['intensity'].data = flux_dn2[i,:,:].value
    filename = "drifted_illum_model_"+"{:04d}".format(new_index_int[i])+".fits"
    filename = os.path.join(out_dir, filename)
    hdul.writeto(filename)
    if (verbose):print(filename)
    hdul.close()


