#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
FILL the files rateints from the output of mirisim lrs :
1) read the cube output of script_process_illum_models.py
2) add a subcube in the sci
3) compute the date, the time, mjd and add to the header
4) write result

input :
cooked cube + wavelength in a fits file
set of rateints fits files

output :
set of rateints fits files updated and copied

Rene Gastaud CEA-SACLAY IRFU  March, 14th 2022
"""

import os
import glob
import astropy.units as u
from astropy.io import fits
from astropy.io import ascii
from astropy.time import Time
import matplotlib.pyplot as plt
import numpy as np
#
from read_cube_wavelength import read_cube_wavelength
from plot_cube import plot_cube

####################### INITIALISATION PARAMETERS #######################
gain = 5.5*u.electron/u.adu
times_file = '/Volumes/KINGSTON/ERS_NGTS10_2022/sed/sed_cst_6/times.dat'
#times_file = '/Volumes/KINGSTON/ERS_NGTS10_2022/stage_2_with_drift/my_spectra/times.dat'
input_cube_file='result.fits'
output_dir='/Volumes/KINGSTON/ERS_NGTS10_2022_bis/rateints'
input_dir = '/Volumes/KINGSTON/ERS_NGTS10_2022/stage_1_no_drift/'
input_pattern = 'ERS_NGTS10_2022_nodrift_seg_*_rateints.fits'
output_pattern = 'ERS_NGTS10_2022_noisy_drift_seg_{:02d}_rateints.fits'
epoch_jd  = 2459808.7904368257 # julian day
epoch_mjd   = 59808.290436825715
Period = 0.7668944
#
#######################  DO IT  #######################
###  read  index, phase, time
#file_index, phase, time = read_exonoodle_time(times_file)
table = ascii.read(times_file)
file_index = table['file_index']
phase = table['phase']
time = table['time']*u.second
#
###  read cube
cube, wavelength = read_cube_wavelength(input_cube_file)
cube_dn = np.float32(cube/gain)
#
in_files = sorted(glob.glob(os.path.join(input_dir, input_pattern )))
nf = len(in_files)
j = 0
for i in np.arange(nf):
    filename = in_files[i]
    print(i, j, os.path.basename(filename))
    hdul = fits.open(filename)
    hdul[0].header['PHASE'] = phase[i]
    hdul[0].header['TIME_0'] = time[i].value
    t_mjd = epoch_mjd + Period*phase[i]
    hdul[0].header['MJD-OBS'] = t_mjd
    t_obj = Time(t_mjd, format='mjd')
    time_obs = t_obj.to_value('iso', 'date_hms')[11:]
    hdul[0].header['TIME-OBS'] = time_obs
    date_obs = t_obj.to_value('iso', 'date')
    hdul[0].header['DATE-OBS'] = date_obs
    #hdul.info()
    sci = hdul['sci'].data
    nz, ny, nx = sci.shape
    hdul['sci'].data = cube_dn[j:j+nz, :,:].value
    j = j+nz
    out_filename = output_pattern.format(i)
    hdul.writeto(os.path.join(output_dir, out_filename))
    hdul.close()
    
