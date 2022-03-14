#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Process the illumination models :
1) read the the raw illumination models
2) add background
3) add companion B
4) add drift
5) add Poisson noise
6) add read out noise
7) write result

input :
set of illumination model fits files

output :
cooked cube + wavelength in a fits file

Rene Gastaud CEA-SACLAY IRFU  March, 14th 2022
"""
import astropy.units as u
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

from mirisim_tso.misc.read_illum_models_dir import read_illum_models_dir
from mirisim_tso.misc.read_illum_models_dir import write_illum_model

###############
from utils_process_illum_models import response_drift
from utils_process_illum_models import add_poisson_noise
from plot_cube import plot_cube


####################### INITIALISATION PARAMETERS #######################
gain = 5.5*u.electron/u.adu
file_B = 'compagnon/illum_models/illum_model_seq1_MIRIMAGE_P750L.fits'
background = 110*u.electron/u.second
pipe_dir='/Volumes/KINGSTON/ERS_NGTS10_2022_bis/mirisim/'
tint = 47.7*u.second
ron = 5*u.adu*np.sqrt(300)/tint
pattern='illum_model_*.fits'
output_filename = 'result.fits'
PLOT = True

#######################  DO IT  #######################
###  read cube
cube, wavelength = read_illum_models_dir(pipe_dir, pattern=pattern)
cube = np.float32(cube)
if (PLOT):plot_cube(cube, 'original cube')

#####  add background
cube = cube + background

#####  add companion B
image_B = fits.getdata(file_B)
nz, ny, nx = cube.shape
cube = cube + image_B.reshape([1,ny,nx])*u.electron/u.second

#####  add drift
signal = cube/gain
print('signal ', signal.shape, signal.min(), signal.max())
signal_difference = response_drift( signal.value)*u.adu/u.second
cube = cube+signal_difference*gain
plot_cube(cube, 'drift')
if (PLOT):plot_cube(signal_difference, 'signal_difference')

## add Poisson Noise
noised_cube = add_poisson_noise(cube, tint)
if (PLOT):plot_cube(noised_cube, 'Poisson Noise')

## add Read Out Noise
ron_noise = np.random.normal(0,1, size=cube.shape)*ron
print('ron_noise', ron_noise.shape, ron_noise.min(), ron_noise.max())
noised2_cube = noised_cube + ron_noise*gain
if (PLOT):plot_cube(ron_noise, 'RON')

write_illum_model(output_filename, noised2_cube, wavelength, overwrite=False)
