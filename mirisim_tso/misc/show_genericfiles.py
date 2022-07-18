# Read the generic spectra files

import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.io import ascii
from readwrite_generic import read_generic_spectra
#####
in_dir='rebin_generic/'
pattern='wasp80_transit*.fits'
period=3.06785234*u.day
my_title = "wasp 80 transit rebin"
######
waves, flux, ferror, mask, times = read_generic_spectra(in_dir, pattern=pattern)#, verbose=True)

waves = waves*u.micron
flux = flux*u.Jy

print(waves.shape, flux.shape, ferror.shape, mask.shape, times.shape)
# ((387,), (387, 450), (387, 450), (387, 450), (450,))

from show_spectra_wave_phase import show_spectra_wave_phase
phase = times/period
show_spectra_wave_phase(flux, waves, phase, my_title)

save_file = 'plot_wasp80b_transit_rebin'
show_spectra_wave_phase(flux, waves, phase, my_title, save_file=save_file)
