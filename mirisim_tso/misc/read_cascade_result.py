# IPython log file

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u

def read_cascade_result(filename):
    hdul = fits.open(filename)
    wave = hdul[1].data['Wavelength']
    depth = hdul[1].data['Depth']
    err_depth = hdul[1].data['Error Depth']
    hdul.close()
    return wave, depth, err_depth
