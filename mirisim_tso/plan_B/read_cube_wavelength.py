
#
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from astropy.io import fits
import astropy.units as u
import datetime

def read_cube_wavelength(filename, verbose=False):
    hdulist = fits.open(filename)
    if (verbose):
        print("read", filename)
        hdulist.info()
    #
    cube = hdulist[0].data
    bunit0 = hdulist[0].header["BUNIT"]
    cube = cube*u.Unit(bunit0)

    wavelength = hdulist['WAVELENGTH'].data
    bunit1 = hdulist['WAVELENGTH'].header["BUNIT"]
    wavelength = wavelength*u.Unit(bunit1)
    #
    hdulist.close()
    return cube, wavelength
