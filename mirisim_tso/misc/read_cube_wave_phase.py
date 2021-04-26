import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from astropy.io import fits
import astropy.units as u
import datetime

################################################################################################################
##########################  READ/WRITE ILLUMINATION MODELS ########################################################
def read_cube_wave_phase(filename, verbose=False):
    '''
    Stacks MIRSim results into a single matrix. If save==True, puts every illum_models simulated in a cube, so as to reduce
    the total weight.
    BEWARE dimension (nt, ny, nx) instead of (ny, nx, nt)
  
    :param file_name: Name of the fits file where results are saved
    :param save: If True: save the results in a fits file
    :type in_dir: str
    
    :type verbose: bool, default False
    
    :return: Matrix of light curves, wavelength
    :rtype: ndarray; dim: (nt, ny, nx); value in electron/s
    '''
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
    orbital_phase = hdulist['PHASE'].data
    #
    hdulist.close()
    return cube, wavelength, orbital_phase

def write_cube_wave_phase(filename, cube, wavelength, orbital_phase, overwrite=False):
    hdu1 = fits.PrimaryHDU(cube.value)
    hdu1.header["BUNIT"] = cube.unit.to_string()
    hdu2 = fits.ImageHDU(wavelength.value)
    hdu2.header["BUNIT"] = wavelength.unit.to_string()
    hdu2.header["EXTNAME"] = 'WAVELENGTH'
    #
    hdu3 = fits.ImageHDU(orbital_phase)
    hdu3.header["EXTNAME"] = 'PHASE'
    #
    hdul = fits.HDUList([hdu1, hdu2, hdu3])
    hdul.writeto(filename, overwrite=overwrite)
    return

