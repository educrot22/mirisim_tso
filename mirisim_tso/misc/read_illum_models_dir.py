import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from astropy.io import fits
import astropy.units as u
import datetime

################################################################################################################
##########################  READ/WRITE ILLUMINATION MODELS ########################################################
def read_illum_models_dir(in_dir, out_filename=None, pattern="simulation*/det_images*/det_image*rate.fits", verbose=False):
    '''
    Stacks MIRSim results into a single matrix. If save==True, puts every illum_models simulated in a cube, so as to reduce
    the total weight.
    BEWARE dimension (nt, ny, nx) instead of (ny, nx, nt)
    
    :param in_dir: Directory where MIRISim results are stored, containing simulation*/ directories.
    :param cube_name: Name of the fits file where results are saved if save==True.
    :param save: If True: save the results in a fits file
    :type in_dir: str
    :type out_filename: str; fits file name, default None
    :type verbose: bool, default False
    
    :return: Matrix of light curves, wavelength
    :rtype: ndarray; dim: (nt, ny, nx); value in electron/s
    '''
    
    data_files = sorted(glob.glob(os.path.join(in_dir, pattern)))
    n_file = data_files.__len__()
    # rene gastaud 17 december 2019 for output of the official MIRI pipeline , LEVEL 2 (slope = rate = DN/s)
    if (n_file == 0):
        data_files = glob.glob(os.path.join(in_dir, "simulation*/det_images*/det_image*rate.fits"))
        n_file = data_files.__len__()
    #print("patch", n_file)

    filename = data_files[0]
    if (verbose):print("read", filename)
    hdulist = fits.open(filename)
    intensity = hdulist[1].data
    bunit1 = hdulist[1].header["BUNIT"]
    wavelength = hdulist[2].data
    bunit2 = hdulist[2].header["BUNIT"]
    hdulist.close()

    ny, nx = intensity.shape
    ##cube = np.zeros([ny, nx, n_file])
    cube = np.zeros([n_file, ny, nx])
    cube[0, :, :] = intensity
    
    i = 0
    for filename in data_files:
        if (verbose):print("read", filename)
        cube[i, :, :] = fits.getdata(filename, ext=1)
        i = i + 1
    #
    cube = cube*u.Unit(bunit1)
    wavelength = wavelength*u.Unit(bunit2)

    if (out_filename is not None):
        if (verbose):print("write", ou_filename)
        write_illum_model(out_filename, cube, wavelength)
    #
    return cube, wavelength

def write_illum_model(filename, cube, wavelength, overwrite=False):
    hdu1 = fits.PrimaryHDU(cube.value)
    hdu1.header["BUNIT"] = cube.unit.to_string()
    hdu2 = fits.ImageHDU(wavelength.value)
    hdu2.header["BUNIT"] = wavelength.unit.to_string()
    hdu2.header["EXTNAME"] = 'WAVELENGTH'
    hdul = fits.HDUList([hdu1, hdu2])
    hdul.writeto(filename, overwrite=overwrite)
    return

