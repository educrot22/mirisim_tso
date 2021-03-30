###  from Cascade with Miri Library
import math
import numpy as np
import os
import glob
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table, Column
import astropy.units as u
import datetime

__version = 'today'


def read_exonoodle_time(time_file, verbose=False):
    '''
    Read the file containing time or orbital phase of each spectrum.

    :param time_file: String; .dat file
    :type time_file:  ascii file, format ecsv with columns : file_index, phase, 'time'

    :return: file_index, phase, time
    :rtype: ndarray: dim: (nt)
    '''

    time_table = ascii.read(time_file)
    if(verbose):print(time_file)
    file_index = time_table['file_index'].data.astype(int)
    phase = time_table['phase'].data
    time = time_table['time'].data*time_table['time'].unit
    
    return file_index, phase, time
    
def read_exonoodle_sed_ecsv(sed_file, verbose=False):
    '''
    Read one spectrum,  output of exonoodle, format ecsv

    :param sed_file: String; .dat file
    :type sed_file: scii file, format ecsv with columns : wavelenvgth, flux

    :return: wavelength, flux
    :rtype: quantity ndarray: dim: (nw)
    '''
    if(verbose):print(sed_file)
    table = ascii.read(sed_file)
    #if(verbose):print(table)
    wavelength = table['wavelength'].data*table['wavelength'].unit
    flux  = table['flux'].data*table['flux'].unit

    return wavelength, flux

def read_exonoodle_sed_fits(sed_file, verbose=False):
    '''
    Read one spectrum,  output of exonoodle, format FITS

    :param sed_file: String; fits file
    :type sed_file: filename,  format FITS with columns : wavelenvgth, flux

    :return: wavelength, flux
    :rtype: quantity ndarray: dim: (nw)
    '''
    if(verbose):print(sed_file)
    hdulist = fits.open(sed_file)
    parameters = hdulist[1].data.names
    wavelength = hdulist[1].data['WAVELENGTH']
    wavelengthUnitString = hdulist[1].columns['WAVELENGTH'].unit
    wavelength = wavelength * u.Unit(wavelengthUnitString)
    #
    flux = hdulist[1].data['FLUX']
    fluxUnitString = hdulist[1].columns['FLUX'].unit
    flux = flux * u.Unit(fluxUnitString)
    hdulist.close()
    return wavelength, flux

def read_absorption(filename, verbose=False):
    '''
    Read absorption, format ecsv

    :param filename: String; .dat file
    :type filename: scii file, format ecsv with columns : wavelength, absorption

    :return: wavelength, absorption
    :rtype: quantity ndarray: dim: (nw)
    '''
    if(verbose):print(filename)
    table = ascii.read(filename)
    #if(verbose):print(table)
    wavelength = table['wavelength'].data*table['wavelength'].unit
    absorption  = table['data'].data
    #absorption  = table['absorption'].data
    return wavelength, absorption


#########################  EXONOODLE OUTPUT SED  ########################################################

def read_exonoodle_sed_dir(dir_path, verbose=False):  # dir_path : "dir/"  FIXED  no need of the final separator "/"
    '''
    Read all spectra simulated with ExoNoodle, and stored in the same directory.
    (see read_exonoodle_sed and write_exonoodle_sed)
    :param Directory path, contains all the SEDXXX.dat for every times. Time file is also in the directory and is
    called times.dat:
    :return: Return the matrix of light curves, wavelength and time arrays
    '''

    # read time file
    time_file = os.path.join(dir_path, "times.dat")
    file_index, phase, time = read_exonoodle_time(time_file, verbose=verbose)

    # read all the seds
    nd = math.ceil(math.log10(file_index.max()))
    if(verbose): print('number of decimals', nd)
    my_format = 'SED_{:0'+str(nd)+'d}.dat'
    sed_file_name =  os.path.join(dir_path, my_format.format(file_index[0]))
    wavelength0, flux0 = read_exonoodle_sed_ecsv(sed_file_name)
    nw = wavelength0.size
    nt = file_index.max() + 1 ### count from zero 
    flux2d = np.zeros([nw, nt])
    flux2d[:,:] = np.nan
    for i in file_index:
        sed_file_name =  os.path.join(dir_path, my_format.format(i))
        wavelength, flux = read_exonoodle_sed_ecsv(sed_file_name, verbose=verbose)
        flux2d[:,i] = flux.value
        if (wavelength0 != wavelength).any():
            print('bad')
            return
    flux2d = flux2d*flux.unit
    return flux2d, wavelength0, time, phase


def read_exonoodle_sed2d(filename, verbose=False):
    """
    It Reads the result of EXONOODLE grouped in one file
    (see read_exonoodle_sed_dir and write_exonoodle_sed)
    
    :param str filename: filename of the fits file.
    
    return:
    :param mat_flux : flux
    :rtype mat_flux: np.array(n_w, n_time)
    
    :param array_wl: wavelength
    :rtype array_wl: np.array(n_w)
    
    :param array_time: time
    :rtype array_time: np.array(n_t)
    
    """

    hdulist = fits.open(filename)
    flux   = hdulist['FLUX'].data
    bunit1 = hdulist['FLUX'].header["BUNIT"]
    flux   = flux*u.Unit(bunit1)
    
    waves  = hdulist['wavelength'].data
    bunit1 = hdulist['wavelength'].header["BUNIT"]
    waves  = waves*u.Unit(bunit1)
    
    times  = hdulist['time'].data
    bunit1 = hdulist['time'].header["BUNIT"]
    times  = times*u.Unit(bunit1)
    
    return flux, waves , times

def write_exonoodle_sed2d(filename, mat_flux, array_wl, array_time, history=None, overwrite=False):
    """
    It writes the result of EXONOODLE
    :param str filename:
    :param mat_flux : flux
    :type mat_flux: np.array(n_w, n_time)
    
    :param array_wl: wavelength
    :type array_wl: np.array(n_w)
    
    :param array_time: time
    :type array_time: np.array(n_t)
    
    :return:
    """
    primary_hdu = fits.PrimaryHDU()
    
    image_hdu1 = fits.ImageHDU(mat_flux)
    image_hdu1.header['EXTNAME'] = 'FLUX'
    image_hdu1.header['BUNIT'] = 'microJansky'
    
    image_hdu2 = fits.ImageHDU(array_wl)
    image_hdu2.header['EXTNAME'] = 'wavelength'
    image_hdu2.header['BUNIT'] = 'micron'
    
    image_hdu3 = fits.ImageHDU(array_time)
    image_hdu3.header['EXTNAME'] = 'time'
    image_hdu3.header['BUNIT'] = 'second'
    
    hdul = fits.HDUList([primary_hdu, image_hdu1, image_hdu2, image_hdu3])
    hdul.writeto(filename, overwrite=overwrite)
    return

def compare_exonoodle_sed2d(flux1, wavelength1, times1, flux2, wavelength2, times2, rtol=5e-8, verbose=False):
    '''
    Compare two 1D spectra (each spectrum contains : a flux 1D, wavelength 1D, time 1D)
    
    :param flux: ndarray: dim: (nw, nt)  quantity (with unit)
    :param wavelength: ndarray: dim: (nw)  quantity (with unit)
    :param time: ndarray: dim: (nt)   numpy array, unit day
    :param verbose: If True: print what is doing
    '''
    result = True
    if not(np.allclose(wavelength1, wavelength2, rtol=rtol)):
        print("error on wavelength")
        result = False
   
    if not(np.allclose(flux1, flux2, rtol=rtol)):
        print("error on flux ")
        result = False

    if not(np.allclose(times1, times2, rtol=rtol)):
        print("error on times")
        result = False
   
    return result

