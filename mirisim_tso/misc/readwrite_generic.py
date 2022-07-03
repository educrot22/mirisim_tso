# Copy from file file_processing.py
#  in CML Cascade with Miri Lrs

import os
import glob
#import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.io import fits

##########################  READ/WRITE/COMPARE GENERIC SPECTRA ########################################################
def read_generic_spectra(in_dir, pattern='*.fits', verbose=False):
    '''
        Stacks cascade generic 1D spectra  into a single matrix (time, wavelength),
        concatenate bjd, check that the wavelength vector is the same for all spectra.
        
        :param in_dir: Directory where fits files  are stored.
        :param pattern: pattern , e.g.  '*.fits'
        :param verbose: If True: print what is doing
        :type in_dir: str
        :type pattern: str;
        :type verbose: bool
        
        :return: wavelength , micron
        :rtype: ndarray; dim: (nw); value in electron/s
        :return: Matrix of light curves
        :rtype: ndarray; dim: (nw, nt); value in electron/s
        :return: Matrix of error
        :rtype: ndarray; dim: (nw nt); value in electron/s
        :return: Matrix of mask
        :rtype: ndarray; dim: (nw nt); value zero if goog
        '''
    data_files = glob.glob(os.path.join(in_dir, pattern))
    n_file = data_files.__len__()
    if (verbose):print(n_file)
    # read the first to get the shape
    filename = data_files[0]
    hdulist = fits.open(filename)
    waves0 =  hdulist[1].data['lambda']
    nw = waves0.size
    hdulist.close()
    #
    flux = np.zeros([nw, n_file])
    ferror = np.zeros([nw, n_file])
    mask = np.zeros([nw, n_file], dtype=int)
    times = np.zeros([n_file])
    i = 0
    for filename in data_files:
        if (verbose): print(i, filename)
        hdulist = fits.open(filename)
        times[i]=hdulist[0].header['TIME_BJD']
        waves = hdulist[1].data['lambda']
        if not(np.array_equal(waves, waves0)):
            print("waves ne waves0")
            return
        flux[:,i] = hdulist[1].data['FLUX']
        ferror[:,i] = hdulist[1].data['FERROR']
        mask[:,i] = hdulist[1].data['MASK']
        hdulist.close()
        i = i+1
    return waves, flux, ferror, mask, times


def write_generic_spectra(in_dir, pattern, wavelength, flux, ferror, mask, times, verbose=False):
    '''
    Write several cascade generic 1D spectra
    :param in_dir: str, name of the output directory
    :param pattern: str, pattern of the name of the output files
    :param wavelength: ndarray: dim: (nw)  quantity (with unit)
    :param flux: ndarray: dim: (nw, nt)  quantity (with unit)
    :param error: ndarray: dim: (nw, nt) quantity (with unit)
    :param mask: ndarray: dim: (nw, nt)
    :param time: ndarray: dim: (nt)   numpy array, unit day
    :param verbose: If True: print what is doing
    '''
    #
    nw, n_file = flux.shape
    my_names = ['LAMBDA', 'FLUX', 'FERROR', 'MASK']
    # should be UPPER CASE, and be carreful lambda is a python language keyword
    for i in np.arange(n_file):
        filename = os.path.join(in_dir, pattern+'{:03}.fits'.format(i))
        print(i, filename)
        #t = Table([wavelength, flux[:,i], ferror[:,i], mask[:,i]], names = my_names)
        col1 = fits.Column(name='FLUX'   , format='E', array=flux[:,i].value,  unit=flux.unit.to_string())
        col2 = fits.Column(name='FERROR' , format='E', array=ferror[:,i].value, unit=ferror.unit.to_string())
        col3 = fits.Column(name='MASK'   , format='J', array=mask[:,i])
        col4 = fits.Column(name='LAMBDA' , format='E', array=wavelength.value, unit=wavelength.unit.to_string())
        #t.meta['TIME_BJD'] = times[i]
        primary_hdu = fits.PrimaryHDU()
        primary_hdu.header['TIME_BJD'] = times[i]
        table_hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4])
        #t.write(filename, format='fits')
        hdul = fits.HDUList([primary_hdu, table_hdu])
        hdul.writeto(filename)
        i = i+1
    return

def compare_generic_spectra(wavelength1, flux1, ferror1, mask1, times1, wavelength2, flux2, ferror2, mask2, times2, verbose=False):
    '''
    Compare several cascade generic 1D spectra
    :param wavelength: ndarray: dim: (nw)  quantity (with unit)
    :param flux: ndarray: dim: (nw, nt)  quantity (with unit)
    :param error: ndarray: dim: (nw, nt) quantity (with unit)
    :param mask: ndarray: dim: (nw, nt)
    :param time: ndarray: dim: (nt)   numpy array, unit day
    :param verbose: If True: print what is doing
    '''
    result = True
    if not(np.allclose(wavelength1, wavelength2, rtol=5e-8)):
        print("error on wavelength")
        result = False
    
    if not(np.allclose(flux1, flux2, rtol=5e-8)):
        print("error on flux ")
        result = False
    
    
    if not(np.allclose(ferror1, ferror2, rtol=5e-8)):
        print("error on ferror")
        result = False
    
    
    if not(np.array_equal(mask1, mask2)):
        print("mask")
        result = False
    
    
    if not(np.allclose(times1, times2, rtol=5e-8)):
        print("error on times")
        result = False
    
    return result
