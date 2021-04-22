#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This is an utilty file for easy testing :
#     from read_kiss_spectra
#
import numpy as np
import os
import glob
from astropy.io import fits
from astropy.io import ascii
#from astropy.table import Table, Column
import astropy.units as u
import sys
#import datetime
#import matplotlib.pyplot as plt


##########################  READ/WRITE/COMPARE GENERIC SPECTRA #########
def read_tso_spectra(in_dir, pattern='*.fits', verbose=False):
    '''
    Stacks cascade generic 1D spectra  into a single matrix (time, wavelength),
    concatenate bjd, check that the wavelength vector is the same for all spectra.
    
    :param in_dir: Directory where fits files  are stored.
    :param pattern: pattern , e.g.  '*.fits'
    :param verbose: If True: print what is doing
    :type in_dir: str
    :type pattern: str;
    :type verbose: bool
    
    :return: Matrix of light curves
    :rtype: ndarray; dim: (nw, nt); value in electron/s
    
    :return: wavelength , micron
    :rtype: ndarray; dim: (nw); value in micron

    '''
    unit_dn = u.def_unit('DN', u.adu)
    #unit_dns = u.def_unit('DN/s', u.adu/u.second)
    u.add_enabled_units(unit_dn)
    print(unit_dn)
    #####
    data_files = sorted(glob.glob(os.path.join(in_dir, pattern)))
    n_file = len(data_files)
    if (verbose):print(n_file)
    if (n_file == 0):
        print('files not found for', os.path.join(in_dir, pattern))
        sys.exit()
    # read the first to get the shape
    filename = data_files[0]
    hdulist = fits.open(filename)
    #
    waves0 =  hdulist['WAVELENGTH'].data
    nw = waves0.size
    unit_lambda = u.Unit(hdulist['WAVELENGTH'].header['BUNIT'])
    #
    flux0 = hdulist['INTENSITY'].data
    ny, nx = flux0.shape
    flux = np.zeros([n_file, ny, nx])
    orbital_phase = np.zeros(n_file)
    unit_flux = u.Unit(hdulist['INTENSITY'].header['BUNIT'])
    i = 0
    for filename in data_files:
        if (verbose): print('x',i, filename)
        hdulist = fits.open(filename)
        waves =  hdulist['WAVELENGTH'].data
        if not(np.isclose(waves, waves0, atol=1e-15).all()):
            diff = np.abs(waves-waves0)
            raise Exception('waves ne waves0: {}'.format(diff.max()))
        flux[i,:, :] = hdulist['INTENSITY'].data
        orbital_phase[i] = hdulist[0].header['phase']
        hdulist.close()
        i = i+1

    # now put unit
    flux = flux*unit_flux
    waves = waves*unit_lambda
    return flux, waves, orbital_phase
