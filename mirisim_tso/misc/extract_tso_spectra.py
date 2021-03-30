#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This is an utilty file for easy testing :
#
#
import numpy as np
import astropy.units as u

#import matplotlib.pyplot as plt


##########################  READ/WRITE/COMPARE GENERIC SPECTRA #########
def extract_tso_spectra(flux3d, wave2d, i0=148, i1=383, verbose=False):
    '''
    Extract spectra

    '''
    n_file, ny, nx = flux3d.shape
    wave1d = wave2d[i0:i1, nx//2]
    cubette = flux3d[:,i0:i1, : ]
    flux2d = cubette.sum(axis=2)
    return flux2d, wave1d
