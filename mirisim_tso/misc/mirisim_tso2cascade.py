#!/usr/bin/env python3
# -*- coding: utf-8 -*-
## TRANSLATE OUTPUT OF MIRISIM-TSO INTO INPUT OF CASCADE
#  correction of orbital phase
#
# see /Users/gastaud/test_dap/cascade/simulation_LRS/drift/plot_eclipse_v0.py
# see kiss_images2spectra.py
#
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from astropy import units as u
import os
import time
#
from cascade.utilities.readwrite_spectra import plot_spectrum_2d
from cascade.utilities.readwrite_spectra import write_spectra_dir
#
from mirisim_tso.misc.mirim_rebin2d import mirim_rebin2d
####
from mirisim_tso.misc.read_tso_spectra import read_tso_spectra
from mirisim_tso.misc.extract_tso_spectra import extract_tso_spectra


def mirisim_tso2cascade(in_dir, in_pattern, object_ephemeris, object_period, zoom_t=1, verbose=False, out_dir=None, out_pattern=None, dphase=0, n_transit=1):

    gain = 5.5
    ####  Read all the spectral images
    flux3d, wave2d, orbital_phase = read_tso_spectra(in_dir, verbose=True)
    if (verbose):
        print(flux3d.shape, wave2d.shape, orbital_phase.shape )
        print('flux3d', flux3d.min(), flux3d.max(), flux3d.shape)
    
    #### compute the spectra from the spectral images, removing pixelS unused by the prism
    #  usefull spectra flux1.shape   (390, 151) wavelength, time

    flux2d, wave1d = extract_tso_spectra(flux3d, wave2d)
    flux2d = flux2d.T
    if (verbose):
        print('flux2d', flux2d.min(), flux2d.max(), flux2d.shape)
        print('wave1d', wave1d.min(), wave1d.max(), wave1d.shape)
    # flux2d 374.52610282798963 DN / s 7785.882042485816 DN / s (1500, 235)
    # wave1d 5.022634983062744 micron 12.008528709411621 micron (235,)

    mjd_times = (orbital_phase-dphase)*object_period + object_ephemeris
    
    if verbose:
        nw, nt = flux2d.shape
        iw, it = np.unravel_index(flux2d.argmax(), flux2d.shape)
        i6 = np.argmin(np.abs(wave1d.value-6))
        i7 = np.argmin(np.abs(wave1d.value-7))
        plt.figure()
        plt.title('spectrum')
        plt.plot(wave1d, flux2d[:,0])
        plt.xlabel('wavelength micron')
        plt.ylabel('flux '+str(flux2d.unit))

        plt.figure()
        plt.title('light curve')
        #plt.plot(orbital_phase, flux2d[0,:], label='0')
        plt.plot(orbital_phase, flux2d[iw,:], label=str(iw))
        plt.plot(orbital_phase, flux2d[iw+1,:], label=str(iw+1))
        plt.xlabel('phase')
        plt.ylabel('flux '+str(flux2d.unit))
        plt.legend()
        #
        plt.figure()
        plt.title('light curve 6 micron')
        plt.plot(orbital_phase, flux2d[i6,:])
        plt.xlabel('phase')
        plt.ylabel('flux '+str(flux2d.unit))
    ######################
    ## Temporal Average of zoom_t
    small_flux = mirim_rebin2d(flux2d, [1,zoom_t])
    small_flux = small_flux/zoom_t ### average
    small_flux = small_flux.value*u.adu/u.second
    if (verbose):print('small_flux ', small_flux.min(), small_flux.max(), small_flux.shape)
    small_nt = nt//zoom_t
    index_t = np.arange(small_nt)*zoom_t + zoom_t//2
    # time origin at the middle of the bin, not at the beginning
    #  pb with cascade ?

    small_mask = np.full( small_flux.shape, False)
    small_ferror = np.zeros(small_flux.shape)*small_flux.unit
    small_mjd_times = mjd_times[index_t]
    small_wavelength = wave1d # does not change

    dtimes = small_mjd_times[1:]-small_mjd_times[:-1]
    if (verbose):print('dtimes', dtimes.min(), dtimes.mean(), dtimes.max(), dtimes.shape)
    sampling_time = dtimes.mean().to(u.second)
    if verbose: print('sampling_time', sampling_time)
    # sampling_time 206.69999994538927 s
    signal_electron = small_flux*sampling_time.decompose()*n_transit
    signal_electron = signal_electron.value*gain
    small_ferror = np.sqrt(signal_electron)/gain/sampling_time.value/n_transit
    small_ferror =  small_ferror*small_flux.unit
   
    
    if verbose:
        plt.figure()
        plt.title('plot noisy flux '+str(n_transit))
        plt.plot(orbital_phase, flux2d[i6,:])
        plt.plot(orbital_phase[index_t], small_flux[i6,:]+small_ferror[i6,:])
        plt.plot(orbital_phase[index_t], small_flux[i6,:]-small_ferror[i6,:])
        #plt.plot(orbital_phase[index_t], flux2d[i6,index_t]-small_ferror[i6,:])

    ###  NORMALISE ####
    last_spectrum = flux2d[:,-1]
    factor = last_spectrum.value.reshape([nw, 1])
    small_flux = small_flux/factor
    small_ferror = small_ferror/factor
    
    write_spectra_dir(out_dir, out_pattern, small_wavelength, small_flux, small_ferror, small_mask, small_mjd_times, verbose=True, overwrite=True)
    #
    return

############
if __name__ == '__main__':
    #
    start_time = time.time()
    in_dir='/Volumes/Transcend/exoplanetA/Wasp43b_version0/data1/'
    in_pattern='drifted_illum_model_*.fits'
    #
    object_ephemeris = 2455726.54186*u.day
    object_period = 0.81347753*u.day
    out_pattern = 'Wasp43B_delivery0'
    out_dir = 'bidon_152'
    zoom_t=10
    verbose=True
    n_transit=1
    mirisim_tso2cascade(in_dir, in_pattern, object_ephemeris, object_period, zoom_t=zoom_t, verbose=verbose, out_dir=out_dir, out_pattern=out_pattern, n_transit=n_transit)
    elapsed_time = time.time() - start_time
    print('elapsed time:', elapsed_time)
#elapsed time: 6.282442092895508
