#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Inspect the output of exonoodle
There is 2 inputs:
    the directory of the output sed files
    the name of the configuration file
    
Created on Frebuary the 3rd, 2021
@author: Rene Gastaud

history :
RG 10 february  2021 in check_input_file if a directory is specified, do not add source_dir

"""

import configparser
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import glob
import astropy.units as u
from astropy.io import ascii
#
from .readwrite_exo import read_exonoodle_sed_dir
from .readwrite_exo import read_absorption
from .readwrite_exo import read_exonoodle_sed_fits
#
##  BUGGGG
#from .readwrite import plot_spectrum_2d
import exonoodle  # just for the path
####
import batman
########
##  to be improved replace 'none' by None
def check_input_file(input_file, source_dir, verbose=True):
    if (input_file != 'none'):
        if (len(os.path.dirname(input_file)) == 0):
            full_input_file = os.path.join(source_dir, input_file)
        else:
            full_input_file = input_file
        if not os.path.exists(full_input_file):
            raise OSError('File does not exist: %s' % full_input_file)
    else:
        full_input_file = 'none'
    #
    if (verbose): print(full_input_file, os.path.exists(full_input_file))
    return full_input_file
#################
def compute_lightcurve(config_file, absorption=None, phase=None, limb_dark = "uniform", ldu=None, verbose=False):
    """
    Compute the light curve with Batman
    Beware : no limb darkening !

    :param str config_file: name of the configuration file for exonoodle
    :param float absorption : absorption
    :param nd.array(nt) phase : orbital phase array

    :return: light curve normalised to one
    :rtype: nd.array(nt)

    """
    config = configparser.ConfigParser()
    config.read(config_file)
    planet_radius = float(config['planet']['radius'])*u.Rjup
    star_radius = float(config['star']['radius'])*u.Rsun
    rp = (planet_radius/star_radius).decompose().value
    if(verbose):print('rp*rp', rp*rp, absorption)
    if (absorption is not None): rp = np.sqrt(absorption)
    semi_major_axis  =  float(config['orbit']['semi_major_axis'])*u.AU
    a = (semi_major_axis/star_radius).decompose()
    inclination = float(config['orbit']['inclination'])*u.degree
    period = float(config['orbit']['orbital_period'])*u.day
    time_sampling  = float(config['exonoodle']['time_sampling'])*u.second
    phase_min = float(config['exonoodle']['phase_min'])
    phase_max = float(config['exonoodle']['phase_max'])
    phase_step = (time_sampling/period).decompose()
    #
    params = batman.TransitParams()      #object to store transit parameters
    params.t0 = 0.                       #time of inferior conjunction
    params.per = 1.                      #orbital period
    params.rp = rp               #planet radius (in units of stellar radii)
    params.a = a.value                #semi-major axis (in units of stellar radii)
    params.inc = inclination.value    #orbital inclination (in degrees)
    params.ecc = 0.                       #eccentricity
    params.w = 90.                        #longitude of periastron (in degrees)
    if (limb_dark == "uniform"):
        params.limb_dark = "uniform"        #limb darkening model
        params.u = []      #limb darkening coefficients for uniform LD
    else:
        params.limb_dark = limb_dark
        params.u = ldu
    #
    orbital_phases = np.arange(start=phase_min, stop=phase_max, step=phase_step)
    if (verbose):print('orbital_phases', orbital_phases.min(), orbital_phases.max(), orbital_phases.shape)
    if (phase is not None):orbital_phases=phase
    m = batman.TransitModel(params, orbital_phases)    #initializes model
    flux = m.light_curve(params)
    return flux


def inspect_exonoodle_output(sed_dir, output_dir='plots', suffix='.pdf'):
    """
    Inspect the output of exonoodle = a set of sed (spectra) files and a file for the time and orbital phase

    :param str sed_dir : name of the directory for exonoodle outputs
    """
    
    # get just the name of the configuration for the title of the plots
    ini_file = glob.glob(os.path.join(sed_dir, "*.ini"))
    if (len(ini_file) != 1):
        print('there must be one and one only ini file', ini_file)
        return
    config_file = ini_file[0]
    my_base   = os.path.basename(config_file)
    my_config = os.path.splitext(my_base)[0]
    
    ######  first read ########################
    # read the ouput of exoonoodle
    flux2d, wavelength_exo, time, phase = read_exonoodle_sed_dir(sed_dir)
    plot_spectrum_2d(flux2d.value, wavelength_exo, time, my_config, kw=100, normalise=False)
    plt.savefig(os.path.join(output_dir, 'plot_'+my_config+'spectrum2d'+suffix))
    print("flux2d, wavelength, time, phase ", flux2d.shape, wavelength_exo.shape, time.shape, phase.shape, "shape")
    #((1335, 151), (1335,), (151,))
    nw, nt = flux2d.shape
    if(phase.size != nt):print('error on phase size')
    if(wavelength_exo.size != nw):print('error on wavelength_exo size')

    # read the configuration
    config = configparser.ConfigParser()
    configuration = config.read(config_file)
    root_dir, sub_dir = os.path.split(exonoodle.__path__[0])
    source_dir = os.path.join(root_dir, 'source')
    spectral_file   = check_input_file(config['star']['spectral_file'], source_dir,)
    try:
        LD_file = check_input_file(config['star']['LD_file'], source_dir )
        #read the limb darkening
        data_LD = ascii.read(LD_file)
        wavelength_LD = data_LD['wavelength']
        nc = len(data_LD.colnames)
    except:
        LD_file = None
        nc = 0
    absorption_file = check_input_file(config['planet']['absorption_file'], source_dir)

    # read the absorption
    wavelength_a, absorption = read_absorption(absorption_file)
    print('absorption shape', wavelength_a.shape, absorption.shape)
    # ((230,), (230,)

    #read the input spectrum of the star
    wavelength_star, flux_star = read_exonoodle_sed_fits(spectral_file)
    print('input star spectrum shape', wavelength_star.shape, flux_star.shape)
    # ((407002,), (407002,))

    
    ######  check the first spectrum  ########################
    ##  this is cut for time = 0, we check the wavelength dependancy
    exo_flux = flux2d[:, 0]
    flux_interpolated = np.interp(wavelength_exo.value, wavelength_star.value, flux_star.value)
    rdiff = (flux2d[:, 0].value-flux_interpolated)/flux2d[:, 0].value*100
    #
    plt.figure()
    plt.title(my_config+' exonoodle flux for the first sed')
    plt.plot(wavelength_star, flux_star, label='input star flux')
    plt.plot(wavelength_exo, exo_flux, label='output of exonoodle')
    plt.xlabel('wavelength micron')
    plt.ylabel('flux')
    #plt.xlim([wavelength_exo.value.min(), wavelength_exo.value.max() ])
    plt.show()
    plt.savefig(os.path.join(output_dir, 'plot_'+my_config+'star_flux'+suffix))
    #
    plt.figure()
    plt.title(my_config+' exonoodle relative error on flux for the first sed')
    plt.plot(wavelength_exo, rdiff)
    plt.xlabel('wavelength micron')
    plt.ylabel('relative flux error %')
    plt.show()
    plt.savefig(os.path.join(output_dir, 'plot_'+my_config+'relative_flux_error_exonoodle'+suffix))

    ######  check one lightcurve  ########################
    ##  this is cut for middle wavelength, we check the time dependancy
    light_curve = flux2d[nw//2, :]
    wavelength_lc = wavelength_exo[nw//2]
    absorption_lc = np.interp(wavelength_lc.value, wavelength_a.value, absorption)
    # limb darkening
    LD_u = []
    if (nc ==5):
        limb_dark="nonlinear"
    elif(nc ==3):
        limb_dark="quadratic"
    else:
        limb_dark="uniform"
    if (nc > 2):
        i_LD = np.argmin(np.abs(wavelength_LD - wavelength_lc))
        for i in np.arange(nc-1):LD_u.append(data_LD[data_LD.colnames[i+1]][i_LD])
    #
    print('limb_darkening', limb_dark, nc)
    #
    light_curve_b = compute_lightcurve(config_file, absorption=absorption_lc, phase=phase,limb_dark=limb_dark, ldu=LD_u, verbose=True)
    #light_curve_c = compute_lightcurve(config_file, absorption=absorption_lc, verbose=True)
    plt.figure()
    plt.title(my_config+' difference of lightcurve with Batman model')
    plt.plot(phase, light_curve, label='exonoodle')
    plt.plot(phase, light_curve_b*light_curve[-1], label='batman')
    #plt.plot(phase, light_curve_b*light_curve.max(), label='batman')
    plt.legend()
    plt.xlabel("orbital phase")
    plt.ylabel(" flux")
    plt.show()
    plt.savefig(os.path.join(output_dir,'plot_'+my_config+'lighcurve_exonoodle_batman'+suffix))
    
    #
    rdiff = (light_curve/light_curve.max()-light_curve_b)*1e6
    plt.figure()
    plt.title(my_config+' relative difference of lightcurve with Batman')
    plt.plot(phase, light_curve-light_curve_b*light_curve.max())
    plt.xlabel("orbital phase")
    plt.ylabel("relative difference of flux ppm")
    plt.show()
    plt.savefig(os.path.join(output_dir,'plot_'+my_config+'rdiff_lighcurve_exonoodle_batman'+suffix))
    
    ######  check the absorption = depth ############
    # there is no noise, no systematics,
    #  the depth for each wavelength can be computed with the minimum and the maximum
    flux_in  =  np.min(flux2d,1)
    flux_out =  np.max(flux2d,1)
    depth = (flux_out - flux_in)/flux_out*100
    #
    plt.figure()
    plt.plot(wavelength_a, absorption*100, label='absorption')
    plt.plot(wavelength_exo, depth, '.', label='exonoodle output minmax')
    plt.legend()
    plt.xlabel('wavelength micron')
    plt.ylabel('depth in %')
    plt.title(my_config+' absorption and exoNoodle depth (min/max method)')
    plt.xlim([3,11])
    plt.savefig(os.path.join(output_dir,'plot_'+my_config+'exonoodle_depth_and_absorption'+suffix))

############
if __name__ == '__main__':
  if len(sys.argv) == 2:
      inspect_exonoodle_output(sys.argv[1], suffix='.png')
  else:
      print('syntax: inspect_exonoodle_output.py sed_dir')

# inspect_exonoodle_output output/Noodles_2021-02-01 configurationWasp43bDark.ini
