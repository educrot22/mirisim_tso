#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Rene Gastaud 11 March 2021 for the data challenge
# compute the duration of the transit
#  Source:
# final thesis of Giuseppe Morello, 2016, formula 3.4
# https://discovery.ucl.ac.uk/id/eprint/1531000/1/final_thesis_Giuseppe_Morello.pdf
#

from astropy import units as u
#from astropy import constants as c
import numpy as np

def compute_eclipse_phase(rp, a, inc, verbose=False):
    """
    Compute the phase duration of the transit
    :param rp: the ratio of planet radius by star radius
    :type rp: float
    
    :param a: the ratio of semi grand axis of the planet  by star radius
    :type a: float
    
    :param inc: the inclination of the planet orbit, degree
    :type inc: float
    

    :return: duration, same unit than period
    :rtype: float
    """
    A = (1+rp) ** 2
    B = (a * np.cos(np.radians(inc))) ** 2
    C =  a * np.sin(np.radians(inc))
    #transit_time_14 = period/np.pi * np.arcsin(np.sqrt(A-B)/C)
    transit_phase_14 = np.arcsin(np.sqrt(A-B)/C).to(u.radian)
    transit_phase_14 = transit_phase_14.value/np.pi
    if (verbose):
        print('transit_phase_14 ', transit_phase_14 )
        print('ecipse', 0.5-transit_phase_14/2, 0.5+transit_phase_14/2, transit_phase_14/2)
    return transit_phase_14

def read_configuration_file(config_file, verbose=False):
    import configparser
    config = configparser.ConfigParser()
    config.read(config_file)
    planet_radius = float(config['planet']['radius'])*u.Rjup
    star_radius = float(config['star']['radius'])*u.Rsun
    rp = (planet_radius/star_radius).decompose().value
    semi_major_axis  =  float(config['orbit']['semi_major_axis'])*u.AU
    a = (semi_major_axis/star_radius).decompose()
    inclination = float(config['orbit']['inclination'])*u.degree
    period = float(config['orbit']['orbital_period'])*u.day
    if (verbose):
        print("rp={}, a={}, inclination={}, period={}".format(rp, a, inclination, period))
    return rp, a, inclination, period

def compute_eclipse_time( config_file, verbose=False):
    """
    Compute the  duration of the transit
    :param config_file: the configration file name
    :type rp: string

    :return: duration, same unit than period
    :rtype: quantity
    """
    rp, a, inclination, period = read_configuration_file(config_file, verbose=verbose)
    transit_phase_14 = compute_eclipse_phase(rp, a, inclination, verbose=verbose)
    transit_time_14  = transit_phase_14*period
    return transit_time_14
    
############
if __name__ == '__main__':
    import sys
    if len(sys.argv) == 2:
      config_file = sys.argv[1]
      rp, a, inclination, period = read_configuration_file(config_file, verbose=True)
      transit_phase_14 = compute_eclipse_phase(rp, a, inclination, verbose=True)
      print('transit_time_14', transit_phase_14*period)
    else:
      print('syntax: compute_eclipse_phase configuration_file')


