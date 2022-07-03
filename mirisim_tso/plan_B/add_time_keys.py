#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import glob
import astropy.units as u
from astropy.io import fits
from astropy.io import ascii
from astropy.time import Time
import numpy as np
#
def add_time_keys(input_dir, pattern, epoch, period, time_filename, nint=1 ):
    """
    Compute the date, the time, mjd and add them to the first header of each file in input_dir+ pattern
    OverWrite the input files

    :param str input_dir: the input directory
    :param str pattern: the name of each file
    :param float period: the period of the exoplanet in day
    :param float epoch: the date of the transit in barycenter julian day
    :param time_file: String; .dat file
    :type time_file:  ascii file, format ecsv with columns : file_index, phase, 'time'
    :param int nint: the number of integrations (ramps) in each read file

    Rene Gastaud CEA-SACLAY IRFU  March, 18th 2022
    """
    #
    ###  read  index, phase, time
    #file_index, phase, time = read_exonoodle_time(times_file)
    table = ascii.read(times_file)
    file_index = table['file_index']
    phase = table['phase']
    time = table['time']*u.second
    #
    ###  read cube
    cube, wavelength = read_cube_wavelength(input_cube_file)
    cube_dn = np.float32(cube/gain)
    #
    in_files = sorted(glob.glob(os.path.join(input_dir, input_pattern )))
    nf = len(in_files)
    j = 0
    for i in np.arange(nf):
        filename = in_files[i]
        print(i, j, os.path.basename(filename))
        hdul = fits.open(filename)
        hdul[0].header['PHASE'] = (phase[nint*i], 'orbital phase, no unit')
        hdul[0].header['TIME_0'] = (time[nint*i].value, 'elapsed time since the start of the observation, second')
        t_bjd = epoch + period*phase[nint*i]
        hdul[0].header['BJD-OBS'] = (t_bjd, 'baycenter julian day')
        t_obj = Time(t_bjd, format='jd')
        time_obs = t_obj.to_value('iso', 'date_hms')[11:]
        hdul[0].header['TIME-OBS'] = (time_obs, 'barycenter utc')
        date_obs = t_obj.to_value('iso', 'date')
        hdul[0].header['DATE-OBS'] = (date_obs, 'barycenter utc')
        hdul[0].header['MJD-OBS'] = (t_obj.mjd,'baycenter modified julian day')
        #
        hdul.writeto(filename, overwrite=True)
        hdul.close()
        
