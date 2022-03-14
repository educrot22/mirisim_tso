# ADD MJD-OBS in the header
# R Gastaud   9 March 2022

from astropy.io import fits
import numpy as np
import os
import glob
#import matplotlib.pyplot as plt
from astropy.time import Time

#filename='/Volumes/KINGSTON/ERS_NGTS10_2022/stage_1_no_drift_poisson_bis/ERS_NGTS10_2022_nodrift_seg_000_rateints.fits'
input_dir = '/Volumes/KINGSTON/ERS_NGTS10_2022/stage_1_no_drift_poisson_bis/'
pattern = 'ERS_NGTS10_2022_nodrift_seg_*_rateints.fits'
files = sorted(glob.glob(os.path.join(input_dir, pattern)))
for filename in files:
    hdul = fits.open(filename)
    #hdul.info()
    date_obs = hdul[0].header['date-obs']
    time_obs = hdul[0].header['time-obs']
    tt = Time(date_obs+'T'+time_obs, format='isot')
    print(os.path.basename(filename), date_obs, time_obs, tt.mjd)
    hdul[0].header['mjd-obs'] =  tt.mjd, 'from date_obs + time_obs'
    hdul.writeto(filename, overwrite=True)
