# Launch Jeroen's script
# R Gastaud, 21 juin 2022
# email René Gastaud 20/06/2022, 20h26
#  From Jeroen stage1_pipeline_used.ipynb
#
## input calibration files from:
###  Preliminary_RTS_Reference_Files
#   https://stsci.app.box.com/folder/161010409068?s=1bpk38ict288ormth0jospc66fdcdnng
#
##  input data from:
##       https://app.box.com/file/964144868818
####   jw01033005001_04103_00001-seg001_mirimage_uncal.fits
##
"""
export CRDS_SERVER_URL=https://jwst-crds.stsci.edu
export CRDS_PATH=/data/gastaud/crds_cache
"""

from astropy.io import fits # for reading FITS file contents
import matplotlib.pyplot as plt    # to display images
import numpy as np
from matplotlib import colors
import glob
import sys
import time
import os
from pathlib import Path
#
from jwst import datamodels
from jwst.pipeline import Detector1Pipeline
from jwst import __version__
#
filename1='/Volumes/Rene/mast_1033/TSO/jw01033005001_04103_00001-seg001_mirimage_uncal.fits'
calibration_dir = Path('/Volumes/Rene/Preliminary_RTS_Reference_Files/')
output_dir = Path('/Volumes/Rene/mast_1033/mon_resultat/')
####
print('CRDS_PATH', os.environ['CRDS_PATH'])
print('CRDS_SERVER_URL', os.environ['CRDS_SERVER_URL']
##
##
print('jwst version', __version__)
# jwst version 1.5.2
##  read
dm = datamodels.open(filename1)
dm.info()
print(dm.data.shape)
print(type(dm.data))
##
cal_dict = {'saturation': str(calibration_dir /
'MIRIMAGE_SATURATION_09.00.00.fits'),
            'linearity' : str(calibration_dir /
'MIRIMAGE_AVERAGE_LINEARITY_09.00.06.fits'),
            'mask': str(calibration_dir /
'MIRI_MIRIMAGE_MASK_09.00.04.fits'),
            'reset': str(calibration_dir /
'MIRI_MIRIMAGE_SLITLESSPRISM_RESET_09.00.05.fits'),
            'dark': str(calibration_dir /
'MIRI_MIRIMAGE_DARK_FASTR1_SLITLESSPRISM_09.01.01.fits'),
            'rscd': str(calibration_dir /
'MIRI_MIRIMAGE_RSCD_09.03.00.fits')
           }

# MIRI_MIRIMAGE_MASK_09.00.04.fits better than MIRI_MIRIMAGE_MASK_Flight_09.00.05.fits

for mycal in cal_dict:
    print(mycal, cal_dict[mycal], os.path.exists(cal_dict[mycal]))


step1_dict = {
    # Detector 1 step:
        'group_scale': {'skip':False},
        #'dq_init': {'skip':False, 'override_mask': cal_dict['mask']},
        'dq_init': {'skip':False}, # Not using override mask at the moment
        'saturation': {'skip':True, 'override_saturation':
cal_dict['saturation']},
        'firstframe': {'skip':True},
        'lastframe': {'skip':False},
        'linearity': {'skip':False, 'override_linearity':
cal_dict['linearity']},
        'reset': {'skip':False, 'override_reset': cal_dict['reset']},
        'rscd': {'skip':True, 'override_rscd': cal_dict['rscd']},
        'dark_current': {'skip':False, 'override_dark': cal_dict['dark']},
        'refpix': {'skip':True},
        'jump': {'skip':False, 'maximum_cores':'none',
'rejection_threshold': 5.0, 'flag_4_neighbors':False},
        'ramp_fit': {'skip':False, 'maximum_cores':'none'},
        'ipc': {'skip': True},
        'gain_scale': {'skip':False}}

for key in step1_dict:
    print(key,step1_dict[key])
    print('')
    
    
detector1 = Detector1Pipeline()
"""
this does not work
rate = detector1.call(dm, output_dir = str(output_dir),
save_results = True, steps = step1_dict)
"""

file_name=filename1
rate = detector1.call(file_name, output_dir = str(output_dir),
save_results = True, steps = step1_dict,
   )

