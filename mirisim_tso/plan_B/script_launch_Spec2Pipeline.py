#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Launch Official Pipeline to create calibrated images and spectra

input :
set of rateints fits files

output :
set of calints and x1dints fits files

Rene Gastaud CEA-SACLAY IRFU  March, 14th 2022
"""
from astropy.io import fits
from astropy.io import ascii
import matplotlib.pyplot as plt
import time
import glob, os
#
from jwst import datamodels
from jwst.pipeline import  Spec2Pipeline
#
input_dir='/Volumes/KINGSTON/ERS_NGTS10_2022_bis/rateints'
pattern = 'ERS_NGTS10_2022_noisy_drift_seg_*_rateints.fits'
output_dir='/Volumes/KINGSTON/ERS_NGTS10_2022_bis/spectres'
bck_patch = 'jwst_miri_extract1d_slitless_withbgr_bis.json'
###
start_time = time.time()
in_files = sorted(glob.glob(os.path.join(input_dir, pattern)))
print(in_files)
steps = {'extract_1d':{'override_extract1d':bck_patch}}
#
for file in in_files:
    elapsed_time = time.time() - start_time
    print(os.path.basename(file), elapsed_time)
    dm = datamodels.open(file)
    dm.meta.wcsinfo.wcsaxes=2
    x1d2 = Spec2Pipeline.call(dm, save_results=True, output_dir=output_dir, steps=steps)
    
