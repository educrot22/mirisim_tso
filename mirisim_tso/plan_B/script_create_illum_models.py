#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Create one illumination model per SED file.
All simulations will be in the folder (new_dir)
The inputs are :
1) the initialisation files (.ini)
2) the directory to read the sed files and the directory to put the illumination models.

output :
set of illumination model fits files

Rene Gastaud CEA-SACLAY IRFU  March, 14th 2022
"""

# Lanch MIRISIM LRS Illumination Model

import numpy as np
import time
import pickle
from astropy.io import fits
import os
import glob
#
from mirisim.config_parser import SimulatorConfig
from mirisim.config_parser import SimConfig
from mirisim.lrssim.lrs import create_lrs_from_scene
import pickle
from mirisim.config_parser import SceneConfig
from mirisim.skysim import scenemaker
#####################

###################   Initialistion Data #################################
sed_dir = '/Volumes/KINGSTON/ERS_NGTS10_2022/sed/sed_cst_all/'
#sed_dir = '/Volumes/KINGSTON/ERS_NGTS10_2022/sed/sed_cst_1/'
output_dir = '/Volumes/KINGSTON/ERS_NGTS10_2022_bis/mirisim'
pointing_file = 'init_files/pointing.save'  # to be replaced by make_dummy_pointing
scene_file = 'init_files/scene_1.ini'
simulator_config_file = 'init_files/simulator.ini'
"""
# to be replaced
simulator_config = SimulatorConfig.from_default()
simulator_config = SimulatorConfig.makeSimulator(
take_webbPsf=False,
include_refpix=True,
include_poisson=False,
include_readnoise=False,
include_badpix=True,
include_dark=False,
include_flat=True,
include_gain=True,
include_nonlinearity=True,
include_drifts=False,
include_latency=False,
cosmic_ray_mode='NONE')  # SOLAR_MIN, SOLAR_MAX, SOLAR_FLARE, NONE
"""
cfgpath =  'LRS_SLITLESS'
#filter_name='P750L'
#ReadDetect = 'SLITLESSPRISM'

###################   DO IT  #################################
sed_files = sorted(glob.glob(os.path.join(sed_dir, 'sed*.dat')))
ns = len(sed_files)
print('ns', ns)
start_time = time.time()
#
with open(pointing_file, 'rb') as pointing_file:
    pointing = pickle.load(pointing_file)
#
scene_config = SceneConfig(scene_file)
scene = scenemaker(scene_config)
#
simulator_config = SimulatorConfig(simulator_config_file)
#
for sed_file in sed_files:
    scene_config['point_1']['sed']['sedfile'] = sed_file
    scene = scenemaker(scene_config)
    illum_model = create_lrs_from_scene(scene, cfgpath, pointing=pointing, rel_obsdate=0, simulator_config=simulator_config, verbose=False)
    illum_model.meta['sedfile'] = os.path.basename(sed_file)
    index = os.path.basename(sed_file)[4:8]
    out_filename = os.path.join(output_dir, 'illum_model_'+index+'.fits')
    illum_model.write(out_filename)
elapsed_time = time.time() - start_time
print('elapsed time: {}'.format(elapsed_time))
