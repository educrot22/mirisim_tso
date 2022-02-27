"""
Post-treatment of MIRISim simulations for exoplanets time-serie observations
(Martin-Lagarde | CEA-Saclay | 2019)
update data['time'] RG 30 April 2021
update metadata 'TSOVISIT': True instead of 'T' RG 28 nov 2021
add the function add_poisson_noise, version '0.7.62', 23 February 2022, R Gastaud
add noise2 in the configuration dictionary, correct the  function response_drift for LRS SLITLESS
              version '0.7.63', 25 February 2022, R Gastaud
"""
import sys
import os
from . import utils
from . import effects
import logging
import sys
import glob
from astropy.io import ascii
import time

from . import version

##  RG DA 16 nov 2021 change read_mask

LOG = logging.getLogger(__name__)


def single_simulation_post_treatment(simulation_folder, t_0, phase,  conf, mask=None):
    """
    Apply post treatment to a single simulation folder.

    If used directly by the user, t_0 should not be used unless the user knows what he's doing. t_0 is mainly usefull
    in a sequential dataset of multiple observations. One might want to use a non-zero value for t_0, meaning the
    isolated observation will have an extra time offset of t_0 for all systematic effects.

    Parameters
    ----------
    simulation_folder: str
        path (relative or absolute) to the MIRISim simulation folder (the one that contains det_images/illum_models folder)
    t_0: float
        Time in second since beginning of the observation (beginning of observation is considered to be t = 0s)
        
    phase: float
            orbital phase, no unit, zero or integer values for mid-transit

    conf: str or dict
        name of the ConfigObj .ini file, or corresponding dictionnary
    mask:
        np.array(bool) - Array of bad pixels (True if bad, False if good)
                 Needed because the bad pixels have non-additive shapes where computation is not applicable. They need
                 to be excluded from the computation.

    Returns
    -------

    """

    if isinstance(conf, str):
        config_dict = utils.get_config(conf)
    elif isinstance(conf, dict):
        config_dict = conf
    else:
        LOG.error("conf parameter needs to be str or dict")
        sys.exit()

    LOG.debug("Post-Treatment for folder {}".format(simulation_folder))

    det_images_filename = glob.glob(os.path.join(simulation_folder, "det_images", "det_image_*.fits"))[0]
        #"det_image_seq1_MIRIMAGE_P750Lexp1.fits")
    illum_models_filename = glob.glob(os.path.join(simulation_folder, "illum_models", "illum_model_*.fits"))[0]
    #"illum_model_1_MIRIMAGE_P750L.fits")

    signal = utils.read_illum_model(illum_models_filename)
    original_ramp, header = utils.read_det_image(det_images_filename)
    LOG.debug("main() | Value check for the original ramp: min={} / max={}".format(original_ramp.min(), original_ramp.max()))

    frame_time = header["TFRAME"]
    
    metadatas = {'history': ["Post processing with MIRISim TSO v{}".format(version.__version__)], 'time_0' : t_0,
                 'phase': phase, 'TSOVISIT': True} # 'TSOVISIT': 'T' RG 28 nov 2021

    new_ramp = original_ramp.copy()
    if config_dict["response_drift"]["active"]:
        metadatas['history'].append("MIRISim TSO: Add Response drift")
        ramp_difference = effects.response_drift(original_ramp, t_0, signal, frame_time)
        new_ramp += ramp_difference
        
    if config_dict["response_drift_one"]["active"]:
        metadatas['history'].append("MIRISim TSO: Add Response drift")
        ramp_difference = effects.response_drift_one(original_ramp, t_0, signal, frame_time)
        new_ramp += ramp_difference

    if config_dict["idle_recovery"]["active"]:
        metadatas['history'].append("MIRISim TSO: Add Idle Recovery")
        ramp_difference = effects.idle_recovery(original_ramp, t_0, signal, frame_time, config_dict)
        new_ramp += ramp_difference

    if config_dict["anneal_recovery"]["active"]:
        metadatas['history'].append("MIRISim TSO: Add Anneal Recovery")
        ramp_difference = effects.anneal_recovery(original_ramp, t_0, frame_time, config_dict)
        new_ramp += ramp_difference


    # Apply poisson noise after all the other effects are applied
    if config_dict["noise"]["active"]:
        if mask is None:
            mask_file = config_dict["CDP"]["mask_file"]
            mode = config_dict["CDP"]["mode"]
            mask = utils.read_mask(mask_file, mode)  # done 16 nov 2021 RG & AD
        metadatas['history'].append("MIRISim TSO: Add Poisson Noise old way")
        new_ramp = effects.poisson_noise(new_ramp, mask)
        
    # Apply poisson noise after all the other effects are applied
    #import pdb
    if 'noise_bis' in config_dict:
        if config_dict["noise_bis"]["active"]:
            if mask is None:
                mask_file = config_dict["CDP"]["mask_file"]
                mode = config_dict["CDP"]["mode"]
                mask = utils.read_mask(mask_file, mode)  # done 16 nov 2021 RG & AD
            metadatas['history'].append("MIRISim TSO: Add Poisson Noise bis")
            new_ramp = effects.add_poisson_noise(new_ramp, mask)
            
    #pdb.set_trace()
    LOG.debug("main() | Value check for the new ramp: min={} / max={}".format(new_ramp.min(), new_ramp.max()))

    # TODO Add the time-stamp in BJD to the file header.
    
    output_folder = config_dict["simulations"]["output_dir"]
    output_filename = os.path.join(output_folder, os.path.basename(simulation_folder))
    
    # Write fits file
    utils.write_det_image_with_effects(det_images_filename, output_filename,  new_data=new_ramp, extra_metadata=metadatas, config=config_dict,
                                       overwrite=config_dict["simulations"]["overwrite"])


def sequential_lightcurve_post_treatment(conf):
    """
    Add the post treatement, integration per integration
    :param conf: name of the ConfigObj .ini file, or corresponding dictionary
    :type conf: str or Dict

    :return:
    """

    # Read configuration file only once, and pass the dictionary to subsequent calls
    if isinstance(conf, str):
        config_dict = utils.get_config(conf)
    elif isinstance(conf, dict):
        config_dict = conf
    else:
        sys.stdout.write("Configuration file parameter needs to be str or dict")
        sys.exit()

    output_folder = config_dict["simulations"]["output_dir"]
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    utils.init_log(log=os.path.join(output_folder, 'mirisim_tso.log'), stdout_loglevel="INFO", file_loglevel="DEBUG")
    LOG.info('Mirisim_TSO launched')
    LOG.info('"mirisim_tso.log" file created')
    LOG.info("Output directory is : " + output_folder)

    start_time = time.time()

    input_folder = config_dict["simulations"]["input_dir"]
    filtername = config_dict["simulations"]["filtername"]
    mask_file = config_dict["CDP"]["mask_file"]
    mode = config_dict["CDP"]["mode"]
    mask = utils.read_mask(mask_file, mode)  # done 16 nov 2021 RG & AD

    # List simulations
    simulations = glob.glob(os.path.join(input_folder, filtername))
    simulations.sort()

    if len(simulations) == 0:
        raise ValueError("No simulations found in {}.".format(input_folder))

    simulation_index = {}  # key: simulation folder ; value: simulation index (int)
    for sim in simulations:
        idx = int(sim.split("_")[-1])
        simulation_index[sim] = idx

    # Create time array
    timedat = os.path.join(input_folder, "times.dat")
    data = ascii.read(timedat)
    t_start = data["time"].data
    # t_start = data["t_start"].data obsolete since version 1.0.0
    orbital_phase = data["phase"]
    index = data['file_index'].astype(int)
    
    simulation_start_time = {}
    simulation_orbital_phase = {}
    for sim in simulations:
        simulation_start_time[sim] = t_start[simulation_index[sim]]
        simulation_orbital_phase[sim] = orbital_phase[simulation_index[sim]]

    # Run each simulation post treatment, one after the other
    nb_simulations = len(simulations)
    if config_dict["simulations"]["nb_simulations"] is not None:
        nb_simulations = config_dict["simulations"]["nb_simulations"]
    LOG.info("sequential_lightcurve_post_treatment() | number of simulation selected={} / total={}".format(nb_simulations, len(simulations)))
    
    # Run each simulation post treatment, one after the other
    simu_i = 0
    LOG.info('Calculating...')
    for simulation in simulations:
        simu_i += 1
        LOG.debug(' ')
        LOG.debug("Run simulation {}: {:.1f}%".format(simulation, (simu_i*100/nb_simulations)))
        single_simulation_post_treatment(simulation, simulation_start_time[simulation], simulation_orbital_phase[simulation], config_dict, mask=mask)
        if (simu_i == nb_simulations):
            print('yoho')
            break

    LOG.info('Done !')
    elapsed_time = time.time() - start_time
    LOG.info(f'Job done in : {elapsed_time} s')
