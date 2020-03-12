"""
Post-treatment of MIRISim simulations for exoplanets time-serie observations
(Martin-Lagarde | CEA-Saclay | 2019)
"""
import os
from . import utils
from . import effects
import logging
import sys
import glob
from astropy.io import ascii

from . import version

LOG = logging.getLogger(__name__)

def single_simulation_post_treatment(simulation_folder, conf, t_0=0., mask=None):
    """
    Apply post treatment to a single simulation folder.

    If used directly by the user, t_0 should not be used unless the user knows what he's doing. t_0 is mainly usefull
    in a sequential dataset of multiple observations. One might want to use a non-zero value for t_0, meaning the
    isolated observation will have an extra time offset of t_0 for all systematic effects.

    Parameters
    ----------
    simulation_folder: str
        path (relative or absolute) to the MIRISim simulation folder (the one that contains det_images/illum_models folder)
    t_0: float [optional]
        Time in second since beginning of the observation (beginning of observation is considered to be t = 0s)
    config: str or dict
        name of the ConfigObj .ini file, or corresponding dictionnary

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

    det_images_filename = os.path.join(simulation_folder, "det_images", "det_image_seq1_MIRIMAGE_P750Lexp1.fits")
    illum_models_filename = os.path.join(simulation_folder, "illum_models", "illum_model_1_MIRIMAGE_P750L.fits")

    signal = utils.read_illum_model(illum_models_filename)
    original_ramp, header = utils.read_det_image(det_images_filename)
    LOG.debug("main() | Value check for the original ramp: min={} / max={}".format(original_ramp.min(), original_ramp.max()))

    frame_time = header["TFRAME"]
    metadatas = {'history': ["Post processing with MIRISim TSO v{}".format(version.__version__)]}

    new_ramp = original_ramp.copy()
    if config_dict["response_drift"]["active"]:
        metadatas['history'].append("MIRISim TSO: Add Response drift")
        ramp_difference = effects.response_drift(original_ramp, t_0, signal, frame_time)
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
            mask = utils.read_mask()  # TODO include path to CDP MASK in ini file
        metadatas['history'].append("MIRISim TSO: Add Poisson Noise")
        new_ramp = effects.poisson_noise(new_ramp, mask)

    LOG.debug("main() | Value check for the new ramp: min={} / max={}".format(new_ramp.min(), new_ramp.max()))


    # Write fits file
    utils.write_det_image_with_effects(det_images_filename, new_data=new_ramp, extra_metadata=metadatas, config=config_dict,
                                       overwrite=config_dict["simulations"]["overwrite"])

def sequential_lightcurve_post_treatment(conf):
    """

    :param str folder: Folder that contain all the simulation folders for that light curve
    :param conf: name of the ConfigObj .ini file, or corresponding dictionary
    :type conf: str or Dict

    :return:
    """

    utils.init_log(stdout_loglevel="DEBUG", file_loglevel="DEBUG")

    # Read configuration file only once, and pass the dictionary to subsequent calls
    if isinstance(conf, str):
        config_dict = utils.get_config(conf)
    elif isinstance(conf, dict):
        config_dict = conf
    else:
        LOG.error("conf parameter needs to be str or dict")
        sys.exit()

    folder = config_dict["simulations"]["dir"]

    mask = utils.read_mask()

    # List simulations
    simulations = glob.glob(os.path.join(folder, "*/"))
    simulations.sort()

    if len(simulations) == 0:
        raise ValueError("No simulations found in {}.".format(folder))

    simulation_index = {}  # key: simulation folder ; value: simulation index (int)
    for sim in simulations:
        idx = int(sim.split("_")[-1][:-1])
        simulation_index[sim] = idx

    # Create time array
    timedat = os.path.join(folder, "times.dat")
    data = ascii.read(timedat)

    t_start = data["t_start"].data
    simulation_start_time = {}
    for sim in simulations:
        simulation_start_time[sim] = t_start[simulation_index[sim]]

    # Run each simulation post treatment, one after the other
    nb_simulations = len(simulations)
    simu_i = 0
    for simulation in simulations:
        simu_i += 1
        LOG.info("Run simulation {}: {:.1f}%".format(simulation, (simu_i*100/nb_simulations)))
        single_simulation_post_treatment(simulation, simulation_start_time[simulation], config_dict, mask=mask)
