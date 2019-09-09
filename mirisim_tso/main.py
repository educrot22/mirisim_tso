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
from astropy.io import fits

from . import version

LOG = logging.getLogger(__name__)

def single_simulation_post_treatment(simulation_folder, t_0, conf):
    """
    Apply post treatment to a single simulation folder

    Parameters
    ----------
    simulation_folder: str
        path (relative or absolute) to the MIRISim simulation folder
    t_0: float
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

    det_images_filename = os.path.join(simulation_folder, "det_images", "det_image_seq1_MIRIMAGE_P750Lexp1.fits")
    illum_models_filename = os.path.join(simulation_folder, "illum_models", "illum_model_1_MIRIMAGE_P750L.fits")

    signal = utils.read_illum_model(illum_models_filename)
    original_ramp, header = utils.read_det_image(det_images_filename)
    LOG.debug("main() | Value check for the original ramp: min={} / max={}".format(original_ramp.min(), original_ramp.max()))

    frame_time = header["TFRAME"]

    new_ramp = original_ramp.copy()
    if config_dict["response_drift"]["active"]:
        ramp_difference = effects.response_drift(original_ramp, t_0, signal, frame_time)
        new_ramp += ramp_difference

    if config_dict["idle_recovery"]["active"]:
        ramp_difference = effects.idle_recovery(original_ramp, t_0, signal, frame_time, config_dict)
        new_ramp += ramp_difference

    if config_dict["anneal_recovery"]["active"]:
        ramp_difference = effects.anneal_recovery(original_ramp, t_0, frame_time, config_dict)
        new_ramp += ramp_difference


    # Apply poisson noise after all the other effects are applied
    if config_dict["noise"]["active"]:
        new_ramp = effects.poisson_noise(new_ramp)

    LOG.debug("main() | Value check for the new ramp: min={} / max={}".format(new_ramp.min(), new_ramp.max()))
    metadatas = {'history':"Post processing with MIRISim TSO v{}".format(version.__version__)}

    # Write fits file
    utils.write_det_image_with_effects(det_images_filename, new_data=new_ramp, extra_metadata=metadatas,
                                       overwrite=config_dict["simulations"]["overwrite"])

def sequential_lightcurve_post_treatment(folder, conf):
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

    # List simulations
    simulations = glob.glob(os.path.join(folder, "*/"))
    simulations.sort()

    if len(simulations) == 0:
        raise ValueError("No simulations found in {}.".format(folder))

    simulation_index = {}  # key: simulation folder ; value: simulation index (int)
    for sim in simulations:
        idx = int(sim.split("_")[-1][:-1])
        simulation_index[sim] = idx

    # Create time arrays
    example_file = glob.glob(os.path.join(simulations[0], "det_images", "*.fits"))[0]
    simulation_start_time = {}
    with fits.open(example_file) as example_hdu:
        ref_meta = example_hdu[0].header

        file_exposure_time = ref_meta["DURATION"]  # in seconds

        for sim in simulations:
            simulation_start_time[sim] = file_exposure_time * simulation_index[sim]  # in seconds

    # Run each simulation post treatment, one after the other
    for simulation in simulations:
        single_simulation_post_treatment(simulation, simulation_start_time[simulation], config_dict)