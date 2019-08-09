"""
Post-treatment of MIRISim simulations for exoplanets time-serie observations
(Martin-Lagarde | CEA-Saclay | 2019)
"""
import os
from . import utils
from . import effects
import logging
import sys

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


    utils.init_log(stdout_loglevel="INFO", file_loglevel="DEBUG")

    det_images_filename = os.path.join(simulation_folder, "det_images", "det_image_seq1_MIRIMAGE_P750Lexp1.fits")
    illum_models_filename = os.path.join(simulation_folder, "illum_models", "illum_model_1_MIRIMAGE_P750L.fits")

    signal = utils.read_illum_model(illum_models_filename)
    original_ramp, header = utils.read_det_image(det_images_filename)

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

    metadatas = {'history':"Post processing with MIRISim TSO v{}".format(version.__version__)}

    # Write fits file
    utils.write_det_image_with_effects(det_images_filename, new_data=new_ramp, extra_metadata=metadatas,
                                       overwrite=config_dict["simulations"]["overwrite"])