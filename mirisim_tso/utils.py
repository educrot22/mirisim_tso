import os
from astropy.io import fits
import numpy as np
import sys

import logging
import collections

import validate
import pkg_resources
import configobj

from . import constants as c

LOG = logging.getLogger(__name__)

def read_det_image(filename):
    """
    . Function to read the variable if it comes from a FITS file. Only takes |.fits containing detector image from
    MIRISim

    . Adapted from exoNoodle project (exonoodle.utils)

    >>>> PARAMETERS
    • illum_model_filename     : String - Name of the entry illum_model_filename containing the det_images to read

    >>>> RETURNS
    • original_ramp : np.array - values in DN of a ramp integration
    """

    try:
        hdulist = fits.open(filename)
        hypercube_image=hdulist[1].data
        header = hdulist[0].header
        hdulist.close()

        return hypercube_image, header
    except OSError:
        LOG.error("The illum_model file shall be a .fits file")
        raise


def read_illum_model(illum_model_filename):
    """
    . Function to read the variable if it comes from a FITS file. Only takes |.fits containing illumination model
    from MIRISim

    . Adapted from exoNoodle project (exonoodle.utils)

    >>>> PARAMETERS
    • illum_model_filename     : String - Name of the entry illum_model_filename containing the illum model to read

    >>>> RETURNS
    • slope_array   : np.array - values in DN/s of the illumination
    """

    (name, ext) = os.path.splitext(illum_model_filename)


    if ext==".fits":
        hdulist      = fits.open(illum_model_filename)
        illumination = hdulist['INTENSITY'].data
        hdulist.close()

        # Add reference pixel to illum models to match size of det_images
        (ny, nx) = illumination.shape
        slope_array = np.zeros((ny, nx+4))
        slope_array[:, 4:] = illumination/c.gain

        return slope_array

    else:
        print("ERROR :: The illum_model file shall be a .fits file !")
        raise ValueError


def extract_signal(illum_models_directory, x, y):
    """
    In order to calculate the effects evolution parameters, one need the theoretical slope. It is available in MIRISim illum_model.
    The file contains the pixel map of theoretical illumination in electrons/second. Therefore it needs conversion in DN/s.


    PARAMETERS:
    • illum_model_filename : [string] file name .fits

    RETURNS:
    • ideal_slope          : [float] flux in DN/s
    """
    illum_model=read_illum_model(illum_models_directory)

    ideal_slope = illum_model[x,y]

    return ideal_slope


def init_log(log="mirisim_tso.log", stdout_loglevel="INFO", file_loglevel="DEBUG", extra_config=None):
    """
    Init logging configuration file. Must be used before any other package you could think of and that might use logging as well.
    If it doesn't_0 work the way you expect, try inverting the imports to see if it changes.

    :param str log: filename where to store logs. By default "mirisim_tso.log"
    :param str stdout_loglevel: log level for standard output (ERROR, WARNING, INFO, DEBUG)
    :param str file_loglevel: log level for log file (ERROR, WARNING, INFO, DEBUG)
    :param dict extra_config: [optional] Set of extra properties to be added to the dict_config for logging
    :return:
    :rtype:
    """

    import logging.config

    log_config = {
        "version": 1,
        "formatters":
            {
                "form01":
                    {
                        "format": "%(asctime)s %(levelname)-8s %(message)s",
                        "datefmt": "%H:%M:%S"
                    },
                "form02":
                    {
                        "format": "%(asctime)s [%(processName)s/%(name)s] %(levelname)s - %(message)s",
                        "datefmt": "%H:%M:%S"
                    },
            },
        "handlers":
            {
                "console":
                    {
                        "class": "logging.StreamHandler",
                        "formatter": "form01",
                        "level": stdout_loglevel,
                        "stream": "ext://sys.stdout",
                    },
                "file":
                    {
                        "class": "logging.FileHandler",
                        "formatter": "form02",
                        "level": file_loglevel,
                        "filename": log,
                        "mode": "w",  # Overwrite file if it exists
                    },
            },
        "loggers":
            {
                "":
                    {
                        "level": "NOTSET",
                        "handlers": ["console", "file"],
                    },
            },
        "disable_existing_loggers": False,
    }

    if extra_config is not None:
        log_config = update_dict(log_config, extra_config)

    logging.config.dictConfig(log_config)


def update_dict(d, u):
    """
    Recursively merge or update dict-like objects.
    i.e, change a value to a key that already exists or
    add a (key, value) that did not previously existed

    source: https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth

    :param dict d: Original dictionnary
    :param dict u: dictionnary of updates to apply to 'd'
    :return dict d: Return updated version of 'd'
    """

    for k, v in u.items():
        if isinstance(v, collections.Mapping):
            d[k] = update_dict(d.get(k, {}), v)
        else:
            d[k] = v
    return d


def write_det_image_with_effects(original_path, new_data, extra_metadata, config, overwrite=True):
    """
    Based on the original .fits file, will overwrite the data cube with the one given
    in parameter.
    This .fits file will be stored in the same simulation folder, but in a new subfolder called
    "det_images_post_processed".


    Parameters
    ----------
    original_path: str
        Full path (relative or absolute) with filename of the input det_image
    new_data: np.ndarray
        new data cube (nb_integrations, nb_groups, nb_y, nb_x)
    extra_metadata: dict
        New metadata you want to add in the future .fits file (as opposed to the original one)
    overwrite: bool
        [optional] By defaut, True. Do you want to overwrite the post_processed .fits file if it already exists?

    Returns
    -------
    Write a new .fits file

    """

    # Construct name depending on effect activated
    final_dir = "det_images"

    if config["response_drift"]["active"]:
        final_dir += "_drift"

    if config["idle_recovery"]["active"]:
        final_dir += "_idle"

    if config["anneal_recovery"]["active"]:
        final_dir += "_anneal"

    if config["noise"]["active"]:
        final_dir += "_noise"

    # If all effects are deactivated, don't write any outputs
    if final_dir == "det_images":
        LOG.warning("All effects are deactivated, not writing any output")
        return

    original_name = os.path.basename(original_path)
    original_dir = os.path.dirname(original_path)

    hdulist = fits.open(original_path)

    # Replace existing data
    hdulist[1].data = new_data

    metadata = hdulist[0].header

    # Treat history separately
    if "history" in extra_metadata:
        history_list = extra_metadata["history"]

        for line in history_list:
            metadata["history"] = line

        del(extra_metadata["history"])

    # Update the rest of the parameters
    metadata.update(extra_metadata)

    new_path = os.path.abspath(os.path.join(original_dir, os.path.pardir, final_dir))

    if not os.path.isdir(new_path):
        os.mkdir(new_path)

    new_path = os.path.join(new_path, original_name)
    hdulist.writeto(new_path, overwrite=overwrite)


def get_nested(data, args):
    """
    Allow to get value in dictionnary tree

    Used in ConfigObj validator

    If ones want to get toto["section1"]["s2]["key"]
    call:
    value = get_nested(toto, ["section1", "s2", "key"])

    Parameter:
    :param dict data: input dict to get data on
    :param list(str) args: list of keys to use recursively

    :return: value corresponding to the list of keys
    """
    value = data.copy()
    for key in args:
        value = value.get(key)

    return value


def get_config(filename):
    """
    Read then validate and convert the input file
    BEWARE there must be a file named 'configspec.ini' in the directory of kiss

    :param str filename: .ini filename

    :return ConfigObj:  config file
    """

    if not os.path.isfile(filename):
        raise ValueError("The file '{}' can't be found".format(filename))

    # Prepare to convert values in the config file
    val = validate.Validator()
    specfile = pkg_resources.resource_filename('mirisim_tso', 'configspec.ini')
    configspec = configobj.ConfigObj(specfile, list_values=False)

    config = configobj.ConfigObj(filename, configspec=configspec, raise_errors=True)

    # Check and convert values in config.ini (i.e str to float or integer/bool)
    results = config.validate(val, preserve_errors=True)

    for entry in configobj.flatten_errors(config, results):

        [section_list, key, error] = entry
        section_list.append(key)

        if not error:
            msg = "The parameter %s was not in the config file\n" % key
            msg += "Please check to make sure this parameter is present and there are no mis-spellings."
            raise ValueError(msg)

        if key is not None and isinstance(error, validate.VdtValueError):
            option_string = get_nested(configspec, section_list)
            msg = "The parameter {} was set to {} which is not one of the allowed values\n".format(
                key, get_nested(config, section_list))
            msg += "Please set the value to be in {}".format(option_string)
            raise ValueError(msg)

    return config


def progressBar(value, endvalue, message, bar_length=20):
    percent = float(value) / endvalue
    arrow = '-' * int(round(percent * bar_length) - 1) + 'x'
    spaces = ' ' * (bar_length - len(arrow))

    sys.stdout.write("\r{0} [{1}] {2}%".format(message, arrow + spaces, int(round(percent * 100))))
    sys.stdout.flush()


def read_mask():
    """
    Read MIRI CDP MASK and return mask cropped for a SLITLESS observation

    :param str filename: path to CDP file

    :return: cropped mask (False if good, True if rejected)
    :rtype: ndarray(bool)
    """

    filename = pkg_resources.resource_filename("mirisim_tso", "data/MIRI_FM_MIRIMAGE_MASK_07.02.01.fits.gz")

    # Fake metadata dictionnary of a SLITLESS LRS det_images
    metadata = {}
    metadata["SUBSTRT1"] = 1  # Starting pixel in axis 1 direction
    metadata["SUBSTRT2"] = 529  # Starting pixel in axis 2 direction
    metadata["SUBSIZE1"] = 72  # Number of pixels in axis 1 direction
    metadata["SUBSIZE2"] = 416  # Number of pixels in axis 2 direction

    x_start = metadata["SUBSTRT1"] - 1
    y_start = metadata["SUBSTRT2"] - 1
    x_stop = x_start + metadata["SUBSIZE1"]
    y_stop = y_start + metadata["SUBSIZE2"]

    with fits.open(filename) as hdulist:
        mask = hdulist[1].data.astype(bool)

    return mask[y_start:y_stop, x_start:x_stop]
