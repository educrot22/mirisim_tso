import os
from astropy.io import fits

import logging
import collections

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
        hdulist.close()

        return hypercube_image
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
    • slope_array   : np.array - values in electrons/s of the illumination
    """
    filename=os.path.join(illum_model_filename, "illum_model_1_MIRIMAGE_P750L.fits")

    (name, ext) = os.path.splitext(illum_model_filename)
    gain=5.75

    if ext==".fits":
        hdulist      = fits.open(illum_model_filename)
        illumination = hdulist['INTENSITY'].data
        hdulist.close()
        slope_array  = illumination/gain
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


def write_det_image_with_effects(original_path, new_data, extra_metadata, overwrite=True):
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

    original_name = os.path.basename(original_path)
    original_dir = os.path.dirname(original_path)

    hdulist = fits.open(original_path)

    # Replace existing data
    hdulist[1].data = new_data

    metadata = hdulist[0].header

    metadata.update(extra_metadata)

    new_path = os.path.abspath(os.path.join(original_dir, os.path.pardir, "det_images_post_processed"))

    if not os.path.isdir(new_path):
        os.mkdir(new_path)

    new_path = os.path.join(new_path, original_name)
    hdulist.writeto(new_path, overwrite=overwrite)