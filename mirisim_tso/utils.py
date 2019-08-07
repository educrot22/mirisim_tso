import os
from astropy.io import fits

def read_det_image(directory, x, y):
    """
    . Function to read the variable if it comes from a FITS file. Only takes |.fits containing detector image from
    MIRISim

    . Adapted from exoNoodle project (exonoodle.utils)

    >>>> PARAMETERS
    • directory     : String - Name of the entry directory containing the det_images to read

    >>>> RETURNS
    • original_ramp : np.array - values in DN of a ramp integration
    """
    filename=os.path.join(directory, "det_image_seq1_MIRIMAGE_P750Lexp1.fits")
    (name, ext) = os.path.splitext(filename)

    if ext==".fits":
        hdulist      = fits.open(filename)
        hypercube_image=hdulist[1].data
        hdulist.close()

        original_ramp = hypercube_image[5, :, x, y]
        original_ramp = original_ramp - hypercube_image[5, :, 4, 4]
        return original_ramp

    else:
        print("ERROR :: The illum_model file shall be a .fits file !")
        raise ValueError


def read_illum_model(directory):
    """
    . Function to read the variable if it comes from a FITS file. Only takes |.fits containing illumination model
    from MIRISim

    . Adapted from exoNoodle project (exonoodle.utils)

    >>>> PARAMETERS
    • directory     : String - Name of the entry directory containing the illum model to read

    >>>> RETURNS
    • slope_array   : np.array - values in electrons/s of the illumination
    """
    filename=os.path.join(directory, "illum_model_1_MIRIMAGE_P750L.fits")

    (name, ext) = os.path.splitext(filename)
    gain=5.75

    if ext==".fits":
        hdulist      = fits.open(filename)
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