
import numpy as np

from . import utils

import logging
LOG = logging.getLogger('mirisim_tso.utils')

def response_drift(original_ramp, t_0, signal, frame=0.19):
    """

    Parameters
    ----------
    original_ramp
                 Original ramp from MIRISim det_image in DN
    t_0
        time since beginning of the observation in seconds
    signal
          illum_image value (either full image or value for one pixel) in DN/s
    frame
         duration of a frame in seconds

    Returns
    -------

    """

    LOG.debug("response_drift() | signal={}".format(signal))

    nb_frames = np.size(original_ramp)

    alpha1 = -0.251006324468 * signal  + 256.617781046
    amp1   = -0.000470839463392 * np.exp(0.0138190609414 * signal) + -9.39756705304
    alpha2 = (11.5503579674 * 236.879448705**2 * signal + 165016897.075) / (236.879448705**2 + (signal - 241.72012081)**2)
    amp2   = -0.0387088940065 * signal + -0.517555367969

    t = t_0 + np.arange(0, nb_frames) * frame  # Time sampling in seconds

    ramp_difference_t_0 = (amp1 * alpha1) * np.exp(-t_0 / alpha1) + (amp2 * alpha2) * np.exp(-t_0 / alpha2)

    # We integrate from t_0, start of the ramp to t_i
    ramp_difference = ramp_difference_t_0 - (amp1 * alpha1) * np.exp(-t / alpha1) - (amp2 * alpha2) * np.exp(-t / alpha2)

    return ramp_difference


def anneal_recovery(original_ramp, t, x, y):

    ramp_difference=0

    return ramp_difference


def idle_recovery(original_ramp, t, x, y):

    ramp_difference=0

    return ramp_difference
