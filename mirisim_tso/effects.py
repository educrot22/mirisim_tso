
import numpy as np

from . import utils

import logging
LOG = logging.getLogger('mirisim_tso.utils')

def response_drift(original_ramp, t_0, signal, frame=0.19):
    """
    Computes the response drift effect on the ramp. Coefficients values come from JPL tests of 2019.

    There are 3 modes:
        - 0 -> 1000    : 2 exponentials
        - 1000 -> 5000 : 1 exponential
        - over 5000    : off (no influence)

    Parameters
    ----------
    original_ramp
                 Original ramp from MIRISim det_image in DN
    t_0
        time since beginning of the observation in seconds
    signal
          illum_image value (either full image or value for one pixel) in DN/s
    frame
         duration of a frame_time in seconds [default for MIRI-LRS = 0.19s]

    Returns
    -------
    ramp_difference
                   Loss of signal in DN on the ramp.

    """
    # Add a test for the original ramp, if it is a 3D array. (x, y, t_0)
    # - If yes, we have a det_image, and the dispatch on the different formulas needs to be done with np.where()
    # - If not, it needs to be turned into a 3D array, before np.where()
    # Can it be done outside this function ?

    (nb_integrations, nb_frames, nb_y, nb_x) = original_ramp.shape

    # For signal < 1000
    lt1000_selection = signal < 1000

    alpha1 = np.ones_like(signal)
    alpha2 = np.ones_like(signal)
    amp1 = np.zeros_like(signal)
    amp2 = np.zeros_like(signal)

    alpha1[lt1000_selection] = -0.251006324468 * signal[lt1000_selection]  + 256.617781046
    amp1[lt1000_selection]   = -0.000470839463392 * np.exp(0.0138190609414 * signal[lt1000_selection]) + -9.39756705304
    alpha2[lt1000_selection] = (11.5503579674 * 236.879448705**2 * signal[lt1000_selection] + 165016897.075) / (236.879448705**2 + (signal[lt1000_selection] - 241.72012081)**2)
    amp2[lt1000_selection]   = -0.0387088940065 * signal[lt1000_selection] + -0.517555367969

    t = t_0 + np.arange(0, nb_frames) * frame  # Time sampling in seconds
    t = t[:,np.newaxis, np.newaxis]  # Prepare broadcasting

    prefactor1 = amp1 * alpha1
    prefactor2 = amp2 * alpha2

    ramp_difference_t_0 = prefactor1 * np.exp(-t_0 / alpha1) + prefactor2 * np.exp(-t_0 / alpha2)

    # We integrate from t_0, need to remove evolution between t_0 and t_i (ramp_difference_t_0)
    ramp_difference = ramp_difference_t_0 - prefactor1 * np.exp(-t / alpha1) - prefactor2 * np.exp(-t / alpha2)

    LOG.debug("response_drift() | computed with success")

    return ramp_difference


def anneal_recovery(original_ramp, anneal_time, t_0, frame=0.19):
    """
    Computes the anneal recovery effect on the ramp. Coefficients values come from JPL tests of 2019.

    Parameters
    ----------
    original_ramp
                 Original ramp from MIRISim det_image in DN
    anneal_time
               time between the end of the anneal phase and the start of the observation in seconds
    t_0
        time since beginning of the observation in seconds
    frame
         duration of a frame_time in seconds [default for MIRI-LRS = 0.19s]

    Returns
    -------
    ramp_difference
                   Loss of signal in DN on the ramp.

    """
    nb_frames = np.size(original_ramp)

    alpha1 = 0.005051236714489755
    amp1   = 11.600852
    alpha2 = 0.0010083051778396745
    amp2   = 0.86786327


    t = t_0 + np.arange(0, nb_frames) * frame  # Time sampling in seconds
    ramp_difference_t_0 = (amp1 * alpha1) * np.exp(-(t_0 + anneal_time) / alpha1) + (amp2 * alpha2) * np.exp(-(t_0 + anneal_time) / alpha2)

    # We integrate from t_0, need to remove evolution between t_0 and t_i (ramp_difference_t_0)
    ramp_difference = ramp_difference_t_0 - (amp1 * alpha1) * np.exp(-(t + anneal_time) * alpha1) - (amp2 * alpha2) * np.exp(-(t + anneal_time) * alpha2)

    LOG.debug("anneal_recovery() | computed with success")
    return ramp_difference


def idle_recovery(original_ramp, idle_time, t_0, signal, frame=0.19):
    """
    Computes the idle recovery effect on the ramp. Coefficients values come from JPL tests of 2019.

    Parameters
    ----------
    original_ramp
                 Original ramp from MIRISim det_image in DN
    idle_time
             duration onf the idle phase before the observation began in seconds
    t_0
        time since beginning of the observation in seconds
    signal
          illum_image value (either full image or value for one pixel) in DN/s
    frame
         duration of a frame_time in seconds [default for MIRI-LRS = 0.19s]

    Returns
    -------
    ramp_difference
                   Loss of signal in DN on the ramp.

    """
    nb_frames = np.size(original_ramp)

    alpha1 = 2975790.21677 * np.exp( -0.0173557837941 * signal ) + 1133.82032361
    amp1   = ((2.04271769088e-05 * signal**2 + -0.0166816654355 * signal + 4.08114118945)/2011.9)*idle_time

    t = t_0 + np.arange(0, nb_frames) * frame  # Time sampling in seconds
    ramp_difference_t_0 = (amp1 * alpha1) * np.exp(-t_0 / alpha1)

    # We integrate from t_0, need to remove evolution between t_0 and t_i (ramp_difference_t_0)
    ramp_difference = ramp_difference_t_0 - (amp1 * alpha1) * np.exp(-t / alpha1)

    LOG.debug("idle_recovery() | computed with success")
    return ramp_difference
