
import numpy as np

from . import utils

import logging
LOG = logging.getLogger(__name__)

def response_drift(original_ramp, t_0, signal, frame=0.19):
    """
    Computes the response drift effect on the ramp. Coefficients values come from JPL tests of 2019.

    There are 3 modes:
        - 0 -> 1000    : 2 exponentials
        - 1000 -> 5000 : 1 exponential
        - over 5000    : off (no influence)

    /!\ Warning: Multiple integrations per file are not supported right now. The problem will come
                 from the t array, because for each integration you need to add the duration of
                 all the previous integrations in the file

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
    bt1000_5000_selection = (signal > 1000) & (signal < 5000)

    alpha1 = np.ones_like(signal)
    alpha2 = np.ones_like(signal)
    amp1 = np.zeros_like(signal)
    amp2 = np.zeros_like(signal)

    alpha1[lt1000_selection] = -0.251006324468 * signal[lt1000_selection]  + 256.617781046
    amp1[lt1000_selection]   = -0.000470839463392 * np.exp(0.0138190609414 * signal[lt1000_selection]) + -9.39756705304
    alpha2[lt1000_selection] = (11.5503579674 * 236.879448705**2 * signal[lt1000_selection] + 165016897.075) / (236.879448705**2 + (signal[lt1000_selection] - 241.72012081)**2)
    amp2[lt1000_selection]   = -0.0387088940065 * signal[lt1000_selection] + -0.517555367969

    alpha1[bt1000_5000_selection] = 30863.54844681253 * np.exp(-0.0005738567944909932 * signal[bt1000_5000_selection])  + 527.584922874043
    amp1[bt1000_5000_selection]   = 2.577897380949174e-06 * signal[bt1000_5000_selection]**2 + -0.029471518433861543 * signal[bt1000_5000_selection] - 10.880600687206693 


    t = t_0 + np.arange(0, nb_frames) * frame  # Time sampling in seconds
    t = t[:,np.newaxis, np.newaxis]  # Prepare broadcasting

    prefactor1 = amp1 * alpha1
    prefactor2 = amp2 * alpha2

    ramp_difference_t_0 = prefactor1 * np.exp(-t_0 / alpha1) + prefactor2 * np.exp(-t_0 / alpha2)

    # We integrate from t_0, need to remove evolution between t_0 and t_i (ramp_difference_t_0)
    ramp_difference = ramp_difference_t_0 - prefactor1 * np.exp(-t / alpha1) - prefactor2 * np.exp(-t / alpha2)

    LOG.debug("response_drift() | ramp shape : {}".format(ramp_difference.shape))

    return ramp_difference


def anneal_recovery(original_ramp, t_0, frame, config):
    """
    Computes the anneal recovery effect on the ramp. Coefficients values come from JPL tests of 2019.

    Parameters
    ----------
    original_ramp
                 Original ramp from MIRISim det_image in DN
    t_0: float
        time since beginning of the observation in seconds
    frame: float
         duration of a frame_time in seconds
    config: dict
        Configuration file for the code

    Returns
    -------
    ramp_difference
                   Loss of signal in DN on the ramp.

    """

    # time between the end of the anneal phase and the start of the observation in seconds
    anneal_time = config["anneal_recovery"]["time"]
    (nb_integrations, nb_frames, nb_y, nb_x) = original_ramp.shape


    alpha1 = 0.005051236714489755 * np.ones((nb_y, nb_x))
    amp1   = 11.600852 * np.ones((nb_y, nb_x))
    alpha2 = 0.0010083051778396745 * np.ones((nb_y, nb_x))
    amp2   = 0.86786327 * np.ones((nb_y, nb_x))

    t = t_0 + np.arange(0, nb_frames) * frame  # Time sampling in seconds
    t = t[:,np.newaxis, np.newaxis]  # Prepare broadcasting

    prefactor1 = amp1 * alpha1
    prefactor2 = amp2 * alpha2

    ramp_difference_t_0 = prefactor1 * np.exp(-(t_0 + anneal_time) / alpha1) + prefactor2 * np.exp(-(t_0 + anneal_time) / alpha2)

    # We integrate from t_0, need to remove evolution between t_0 and t_i (ramp_difference_t_0)
    ramp_difference = ramp_difference_t_0 - prefactor1 * np.exp(-(t + anneal_time) * alpha1) - prefactor2 * np.exp(-(t + anneal_time) * alpha2)

    LOG.debug("anneal_recovery() | ramp shape : {}".format(ramp_difference.shape))
    return ramp_difference


def idle_recovery(original_ramp, t_0, signal, frame, config):
    """
    Computes the idle recovery effect on the ramp. Coefficients values come from JPL tests of 2019.

    Parameters
    ----------
    original_ramp
                 Original ramp from MIRISim det_image in DN
    t_0: float
        time since beginning of the observation in seconds
    signal: np.ndarray
          illum_image value (either full image or value for one pixel) in DN/s
    frame: float
         duration of a frame_time in seconds
    config: dict
        Configuration file for the code

    Returns
    -------
    ramp_difference
                   Loss of signal in DN on the ramp.

    """

    # duration onf the idle phase before the observation began in seconds
    idle_time = config["idle_recovery"]["duration"]
    (nb_integrations, nb_frames, nb_y, nb_x) = original_ramp.shape

    alpha1 = 2975790.21677 * np.exp( -0.0173557837941 * signal ) + 1133.82032361
    amp1   = ((2.04271769088e-05 * signal**2 + -0.0166816654355 * signal + 4.08114118945)/2011.9) * idle_time

    t = t_0 + np.arange(0, nb_frames) * frame  # Time sampling in seconds
    t = t[:,np.newaxis, np.newaxis]  # Prepare broadcasting

    ramp_difference_t_0 = (amp1 * alpha1) * np.exp(-t_0 / alpha1)

    # We integrate from t_0, need to remove evolution between t_0 and t_i (ramp_difference_t_0)
    ramp_difference = ramp_difference_t_0 - (amp1 * alpha1) * np.exp(-t / alpha1)

    LOG.debug("idle_recovery() | ramp shape : {}".format(ramp_difference.shape))
    return ramp_difference


def poisson_noise(original_ramp):
    """
    Compute Poisson noise on all integration of a det_image data cube.

    Parameters
    ----------
    original_ramp
                 Original ramp from MIRISim det_image in DN. Dimensions: (nb_integrations, nb_frames, nb_y, nb_x)

    Returns
    -------
    data cube with Poisson noise added

    """
    #TODO check negative values and decide what to do.
    # Not working at the moment
    frame_differences = np.diff(original_ramp, axis=1)

    first_frame = original_ramp[:,0,:,:]
    first_frame = first_frame[:,np.newaxis,:,:]
    frame_differences = np.append(first_frame, frame_differences, axis=1)

    single_frame_noise = np.random.poisson(abs(frame_differences))

    noised_ramp = np.cumsum(single_frame_noise, axis=1)

    LOG.debug("poisson_noise() |  noised ramp shape : {}".format(noised_ramp.shape))
    return noised_ramp
