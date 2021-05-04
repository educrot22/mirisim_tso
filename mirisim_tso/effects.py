
import numpy as np

from . import constants as c

import logging
LOG = logging.getLogger(__name__)


def response_drift_one(original_ramp, t_0, signal, t_frame=0.159):
    """
    This is a simpiflied version of response_drift, with only one exponential model.
    It computes the response drift effect on the ramp. Coefficients values come from
    JPL tests of 2019.
    For pixels wich values are more than 5000, bright pixel there is no drift.
    For pixels wich values are less than 75, faint pixel there is no drift.
    The sampling time for MIRI LRS is 0.159
https://jwst-docs.stsci.edu/mid-infrared-instrument/miri-instrumentation/miri-detector-overview/miri-detector-subarrays

    /!\ Warning: Multiple integrations per file are not supported right now. The problem will come
                 from the t array, because for each integration you need to add the duration of
                 all the previous integrations in the file

    Parameters
    ----------
    :param original_ramp: hyper-cube of ramps from MIRISim det_image in DN
                 only the shape is used !
                 nb_y=416 and nb_x=72  without the reference columns ???
    :type original_ramp : np.array(nb_integrations, nb_frames, nb_y, nb_x)
    
    :param t_0: time since beginning of the observation in seconds
    :type t_0: float
    
    :param signal: illum_image value (either full image or value for one pixel) in DN/s
         (nb_y=416, nb_x=76)   nb_x includes 4 reference columns
    :type signal : np.array(nb_y, nb_x) or float


    :param t_frame: duration of a frame in seconds [default for MIRI-LRS = 0.159s]
    :type t_frame: float

    -------
    :return: ramp_difference Loss of signal in DN on the ramp, same shape than original_ram
              the dimensions are chosen to be  hypercube to be compatible with original_ramp
    :rtype:np.array(1, nb_integrations, nb_frames, nb_y, nb_x)

    """
    # Add a test for the original ramp, if it is a 3D array. (x, y, t_0)
    # - If yes, we have a det_image, and the dispatch on the different formulas needs to be done with np.where()
    # - If not, it needs to be turned into a 3D array, before np.where()
    # Can it be done outside this function ?
    low_threshold =  75  # to avoid negative pixel
    fading_threshold = 5000

    (nb_integrations, nb_frames, nb_y, nb_x) = original_ramp.shape
    
    (nbs_y, nbs_x) = signal.shape
    if( nbs_x == (nb_x+4)):
        LOG.warning('problem of reference pixels,{} remove the 4 first columns from illumination model{}'.format(original_ramp.shape, signal.shape))
        signal = signal[:,4:]

    # Creating two pixels masks corresponding to two different fitting regimes
    index = np.where( (signal > low_threshold) & (signal < fading_threshold) )

    alpha1 = np.ones_like(signal)
    amp1 = np.zeros_like(signal)

    # Values fitted from JPL test data
     
    alpha1[index] = 2.59295558e+03 * np.exp(-8.57428099e-04 * signal[index]) + 1.20593193e+02

    amp1[index] = 2.94507350e-06 * signal[index]**2 + -3.27886892e-02 * signal[index]\
    + -4.70669170e+00
    

    t = t_0 + np.arange(0, nb_frames) * t_frame  # Time sampling in seconds
    t = t[:, np.newaxis, np.newaxis]  # Prepare broadcasting

    prefactor1 = amp1 * alpha1
   
    ramp_difference_t_0 = prefactor1 * np.exp(-t_0 / alpha1)

    # We integrate from t_0, need to remove evolution between t_0 and t_i (ramp_difference_t_0)
    ramp_difference = ramp_difference_t_0 - prefactor1 * np.exp(-t / alpha1)
    #import pdb pdb.set_trace()
    ramp_difference = ramp_difference.reshape([1, nb_frames, nb_y, nb_x]) # RG convention hypercube
    ramp_difference = np.float32(ramp_difference)
    LOG.debug("response_drift_one() | ramp shape={:}, dtype={:} ".format(ramp_difference.shape, ramp_difference.dtype))

    return ramp_difference


def response_drift(original_ramp, t_0, signal, frame=0.159):
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
         duration of a frame_time in seconds [default for MIRI-LRS = 0.159s]

    Returns
    -------
    ramp_difference
                   Loss of signal in DN on the ramp.

    """
    # Add a test for the original ramp, if it is a 3D array. (x, y, t_0)
    # - If yes, we have a det_image, and the dispatch on the different formulas needs to be done with np.where()
    # - If not, it needs to be turned into a 3D array, before np.where()
    # Can it be done outside this function ?
    transition_threshold =  750
    fading_threshold = 5000

    (nb_integrations, nb_frames, nb_y, nb_x) = original_ramp.shape

    # Creating two pixels masks corresponding to two different fitting regimes
    faint_pixels = signal < transition_threshold
    medium_pixels = (~faint_pixels) & (signal < fading_threshold)

    alpha1 = np.ones_like(signal)
    alpha2 = np.ones_like(signal)
    amp1 = np.zeros_like(signal)
    amp2 = np.zeros_like(signal)

    # Values fitted from JPL test data, obtained with notebook extract_responsedrift_expressions.ipynb
    alpha1[faint_pixels] = -0.24824707509094024 * signal[faint_pixels] + 256.0307052443123
    amp1[faint_pixels] = -3.309038643901135 * np.exp(0.002099958184394631 * signal[faint_pixels]) + -3.1679206249778638
    alpha2[faint_pixels] = (11.128462112080102 * 55332.60403145965 * signal[faint_pixels] + 168405550.01410204) / (55332.60403145965 + (signal[faint_pixels] - 245.45591611338756)**2)
    amp2[faint_pixels] = -0.041975980807126514 * signal[faint_pixels] + 0.28375014129499576

    alpha1[medium_pixels] = 2643.2695863041595 * np.exp(-0.0008722598069934617 * signal[medium_pixels]) + 125.05736916420666
    amp1[medium_pixels] = 3.311621436854758e-06 * signal[medium_pixels]**2 + -0.034646597553672186 * signal[medium_pixels] + -2.551106144746854

    t = t_0 + np.arange(0, nb_frames) * frame  # Time sampling in seconds
    t = t[:, np.newaxis, np.newaxis]  # Prepare broadcasting

    prefactor1 = amp1 * alpha1
    prefactor2 = amp2 * alpha2

    ramp_difference_t_0 = prefactor1 * np.exp(-t_0 / alpha1) + prefactor2 * np.exp(-t_0 / alpha2)

    # We integrate from t_0, need to remove evolution between t_0 and t_i (ramp_difference_t_0)
    ramp_difference = ramp_difference_t_0 - prefactor1 * np.exp(-t / alpha1) - prefactor2 * np.exp(-t / alpha2)
    ramp_difference = np.float32(ramp_difference)
    LOG.debug("response_drift() | ramp shape={:}, dtype={:} ".format(ramp_difference.shape, ramp_difference.dtype))

    return np.float32(ramp_difference)


def anneal_recovery(original_ramp, t_0, frame, config):
    """
    Computes the anneal recovery effect on the ramp. Coefficients values come from JPL tests of 2019.

    In config, will read config["anneal_recovery"]["time"]. This number is the number of seconds between the anneal and
    the beginning of the observation. if 0, this means the observation start right after the anneal.
    The default value, 600, means the observation starts 10 minutes after the anneal, tought to be the minimum time possible.

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

    beta1 = 197.97132 * np.ones((nb_y, nb_x))
    amp1   = 11.600852 * np.ones((nb_y, nb_x))
    beta2 = 991.76323 * np.ones((nb_y, nb_x))
    amp2   = 0.86786327 * np.ones((nb_y, nb_x))

    t = t_0 + np.arange(0, nb_frames) * frame  # Time sampling in seconds
    t = t[:, np.newaxis, np.newaxis]  # Prepare broadcasting

    prefactor1 = amp1 * beta1
    prefactor2 = amp2 * beta2

    ramp_difference_t_0 = prefactor1 * np.exp(-(t_0 + anneal_time) / beta1) + prefactor2 * np.exp(-(t_0 + anneal_time) / beta2)

    # We integrate from t_0, need to remove evolution between t_0 and t_i (ramp_difference_t_0)
    ramp_difference = ramp_difference_t_0 - prefactor1 * np.exp(-(t + anneal_time) / beta1) - prefactor2 * np.exp(-(t + anneal_time) / beta2)
    ramp_difference = np.float32(ramp_difference)
    LOG.debug("anneal_recovery() | ramp shape={:}, dtype={:} ".format(ramp_difference.shape, ramp_difference.dtype))
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

    alpha1 = 2975790.21677 * np.exp(-0.0173557837941 * signal) + 1133.82032361
    amp1 = ((2.04271769088e-05 * signal**2 + -0.0166816654355 * signal + 4.08114118945)/2011.9) * idle_time

    t = t_0 + np.arange(0, nb_frames) * frame  # Time sampling in seconds
    t = t[:, np.newaxis, np.newaxis]  # Prepare broadcasting

    ramp_difference_t_0 = (amp1 * alpha1) * np.exp(-t_0 / alpha1)

    # We integrate from t_0, need to remove evolution between t_0 and t_i (ramp_difference_t_0)
    ramp_difference = ramp_difference_t_0 - (amp1 * alpha1) * np.exp(-t / alpha1)
    ramp_difference = np.float32(ramp_difference)
    LOG.debug("idle_recovery() | ramp shape={:}, dtype={:} ".format(ramp_difference.shape, ramp_difference.dtype))
    return ramp_difference


def poisson_noise(original_ramp, mask):
    """
    Compute Poisson noise on all integration of a det_image data cube.
    We follow the technical note of Massimo Roberto WFC3-2007-12.pdf, paragraph 2.4, equation 1.40
    y_i = y_i-1 + p_i
    where the various p_i are statistically independent packets of electrons.
    So the Poisson distribution is applied to p_i, not y_i
    http://web.ipac.caltech.edu/staff/fmasci/home/astro_refs/SUR_vs_CDS.pdf
    Parameters
    ----------
    original_ramp:
                 Original ramp from MIRISim det_image in DN. Dimensions: (nb_integrations, nb_frames, nb_y, nb_x)
    mask:
        np.array(bool) - Array of bad pixels (True if bad, False if good)
                 Needed because the bad pixels have non-additive shapes where computation is not applicable. They need
                 to be excluded from the computation.

    Returns
    -------
    data cube with Poisson noise added (this is not a ramp difference, this is the full ramp)

    """
    # get the shape and check that the cube is 4D
    (nb_integrations, nb_frames, nb_y, nb_x) =  original_ramp.shape
    # create the hypercube of differences in electron
    #  I duplicate the  first image of difference
    frame_differences = np.diff(original_ramp, axis=1)
    first_difference = frame_differences[0,0,:,:]
    first_difference = first_difference.reshape([1, 1, nb_y, nb_x])
    frame_differences = np.append(first_difference, frame_differences, axis=1)
    # now frame_differences has the same shpae than original_ramp

    # add the noise, comput in electron
    ## BUG PATCH NP.ABS
    frame_noise = np.float32(np.random.poisson(np.abs(frame_differences*c.gain)))/c.gain
    
    # add the first image of the ramp
    frame_noise[0,0,:,:] += (original_ramp[0,0,:,:] - first_difference[0,0,:,:])

    # cumulative sum
    noised_ramp = np.cumsum(frame_noise, axis=1)

    # This works only for the good pixels (which accumulate signal). We use the bad pixels CDP from MIRISim
    # Overwrite bad pixels with the original value.
    bad_pixels = np.broadcast_to(mask, original_ramp.shape)
    noised_ramp[bad_pixels] = original_ramp[bad_pixels]

    LOG.debug("poisson_noise() |  noised ramp shape={:}, dtype={:} ".format(ramp_difference.shape, ramp_difference.dtype))
    return noised_ramp
