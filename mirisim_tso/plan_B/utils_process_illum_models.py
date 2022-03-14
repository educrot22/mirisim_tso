# ADD EFFECTS

import astropy.units as u
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np


###############################################

def add_poisson_noise(cube, sampling_time):
    """
    Compute Poisson noise on all integration of a det_image data cube.

    Parameters
    ----------
    cube  quantity in electron/second
    gain  quantity in electron/adu
    sampling_time  quantity in second

    Returns
    -------
    data cube with Poisson noise added

    """
    cube_elec = (cube*sampling_time).decompose()
    ii  = np.where(cube_elec.value < 0)
    if (len(ii[0]) > 0):
        print("warning "+str(len(ii[0]))+ "negative values" )
    ii  = np.where(cube_elec.value >= 0)
    noised_cube = np.zeros(cube.shape)
    noised_cube[ii] = np.random.poisson(cube_elec[ii].value)
    noised_cube = noised_cube*cube_elec.unit
    noised_cube = noised_cube/sampling_time
    return np.float32(noised_cube)

###############################################
def response_drift( signal, tint=47.7):
    """
    Computes the response drift effect on the signl. Coefficients values come from JPL tests of 2019.

        - 0 -> 1000    : 2 exponentials
       
    Parameters
    ----------

    signal
          illum_image value (either full image or value for one pixel) in DN/s


    Returns
    -------
    signal_difference
                   Loss of signal in DN/s on the ramp.

    """
    # Add a test for the original ramp, if it is a 3D array. (x, y, t_0)
    # - If yes, we have a det_image, and the dispatch on the different formulas needs to be done with np.where()
    # - If not, it needs to be turned into a 3D array, before np.where()
    # Can it be done outside this function ?
    transition_threshold =  750
    fading_threshold = 5000

    (nb_frames, nb_y, nb_x) = signal.shape

    # Creating two pixels masks corresponding to two different fitting regimes
    if (signal.min() < 0):
        print('EROR SIGNAL < 0')
        return None
    if (signal.max() > transition_threshold):
        print('EROR SIGNAL > ', transition_threshold)
        return None

    alpha1 = np.ones_like(signal)
    alpha2 = np.ones_like(signal)
    amp1 = np.zeros_like(signal)
    amp2 = np.zeros_like(signal)

    # Values fitted from JPL test data, obtained with notebook extract_responsedrift_expressions.ipynb
    alpha1 = -0.24824707509094024 * signal + 256.0307052443123
    amp1   = -3.309038643901135 * np.exp(0.002099958184394631 * signal) + -3.1679206249778638

    alpha2 = (11.128462112080102 * 55332.60403145965 * signal + 168405550.01410204) / (55332.60403145965 + (signal - 245.45591611338756)**2)

    amp2 = -0.041975980807126514 * signal + 0.28375014129499576


    #
    t = np.arange(nb_frames)*tint
    t = t.reshape([nb_frames, 1, 1])
    signal_difference =  amp1*np.exp(-t / alpha1) + amp2*np.exp(-t / alpha2)
    signal_difference = np.float32(signal_difference)

    return signal_difference
