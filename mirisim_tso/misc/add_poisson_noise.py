
# https://gitlab.com/mmartin-lagarde/mirisim_tso/-/tree/corrected_rd/mirisim_tso
#  file effects.py function  poisson_noise

import numpy as np
from astropy import units as u


def add_poisson_noise(cube, gain, sampling_time):
    """
    Compute Poisson noise on all integration of a det_image data cube.

    Parameters
    ----------
    cube  quantity in adu/second
    gain  quantity in electron/adu
    sampling_time  quantity in second

    Returns
    -------
    data cube with Poisson noise added

    """
    cube_elec = (cube*gain*sampling_time).decompose()
    ii  = np.where(cube_elec.value > 0)
    if (len(ii[0]) > 0):
        print("warning "+str(len(ii[0]))+ "negative values" )
    noised_cube = np.zeros(cube.shape)
    noised_cube[ii] = np.random.poisson(cube_elec[ii].value)*cube_elec.unit
    noised_cube = noised_cube/gain/sampling_time
    return noised_cube
