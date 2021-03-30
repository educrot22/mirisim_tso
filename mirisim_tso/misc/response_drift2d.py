#!/usr/bin/env python
# -*- coding: utf-8 -*-

# https://gitlab.com/mmartin-lagarde/mirisim_tso/-/tree/corrected_rd/mirisim_tso
#  file effects.py function response_drift


# history add variable starting with underscore  _time, _input_flux, _marine_flux
#  for minimisation purpose
#  flux = slope DN/s
import numpy as np
import matplotlib.pyplot as plt


def response_drift2d( image, t ):
    """
    Computes the response drift effect on the ramp. Coefficients values come from JPL tests of 2019.
    The response drift has been simplified with 2 modes only : low fluxes, high fluxes.

    There are 2 modes:
        - 0 -> 5000 : one exponential
        - over 5000    : off (no influence)

    Parameters
    ----------
    
    image
          illum_images value (either full image or value for one pixel) in DN/s
          (np.array dimension ny, nx)
    t
         time grid second  (np.array dimension nb_frames)

    Returns
    -------
    flux_difference Loss of flux in DN/s  (np.array dimension nb_frames, nb_y, nb_x)

    """
    # Add a test for the original ramp, if it is a 3D array. (x, y, t_0)
    
    fading_threshold = 5000   # *u.adu/u.second  # DN/s

    (nb_y, nb_x) = image.shape
    nb_frames = t.shape[0] ## I hate python

    # Creating two pixels masks corresponding to two different fitting regimes
    low_pixels = image < fading_threshold

    alpha1 = np.ones_like(image)
    amp1 = np.zeros_like(image)

    # Values fitted from Marine Martin_Lagarde Mirisim_tso modelled data,
    # obtained by fitting the data from nnotebook extract_responsedrift_expressions.ipynb
    alpha1[low_pixels] = 2.59295558e+03 * np.exp(-8.57428099e-04 * image[low_pixels]) + 1.20593193e+02
    amp1[low_pixels] = 2.94507350e-06 * image[low_pixels]**2 + -3.27886892e-02 * image[low_pixels]\
                          + -4.70669170e+00

    # (nb_frames, nb_y, nb_x)
    # Prepare broadcasting
    tt = t.reshape([nb_frames,1, 1])
    amp1 = amp1.reshape([1, nb_y, nb_x])
    alpha1 = alpha1.reshape([1,  nb_y, nb_x])
    
    flux_difference = amp1 * np.exp(-tt / alpha1)

    return flux_difference

