# Extract spectra for a cube of images
# very simple algorithm :
#    sum all the good pixels in a line (good defined by a mask)
#    compute the background
#    subtract the backtround (bad pixels, as given by 1-mask)
# This does not take into account the curvature
# An optimal spectrum extraction procedure is given by Horne, 1986, and implemented in Cascade
#  So use cascade for optimal result !
# Author Rene Gastaud, CEA-Saclay
# Date June 2022

from astropy.io import fits
from astropy.io import ascii
import matplotlib.pyplot as plt
import numpy as np

def extract_spectra_with_mask(cube, mask, plot=False):
    #ny, nx = image.shape
    # (44, 416, 72)  nz, ny, nx
    nb_pixel_in = mask.sum(axis=1)
    background2d = np.median((cube*(1-mask)), axis=2)
    flux2d = (cube*mask).sum(axis=2)
    flux2d_bkg = (flux2d - background2d*nb_pixel_in)
    #
    return flux2d_bkg
