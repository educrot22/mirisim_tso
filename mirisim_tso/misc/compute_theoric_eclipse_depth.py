# IPython log file

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
#

from astropy.modeling import models

def compute_theoric_eclipse_depth(waves,planet_temperature=1600.*u.Kelvin,star_temperature = 4400.*u.Kelvin ):

    bb_planet = models.BlackBody(temperature=planet_temperature)
    bb_star = models.BlackBody(temperature=star_temperature)
    depth = bb_planet(waves)/bb_star(waves)* 0.15961341**2
    return depth
