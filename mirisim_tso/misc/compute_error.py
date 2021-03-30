
from astropy import units as u
import numpy as np

def compute_error(wavelength, flux, depth, tint=103.5*u.second):
    print(flux.shape, wavelength.shape)
    # converts the flux into u.photon/u.m**2/u.micron/u.s/u.steradian using the wavelength equivalency
    flux2 = flux.to(u.photon / u.m ** 2 / u.micron / u.s, equivalencies=u.spectral_density(wavelength))
    surface = 25*u.meter*u.meter
    dwave=0.06*u.micron
    flux3 = (flux2*dwave*surface*tint).decompose()
    print('flux3.unit', flux3.unit)
    error = (depth + np.sqrt(2))/np.sqrt(flux3.value)
    return error
