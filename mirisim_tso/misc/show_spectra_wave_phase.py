# Plot Spectra function of phase and wavelength

import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u

"""
my_title = 'WASP-80 transit'
suffix = '.png'
"""
########################
def show_spectra_wave_phase(flux, wavelength, phase, my_title, suffix = '.png', save_file=None):
    ref_waves = np.array([5.152243, 7.014417, 8.881744, 12.])#*u.micron
    # x1dints 12.008793404814282 5.022831579805857
    ##
    indexes = np.zeros(ref_waves.size, dtype='int')
    for i in np.arange(ref_waves.size):
        indexes[i] = np.argmin(np.abs(wavelength.value-ref_waves[i]))

    for i in np.arange(ref_waves.size):
        print('{} wavelength={:4.3f},'.format(i, wavelength[indexes[i]]))
    """
    0 sed_wavelength=5.154 micron, out_wavelength=5.153
    1 sed_wavelength=7.014 micron, out_wavelength=7.015
    2 sed_wavelength=8.880 micron, out_wavelength=8.880
    3 sed_wavelength=12.000 micron, out_wavelength=12.009
    """
    #########################################################
    ################## last spectrum ##########################
    plt.figure()
    plt.title(my_title+' last spectrum')
    plt.plot(wavelength, flux[:,-1].to(u.milliJansky))
    plt.xlabel('wavelength micron')
    plt.ylabel('flux milliJansky')
    plt.legend()
    if (save_file is not None):
        plt.savefig(save_file+"last_spectrum"+suffix)
        
    #########################################################
    ################## light curve ##########################

    for i in np.arange(indexes.size):
        index = indexes[i]
        plt.figure()
        plt.title(my_title+' Light Curve for wavelength={:4.3f}'.format(wavelength[index]))
        plt.plot(phase, flux[index,:].to(u.milliJansky))
        #, label='wavelength {:4.3f}'.format(wavelength[index]) )
        plt.ylabel('flux milliJansky')
        plt.xlabel('phase')
        
        if (save_file is not None):
            plot_name = save_file+'lightCurve_wavelength_{:4.3f}'.format(wavelength[index])+suffix
            plt.savefig(plot_name)
    return
