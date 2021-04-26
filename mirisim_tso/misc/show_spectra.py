# IPython log file

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u

from mirisim_tso.misc.compute_theoric_eclipse_depth import compute_theoric_eclipse_depth
from mirisim_tso.misc.readwrite_spectra import write_depth_ecsv

def show_spectra(flux2d, wave1d, orbital_phase, period, SAVE=False, model="", pattern=""):

    ttn = (orbital_phase - orbital_phase.min())*period
    ttn = ttn.to(u.second)

    nw, nt = flux2d.shape
    i6 = np.argmin(np.abs(wave1d -6*u.micron))
    i61 = np.argmin(np.abs(wave1d -6.1*u.micron))
    i7 = np.argmin(np.abs(wave1d-7*u.micron))
    i10 = np.argmin(np.abs(wave1d-10*u.micron))
    i11 = np.argmin(np.abs(wave1d-11*u.micron))
    print(nw, nt)

    ####  SPECTRA  ####
    ip = np.argmin(np.abs(orbital_phase+0.5))
    plt.figure()
    plt.title(model+'spectrum')
    plt.plot(wave1d, flux2d[:, 0], label="phase={:.2f}".format(orbital_phase[0]))
    plt.plot(wave1d, flux2d[:, ip], label="phase={:.2f}".format(orbital_phase[ip]))
    plt.legend()
    plt.xlabel('wavelength '+str(wave1d.unit))
    plt.ylabel('flux '+str(flux2d.unit))
    if (SAVE):plt.savefig('plot_spectra_'+pattern+'.png')
   
    ####  LIGHT CURVES
    plt.figure()
    plt.title(model+'light curve')
    plt.plot(orbital_phase, flux2d[i6,:], label="wavelength={:.2f}".format(wave1d[i6]))
    #plt.plot(orbital_phase, flux2d[i61,:], label="wavelength={:.2f}".format(wave1d[i61]))
    plt.legend()
    plt.xlabel('phase')
    plt.ylabel('flux '+str(flux2d.unit))
    if (SAVE):plt.savefig('plot_eclipse_'+pattern+'_lightcurve_6_6p1.png')

    ####  NORMALISED LIGHT CURVES
    plt.figure()
    plt.title(model+'normalised light curve')
    plt.plot(orbital_phase, flux2d[i6,:]/flux2d[i6,:].max(), label='6 micron depth=0.57%')
    plt.plot(orbital_phase, flux2d[i7,:]/flux2d[i7,:].max(), label='7 micron depth=0.62%')
    plt.legend()
    plt.xlabel('orbital_phase')
    plt.ylabel('normalised flux')
    if (SAVE):plt.savefig('plot_eclipse_'+pattern+'_normalised_lightcurve_6_7.png')


    ####  ECLIPSED DEPTH
    flux_min = flux2d.min(1)
    flux_max = flux2d.max(1)
    depth = (flux_max-flux_min)/flux_max*100
    plt.figure()
    plt.title(model+'measured depth minmax method')
    plt.plot(wave1d, depth)
    plt.xlabel('wavelength micron')
    plt.ylabel('depth %')
    plt.plot([6,6], [0.3,0.568])
    plt.plot([7,7], [0.3,0.626])
    if (SAVE):
        plt.savefig('plot_measured_depth_'+pattern+'.png')
        write_depth_ecsv('table_depth_'+pattern+'.dat', depth, wave1d.to(u.micron).value)
    return
