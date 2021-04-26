# IPython log file

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from scipy.ndimage import gaussian_filter1d

from mirisim_tso.misc.compute_theoric_eclipse_depth import compute_theoric_eclipse_depth
from mirisim_tso.misc.readwrite_spectra import write_depth_ecsv

def show_spectra_tso(flux2d, wave1d, orbital_phase, period, SAVE=False, model="", pattern=""):

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
    plt.title(model+'light curve wavelength={:.2f}'.format(wave1d[i6]))
    plt.plot(orbital_phase, flux2d[i6,:])
    #plt.title(model+'light curve')
    #plt.plot(orbital_phase, flux2d[i6,:],  label="wavelength={:.2f}".format(wave1d[i6]))
    #plt.plot(orbital_phase, flux2d[i61,:], label="wavelength={:.2f}".format(wave1d[i61]))
    #plt.legend()
    plt.xlabel('phase')
    plt.ylabel('flux '+str(flux2d.unit))
    if (SAVE):plt.savefig('plot_eclipse_'+pattern+'_lightcurve_6_6p1.png')

    ####  NORMALISED LIGHT CURVES
    ligne6 = gaussian_filter1d(flux2d[i6,:].value, 3)
    ligne7 = gaussian_filter1d(flux2d[i7,:].value, 3)
    plt.figure()
    plt.title(model+'normalised light curve')
    plt.plot(ttn, ligne6/ligne6.max(), 'g', label='wavelength 6. micron depth=0.57%')
    plt.plot(ttn, flux2d[i6,:].value/ligne6.max(),'+g')
    plt.plot(ttn, ligne7/ligne7.max(), 'b', label='wavelength 7. micron depth=0.62%')
    plt.plot(ttn, flux2d[i7,:].value/ligne7.max(),'+b')
    plt.legend()
    plt.xlabel('orbital_phase')
    plt.ylabel('normalised flux')
    if (SAVE):plt.savefig('plot_eclipse_'+pattern+'_normalised_lightcurve_6_7.png')


    ####  ECLIPSED DEPTH
    flux_min = flux2d.min(1)
    flux_max = flux2d.max(1)
    depth = (flux_max-flux_min)/flux_max*100
    #
    t = np.arange(nt)
    dans = np.where((t > 57) & (t < 80))
    dans = dans[0]
    dehors = np.where((t < 47) | (t > 47+44))
    dehors = dehors[0]
    deg=5
    depth = np.zeros(nw)
    for i in np.arange(0, nw):
        coeff = np.polyfit(ttn[dehors].value, flux2d[i,dehors].value, deg)
        yfit = np.poly1d(coeff)(ttn.value)
        diff = flux2d[i,:].value-yfit
        depth[i] = -diff[dans].mean()/coeff[-1]*100
    #
    plt.figure()
    plt.title(model+'measured depth polyfit method')
    plt.plot(wave1d, depth)
    plt.xlabel('wavelength micron')
    plt.ylabel('depth %')
    plt.plot([6,6], [0.3,0.568])
    plt.plot([7,7], [0.3,0.626])
    if (SAVE):
        plt.savefig('plot_measured_depth_'+pattern+'.png')
        write_depth_ecsv('table_measured_depth_'+pattern+'.dat', depth, wave1d.to(u.micron).value)
    return
