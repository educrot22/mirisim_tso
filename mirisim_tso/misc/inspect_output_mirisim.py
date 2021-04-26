# IPython log file

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u

from mirisim_tso.misc.extract_tso_spectra import extract_tso_spectra
from read_illum_models_dir import read_illum_models_dir

in_dir='data/'
bckg_file = 'data/illum_model_background.fits'
bckg = fits.getdata(bckg_file )*u.electron/u.second
plt.figure()
plt.imshow(bckg.value)
#
flux3d0,  wave2d = read_illum_models_dir(in_dir, pattern="SED*.fits")
ny, nx, nt = flux3d0.shape
print(flux3d0.shape)
flux3d  = flux3d0 - bckg.reshape([ny, nx, 1])
i0=140
i1=375
wave1d = wave2d[i0:i1, nx//2]
cubette = flux3d[i0:i1,:, : ]
print(cubette.shape)
flux2d = cubette.sum(axis=1)


nw, nt = flux2d.shape
i6 = np.argmin(np.abs(wave1d.value-6))
i7 = np.argmin(np.abs(wave1d.value-7))
i10 = np.argmin(np.abs(wave1d.value-10))
i11 = np.argmin(np.abs(wave1d.value-11))
print(nw, nt)
ttn= np.arange(nt)*100.
plt.figure()
plt.title('spectrum')
plt.plot(wave1d, flux2d[:,0])
plt.xlabel('wavelength micron')
plt.ylabel('flux '+str(flux2d.unit))

plt.figure()
plt.plot(ttn, flux2d[i6,:], label='wavelength 6. micron')
plt.plot(ttn, flux2d[i7,:], label='wavelength 7. micron')
plt.plot(ttn, flux2d[i10,:], label='wavelength 10. micron')
plt.plot(ttn, flux2d[i11,:], label='wavelength 11. micron')
plt.legend()
plt.xlabel('time in second')
plt.ylabel(' flux'+str(flux2d.unit))
plt.title('wasp43-b light curve, output of mirisim')

plt.figure()
plt.plot(ttn, flux2d[i6,:]/flux2d[i6,:].max(), label='wavelength 6. micron')
plt.plot(ttn, flux2d[i7,:]/flux2d[i7,:].max(), label='wavelength 7. micron')
plt.plot(ttn, flux2d[i10,:]/flux2d[i10,:].max(), label='wavelength 10. micron')
plt.plot(ttn, flux2d[i11,:]/flux2d[i11,:].max(), label='wavelength 11. micron')
plt.legend()
plt.xlabel('time in second')
plt.ylabel('normalised flux')
plt.title('wasp43-b normalised light curve, output of mirisim')
plt.savefig('plot_normalised_lightcurve6711_mirisim.png')

flux_min = flux2d.min(1)
flux_min.shape
flux_max = flux2d.max(1)
depth = (flux_max-flux_min)/flux_max*100
plt.figure()
plt.title('mirisim output, measured depth minmax method')
plt.plot(wave1d, depth)
plt.xlabel('wavelength micron')
plt.ylabel('depth %')
plt.plot([6,6], [0.3,0.55])
plt.plot([7,7], [0.3,0.624])
plt.savefig('plot_mirisim_measured_depth.png')


from mirisim_tso.misc.compute_theoric_eclipse_depth import compute_theoric_eclipse_depth
theoric_depth = compute_theoric_eclipse_depth(wave1d)
plt.figure()
plt.title('mirisim output measured depth minmax method')
plt.plot(wave1d, depth, label='measured')
plt.plot(wave1d, theoric_depth*100, label='theoric')
plt.legend()
plt.xlabel('wavelength micron')
plt.ylabel('depth %')
plt.plot([6,6], [0.3,0.55])
plt.plot([7,7], [0.3,0.624])
plt.savefig('plot_mirisim_measured_theoric_depth.png')
