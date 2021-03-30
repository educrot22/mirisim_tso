# IPython log file

from astropy import units as u
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table, Column
import datetime

# format='fits'
def write_sed(fileName, wavelength, flux, mask=None, meta=None, format='ascii.ecsv'):
    col1 = Column(wavelength.value, name='wavelength', unit=str(wavelength.unit), dtype='float32')
    col2 = Column(flux.value, name='flux', unit=str(flux.unit), dtype='float32')
    t = Table([col1, col2])
    if mask is not None:
        col3 = Column(mask, name='mask', dtype='unit32')
        t = Table([col1, col2, col3])
    if meta is not None:
        t.meta=meta
    t.meta['DATE'] = str(datetime.datetime.now())
    #
    t.write(fileName, format=format)
    return
