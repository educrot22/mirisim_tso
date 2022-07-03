# read_uncal_dir : read uncal  from jwst official pipeline
# and compute the slope with simple algorithm
##  from outils_bis
##   July, 2nd 2022    René Gastaud, CEA-Saclay
##  BEWARE RETURN IN DN/S

from astropy.io import fits
import astropy.units as u
import numpy as np
import glob
import os
#
#

def read_uncal_dir(input_dir, pattern='*_uncal.fits', verbose=False, out_file=None):
    files = sorted(glob.glob(os.path.join(input_dir, pattern)))
    nn = len(files)
    print('number of files : {} '.format(nn))
    if( nn ==0):return (0,0)
    if verbose:
        print('((before last frame) - (first frame))/tframe/(ngroups-2)')
        print(' version 1')
        
    header0 = fits.getheader(files[0])
    sci = fits.getdata(files[0])
    nt, nz, ny, nx = sci.shape
    nints = header0['nints']
    tint = header0['EFFINTTM']
    tframe= header0['TFRAME']
    ngroups =  header0['ngroups']
    tint2 = tframe*(ngroups-2)
    keys=['TARGPROP', 'TARGNAME', 'NINTS', 'INTSTART', 'INTEND', 'NGROUPS', 'TFRAME', 'TGROUP', 'EFFINTTM', 'DURATION']
    for key in keys:print(key, header0[key])
    cube = np.zeros([nints, ny, nx])
    if(verbose):
        print(nints, ny, nx)
        print('cube.shape', cube.shape, 'sci.shape', sci.shape)
    cube[:]=np.nan
    for input_file in files:
        hdul = fits.open(input_file)
        nints = hdul[0].header['nints']
        int_start = hdul[0].header['INTSTART']
        int_end = hdul[0].header['INTEND']
        sci = hdul['sci'].data
        if(verbose):print(input_file, int_start, int_end, int_end-int_start, nints)
        ii = np.where(sci == np.nan)
        if (len(ii[0]) > 0):
            print('warning nan in sci')
        # beware sci unsigned integer !
        slopes = sci[:,-2,:,:]/tint2
        slopes = slopes - sci[:,0,:,:]/tint2
        cube[int_start-1:int_end,:,:] = slopes
        hdul.close()
        if(verbose):print(" ")
    #
    if (out_file is not None):
        hdu = fits.PrimaryHDU(cube)
        hdu.header = header0
        hdu.header['in_dir'] = input_dir
        hdu.header['pattern'] = pattern
        hdu.header['bunit'] = 'DN/second'
        hdu.header['history'] = 'read_uncal version 1'
        hdu.header['history'] = '[(before last frame) - (first frame)]/tframe/(ngroups-2)'
        hdu.header['tint2'] = tint2
        hdu.writeto(out_file)
    #
    return cube

