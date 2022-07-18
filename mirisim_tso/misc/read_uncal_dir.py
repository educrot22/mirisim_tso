# read_uncal_dir : read uncal  from jwst official pipeline
# and compute the slope with simple algorithm
##  from outils_bis
##   July, 2nd 2022    René Gastaud, CEA-Saclay
###  July, 4th  defensive coding in print key
###  July, 4th  defensive coding remove int_start and int_end, see read_rateints_dir
##
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
    nints0 = header0['nints']
    tint = header0['EFFINTTM']
    tframe= header0['TFRAME']
    ngroups =  header0['ngroups']
    tint2 = tframe*(ngroups-2)
    keys=['TARGPROP', 'TARGNAME', 'NINTS', 'INTSTART', 'INTEND', 'NGROUPS', 'TFRAME', 'TGROUP', 'EFFINTTM', 'DURATION', 'TINT']
    #
    ###########  bad definition of nints in merge_sim_files.py
    # nints is the total number of integration, not the number of integration in this file
    nints = 0
    my_ends = []
    for input_file in files:
        header1 = fits.getheader(input_file, 1)
        naxis3 = header1['naxis4']
        nints = nints+naxis3
        my_ends.append(nints)

    cube = np.zeros([nints, ny, nx])
    if(verbose):
        for key in keys:
            if key in header0:
                print(key, header0[key])
            else:
                print(key, 'not found')
        print(nints, ny, nx)
        print('cube.shape', cube.shape, 'sci.shape', sci.shape)
    cube[:]=np.nan
    i = 0
    for input_file in files:
        hdul = fits.open(input_file)
        #nints = hdul[0].header['nints']
        #int_start = hdul[0].header['INTSTART']
        #int_end = hdul[0].header['INTEND']
        sci = hdul['sci'].data
        nt,nz, ny, nx = sci.shape
        if(verbose):print(input_file, nt, my_ends[i]-nt)
        ii = np.where(sci == np.nan)
        if (len(ii[0]) > 0):
            print('warning nan in sci')
        # beware sci unsigned integer !
        slopes = sci[:,-2,:,:]/tint2
        slopes = slopes - sci[:,0,:,:]/tint2
        cube[my_ends[i]-nt:my_ends[i],:,:] = slopes
        hdul.close()
        i = i+1
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

