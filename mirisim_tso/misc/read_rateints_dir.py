# readwrite_rateints_dir : read write  cube of slopes  from jwst official pipeline
## this version does not use the keywords nints, intstart, intend
###  defensive coding !
###  date 28 jui 2022 author R Gastaud
##  from outils_bis

from astropy.io import fits
import astropy.units as u
import numpy as np
import glob
import os
#
#

def read_rateints_dir(input_dir, pattern='*_rateints.fits', verbose=False, out_file=None):
    files = sorted(glob.glob(os.path.join(input_dir, pattern)))
    nn = len(files)
    print('number of files : {} '.format(nn))
    if( nn ==0):return (0,0)
    hdr0 = fits.getheader(files[0])
    sci = fits.getdata(files[0])
    nz, ny, nx = sci.shape
    nints0 = hdr0['nints']

    ###########  bad definitoin of nints in merge_sim_files.py
    # nints is the total number of integration, not the number of integration in this file
    nints = 0
    my_ends = []
    for input_file in files:
        header1 = fits.getheader(input_file, 1)
        naxis3 = header1['naxis3']
        nints = nints+naxis3
        my_ends.append(nints)
    print('nints', nints0, nints)
    print(my_ends)
    ###  end of the patch
    #
    cube = np.zeros([nints, ny, nx])
    if(verbose):
        print(nints, ny, nx)
        print('cube.shape', cube.shape, 'sci.shape', sci.shape)
    cube[:]=np.nan
    i = 0
    for input_file in files:
        hdul = fits.open(input_file)
        """
        nints = hdul[0].header['nints']
        int_start = hdul[0].header['INTSTART']
        int_end = hdul[0].header['INTEND']
        """
        sci = hdul['sci'].data
        nz, ny, nx = sci.shape
        #if(verbose):print(input_file, int_start, int_end, int_end-int_start, nints)
        ii = np.where(sci == np.nan)
        if (len(ii[0]) > 0):
            print('warning nan in sci')
        cube[my_ends[i]-nz:my_ends[i],:,:] = sci
        hdul.close()
        i = i+1
        if(verbose):print(input_file, "read")
    #
    #
    return cube

