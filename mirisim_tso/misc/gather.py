# GATHER
# gathers all  the files which contains only one integration per file of 2 GB
# Author Rene Gastaud, rene.gastaud@cea.fr
#    11 May 2021

from astropy.io import fits
import glob
import os
import math
import numpy as np
import sys
#
# pattern='/Users/gastaud/miri/coulais/data_challenge/DataChallenge2102_noisy/simulation_*det_images_noise.fits'

def gather(pattern):
    ###  hard-coded  ###
    file_size_max=2048 # megabytes
    ##file_size_max=100
    tframe = 0.15904
    t_file = 65*tframe/3600/24   # NGROUPS = 65  Number of groups in integration
    # time in day, because EXPSTART, EXPEND, EXPMI are in MJD
    out_file = 'big_file_'
    
    ##
    files = sorted(glob.glob(pattern))
    nints = len(files)

    #######  measure total size
    total_size = 0
    for file in files:
        total_size += os.stat(file).st_size/1024/1024 # megabytes
    mean_size = total_size/len(files)

    ####### compute  exsegtot

    exsegtot = math.ceil(total_size/file_size_max)
    number_integrations_per_file = math.floor(file_size_max/mean_size)

    ##############################
    last_file = os.path.basename(files[-1])
    current_size=0
    exsegnum = 1
    k = 0
    intstart = k+1
    sci = []
    refout = []

    for file in files:
        hdul = fits.open(file)
        sci.append(hdul['SCI'].data)
        refout.append(hdul['REFOUT'].data)
        current_file = os.path.basename(file)
        current_size += os.stat(file).st_size/1024/1024
        k = k +1
        if (current_size >= (file_size_max-mean_size)) or (current_file == last_file):
                print(k, exsegnum, current_size, current_file)
                intend = k
                ## write hdu
                hdul_out = hdul.copy()
                sci = np.array(sci).squeeze()
                refout = np.array(refout).squeeze()
                hdul_out['SCI'].data = sci
                hdul_out['REFOUT'].data = refout
                #
                hdul_out[0].header['EXSEGNUM'] = exsegnum
                hdul_out[0].header['EXSEGTOT'] = exsegtot
                hdul_out[0].header['NINTS'] = nints
                hdul_out[0].header['INTSTART'] = intstart
                hdul_out[0].header['INTEND'] = intend
                #
                hdul_out[0].header['EXPSTART'] = intstart*t_file
                hdul_out[0].header['EXPEND']   = intend*t_file
                hdul_out[0].header['EXPMID']   = (intstart+intend)*t_file/2
                hdul_out.writeto(out_file+str(exsegnum)+'.fits')
                #  reset
                intstart = intend+1
                current_size=0
                exsegnum = exsegnum +1
                sci = []
                refout = []
    return
#
if __name__ == '__main__':
    print(sys.argv)
    if len(sys.argv) == 2:
      pattern = sys.argv[1]
      gather(pattern)
    else:
      print('syntax: gather pattern')


         
            
