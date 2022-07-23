"""
Given MIRISim simulations for TSO, will merge integrations into bigger files. Ultimately, files will be split into
smaller files to make up for JWST pipeline rules w.r.t. TSOs

source: https://jwst-pipeline.readthedocs.io/en/latest/jwst/data_products/file_naming.html#segmented-files

Authors: Christophe Cossou, christophe.coussou@cea.fr, Rene Gastaud

History :
20 July RG intstart, intend, nints corrected, and add overwrite option
4 March 2022 correct header do not use the first header for all the output files !
16 March 2022 starts counting from zero for the output files
create hdri
"""
import glob
import os
import numpy as np
from astropy.io import fits
import math

# I don't use mirage.utils.constants.FILE_SPLITING_LIMIT) which gives files way too small.
# this is for Nircam which has 10 detectors

#FILE_SPLITTING_LIMIT = 200E6 # 90167050 # 45083520*2+10  for testing
FILE_SPLITTING_LIMIT =  2000*1E6 # 2 Giga Bytes

ma_version='deuxieme'

def merge_files(input_dir, output_dir, pattern='simulation*.fits', verbose=False, overwrite=False):
    """
    We assume that all files have the same characteristics (in particular same number of integrations, one integration)

    :param str input_dir: the input directory
    :param str output_dir: the output directory
    :param str pattern: the name of each file
    :return:
    """

    input_files = sorted(glob.glob(os.path.join(input_dir, pattern)))
    nb_input_files = len(input_files)

    if not input_files:
        raise ValueError(f"No files found with input_dir: '{input_dir}' and pattern: '{pattern}'")

    # Get parameters from first file
    file_name = input_files[0]
    file_size = os.path.getsize(file_name)
    #  total number of integrations
    nint = int(math.floor(FILE_SPLITTING_LIMIT/file_size))
    #RG 16 March 2022  just a trick for our case, to be removed
    ## replace floor by ceil
    ## nint = int(math.ceil(FILE_SPLITTING_LIMIT/file_size))
    if((nb_input_files%nint) ==1):
        print('warning last segment one file')
        one_flag = True
    else:
        one_flag = False
        
    # total number of output files  EXSEGTOT
    nb_segments = math.ceil(nb_input_files/nint)
    
    ## read the first file for the headers and the shape
    hdulist = fits.open(file_name)
    if(verbose): hdulist.info()
    hdr0 = hdulist[0].header
    sci = hdulist['SCI'].data
    hdr_sci = hdulist['SCI'].header
    pixel_dq = hdulist['PIXELDQ'].data
    hdr_pixel_dq = hdulist['PIXELDQ'].header
    refout = hdulist['REFOUT'].data
    hdr_refout = hdulist['REFOUT'].header
    asdf = hdulist['ASDF'].data
    hdr_asdf = hdulist['ASDF'].header

    # create  cube
    nt, nz, ny, nx = sci.shape
    cube = np.zeros([nint, nz, ny, nx], dtype='float32')
    #
    nt, nz, ny, nx = refout.shape
    cube_ref = np.zeros([nint, nz, ny, nx], dtype='float32')
    #
    # read cube
    i = 0             ## count the number of read files  (small files), reset to zero each nint <==> modulo nints
    i_big_file = 0    ## count the number of written files  (big files)  RG 16 March 2022
    i_small_file = 0  ## count the number of read files  (small files)
    for filename in input_files:
        if(verbose):print('reading',i, filename)
        if (i ==0): hdri = fits.getheader(filename)  # 4 March 2022
        cube[i, :, :, :]     = fits.getdata(filename, extname='SCI')
        cube_ref[i, :, :, :] = fits.getdata(filename, extname='REFOUT')
        i = i + 1
        i_small_file = i_small_file + 1
        if (i == nint):
            i=0
            if(verbose):print('write cube', i_big_file)
            if ((i_big_file == (nb_segments-2)) and one_flag):
                print('special processing add an extra integration')
            i_start = i_small_file - nint + 1
            i_end   = i_small_file # + 1  RG 11 jul 2022
            output_filename = fill_header( hdri,  nb_segments, nb_input_files, i_start, i_end, i_big_file)
            output_filename = os.path.join(output_dir, output_filename)
            write_detector_model(output_filename, hdri,  cube, hdr_sci, pixel_dq, hdr_pixel_dq, cube_ref, hdr_refout, asdf, hdr_asdf, overwrite=overwrite )
            i_big_file = i_big_file+1
    if(i>0):
        if(verbose):print('write last big file', i_big_file, i)
        i_start = i_small_file - i + 1
        i_end   = i_small_file # + 1  RG 11 jul 2022
        output_filename = fill_header( hdri,  nb_segments, nb_input_files, i_start, i_end, i_big_file)
        output_filename = os.path.join(output_dir, output_filename)
        cube = cube[0:i, :, :, :]
        cube_ref = cube_ref[0:i, :, :, :]
        write_detector_model(output_filename, hdri,  cube, hdr_sci, pixel_dq, hdr_pixel_dq, cube_ref, hdr_refout, asdf, hdr_asdf, overwrite=overwrite )

        

  
def write_detector_model(filename, hdr0,  cube, hdr_sci, pixel_dq, hdr_pixel_dq, cube_ref, hdr_refout, asdf, hdr_asdf, overwrite ):
    '''
    Write a cube of data, cube of error, cube of flags
    :param filename: str, name of the output file
    :param hdr0: str, header
    :param cube: ndarray: dim: (nt, nz, ny, nx)
    :param hdr_sci: str, header
    :param pixel_dq: ndarray: dim: (nt, nz, ny, nx)
    :param hdr_pixel_dq: str, header
    :param cube_ref: ndarray: dim: (nt, nz, ny, nx)
    :param hdr_refout: str, header
    :param asdf: binary table
    :param hdr_asdf: str, header
    '''
    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header = hdr0

    image_hdu1 = fits.ImageHDU(cube)
    image_hdu1.header = hdr_sci

    image_hdu2 = fits.ImageHDU(pixel_dq)
    image_hdu2.header = hdr_pixel_dq
    
    image_hdu3 = fits.ImageHDU(cube_ref)
    image_hdu3.header = hdr_refout

    image_hdu4 = fits.BinTableHDU(asdf)
    image_hdu4.header = hdr_asdf
    
    hdul = fits.HDUList([primary_hdu, image_hdu1, image_hdu2, image_hdu3, image_hdu4])
    #hdul = fits.HDUList([primary_hdu])
    hdul.writeto(filename, overwrite=overwrite)
    return

##merge_files(XX, XX, pattern='simulation*.fits', verbose=False)


def fill_header( hdri,  nb_segments, nb_input_files, i_start, i_end, i_big_file):
    """
    Fill the header

    :param ?? hdri: header
    :param int nb_segments: Number of segments, id est big files
    :param int nb_input_files: Number of input files (small)
    :param int i_start: start index for small input files. i_small_file-nint+1
    :param int i_end: end index for small input files. i_small_file+1
    :param int i_big_file : end index for small input files. i_small_file+1
    :return: name of the output file
    :rtype: str
    """

    hdri["TSOVISIT"] = True
    hdri["EXSEGTOT"] = (nb_segments, "The total number of segments")
    hdri["DURATION"] = (nb_input_files*hdri["EFFINTTM"], '[s] Effective exposure time')
    hdri["EFFEXPTM"] = (nb_input_files*hdri["EFFINTTM"], '[s] Total duration of exposure')
    #
    hdri["INTSTART"] = (i_start, "The starting integration number of the data in this segment, zero-counting")
    hdri["INTEND"]   = (i_end, "The ending integration number of the data in this segment")
    hdri["EXSEGNUM"] = (i_big_file, "The segment number of the current product")
    hdri["NINTS"] = nb_input_files # i_end - i_start
    output_filename = "miri_big_cube_{:03d}.fits".format(i_big_file)
    return output_filename


