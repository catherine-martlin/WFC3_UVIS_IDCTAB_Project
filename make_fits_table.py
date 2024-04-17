#! /usr/bin/env python

"""Generates the IDCTAB. 

This script takes the output textfile of readwf3poly_onefilter and
creates a FITs table from it that is the final IDCTAB FITs table.

Author: 
	Catherine Martlin, October 2016, 2024.

Use: 

python make_fits_table.py --path='/grp/hst/wfc3v/martlin/idctab_creation/current_idctab_results/feb_06_2017_results/UVIS_F606W_CM_2017_02_06_idc.textfile' --output_filename='UVIS_F606W_CM_2017_02_10_idc.fits'


Outputs:

Dependencies:

Notes:
"""

import argparse
import numpy as np
import os
from astropy.io import fits
import csv
# ---------------------------------------------------------------------
def readfile(path, delim):
    
    lines =[]
    
    with open(path, 'rb') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            for line in row:
                lines.append(line)
    
    return lines
# ---------------------------------------------------------------------
# The main controller
# ---------------------------------------------------------------------

def make_fits_table_main(path, output_filename):
    """ The main controller of the make_idctab_main

        Parameters:

        Returns:

        Outputs:

        """
    #Read in the file: 
    idctab_all_vals = readfile(path, ' ')

    #Split the file into chip 1 and 2: 
    chip1 = idctab_all_vals[0:39]
    chip2 = idctab_all_vals[39:79]
    chip2.remove('')

    # Now set the various values as the type they need to be: 
    chip1[0] = np.int16(chip1[0])
    chip1[1] = np.str(chip1[1])
    chip1[2] = np.str(chip1[2])
    chip1[3] = np.int32(chip1[3])
    chip1[4] = np.int32(chip1[4])
    chip1[5:39] = np.float32(chip1[5:39])

    chip2[0] = np.int16(chip2[0])
    chip2[1] = np.str(chip2[1])
    chip2[2] = np.str(chip2[2])
    chip2[3] = np.int32(chip2[3])
    chip2[4] = np.int32(chip2[4])
    chip2[5:39] = np.float32(chip2[5:39])

    # Make all the arrays: 
    a1 = np.array([chip1[0],chip2[0]])
    a2 = np.array([chip1[1],chip2[1]])
    a3 = np.array([chip1[2],chip2[2]]) 
    a4 = np.array([chip1[3],chip2[3]])
    a5 = np.array([chip1[4],chip2[4]])
    a6 = np.array([chip1[5],chip2[5]])
    a7 = np.array([chip1[6],chip2[6]])
    a8 = np.array([chip1[7],chip2[7]])
    a9 = np.array([chip1[8],chip2[8]])
    a10 = np.array([chip1[9],chip2[9]])
    a11 = np.array([chip1[10],chip2[10]])
    a12 = np.array([chip1[11],chip2[11]])
    a13 = np.array([chip1[12],chip2[12]])
    a14 = np.array([chip1[13],chip2[13]])
    a15 = np.array([chip1[14],chip2[14]])
    a16 = np.array([chip1[15],chip2[15]])
    a17 = np.array([chip1[16],chip2[16]])
    a18 = np.array([chip1[17],chip2[17]])
    a19 = np.array([chip1[18],chip2[18]])
    a20 = np.array([chip1[19],chip2[19]])
    a21 = np.array([chip1[20],chip2[20]])
    a22 = np.array([chip1[21],chip2[21]])
    a23 = np.array([chip1[22],chip2[22]])
    a24 = np.array([chip1[23],chip2[23]])
    a25 = np.array([chip1[24],chip2[24]])
    a26 = np.array([chip1[25],chip2[25]])
    a27 = np.array([chip1[26],chip2[26]])
    a28 = np.array([chip1[27],chip2[27]])
    a29 = np.array([chip1[28],chip2[28]])
    a30 = np.array([chip1[29],chip2[29]])
    a31 = np.array([chip1[30],chip2[30]])
    a32 = np.array([chip1[31],chip2[31]])
    a33 = np.array([chip1[32],chip2[32]])
    a34 = np.array([chip1[33],chip2[33]])
    a35 = np.array([chip1[34],chip2[34]])
    a36 = np.array([chip1[35],chip2[35]])
    a37 = np.array([chip1[36],chip2[36]])
    a38 = np.array([chip1[37],chip2[37]])
    a39 = np.array([chip1[38],chip2[38]])

    # Make all the columns: 
    col1 = fits.Column(name = 'DETCHIP', format='1I', unit ='chip', null=-32767, disp='I4',array=a1)
    col2 = fits.Column(name = 'DIRECTION', format='8A', disp='A8',array=a2 )
    col3 = fits.Column(name = 'FILTER', format='8A', disp='A8',array=a3 )
    col4 = fits.Column(name = 'XSIZE', format='1J', disp='I8', null=-2147483647,array=a4 ) # Fixed 02/08/2017
    col5 = fits.Column(name = 'YSIZE', format='1J', disp='I8', null=-2147483647,array=a5 ) # Fixed 02/08/2017
    col6 = fits.Column(name = 'XREF', format='1E', disp='F10.2',array=a6 )
    col7 = fits.Column(name = 'YREF', format='1E', disp='F10.2',array=a7 )
    col8 = fits.Column(name = 'THETA', format='1E', disp='F12.4',array=a8 )
    col9 = fits.Column(name = 'SCALE', format='1E', disp='F12.4',unit ='arcsec',array=a9 )
    col10 = fits.Column(name = 'V2REF', format='1E', disp='F12.4',array=a10 )
    col11 = fits.Column(name = 'V3REF', format='1E', disp='F12.4',array=a11 )
    col12 = fits.Column(name = 'CX10', format='1E', disp='E20.7',array=a12 )
    col13 = fits.Column(name = 'CX11', format='1E', disp='E20.7',array=a13 )
    col14 = fits.Column(name = 'CX20', format='1E', disp='E20.7',array=a14 )
    col15 = fits.Column(name = 'CX21', format='1E', disp='E20.7',array=a15 )
    col16 = fits.Column(name = 'CX22', format='1E', disp='E20.7',array=a16 )
    col17 = fits.Column(name = 'CX30', format='1E', disp='E20.7',array=a17 )
    col18 = fits.Column(name = 'CX31', format='1E', disp='E20.7',array=a18 )
    col19 = fits.Column(name = 'CX32', format='1E', disp='E20.7',array=a19 )
    col20 = fits.Column(name = 'CX33', format='1E', disp='E20.7',array=a20 )
    col21 = fits.Column(name = 'CX40', format='1E', disp='E20.7',array=a21 )
    col22 = fits.Column(name = 'CX41', format='1E', disp='E20.7',array=a22 )
    col23 = fits.Column(name = 'CX42', format='1E', disp='E20.7',array=a23 )
    col24 = fits.Column(name = 'CX43', format='1E', disp='E20.7',array=a24 )
    col25 = fits.Column(name = 'CX44', format='1E', disp='E20.7',array=a25 )
    col26 = fits.Column(name = 'CY10', format='1E', disp='E20.7',array=a26 )
    col27 = fits.Column(name = 'CY11', format='1E', disp='E20.7',array=a27 )
    col28 = fits.Column(name = 'CY20', format='1E', disp='E20.7',array=a28 )
    col29 = fits.Column(name = 'CY21', format='1E', disp='E20.7',array=a29 )
    col30 = fits.Column(name = 'CY22', format='1E', disp='E20.7',array=a30 )
    col31 = fits.Column(name = 'CY30', format='1E', disp='E20.7',array=a31 )
    col32 = fits.Column(name = 'CY31', format='1E', disp='E20.7',array=a32 )
    col33 = fits.Column(name = 'CY32', format='1E', disp='E20.7',array=a33 )
    col34 = fits.Column(name = 'CY33', format='1E', disp='E20.7',array=a34 )
    col35 = fits.Column(name = 'CY40', format='1E', disp='E20.7',array=a35 )
    col36 = fits.Column(name = 'CY41', format='1E', disp='E20.7',array=a36 )
    col37 = fits.Column(name = 'CY42', format='1E', disp='E20.7',array=a37 )
    col38 = fits.Column(name = 'CY43', format='1E', disp='E20.7',array=a38 )
    col39 = fits.Column(name = 'CY44', format='1E', disp='E20.7',array=a39 )

    # Make what will become the fits file: 
    cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,
                     col18,col19,col20,col21,col22,col23,col24,col25,col26,col27,col28,col29,col30,col31,col32,
                     col33,col34,col35,col36,col37,col38,col39])

    # Make the fits file: 
    tbhdu = fits.BinTableHDU.from_columns(cols)
    temp_output = 'temp_fits_file_' + output_filename + '.fits'
    tbhdu.writeto(temp_output)

    # Now open the file to add all the comments/Primary Header: 
    my_test_image = fits.open(temp_output)

    # Add the comments to all the values: 
    my_test_image[1].header.comments['TTYPE1'] = 'label for field   1  '
    my_test_image[1].header.comments['TFORM1'] = 'data format of field: 2-byte INTEGER '
    my_test_image[1].header.comments['TUNIT1'] = 'physical unit of field '
    my_test_image[1].header.comments['TTYPE2'] = 'label for field   2  '
    my_test_image[1].header.comments['TFORM2'] = 'data format of field: ASCII Character '
    my_test_image[1].header.comments['TTYPE3'] = 'label for field   3  '
    my_test_image[1].header.comments['TFORM3'] = 'data format of field: ASCII Character '
    my_test_image[1].header.comments['TTYPE4'] = 'label for field   4  '
    my_test_image[1].header.comments['TFORM4'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE5'] = 'label for field   5  '
    my_test_image[1].header.comments['TFORM5'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE6'] = 'label for field   6  '
    my_test_image[1].header.comments['TFORM6'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE7'] = 'label for field   7  '
    my_test_image[1].header.comments['TFORM7'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE8'] = 'label for field   8  '
    my_test_image[1].header.comments['TFORM8'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE9'] = 'label for field   9  '
    my_test_image[1].header.comments['TFORM9'] = 'data format of field: 2-byte INTEGER '
    my_test_image[1].header.comments['TUNIT9'] = 'physical unit of field '
    my_test_image[1].header.comments['TTYPE10'] = 'label for field   10  '
    my_test_image[1].header.comments['TFORM10'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE11'] = 'label for field   11  '
    my_test_image[1].header.comments['TFORM11'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE12'] = 'label for field   12  '
    my_test_image[1].header.comments['TFORM12'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE13'] = 'label for field   13  '
    my_test_image[1].header.comments['TFORM13'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE14'] = 'label for field   14  '
    my_test_image[1].header.comments['TFORM14'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE15'] = 'label for field   15  '
    my_test_image[1].header.comments['TFORM15'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE16'] = 'label for field   16  '
    my_test_image[1].header.comments['TFORM16'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE17'] = 'label for field   17  '
    my_test_image[1].header.comments['TFORM17'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE18'] = 'label for field   18  '
    my_test_image[1].header.comments['TFORM18'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE19'] = 'label for field   19  '
    my_test_image[1].header.comments['TFORM19'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE20'] = 'label for field   20  '
    my_test_image[1].header.comments['TFORM20'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE21'] = 'label for field   21  '
    my_test_image[1].header.comments['TFORM21'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE22'] = 'label for field   22  '
    my_test_image[1].header.comments['TFORM22'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE23'] = 'label for field   23  '
    my_test_image[1].header.comments['TFORM23'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE24'] = 'label for field   24  '
    my_test_image[1].header.comments['TFORM24'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE25'] = 'label for field   25  '
    my_test_image[1].header.comments['TFORM25'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE26'] = 'label for field   26  '
    my_test_image[1].header.comments['TFORM26'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE27'] = 'label for field   27  '
    my_test_image[1].header.comments['TFORM27'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE28'] = 'label for field   28  '
    my_test_image[1].header.comments['TFORM28'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE29'] = 'label for field   29  '
    my_test_image[1].header.comments['TFORM29'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE30'] = 'label for field   30  '
    my_test_image[1].header.comments['TFORM30'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE31'] = 'label for field   31  '
    my_test_image[1].header.comments['TFORM31'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE32'] = 'label for field   32  '
    my_test_image[1].header.comments['TFORM32'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE33'] = 'label for field   33  '
    my_test_image[1].header.comments['TFORM33'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE34'] = 'label for field   34  '
    my_test_image[1].header.comments['TFORM34'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE35'] = 'label for field   35  '
    my_test_image[1].header.comments['TFORM35'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE36'] = 'label for field   36  '
    my_test_image[1].header.comments['TFORM36'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE37'] = 'label for field   37  '
    my_test_image[1].header.comments['TFORM37'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE38'] = 'label for field   38  '
    my_test_image[1].header.comments['TFORM38'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TTYPE39'] = 'label for field   39  '
    my_test_image[1].header.comments['TFORM39'] = 'data format of field: 4-byte INTEGER '
    my_test_image[1].header.comments['TDISP1'] = 'display format  '
    my_test_image[1].header.comments['TDISP2'] = 'display format  '
    my_test_image[1].header.comments['TDISP3'] = 'display format  '
    my_test_image[1].header.comments['TDISP4'] = 'display format  '
    my_test_image[1].header.comments['TDISP5'] = 'display format  '
    my_test_image[1].header.comments['TDISP6'] = 'display format  '
    my_test_image[1].header.comments['TDISP7'] = 'display format  '
    my_test_image[1].header.comments['TDISP8'] = 'display format  '
    my_test_image[1].header.comments['TDISP9'] = 'display format  '
    my_test_image[1].header.comments['TDISP10'] = 'display format  '
    my_test_image[1].header.comments['TDISP11'] = 'display format  '
    my_test_image[1].header.comments['TDISP12'] = 'display format  '
    my_test_image[1].header.comments['TDISP13'] = 'display format  '
    my_test_image[1].header.comments['TDISP14'] = 'display format  '
    my_test_image[1].header.comments['TDISP15'] = 'display format  '
    my_test_image[1].header.comments['TDISP16'] = 'display format  '
    my_test_image[1].header.comments['TDISP17'] = 'display format  '
    my_test_image[1].header.comments['TDISP18'] = 'display format  '
    my_test_image[1].header.comments['TDISP19'] = 'display format  '
    my_test_image[1].header.comments['TDISP20'] = 'display format  '
    my_test_image[1].header.comments['TDISP21'] = 'display format  '
    my_test_image[1].header.comments['TDISP22'] = 'display format  '
    my_test_image[1].header.comments['TDISP23'] = 'display format  '
    my_test_image[1].header.comments['TDISP24'] = 'display format  '
    my_test_image[1].header.comments['TDISP25'] = 'display format  '
    my_test_image[1].header.comments['TDISP26'] = 'display format  '
    my_test_image[1].header.comments['TDISP27'] = 'display format  '
    my_test_image[1].header.comments['TDISP28'] = 'display format  '
    my_test_image[1].header.comments['TDISP29'] = 'display format  '
    my_test_image[1].header.comments['TDISP30'] = 'display format  '
    my_test_image[1].header.comments['TDISP31'] = 'display format  '
    my_test_image[1].header.comments['TDISP32'] = 'display format  '
    my_test_image[1].header.comments['TDISP33'] = 'display format  '
    my_test_image[1].header.comments['TDISP34'] = 'display format  '
    my_test_image[1].header.comments['TDISP35'] = 'display format  '
    my_test_image[1].header.comments['TDISP36'] = 'display format  '
    my_test_image[1].header.comments['TDISP37'] = 'display format  '
    my_test_image[1].header.comments['TDISP38'] = 'display format  '
    my_test_image[1].header.comments['TDISP39'] = 'display format  '
    my_test_image[1].header.comments['TNULL1'] = 'undefined value for column '
    my_test_image[1].header.comments['TNULL4'] = 'undefined value for column '
    my_test_image[1].header.comments['TNULL5'] = 'undefined value for column '

    # Add header information: 
    prihdr = fits.Header()
    prihdr['NORDER'] = 4
    prihdr['TELESCOP'] = 'HST ' # Added 02/08/2017
    prihdr['INSTRUME'] = 'WFC3 '
    prihdr['DETECTOR'] = 'UVIS '
    prihdr['FILETYPE'] = 'Distortion Coefficients'
    prihdr['COMMENT'] = 'Created by C. Martlin, February, 2017'  
    prihdr['DESCRIP'] = 'Geometric Distortion Coefficients'
    prihdu = fits.PrimaryHDU(header=prihdr)

    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist[0]._bitpix=16 # Added 02/10/2017
    thdulist.writeto(output_filename)

# -----------------------------------------------------------------------------
# For command line execution
# -----------------------------------------------------------------------------

def parse_args():
    """Parses command line arguments.

    Parameters:
        nothing

    Returns:
        args : argparse.Namespace object
            An argparse object containing all of the added arguments.

    Outputs:
        nothing
    """

    #Create help string:
    path_help = 'List of path to the textfile that needs to become the FITs table.'
    output_filename_help = 'The name of the final FITs file.'

    # Add arguments:
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', '-p', dest = 'path', action = 'store',
                        type = str, required = True, help = path_help)
    parser.add_argument('--output_filename', '-ofile', dest = 'output_filename', action = 'store',
                        type = str, required = True, help = output_filename_help)
    
    # Parse args:
    args = parser.parse_args()

    return args
 # -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_args()

    make_fits_table_main(args.path, args.output_filename)
