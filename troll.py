#! /usr/bin/env python
'''
ABOUT:


--filename_list - A list of the paths to the files in question.

--path_to_data (Optional)- The directory path to the files.

-- outfile - Name of file the products are printed to.

-- extension - UVIS chip to run on. Either 1 or 4 for UVIS 2 or 1.

EXAMPLE RUN:

python troll.py --filename_list --outfile_date --extension

python troll.py -fl=all_filenames.txt -out_date=mar_31_2016 -ext=4
'''

#import all the python things here!
import math
from astropy.io import fits
import os
import argparse
import numpy as np
import sys

# -----------------------------------------------------------------------------
def upload_header_info(flt_file, extension, path_to_data):
    '''
    We want to get these values from the header instead of looking
    them up and plugging them in when we run the code.

    dec0 - The DEC_TARG or delcination of the
            target that is found in the header.

    roll0 - The PA_V3 value, or the angle between
            North and the V3 vector, that is found in
            the header.

    '''

    flt_path = path_to_data + flt_file

    hdu = fits.open(flt_path)
    dec0 = hdu[0].header['DEC_TARG']
    roll0 = hdu[0].header['PA_V3']
    ra0 = hdu[0].header['RA_TARG']
    if extension == 4:
        orientat_header = hdu[4].header['ORIENTAT']
        crval1 = hdu[4].header['crval1']
        crval2 = hdu[4].header['crval2']
        v_factor0 = hdu[4].header['VAFACTOR']
    elif extension == 1:
        orientat_header = hdu[1].header['ORIENTAT']
        crval1 = hdu[1].header['crval1']
        crval2 = hdu[1].header['crval2']
        v_factor0 = hdu[1].header['VAFACTOR']
    hdu.close()

    return dec0,roll0,ra0,crval1,crval2,orientat_header,v_factor0
# -----------------------------------------------------------------------------
def convert_to_rads(item):
    '''
    The trigonometric functions take in radians so here we convert
    from arcseconds into radians.
    '''

    item_degree = item * (1./3600.)
    item_rad = math.radians(item_degree)

    return item_rad
# -----------------------------------------------------------------------------
def read_in_angle_stat(path, path_to_angle, extension):
    '''
    Read in the values from angle.stat and match to the file being
    calculated and provide the calculated orientat and the EPS for
    that file. Make sure the angle.stat file is correct for whichever
    UVIS chip you are using.
    '''
    # Cut out the unique name of the file from the path:

    length = len(path)
    mid = length - 11 - 7
    end = length - 7
    end2 = length - 9
    current_I = path[:end2]

    name_list = []
    orientat_calc_list = []
    EPS_list = []
    angle_file = path_to_angle + "angle.stat"

    global orientat_calc
    global EPS

    # Read in all the angle.stat info and match to the unique name.
    if extension == 1:
        with open(angle_file, 'r') as angle_stats_1:
            for line in angle_stats_1:
                curr_name = line[0:9]
                extension_val = line[10:11]
                curr_orientat = line[21:32]
                curr_EPS = line[38:47]
                name_list.append(curr_name)
                orientat_calc_list.append(curr_orientat)
                EPS_list.append(curr_EPS)

            for x in xrange(len(name_list)):
                if name_list[x] == current_I :
                    orientat_calc = orientat_calc_list[x]
                    EPS = EPS_list[x]

    elif extension == 4:
        with open(angle_file, 'r') as angle_stats_4:
            for line in angle_stats_4:
                curr_name = line[0:9]
                extension_val = line[10:11]
                curr_orientat = line[21:32]
                curr_EPS = line[38:47]
                name_list.append(curr_name)
                orientat_calc_list.append(curr_orientat)
                EPS_list.append(curr_EPS)

            for x in xrange(len(name_list)):
                if name_list[x] == current_I :
                    orientat_calc = orientat_calc_list[x]
                    EPS = EPS_list[x]


    return orientat_calc, EPS
# -----------------------------------------------------------------------------
# The main controller
# -----------------------------------------------------------------------------

def troll_main(filename_list, outfile_date, extension, path_to_data,path_to_coeffs,screen_outputfile):
    '''
    The main controller/conversion tool.
    '''
    orig_stdout = sys.stdout
    out_f = file(screen_outputfile, 'a')
    sys.stdout = out_f

    with open(filename_list, 'r') as f:
        for path in f:
            path = path.rstrip()
            # Get the orientat_calc and EPS from the angle.stat:
            path_to_angle = path_to_coeffs
            #path_to_angle = '/grp/hst/wfc3o/martlin/idctab_vera/make_idctab_codes/'
            orientat_calc, EPS = read_in_angle_stat(path, path_to_angle, extension)

            # Read in the files and the values
            dec0, roll0, ra0, crval1, crval2, orientat_header, v_factor0 = upload_header_info(path, extension, path_to_data)

            # Convert the stuff:
            pi = math.pi
            v2 = -27.4596 # New one from 2012
            v3 = -33.2604
            v2_rad = convert_to_rads(v2)
            v3_rad = convert_to_rads(v3)

            # Now work out the roll:
            v2_sin = math.sin(v2_rad)**2
            v3_sin = math.sin(v3_rad)**2
            v2_sin1 = math.sin(v2_rad)
            v3_sin1 = math.sin(v3_rad)
            sin_rho_sqr = v2_sin + v3_sin - v2_sin*v3_sin

            if sin_rho_sqr > 0.0:
                sin_rho = math.sqrt(sin_rho_sqr)
                cos_rho = math.sqrt(1.0 - sin_rho_sqr)
                if v2 > 0.0:
                    beta = math.asin(v3_sin1/sin_rho)
                else:
                    beta = pi - math.asin(v3_sin1/sin_rho)
                if v3 > 0.0:
                    gamma = math.asin(v2_sin1/sin_rho)
                else:
                    gamma = pi - math.asin(v2_sin1/sin_rho)

                angle_a = pi/2.0 + math.radians(roll0) - beta
                sin_a = math.sin(angle_a)
                cos_a = math.cos(angle_a)
                dec0_rad = math.radians(dec0)
                cos_dec0 = math.cos(dec0_rad)
                sin_dec0 = math.sin(dec0_rad)

                angle_b = math.atan2(sin_a * cos_dec0, sin_dec0 *sin_rho - cos_dec0*cos_rho*cos_a)

                roll_rad = pi - (gamma + angle_b)
                if roll_rad < 0.0:
                    roll_rad = roll_rad + 2.0*pi # To keep roll between +/- pi
                roll_deg = math.degrees(roll_rad)
            else:
                roll_deg = roll0


            # Then print all the things to a file?
            outfile = 'troll_output_{}_{}.txt'.format(extension,outfile_date)
            f = open(outfile, 'a')
            filename = path + '[' + str(extension) + ']'
            ra0_prec = '{:}'.format(np.float64(ra0)) #{:.6e}
            dec0_prec = '{:}'.format(np.float64(dec0)) #{:.6e}
            crval1_prec = '{:}'.format(np.float64(crval1)) #{:.6e}
            crval2_prec = '{:}'.format(np.float64(crval2)) #{:.6e}
            orientat_header_prec = '{:}'.format(np.float64(orientat_header)) #{:.6e}
            orientat_calc_prec = '{:}'.format(np.float64(orientat_calc)) #{:.6e}
            EPS_prec = '{:}'.format(np.float64(EPS)) #{:.6e}
            roll0_prec = '{:}'.format(np.float64(roll0)) #{:.8e}
            roll_deg_prec = '{:}'.format(np.float64(roll_deg)) #{:.8e}
            v_factor0_prec = '{:}'.format(np.float64(v_factor0)) #{:.6e}
            allangle = filename + ',' + str(ra0) + ',' + str(dec0) + ',' + str(crval1) + ',' + str(crval2) + ',' + \
                        str(orientat_header) + ',' + str(orientat_calc_prec) + ',' + str(EPS_prec) + ',' + str(roll0_prec) + ',' + \
                        str(roll_deg_prec) + ',' + str(v_factor0)
            f.write(allangle + '\n')
            f.close()
    sys.stdout = orig_stdout
    out_f.close()

# -----------------------------------------------------------------------------
#def parse_args():
#    '''
 #   Parse command line arguments, returns args object.
  #  '''

    #Create help string:
 #   filename_list_help = 'List of paths to files you need to find the roll at v2, v3 for.'
  #  outfile_date_help = 'Date file made.'
   # extension_help = 'UVIS chip extension that will specify which angle.stat file and which header (other than [0]) to get values from.'
    #path_to_data_help = 'The pwd to the subdirectory where the data is.'
    # Add arguments:
 #   parser = argparse.ArgumentParser()
 #   parser.add_argument('--filename_list', '-fl', dest = 'filename_list',
 #                       action = 'store', type = str, required = True,
 #                       help = filename_list_help)
 #   parser.add_argument('--outfile_date', '-out_date', dest = 'outfile_date',
 #                       action = 'store', type = str, required = True,
 #                       help = outfile_date_help)
 #   parser.add_argument('--extension', '-ext', dest = 'extension',
 #                       action = 'store', type = int, required = True,
 #                       help = extension_help)
#    parser.add_argument('--path_to_data', '-path', dest = 'path_to_data',
#                        action = 'store', type  = str, required = False,
#                        help = path_to_data_help,
#                        default = '/grp/hst/wfc3o/martlin/distortion_calibration_work_martlin/beginning_data')

    # Parse args:
 #   args = parser.parse_args()

  #  return args
## -----------------------------------------------------------------------------

#if __name__ == '__main__':

 #   args = parse_args()
 #   troll_main(args.filename_list, args.outfile_date, args.extension, args.path_to_data)
