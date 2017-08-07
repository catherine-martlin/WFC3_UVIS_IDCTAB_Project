#! /usr/bin/env python

"""Generates the IDCTAB. 

This script serves as a pipeline to create the IDCTAB from various 
other inputs and is a wrapper around several modules that perform 
subtasks of the IDCTAB creation algorithm.

Author: 
	Catherine Martlin, 2016.

Use: python make_individual_idctab.py --filter_name='F606W' 
--filename_list='f606w_oct_2016_flc_files_removed_mismatched_odd_file_removed.txt' 
--path_to_data='/grp/hst/wfc3h/verap/F606W_FLC/xymc/' 
--path_to_coeff_files='/grp/hst/wfc3h/verap/F606W_FLC/meta_pix.v5_new/'


Outputs:

Dependencies:

Notes:
The angle.stat file is currently needed to be in the same directory as the path_to_coeffs

"""
import argparse
from troll import troll_main
#from coeff2v23_rot2mass import coeff2v23_rot2mass_main
from coeff2v23_rot2mass_vanLeeuween_update_theta import coeff2v23_rot2mass_vanLeeuween_update_theta_main
from veracoeff_scaling import veracoeff_scaling_main
from readwf3poly_onefilter import readwf3poly_onefilter_main
import datetime
import numpy as np
import glob
# ---------------------------------------------------------------------
def prep_files_for_readwf3poly(veracoeff_output, coeff2v23_output):
    """

        """

    list_of_values = []
    veracoeff_output_correct_order = [veracoeff_output[1],veracoeff_output[0]]

    for i in veracoeff_output_correct_order:
        with open(i,'r') as output:
            for line in output:
                line = line.strip('\n')
                list_of_values.append(line)

    print list_of_values
    for j in coeff2v23_output:
        list_of_values.remove(j)
    print list_of_values

    if 'vafactor' in i:
        outfile_name = 'veracoeff_combined_output_vafactor.txt'
    else: 
        outfile_name = 'veracoeff_combined_output.txt'
    f = open(outfile_name, 'w')
    for k in list_of_values:
        print k
        f.write(k + '\n')
    f.close()

    return outfile_name

# ---------------------------------------------------------------------
# The main controller
# ---------------------------------------------------------------------

def make_individual_idctab_main(filter_name, filename_list, path_to_data,path_to_coeff_files):
    """ The main controller of the make_idctab_main

        Parameters:

        Returns:

        Outputs:

        """
    outfile_date = datetime.date.today()
    screen_outputfile = filter_name + '_'+ str(outfile_date) + '_outputs.info'

    individual_files = []
    with open(filename_list,'r') as all_files:
        for line in all_files:
            line = line.strip('\n')
            individual_files.append(line)

    for k in individual_files:
        name = k + '_run_individual_file'
        new_file = open(name, 'w')
        new_file.write(k)
        new_file.close()

    running_files = glob.glob('*_run_individual_file')

    for j in running_files:
        extension_list = [1,4]

        troll_outputs = []
        for i in extension_list:
            print " "
            print "RUNNING troll.py:"
            troll_main(j, outfile_date, i, path_to_data,path_to_coeff_files,screen_outputfile)
            troll_outfile = 'troll_output_{}_{}.txt'.format(i,outfile_date)
            troll_outputs.append(troll_outfile)
            print troll_outfile

        uvis_names = ['uvis2','uvis1']
        uvis_troll = zip(troll_outputs,uvis_names)

        coeff2v23_output = []
        for i in uvis_troll:
            print " "
            print "RUNNING coeff2v23_rot2mass.py:"
            uvis = i[1]
            troll_output = i[0]
            coeff2v23_rot2mass_vanLeeuween_update_theta_main(uvis,filter_name,outfile_date,troll_output,path_to_coeff_files,screen_outputfile)
            coeff_outfile_name = 'coeff2v23_vafactor_output_vL_mean_old_cat_theta_{}_{}_{}.txt'.format(uvis, filter_name, outfile_date)
            coeff2v23_output.append(coeff_outfile_name)

        chip_number = ['2','1']
        chip_coeff2v23 = zip(chip_number,coeff2v23_output)
        veracoeff_output = []
        veracoeff_output_vafactor = []
        for i in chip_coeff2v23:
            print " "
            print "RUNNING veracoeff.py:"
            chipnumber = i[0]
            infile = i[1]
            veracoeff_scaling_main(infile, chipnumber,outfile_date,screen_outputfile)
            detector = "UVIS" + str(chipnumber)
            outfile_veracoeff = "veracoeff_output_{}_{}.txt".format(detector,outfile_date)
            outfile_veracoeff_vafactor = "veracoeff_output_vafactor_{}_{}.txt".format(detector,outfile_date)
            veracoeff_output.append(outfile_veracoeff)
            veracoeff_output_vafactor.append(outfile_veracoeff_vafactor)
            print outfile_veracoeff
            print outfile_veracoeff_vafactor

        #Prep output of veracoeff for readwf3poly:
        combined_veracoeff_output = prep_files_for_readwf3poly(veracoeff_output, coeff2v23_output)
        combined_veracoeff_output_vafactor = prep_files_for_readwf3poly(veracoeff_output_vafactor, coeff2v23_output)

        print " "
        print "RUNNING readwf3poly_perfilter.py:"
        readwf3poly_outfile = 'UVIS_IDC_{}_{}_vL.txt'.format(outfile_date,j)
        readwf3poly_outfile_vafactor = 'UVIS_IDC_vafactor_{}_{}_vL.txt'.format(outfile_date,j)
        readwf3poly_onefilter_main(combined_veracoeff_output,readwf3poly_outfile,screen_outputfile)
        print combined_veracoeff_output
        readwf3poly_onefilter_main(combined_veracoeff_output_vafactor,readwf3poly_outfile_vafactor,screen_outputfile)
        print combined_veracoeff_output_vafactor

    print " "
    print "This IDCTAB was found using the VanLeeuween rotations and the old catalogue theta value."

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
    filename_list_help = 'List of paths to files you need to find the roll at v2, v3 for.'
    filter_name_help = 'Filter name (e.g. f606w).'
    path_to_coeff_files_help = 'Path to where the UVIS .coeff files are stored.'
    path_to_data_help = 'The pwd to the subdirectory where the data is.'
    # Add arguments:
    parser = argparse.ArgumentParser()
    parser.add_argument('--filter_name', '-filter', dest = 'filter_name', action = 'store',
                        type = str, required = False, help = filter_name_help, default = 'f606w')
    parser.add_argument('--path_to_coeff_files', '-path_coeff', dest = 'path_to_coeff_files', 
                        action = 'store', type = str, required = False, help = path_to_coeff_files_help,
                        default = '/grp/hst/wfc3o/martlin/distortion_calibration_work_martlin/coeffs/')
    parser.add_argument('--path_to_data', '-path', dest = 'path_to_data',
                        action = 'store', type  = str, required = False,
                        help = path_to_data_help,
                        default = '/grp/hst/wfc3o/martlin/distortion_calibration_work_martlin/beginning_data/')
    parser.add_argument('--filename_list', '-fl', dest = 'filename_list',
                        action = 'store', type = str, required = True,
                        help = filename_list_help)

    # Parse args:
    args = parser.parse_args()

    return args
 # -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

if __name__ == '__main__':
	args = parse_args()

	make_individual_idctab_main(args.filter_name, args.filename_list, args.path_to_data, args.path_to_coeff_files)

