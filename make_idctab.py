#! /usr/bin/env python

"""Generates the IDCTAB. 

This script serves as a pipeline to create the IDCTAB from various 
other inputs and is a wrapper around several modules that perform 
subtasks of the IDCTAB creation algorithm.

Author(s)
---------
    Catherine Martlin, January 2017, 2024.

Use
--- 

python make_idctab.py --filter_name='F606W' --filename_list='f606w_oct_2016_flc_files_removed_mismatched_odd_file_removed.txt' --path_to_data='/grp/hst/wfc3h/verap/F606W_FLC/xymc/' --path_to_coeff_files='/grp/hst/wfc3h/verap/F606W_FLC/meta_pix.v32_new/'


python make_idctab.py --filter_name='F410M' --filename_list='F410M_oldcat_filenames.txt' 
--path_to_data='/grp/hst/wfc3h/verap/CAL_14393/F410M_oldcat/xymc/' 
--path_to_coeff_files='/grp/hst/wfc3h/verap/CAL_14393/F410M_oldcat/meta_pix.v3/'

Dependencies
------------

Notes
-----
"""


import argparse
from troll import troll_main
#from coeff2v23_rot2mass import coeff2v23_rot2mass_main
from coeff2v23_rot2mass_vanLeeuween_update_theta import coeff2v23_rot2mass_vanLeeuween_update_theta_main
from veracoeff_scaling import veracoeff_scaling_main
from readwf3poly_onefilter import readwf3poly_onefilter_main
import datetime
import numpy as np
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

    print(list_of_values)
    for j in coeff2v23_output:
        list_of_values.remove(j)
    print(list_of_values)

    if 'vafactor' in i:
        outfile_name = 'veracoeff_combined_output_vafactor.txt'
    else: 
        outfile_name = 'veracoeff_combined_output.txt'
    f = open(outfile_name, 'w')
    for k in list_of_values:
        print(k)
        f.write(k + '\n')
    f.close()

    return outfile_name

# ---------------------------------------------------------------------
# The main controller
# ---------------------------------------------------------------------

def make_idctab_main(filter_name, filename_list, path_to_data,path_to_coeff_files):
    """ The main controller of the make_idctab_main

        Parameters
        ----------

        Returns
        -------

        """
    outfile_date = datetime.date.today()
    extension_list = [1,4]
    uvis_list = ['uvis1','uvis2']

    screen_outputfile = filter_name + '_'+ str(outfile_date) + '_outputs.info'

    troll_outputs = []
    print( " ")
    print("RUNNING troll.py:")
    for i in extension_list:
        print(screen_outputfile)
        print(path_to_coeff_files)
        troll_main(filename_list, outfile_date, i, path_to_data,path_to_coeff_files,screen_outputfile)
        troll_outfile = 'troll_output_{}_{}.txt'.format(i,outfile_date)
        troll_outputs.append(troll_outfile)
        print(troll_outfile)

    uvis_names = ['uvis2','uvis1']
    uvis_troll = zip(troll_outputs,uvis_names)

    coeff2v23_output = []
    print(" ")
    print("RUNNING coeff2v23_rot2mass.py:")
    for k in uvis_troll:
        print(k)
        uvis = k[1]
        troll_output = k[0]
        #if vanL == 'True': 
        coeff2v23_rot2mass_vanLeeuween_update_theta_main(uvis,filter_name,outfile_date,troll_output,path_to_coeff_files,screen_outputfile)
        #coeff2v23_rot2mass_vanLeeuween_main(uvis,filter_name,outfile_date,troll_output,path_to_coeff_files,screen_outputfile)
        #coeff_outfile_name = 'coeff2v23_with_vafactor_output_vL_mean_theta_{}_{}_{}.txt'.format(uvis, filter_name, outfile_date)
        
        # If you want to use the new_cat theta then we need
        # to go into coeff2v23... and change the theta value 
        # being used and uncomment line 111 and comment 116. 
        coeff_outfile_name = 'coeff2v23_vafactor_output_vL_mean_old_cat_theta_{}_{}_{}.txt'.format(uvis, filter_name, outfile_date)
        print(coeff_outfile_name)
        coeff2v23_output.append(coeff_outfile_name)
            
    chip_number = ['2','1']
    chip_coeff2v23 = zip(chip_number,coeff2v23_output)
    veracoeff_output = []
    veracoeff_output_vafactor = []
    print(" ")
    print("RUNNING veracoeff.py:")
    for l in chip_coeff2v23:
        print(chip_coeff2v23)
        chipnumber = l[0]
        infile = l[1]
        veracoeff_scaling_main(infile, chipnumber,outfile_date,screen_outputfile)
        detector = "UVIS" + str(chipnumber)
        outfile_veracoeff = "veracoeff_output_{}_{}.txt".format(detector,outfile_date)
        outfile_veracoeff_vafactor = "veracoeff_output_vafactor_{}_{}.txt".format(detector,outfile_date)
        print(outfile_veracoeff)
        print(outfile_veracoeff_vafactor)
        veracoeff_output.append(outfile_veracoeff)
        veracoeff_output_vafactor.append(outfile_veracoeff_vafactor)

    #Prep output of veracoeff for readwf3poly:
    combined_veracoeff_output = prep_files_for_readwf3poly(veracoeff_output, coeff2v23_output)
    combined_veracoeff_output_vafactor = prep_files_for_readwf3poly(veracoeff_output_vafactor, coeff2v23_output)
    print(" ")

    readwf3poly_outfile = 'UVIS_{}_CM_{}_idc.txt'.format(filter_name,outfile_date)
    readwf3poly_outfile_vafactor = 'UVIS_vafactor_{}_CM_{}_idc.txt'.format(filter_name,outfile_date)
    print(" ")
    print("RUNNING readwf3poly_perfilter.py:")
    readwf3poly_onefilter_main(combined_veracoeff_output,readwf3poly_outfile,screen_outputfile)
    print(combined_veracoeff_output)
    readwf3poly_onefilter_main(combined_veracoeff_output_vafactor,readwf3poly_outfile_vafactor,screen_outputfile)
    print(combined_veracoeff_output_vafactor)

    print(" ")
    print("This IDCTAB was found using the VanLeeuween rotations and the mean old catalogue theta value.")


# -----------------------------------------------------------------------------
# For command line execution
# -----------------------------------------------------------------------------

def parse_args():
    """Parses command line arguments.

    Parameters
    ----------
        nothing

    Returns
    -------
        args : argparse.Namespace object
            An argparse object containing all of the added arguments.

    """

    #Create help string:
    filename_list_help = 'List of paths to files you need to find the roll at v2, v3 for.'
    filter_name_help = 'Filter name (e.g. f606w).'
    path_to_coeff_files_help = 'Path to where the UVIS .coeff files are stored.'
    path_to_data_help = 'The pwd to the subdirectory where the data is.'
    #vanL_help = 'Do you want to use the VanLeeuween rotations to solve?'
    # Add arguments:
    parser = argparse.ArgumentParser()
    parser.add_argument('--filter_name', '-filter', dest = 'filter_name', action = 'store',
                        type = str, required = False, help = filter_name_help, default = 'f606w')
    parser.add_argument('--path_to_coeff_files', '-path_coeff', dest = 'path_to_coeff_files', 
                        action = 'store', type = str, required = False, help = path_to_coeff_files_help,
                        default = '/grp/hst/wfc3v/martlin/idctab_creation/coeffs/')
    parser.add_argument('--path_to_data', '-path', dest = 'path_to_data',
                        action = 'store', type  = str, required = False,
                        help = path_to_data_help,
                        default = '/grp/hst/wfc3v/martlin/idctab_creation/beginning_data/')
    parser.add_argument('--filename_list', '-fl', dest = 'filename_list',
                        action = 'store', type = str, required = True,
                        help = filename_list_help)
    #parser.add_argument('--vanL', '-vanL', dest = 'vanL',
    #                    action = 'store', type = str, required = True,
    #                    help = vanL_help)    

    # Parse args:
    args = parser.parse_args()

    return args
 # -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_args()

    make_idctab_main(args.filter_name, args.filename_list, args.path_to_data, args.path_to_coeff_files)

