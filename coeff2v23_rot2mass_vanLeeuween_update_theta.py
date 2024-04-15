'''
ABOUT
-----
Calculates the rotation angles from individual UVIS detector x and y into
the V2V3 coordinate system. It reads in output from troll.py and it creates
the mean of each coefficient which is used in the next step to making the 
IDCTAB. This code uses the Van Leeuween transformations at one point. We will
also have a set that doesn't do that in order to test which may be better. (10/10/16)

NOTE
----
Need a file in directory you are running the code that has name 'uvis1'
or 'uvis2' that has a list of the .coeff files you want to run this on in it.

--uvis - UVIS chip being used and name of .coeff list file. Examples include:
         'uvis1' or 'uvis2'.

--filter_name - The name of the filter being run.

--outfile_date - The date this is being run.

--path_troll_output (Optional) - Path to the troll.py output file.
                                 Default is set.

--path_to_coeff_files (Optional) - Path to directory where the UVIS .coeff files
                                   are stored. Default is set.


EXAMPLE
-------

python coeff2v23_rot2mass.py --uvis --filter_name --outfile_date
                             --path_troll_output --path_to_coeff_files

python coeff2v23_rot2mass.py --uvis='uvis1' -filter='f606w'
                             -out_date='jan_11_2016'
'''
#import all the python things here!
import math
from astropy.io import fits
import os
import argparse
import numpy as np
import argparse
import pandas as pd
import sys
# -----------------------------------------------------------------------------
def convert_to_rads(item):
    '''
    The trigonometric functions take in radians so here we convert
    from arcseconds to degrees and then degrees into radians.
    '''

    item_degree = item * (1./3600.)
    item_rad = math.radians(item_degree)

    return item_rad
# -----------------------------------------------------------------------------
def get_troll_data(path_troll_output):
    '''
    Read in the values from orientat_pav3_allangle.new.stat which was the
    output of the troll.py. Be sure that the path to the file is correct.
    '''

    troll_output = pd.read_csv(path_troll_output, sep=",", header = None)
    troll_output.columns = ['name','ra_t','dec_t','crval1','crval2','orient_h','orient_c','eps','pav3','pav3corr','vafactor']

    return troll_output

# -----------------------------------------------------------------------------
def file_coeff_names(rootnames, uvis):
    '''
    '''

    if uvis == 'uvis2':
        file_end = '1.coeffs'
    if uvis == 'uvis1':
        file_end = '4.coeffs'

    coeff_file_names = [i + file_end for i in rootnames] 


    return coeff_file_names

# -----------------------------------------------------------------------------

def coeff2v23_rot2mass_vanLeeuween_update_theta_main(uvis, filter_name, outfile_date, path_troll_output, path_to_coeff_files,screen_outputfile):

    '''
    The main controller/conversion tool.
    '''
    orig_stdout = sys.stdout
    out_f = open(screen_outputfile, 'a')
    sys.stdout = out_f

    print("  ")
    print("RUNNING coeff2v23_rot2mass_vanLeeuween.py:")

    main_dir = os.getcwd()

    # Define path to troll output file, Read in the data from troll.py output:
    troll_output = get_troll_data(path_troll_output)

    rootnames = [i[0:10] for i in troll_output['name']]
    coeff_file_names = file_coeff_names(rootnames, uvis)
    troll_output['coeff_file_names'] = pd.Series(coeff_file_names)

    os.chdir(path_to_coeff_files)


    # From Van Leeuween:
    cx_3 = 3.5854120e-10
    cx_33 = 2.4239940e-07

    que_rad = cx_3/cx_33
    que_degree = math.degrees(que_rad)
    cos_que = math.cos(que_rad)
    sin_que = math.sin(que_rad)

    # Back to solving:

    coeff_nums = range(0,15)

    ay = []
    ax = []
    for i in troll_output['coeff_file_names']:
        with open(i, 'r') as coeff_file:
            coeff_file_values = [filter(None, line.strip().split(' ')) for line in coeff_file]
            coefficient_numbers = [coeff_file_values[i][0] for i in coeff_nums]
            ax_temp = [coeff_file_values[i][1] for i in coeff_nums]
            ay_temp = [coeff_file_values[i][4] for i in coeff_nums]
            ax_tfloats = [np.float64(i) for i in ax_temp]
            ay_tfloats = [np.float64(i) for i in ay_temp]
            ax.append(ax_tfloats)
            ay.append(ay_tfloats)

    coeff_output = pd.DataFrame({'coeff_names' : pd.Series(coeff_file_names),
                                'ax' : pd.Series(ax),
                                'ay' : pd.Series(ay)})

    all_ax_vafactor = []
    all_ay_vafactor = []
    for i in range(0,len(coeff_file_names)):
        ax_vafactor = (np.array(coeff_output['ax'][i])/np.array(troll_output['vafactor'][i]))
        ay_vafactor = (np.array(coeff_output['ay'][i])/np.array(troll_output['vafactor'][i]))
        all_ax_vafactor.append(ax_vafactor)
        all_ay_vafactor.append(ay_vafactor)
    coeff_output['ax_vafactor'] = pd.Series(all_ax_vafactor)
    coeff_output['ay_vafactor'] = pd.Series(all_ay_vafactor)

    frame_n_ax = [coeff_output['ax'][i] for i in range(0,len(coeff_file_names))]
    frame_n_ay = [coeff_output['ay'][i] for i in range(0,len(coeff_file_names))]

    print(" ")
    #print("Coeffs not divided by vafactor:", ax[10])
    #print("Coeffs divided by vafactor:",all_ax_vafactor[10])
    #print("VAFACTORS:", troll_output['vafactor'][10])
    print(" ")

    cx = []
    cy = []
    for i in range(0, len(frame_n_ax)):
        frame_temp_ax = frame_n_ax[i]
        frame_temp_ay = frame_n_ay[i]
        frame_axcos = [(np.float64(i) * cos_que) for i in frame_temp_ax]
        frame_aysin = [(np.float64(i) * sin_que) for i in frame_temp_ay]
        frame_axsin = [(np.float64(i) * sin_que) for i in frame_temp_ax]
        frame_aycos = [(np.float64(i) * cos_que) for i in frame_temp_ay]
        frame_ax_val = [frame_axcos[i] - frame_aysin[i] for i in range(0, len(frame_axcos))] 
        frame_ay_val = [frame_axsin[i] + frame_aycos[i] for i in range(0, len(frame_axsin))]
        cx.append(frame_ax_val)
        cy.append(frame_ay_val)

    coeff_output['cx'] = pd.Series(cx)
    coeff_output['cy'] = pd.Series(cy)

    #the above for vafactor: 
    frame_n_ax_va = [coeff_output['ax_vafactor'][i] for i in range(0,len(coeff_file_names))]
    frame_n_ay_va = [coeff_output['ay_vafactor'][i] for i in range(0,len(coeff_file_names))]

    cx_va = []
    cy_va = []
    for i in range(0, len(frame_n_ax_va)):
        frame_temp_ax_va = frame_n_ax_va[i]
        frame_temp_ay_va = frame_n_ay_va[i]
        frame_axcos_va = [(float(i) * cos_que) for i in frame_temp_ax_va]
        frame_aysin_va = [(float(i) * sin_que) for i in frame_temp_ay_va]
        frame_axsin_va = [(float(i) * sin_que) for i in frame_temp_ax_va]
        frame_aycos_va = [(float(i) * cos_que) for i in frame_temp_ay_va]
        frame_ax_val_va = [frame_axcos_va[i] - frame_aysin_va[i] for i in range(0, len(frame_axcos_va))] 
        frame_ay_val_va = [frame_axsin_va[i] + frame_aycos_va[i] for i in range(0, len(frame_axsin_va))]
        cx_va.append(frame_ax_val_va)
        cy_va.append(frame_ay_val_va)

    coeff_output['cx_vafactor'] = pd.Series(cx_va)
    coeff_output['cy_vafactor'] = pd.Series(cy_va)


    #define Theta:
    #(Uses new catalogue values from NGC5139_WFC3UV_F606W.fits (/grp/hst/wfc3h/verap/STANDARD_CAT/JAYWORK/99.EXPORT))
    ra_c = 2.016970062256e+02
    dec_c = -4.74794731140e+01
    #dec_c =-4.747222222222e+01 # Old Catalogue.

    epsilon_deg_list = []
    ra_t_all = troll_output['ra_t']
    for x in xrange(len(ra_t_all)):
        ra_t = ra_t_all[x]
        d_ra = (ra_t - ra_c)
        d_ra_rad = math.radians(d_ra)
        dec_c_rad = math.radians(dec_c)
        tan_d_ra = math.tan(d_ra_rad)
        sin_dec_c = math.sin(dec_c_rad)
        atan_all = math.atan(tan_d_ra * sin_dec_c)
        epsilon_rad = atan_all
        epsilon_deg = math.degrees(epsilon_rad)
        epsilon_deg_list.append(epsilon_deg)

    theta_rad_list = []
    #sin_t_list = []
    #cos_t_list = []
    orient_c = troll_output['orient_c']
    eps = troll_output['eps']
    pav3corr = troll_output['pav3corr']
    #print(rootnames[64])
    #print(orient_c[64])
    #print(eps[64])
    #print(pav3corr[64])
  

    for x in xrange(len(orient_c)):
        epsilon_deg_now = epsilon_deg_list[x]
        #print(x, " ", orient_c[x], " ", eps[x], " ",epsilon_deg_now, " ",pav3corr[x])
        theta_deg_now = (orient_c[x] - eps[x]) + epsilon_deg_now - pav3corr[x] # This is in degrees (not arcseconds).
        theta_rad = math.radians(theta_deg_now)
        theta_rad_list.append(theta_rad)
        #print(x, ' ', theta_rad)
        # Define some trig values from theta calculations:
        #sin_t = math.sin(theta_rad)
        #sin_t_list.append(sin_t)
        #cos_t = math.cos(theta_rad)
        #cos_t_list.append(cos_t)

    #mean_theta_rad = np.mean(theta_rad_list)
    #mean_theta_deg = math.degrees(mean_theta_rad)
    #print("Theta list: ", theta_rad_list)
    #print("Mean Theta new catalogue (rad) = ", mean_theta_rad)
    #print("Mean Theta new catalogue(deg) = ", mean_theta_deg)

    # 12/21/2016: We want to try to calculate with the mean theta from the old catalogue:
    # old_dec_theta=44.663478; old_rad_theta=0.77952474
    
    mean_old_theta_rad = 0.77952474
    mean_old_theta_deg = 44.663478

    print("Mean Theta old catalogue (rad) = ", mean_old_theta_rad)
    print("Mean Theta old catalogue(deg) = ", mean_old_theta_deg)

    sin_mean_theta = math.sin(mean_old_theta_rad)
    cos_mean_theta = math.cos(mean_old_theta_rad)

    #coeff_output['cos_t'] = pd.Series(cos_t_list)
    #coeff_output['sin_t'] = pd.Series(sin_t_list)

    r_cx_list = []
    r_cy_list = []
    r_cy_list_va = []
    r_cx_list_va = []


    for i in range(0,len(coeff_file_names)):
        coeff_x_sin = np.array(coeff_output['cx'][i]) * sin_mean_theta
        coeff_y_sin = np.array(coeff_output['cy'][i]) * sin_mean_theta
        coeff_x_cos = np.array(coeff_output['cx'][i]) * cos_mean_theta
        coeff_y_cos = np.array(coeff_output['cy'][i]) * cos_mean_theta
        cx_r_before_mean = (-coeff_x_cos + coeff_y_sin)
        cy_r_before_mean = (coeff_x_sin + coeff_y_cos)
        r_cx_list.append(cx_r_before_mean)
        r_cy_list.append(cy_r_before_mean)


    for i in range(0,len(coeff_file_names)):
        coeff_x_sin_va = np.array(coeff_output['cx_vafactor'][i]) * sin_mean_theta
        coeff_y_sin_va = np.array(coeff_output['cy_vafactor'][i]) * sin_mean_theta
        coeff_x_cos_va = np.array(coeff_output['cx_vafactor'][i]) * cos_mean_theta
        coeff_y_cos_va = np.array(coeff_output['cy_vafactor'][i]) * cos_mean_theta
        cx_r_before_mean_va = (-coeff_x_cos_va + coeff_y_sin_va)
        cy_r_before_mean_va = (coeff_x_sin_va + coeff_y_cos_va)
        r_cx_list_va.append(cx_r_before_mean_va)
        r_cy_list_va.append(cy_r_before_mean_va)
        #print(coeff_output['coeff_names'][i])
        #print(coeff_output['cx'][i])
        #print(' ')
        #print(troll_output['name'][i])
        #print(troll_output['vafactor'][i])
        #print(' ')
        #print((np.array(coeff_output['cx'][i])/np.array(troll_output['vafactor'][i])))


    #print(coeff_output)
    #vafactors = troll_output['vafactor']
    #mean_vafactor = sum(vafactors)/len(vafactors)
    #print(" ")
    #print("mean_vafactor value = ", mean_vafactor)
    #print(" ")

    m_cx = sum(r_cx_list)/len(r_cx_list)
    m_cy = sum(r_cy_list)/len(r_cy_list)
    m_cx_va = sum(r_cx_list_va)/len(r_cx_list_va)
    m_cy_va = sum(r_cy_list_va)/len(r_cy_list_va)
    
    os.chdir(main_dir)
    coeff_num = range(1,16)
    m_cx_prec = []
    m_cy_prec = []
    m_cx_prec_va = []
    m_cy_prec_va = []
    for i in range(0,len(m_cx)):
        m_cx_prec_t = '{:}'.format(np.float64(m_cx[i]))
        m_cx_prec.append(np.float64(m_cx_prec_t))
        m_cy_prec_t = '{:}'.format(np.float64(m_cy[i]))
        m_cy_prec.append(np.float64(m_cy_prec_t))
        m_cx_prec_t_va = '{:}'.format(np.float64(m_cx_va[i]))
        m_cx_prec_va.append(np.float64(m_cx_prec_t_va))
        m_cy_prec_t_va = '{:}'.format(np.float64(m_cy_va[i]))
        m_cy_prec_va.append(np.float64(m_cy_prec_t_va))
    col1,col2,col3,col4,col5,col6 = coeff_num, m_cx_prec,m_cy_prec,coeff_num,m_cx_prec_va,m_cy_prec_va
    output = '\n'.join('\t'.join(map(str,row))for row in zip(col1,col2,col3,col4,col5,col6))
    filter_name = filter_name
    chip = uvis
    outfile_name = 'coeff2v23_vafactor_output_vL_mean_old_cat_theta_{}_{}_{}.txt'.format(chip, filter_name, outfile_date)
    f = open(outfile_name, 'w')
    f.write(output)
    f.close()
    #print((output))
    #print("   ")

    sys.stdout = orig_stdout
    out_f.close()

# -----------------------------------------------------------------------------
#def parse_args():
#    '''
#    Parse command line arguments, returns args object.
#    '''

#    Create help string:
#    filename_list_help = 'List of paths to files you need to find the roll at v2, v3 for.'
#    uvis_help = 'UVIS chip input. (e.g. uvis1 or uvis2).'
#    filter_name_help = 'Filter name (e.g. f606w).'
#    outfile_date_help = 'Date file is created (e.g. jan_11_2016; feb_22_2045; etc.'
#    path_troll_output_help = 'Path to the troll.py file.'
#    path_to_coeff_files_help = 'Path to where the UVIS .coeff files are stored.'
#     Add arguments:
#    parser = argparse.ArgumentParser()
#    parser.add_argument('--uvis', dest = 'uvis', action = 'store', type = str,
#                        required = False, help = uvis_help, default = 'uvis1')
#    parser.add_argument('--filter_name', '-filter', dest = 'filter_name', action = 'store',
#                        type = str, required = False, help = filter_name_help, default = 'f606w')
#    parser.add_argument('--outfile_date', '-out_date', dest = 'outfile_date', action = 'store',
#                        type = str, required = False, help = outfile_date_help, default = '01.01.2016')
#    parser.add_argument('--path_troll_output', '-path_troll', dest = 'path_troll_output', action='store', type = str,
#                        required = True, help = path_troll_output_help)
#    parser.add_argument('--path_to_coeff_files', '-path_coeff', dest = 'path_to_coeff_files', 
#                        action = 'store', type = str, required = True, help = path_to_coeff_files_help)
#    Parse args:
#    args = parser.parse_args()

#    return args
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------

#if __name__ == '__main__':

 #   args = parse_args()
  #  update_coeff2v23_rot2mass_main(args.uvis, args.filter_name, args.outfile_date, args.path_troll_output, args.path_to_coeff_files)
