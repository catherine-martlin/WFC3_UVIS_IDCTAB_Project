#! /usr/bin/env python
'''
ABOUT
-----
Scales the results from previous code. It reads in output from coeff2v23_rot2mass.py
and it creates output that is used in the next step to making the IDCTAB.

--  .

-- infile - input file name.

-- chipnumber - Chipnumber is either 1 or 2 for the UVIS chips or 3 for the IR.

EXAMPLE
-------

python veracoeff_scaling.py --infile --chipnumber

python veracoeff_scaling.py --infile = 'uvis2.f606w.031011.vL.test.idl.comma.meancoeff'
                            --chipnumber = '2'
'''
#import all the python things here!
import math
from astropy.io import fits
import os
import argparse
import numpy as np
import sys
# -----------------------------------------------------------------------------
def veracoeff_scaling_main(chipfile, chipnumber,outfile_date,screen_outputfile):

    '''
    The main controller/conversion tool.
    '''
    orig_stdout = sys.stdout
    out_f = file(screen_outputfile, 'a')
    sys.stdout = out_f

    #pathname = '/grp/hst/wfc3o/martlin/idctab_vera/make_idctab_codes/'
    chipfilepath = chipfile

    order = 4
    terms = (order+1)*(order+2)/2

    # Defining various matrices:
    x = []
    y = []
    x_va = []
    y_va = []

    nim = order + 1
    a = np.zeros((nim,nim),np.float64)
    b = np.zeros((nim,nim),np.float64)
    a_va = np.zeros((nim,nim),np.float64)
    b_va = np.zeros((nim,nim),np.float64)

    # Read in information from previous IDCTAB step:
    print chipfilepath
    with open(chipfilepath, 'r') as chipfile2:
        for line in chipfile2:
            line = line.strip()
            print line
            a_temp = line.split('\t')
            #c,d,e,f,g = a_temp.split(',')
            x.append(np.float64(a_temp[1]))
            y.append(a_temp[2])
            x_va.append(a_temp[4])
            y_va.append(a_temp[5])
    
    test = np.array(x) * 0.04
    # Fill Matrices with values:
    # i_list and j_list are created from the IDL code and made to mimic that
    # method to ensure the results are the same.
    k = 0
    i_list = [0,1,1,2,2,2,3,3,3,3,4,4,4,4,4]
    j_list = [0,1,0,2,1,0,3,2,1,0,4,3,2,1,0]
    for m in xrange(len(i_list)):
        i = i_list[m]
        j = j_list[m]
        a[j,i] = x[k]
        b[j,i] = y[k]
        a_va[j,i] = x_va[k]
        b_va[j,i] = y_va[k]
        k = k+1

    #print x
    #print a[1,1]
    #print a[0,1]
    # Choosing a scaling value based on UVIS or IR chip then scaling values:
    if chipnumber == 3:
        pix = 0.128 #arcsec/pixel
    else:
        pix = 0.04 #arcsec/pixel

    a_scaled = (a*pix)
    b_scaled = (b*pix)

    a_sqr_2 = (a_scaled[1,1])**(2.0)
    a_sqr_3 = (a_scaled[0,1])**(2.0)
    b_sqr_2 = (b_scaled[1,1])**(2.0)
    b_sqr_3 = (b_scaled[0,1])**(2.0)

    #Now for vafactor affected ones: 
    a_va_scaled = (a_va*pix)
    b_va_scaled = (b_va*pix)

    a_sqr_2_va = (a_va_scaled[1,1])**(2.0)
    a_sqr_3_va = (a_va_scaled[0,1])**(2.0)
    b_sqr_2_va = (b_va_scaled[1,1])**(2.0)
    b_sqr_3_va = (b_va_scaled[0,1])**(2.0)

    # Defining the x and y scales and the opening angle:
    xscale = math.sqrt(a_sqr_2 + b_sqr_2)
    yscale = math.sqrt(a_sqr_3 + b_sqr_3)
    print 'xscale ',xscale
    print 'yscale ',yscale

    betax_rad = math.atan2(a[1,1], b[1,1])
    betay_rad = math.atan2(a[0,1], b[0,1])
    betax = math.degrees(betax_rad)
    betay = math.degrees(betay_rad)

    print 'betax', betax
    print 'betay', betay
    print 'Opening Angle', betay-betax

    # Defining the x and y scales and the opening angle for the vafactor ones:
    xscale_va = math.sqrt(a_sqr_2_va + b_sqr_2_va)
    yscale_va = math.sqrt(a_sqr_3_va + b_sqr_3_va)
    print 'xscale_va ',xscale_va
    print 'yscale_va ', yscale_va

    betax_rad_va = math.atan2(a_va[1,1], b_va[1,1])
    betay_rad_va = math.atan2(a_va[0,1], b_va[0,1])
    betax_va = math.degrees(betax_rad_va)
    betay_va = math.degrees(betay_rad_va)

    print 'betax_va', betax_va
    print 'betay_va', betay_va
    print 'Opening Angle_va', betay_va-betax_va

    # Creating a file name based on the chip value:
    if chipnumber == 3:
        detector = "IR"
    else:
        detector = "UVIS" + str(chipnumber)
    print detector

    outfile = "veracoeff_output_{}_{}.txt".format(detector,outfile_date)

    outfile_va = "veracoeff_output_vafactor_{}_{}.txt".format(detector,outfile_date)
    print outfile_va

    name_info = str(chipfile), detector

    # Now reordering the scaled values:
    la = np.zeros(terms,np.float64)
    lb = np.zeros(terms,np.float64)
    k = 0
    j_list_2 = [0,0,1,0,1,2,0,1,2,3,0,1,2,3,4] # This list was made to match the indexing of the IDL code.
    for x in xrange(len(i_list)):
        i = i_list[k]
        j = j_list_2[k]
        la[k] = a_scaled[j,i]
        lb[k] = b_scaled[j,i]
        #print k
        #print la[k], ' ', lb[k]
        k = k+1
    #print la
    #print lb

    # Saving scaled and reordered values to a file:
    with open(outfile, 'w') as outputfile:
        file_info = str(chipfile) + '\n'
        outputfile.write(file_info)
        for x in xrange(len(la)):
            term_la = la[x]
            term_lb = lb[x]
            line = str(term_la) + ',' + str(term_lb) + '\n'
            print line
            outputfile.write(str(line))

     # Now reordering the scaled values for vafactor:
    la_va = np.zeros(terms,np.float64)
    lb_va = np.zeros(terms,np.float64)
    k = 0
    j_list_2 = [0,0,1,0,1,2,0,1,2,3,0,1,2,3,4] # This list was made to match the indexing of the IDL code.
    for x in xrange(len(i_list)):
        i = i_list[k]
        j = j_list_2[k]
        la_va[k] = a_va_scaled[j,i]
        lb_va[k] = b_va_scaled[j,i]
        #print k
        #print la[k], ' ', lb[k]
        k = k+1
    #print la_va
    #print lb_va

    # Saving scaled and reordered values to a file:
    with open(outfile_va, 'w') as outputfile_va:
        file_info_va = str(chipfile) + '\n'
        outputfile_va.write(file_info_va)
        for x in xrange(len(la_va)):
            term_la_va = la_va[x]
            term_lb_va = lb_va[x]
            line_va = str(term_la_va) + ',' + str(term_lb_va) + '\n'
            print line_va
            outputfile_va.write(str(line_va))

    sys.stdout = orig_stdout
    out_f.close()

# -----------------------------------------------------------------------------

#if __name__ == '__main__':

    #args = parse_args()
#    veracoeff_main(chipfile = 'coeff2v23_with_vafactor_output_uvis2_f606w_test_2016_04_26.txt', chipnumber = '2',outfile_date='test_2016_04_26')
