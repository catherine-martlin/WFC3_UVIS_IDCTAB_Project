#! /usr/bin/env python
'''
ABOUT
-----
It reads in output from veracoeff.py and it creates the output IDC file.

-- infile - input file name.

-- other variable?

EXAMPLE
-------

python readwf3poly_onefilter --infile --other variable?

python  readwf3poly_onefilter --infile = 'test_output_from_veracoeff.txt'
'''
#import all the python things here!
import math
from astropy.io import fits
import os
import argparse
import numpy as np
import sys
# -----------------------------------------------------------------------------
def readwf3poly_onefilter_main(filename,outfile_name,screen_outputfile):

    '''
    The main controller/conversion tool.
    '''
    #pathname = '/grp/hst/wfc3o/martlin/idctab_vera/make_idctab_codes/'
    orig_stdout = sys.stdout
    out_f = file(screen_outputfile, 'a')
    sys.stdout = out_f

    print("   ")
    print("RUNNING readwf3poly_onefilter.py:")

    filenamepath = filename
    order = 4
    terms = (order+1)*(order+2)/2
    a = []
    b = []
    a1 = []
    a2 = []
    b1 = []
    b2 = []
    c = 0

    nim = order + 1

    nt = xrange(terms - 1)


    with open(filenamepath, 'r') as chipfile2:
        for x in xrange(len(nt)):
            for line in chipfile2:
                line = line.strip()
                temp_val = line.split(',')
                a.append(temp_val[0])
                b.append(temp_val[1])

    a1 = a[0:15]
    a2 = a[15:31]
    b1 = b[0:15]
    b2 = b[15:31]


    new_a1 = []
    new_a2 = []
    new_b1 = []
    new_b2 = []

    list_1 = [0,1,3,6,10]
    list_2 = [1,3,6,10,15]

    for x in xrange(0,order+1):
        y = list_1[x]
        z = list_2[x]
        new_a1.append(a1[y:z])
        new_a2.append(a2[y:z])
        new_b1.append(b1[y:z])
        new_b2.append(b2[y:z])

    print('A1')
    for x in xrange(0,order+1):
        print(new_a1[x])
    print('B1')
    for x in xrange(0,order+1):
        print(new_b1[x])
    print('A2')
    for x in xrange(0,order+1):
        print(new_a2[x])
    print('B2')
    for x in xrange(0,order+1):
        print(new_b2[x])

    print(' ')
    print('Reference points')
    print(a1[0], b1[0])
    v21 = a1[0]
    v31 = b1[0]

    print(a2[0], b2[0])
    v22 = a2[0]
    v32 = b2[0]

    # Derived quantites:
    for x in xrange(len(a1)):
        a1[x] = np.float64(a1[x])
        b1[x] = np.float64(b1[x])
        a2[x] = np.float64(a2[x])
        b2[x] = np.float64(b2[x])
    a_sqr_12 = a1[2]**(2.0)
    b_sqr_12 = b1[2]**(2.0)
    a_sqr_11 = a1[1]**(2.0)
    b_sqr_11 = b1[1]**(2.0)
    a_sqr_22 = a2[2]**(2.0)
    b_sqr_22 = b2[2]**(2.0)
    a_sqr_21 = a2[1]**(2.0)
    b_sqr_21 = b2[1]**(2.0)

    x1scale = math.sqrt(a_sqr_12 + b_sqr_12)
    y1scale = math.sqrt(a_sqr_11 + b_sqr_11)

    x2scale = math.sqrt(a_sqr_22 + b_sqr_22)
    y2scale = math.sqrt(a_sqr_21 + b_sqr_21)

    betax1_rad = math.atan2(a1[2], b1[2])
    betay1_rad = math.atan2(a1[1], b1[1])
    betax1 = math.degrees(betax1_rad)
    betay1 = math.degrees(betay1_rad)

    betax2_rad = math.atan2(a2[2], b2[2])
    betay2_rad = math.atan2(a2[1], b2[1])
    betax2 = math.degrees(betax2_rad)
    betay2 = math.degrees(betay2_rad)

    print(' ')
    print('X-scale', 'Y-scale', 'Beta-X', 'Beta-y', 'Beta-y - Beta-x')
    print(x1scale, y1scale, betax1, betay1, (betay1-betax1))
    print(x2scale, y2scale, betax2, betay2, (betay2-betax2))
    print(' ')

    #Shift to known UVIS reference from alignment measurement from June 2012:
    uv2v2 = -27.4596
    uv2v3 = -33.2604
    v21 = a1[0]-a2[0] + uv2v2
    v31 = b1[0]-b2[0] + uv2v3
    v22 = uv2v2
    v32 = uv2v3

    print('New V-values', v21, v31, v22, v32)
    print(' ')
    print(a1[0])
    print(a2[0])
    print(b1[0])
    print(b2[0])

    #IDC table output - Drizzle reference frame
    theta = 45.0
    neg_a1 = [-x for x in a1]
    neg_a2 = [-x for x in a2]
    ad1 = []
    ad2 = []
    bd1 = []
    bd2 = []
    for x in xrange(len(a1)):
        ad1.append((neg_a1[x] + b1[x])/math.sqrt(2.0))
        bd1.append((a1[x] + b1[x])/math.sqrt(2.0))
        ad2.append((neg_a2[x] + b2[x])/math.sqrt(2.0))
        bd2.append((a2[x] + b2[x])/math.sqrt(2.0))


    uvisfilters = ['F200LP', 'F300X ', 'F350LP', 'F475X ', 'F600LP', 'F850LP', 'F218W ','F225W',
                    'F275W ', 'F336W ', 'F390W ', 'F438W ', 'F475W ', 'F555W ', 'F606W ',
                    'F625W ', 'F775W ', 'F814W ', 'F390M ', 'F410M ', 'FQ422M', 'F467M ', 'F547M',
                    'F621M ', 'F689M ', 'F763M ', 'F845M ', 'FQ232N', 'FQ243N', 'F280N ',
                    'F343N ', 'F373N ', 'FQ378N', 'FQ387N', 'F395N ', 'FQ436N', 'FQ437N', 'F469N',
                    'F487N ', 'FQ492N', 'F502N ', 'FQ508N', 'FQ575N', 'FQ619N', 'F631N ',
                    'FQ634N', 'F645N ', 'F656N ', 'F657N ', 'F658N ', 'F665N ', 'FQ672N', 'F673N',
                    'FQ674N', 'F680N ', 'FQ727N', 'FQ750N', 'FQ889N', 'FQ906N', 'FQ924N', 'FQ937N', 'F953N ', 'G280 ']

    filtername = 'F606W'
    xsize = 4096
    ysize = 2051
    xref = 2048.0
    yref = 1026.0
    scale = 0.03962
    idcfile = outfile_name

    ad1_prec = []
    bd1_prec = []
    ad2_prec = []
    bd2_prec = []
    for x in xrange(len(ad1)):
        ad1_prec_t = '{:.7e}'.format(np.float64(ad1[x]))
        ad1_prec.append(np.float64(ad1_prec_t))
        bd1_prec_t = '{:.7e}'.format(np.float64(bd1[x]))
        bd1_prec.append(np.float64(bd1_prec_t))
        ad2_prec_t = '{:.7e}'.format(np.float64(ad2[x]))
        ad2_prec.append(np.float64(ad2_prec_t))
        bd2_prec_t = '{:.7e}'.format(np.float64(bd2[x])) 
        bd2_prec.append(np.float64(bd2_prec_t))

    with open(idcfile, 'w') as idctab_file:
        line1 = str(1) + '\t' + str('FORWARD') + '\t' + str(filtername) + '\t' + str(xsize) +'\t' + \
                str(ysize) + '\t' + str(xref) + '\t' + str(yref) +'\t' + str(theta) + '\t' \
                + str(scale) + '\t' + str(v21) + '\t' + str(v31) +'\n'
        line2 = ad1_prec[1:terms]
        line3 = bd1_prec[1:terms]
        line2_hardcode = str(line2[0]) + '\t' + str(line2[1]) + '\t' + str(line2[2]) + '\t' + str(line2[3]) +'\t' + \
                str(line2[4]) + '\t' + str(line2[5]) + '\t' + str(line2[6]) +'\t' + str(line2[7]) + '\t' \
                + str(line2[8]) + '\t' + str(line2[9]) + '\t' + str(line2[10]) + '\t' + str(line2[11]) + '\t' + str(line2[12]) + '\t' + \
                str(line2[13]) + '\n'
        line3_hardcode = str(line3[0]) + '\t' + str(line3[1]) + '\t' + str(line3[2]) + '\t' + str(line3[3]) +'\t' + \
                str(line3[4]) + '\t' + str(line3[5]) + '\t' + str(line3[6]) +'\t' + str(line3[7]) + '\t' \
                + str(line3[8]) + '\t' + str(line3[9]) + '\t' + str(line3[10]) + '\t' + str(line3[11]) + '\t' + str(line3[12]) + '\t' + \
                str(line3[13]) + '\n'
        line4 = str(2) + '\t' + str('FORWARD') + '\t' + str(filtername) + '\t' + str(xsize) +'\t' + \
                str(ysize) + '\t' + str(xref) + '\t' + str(yref) +'\t' + str(theta) + '\t' \
                + str(scale) + '\t' + str(v22) + '\t' + str(v32) +'\n'
        line5 = ad2_prec[1:terms]
        line6 = bd2_prec[1:terms]
        line5_hardcode = str(line5[0]) + '\t' + str(line5[1]) + '\t' + str(line5[2]) + '\t' + str(line5[3]) +'\t' + \
                str(line5[4]) + '\t' + str(line5[5]) + '\t' + str(line5[6]) +'\t' + str(line5[7]) + '\t' \
                + str(line5[8]) + '\t' + str(line5[9]) + '\t' + str(line5[10]) + '\t' + str(line5[11]) + '\t' + str(line5[12]) + '\t' + \
                str(line5[13]) + '\n'
        line6_hardcode = str(line6[0]) + '\t' + str(line6[1]) + '\t' + str(line6[2]) + '\t' + str(line6[3]) +'\t' + \
                str(line6[4]) + '\t' + str(line6[5]) + '\t' + str(line6[6]) +'\t' + str(line6[7]) + '\t' \
                + str(line6[8]) + '\t' + str(line6[9]) + '\t' + str(line6[10]) + '\t' + str(line6[11]) + '\t' + str(line6[12]) + '\t' + '\t' + \
                str(line6[13]) + '\n'
        idctab_file.writelines([line1, line2_hardcode, line3_hardcode, line4, line5_hardcode, line6_hardcode])


    print(' ')
    print('Data appended to', idcfile)

    sys.stdout = orig_stdout
    out_f.close()

# -----------------------------------------------------------------------------
#if __name__ == '__main__':

    #args = parse_args()
#    readwf3poly_onefilter_main(filename = 'veracoeff_combinded_new_vafactor_all_dec.txt')
