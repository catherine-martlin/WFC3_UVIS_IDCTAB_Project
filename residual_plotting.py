"""Runs tweakreg and creates plot

This script reads in *fit.match files from tweakreg output and makes .png 
residual plots for each filter. This code must be ran with python 3 for correct rms. 


Author:

--------- Myles McKay, October 20, 2016

Use:

python residual_plotting.py --path='/grp/hst/wfc3v/martlin/Myles_GD_files/medium_filter_test/F410M/old_files' --filter='F####' --new_or_old_version='old'


Output:

Residual plots in png format

Dependencies:

Notes:
"""


import glob
from stsci.tools import teal
#from drizzlepac import tweakreg
import os
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import argparse

def residual_plotting_main(path,filter,new_or_old_version):

	os.chdir(path)
	os.system('pwd')               
	list=sorted(glob.glob('*fit.match'))
	image=sorted(glob.glob('*_flc.fits'))
	
	for file,im in zip(list, image): 	
		with open(file) as f:
			content= f.readlines(1)
			content=f.readlines(2)
			content=f.readlines(3)
			xRMS=content[0].split()[5]
			yRMS=content[0].split()[6]
			xRMS='%.3f' % float(xRMS)
			yRMS='%.3f' % float(yRMS)	
			hdu= fits.open(im)
			print (filter)	
	
			fig, [(ax0, ax2), (ax1, ax3)] = plt.subplots(2, 2, figsize=(12, 6),sharex='col',sharey='row')
			fig.tight_layout(pad=5, w_pad=0.5, h_pad=.625)
			table=np.loadtxt(file,usecols=(0,1,6,7))
			stars=sum(1 for _ in table)
			fig.suptitle('RMS(X)= {}, RMS(Y)= {}, Filter: {},\n #{}'.format(xRMS,yRMS,filter,stars))
			
			
		#-----------------------------------------------										
		######### X1 VS DX PLOT ##############
		#------------------------------------------------
			
		
			ax0.plot(table[:,0],table[:,2],'o',markersize=1, label='points')
			ax0.set_ylabel('DX (pixels)')
			ax0.set_xlim([np.min(table[:,0])-500, np.max(table[:,0])+500])
			ax0.set_ylim([-0.4,0.4])
			ax0.axhline(xmin=1,xmax=0,linewidth=2, color='red')
			polyfit=np.polyfit(table[:,0],table[:,2],5)
			sort=np.sort(table[:,0])
			polyfit_data=((polyfit[0]*(sort**5))+
			              (polyfit[1]*(sort**4))+
			              (polyfit[2]*(sort**3))+
			              (polyfit[3]*(sort**2))+
			              (polyfit[4]*(sort**1))+
			              (polyfit[5]))
			ax0.plot(sort,polyfit_data,color='green', linewidth=3)		
		
		#-----------------------------------------------										
		######## X1 VS DY PLOT ##############
		#------------------------------------------------
		
			ax1.plot(table[:,0],table[:,3],'o',markersize=1)
			ax1.set_xlabel('X1 (pixels)')
			ax1.set_ylabel('DY (pixels)')
			ax1.set_xlim([-2500, 2500])
			ax1.set_ylim([-0.4,0.4])
			ax1.axhline(xmin=1,xmax=0,linewidth=2, color='red')
			polyfit=np.polyfit(table[:,0],table[:,3],5)
			sort=np.sort(table[:,0])
			polyfit_data=((polyfit[0]*(sort**5))+
			              (polyfit[1]*(sort**4))+
			              (polyfit[2]*(sort**3))+
			              (polyfit[3]*(sort**2))+
			              (polyfit[4]*(sort**1))+
			              (polyfit[5]))
			ax1.plot(sort,polyfit_data,color='green', linewidth=3)
		
		
		#-----------------------------------------------										
		######## Y1 VS DX PLOT ##############
		#------------------------------------------------										
												
												
			ax2.plot(table[:,1],table[:,2],'o',markersize=1)
			ax2.axhline(xmin=1,xmax=0,linewidth=2, color='red')
			ax2.set_xlim([np.min(table[:,1])-500, np.max(table[:,1])+500])
			ax2.set_ylim([-0.4,0.4])
			polyfit=np.polyfit(table[:,1],table[:,2],5)
			sort=np.sort(table[:,1])
			polyfit_data=((polyfit[0]*(sort**5))+
			              (polyfit[1]*(sort**4))+
			              (polyfit[2]*(sort**3))+
			              (polyfit[3]*(sort**2))+
			              (polyfit[4]*(sort**1))+
			              (polyfit[5]))
			ax2.plot(sort,polyfit_data,color='green',markersize=1, linewidth=3)
		
		
		
		#-----------------------------------------------										
		######## Y1 VS DY PLOT ##############
		#------------------------------------------------
												
												
			ax3.plot(table[:,1],table[:,3],'o',markersize=1)
			ax3.set_xlabel('Y1 (pixels)')
			ax3.axhline(xmin=1,xmax=0,linewidth=2, color='red')
			ax3.set_xlim([np.min(table[:,1])-500, np.max(table[:,1])+500])
			ax3.set_ylim([-0.4,0.4])
			polyfit=np.polyfit(table[:,1],table[:,3],5)
			sort=np.sort(table[:,1])
			polyfit_data=((polyfit[0]*(sort**5))+
			              (polyfit[1]*(sort**4))+
			              (polyfit[2]*(sort**3))+
			              (polyfit[3]*(sort**2))+
			              (polyfit[4]*(sort**1))+
			              (polyfit[5]))
			ax3.plot(sort,polyfit_data,color='green', linewidth=3)
			plt.savefig('{}_{}.png'.format(new_or_old_version,file[0:13]))
			plt.clf()
			hdu.close()
		
#	os.system('mkdir *.fits.match.png /grp/hst/wfc3v/martlin/median_filter_tests/')	
#	os.system('cp *flc.png /grp/hst/wfc3v/martlin/Myles_GD_files/medium_filter_test/old_resids_plots')	
		
# -------------------------------------------------------------------
# For command line execution
# -------------------------------------------------------------------

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
    path_help = 'Path to the folder with files to run tweakreg.'
    filter_help = 'The filter that is being used'
    new_or_old_version_help= 'Are you using old or new idctab, npol and d2im files'
    # Add arguments:
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', '-path', dest = 'path', action = 'store',
                        type = str, required = True, help = path_help)
    parser.add_argument('--filter', '-filter', dest = 'filter', action = 'store',
                        type = str, required = True, help = filter_help)
    parser.add_argument('--new_or_old_version', '-new_or_old_version', dest = 'new_or_old_version', action = 'store',
                        type = str, required = True, help = new_or_old_version_help)                    
    # Parse args:
    args = parser.parse_args()

    return args
# -------------------------------------------------------------------
if __name__ == '__main__':
	args = parse_args()

	residual_plotting_main(args.path, args.filter, args.new_or_old_version)

		
		
		
		
		
		