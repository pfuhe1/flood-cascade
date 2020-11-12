# calc_q_returnperiod.py
# Takes input from mizuRoute discharge and calculate yearly maximum values, then a percentile e.g. 50th percentile to use as bankfull discharge
#
# Peter Uhe
# 2020/01/28
#

# Load modules
import os,sys,glob,pickle,shutil,socket
import numpy as np
from netCDF4 import Dataset,num2date
from osgeo import ogr
import csv
import datetime
from write_bdy_bci import write_bdy, write_bci_v2
import subprocess

import warnings
warnings.filterwarnings("error")

###############################################################################################
# Main script: Loop over discharge files
# First calculate yearly maximum discharge using cdo
# Then load that data into a big array
#
host = socket.gethostname()
if host[:3]=='bp1':
	mizuroute_outdir = '/work/pu17449/mizuRoute/output/'
	runname0 = 'GBM-p1deg_MSWEP2-2-ERA5'
	fpattern = os.path.join(mizuroute_outdir,'GBM-p1deg_90?_MSWEP2-2-ERA5-calibrated?_MSWEP2-2-ERA5','q_*.nc')
	# Template file to use for output
	template = os.path.join(mizuroute_outdir,'template.nc')
percentile = 50 # percentile assumed for bankfull flow

yrmax_dir = os.path.join(mizuroute_outdir,'yrmax_data')
files = glob.glob(fpattern)
fulldata = np.zeros([len(files),3527]) # river segment
for i,f_discharge in enumerate(files):
	#f_discharge = '/home/pu17449/data2/mizuRoute/merithydro/q_GBM_MERIT-Hydro_1988-1-1.nc'
	fname = os.path.basename(f_discharge)
	print(fname)
	runname = fname[2:-7]
	outdir  = os.path.join(yrmax_dir,runname)
	maxfile = os.path.join(outdir,fname[:-3]+'_yrmax.nc')
	if os.path.exists(maxfile):
		print('Yrmax Result already exists, skipping calculation')
	else:
		if not os.path.exists(outdir):
			os.makedirs(outdir)
			cdo_cmd = ['cdo','yearmax','-selvar,IRFroutedRunoff',f_discharge,maxfile]
			print(' '.join(cdo_cmd))
			ret = subprocess.call(cdo_cmd)
			if not ret==0:
				raise Exception('Error with cdo command')
	with Dataset(maxfile,'r') as f:
		fulldata[i,:] = f.variables['IRFroutedRunoff'][0,:]

# Now calculate bankfull from 50th percentile (e.g. 2 year return period flow)
bankfull = np.percentile(fulldata,percentile,axis=0)
print('Debug,bankfull:',bankfull.min(),bankfull.mean(),bankfull.max())

# Now write out bankfull to netcdf file, based on template
bankfull_file = os.path.join(mizuroute_outdir,'q_bankfull_'+runname0+'_'+str(percentile)+'ile.nc')
shutil.copy(template,bankfull_file)
with Dataset(bankfull_file,'a') as f:
	f.variables['IRFroutedRunoff'][0,:] = bankfull
print('Written',bankfull_file)
