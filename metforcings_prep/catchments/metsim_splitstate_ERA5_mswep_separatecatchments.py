import pickle,os,glob
import subprocess
import numpy as np

# Peter Uhe Nov 2019:
#
# Loop over catchments
# Reads in data calculated by `metsim_merge_ERA5_mswep_separatecatchments.py`
# Splits 'metsim-forcing' file into a state file (first 90 days) and 'forcing' file (remainder of the timeseries)

#############################################################################
# Input paths
fcatchment_list = '/export/anthropocene/array-01/pu17449/FUSE_inputs/catchments_p1deg/GBM-p1deg_donorsbest3.txt'
timeperiod = '19790101-20171231'
var = 'metsim-forcing'
datadir = '/export/anthropocene/array-01/pu17449/FUSE_inputs/catchments_p1deg_v3'




#############################################################################
# Loop over catchmentsn and calculate forcings
with open(fcatchment_list,'r') as f:
	for line in f:
		grdcstr = line.strip()
		grdcid = int(grdcstr)
		catchment_files = []
		varfile = '_'.join([grdcstr,var,timeperiod+'.nc'])
		fpath_var = os.path.join(datadir,var,varfile)

		statefile = varfile = '_'.join([grdcstr,'state_19790101-19790331.nc'])
		fpath_state = os.path.join(datadir,var,statefile)
		cdo_cmd = ['cdo','seldate,1979-01-01,1979-03-31',fpath_var,fpath_state]
		print(' '.join(cdo_cmd))
		ret = subprocess.call(cdo_cmd)
		if not ret==0:
			print('Error, cdo command failed')

		forcingfile = varfile = '_'.join([grdcstr,'forcing_19790401-20171031.nc'])
		fpath_forcing = os.path.join(datadir,var,forcingfile)
		cdo_cmd = ['cdo','seldate,1979-04-01,2017-10-31',fpath_var,fpath_forcing]
		print(' '.join(cdo_cmd))
		ret = subprocess.call(cdo_cmd)
		if not ret==0:
			print('Error, cdo command failed')
		#raise Exception('Error, cdo command failed')
