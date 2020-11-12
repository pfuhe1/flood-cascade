import pickle,os,glob
import subprocess
import numpy as np

# Peter Uhe Nov 2019:
#
# Merges meterological inputs for catchments into single file
# Run this script after `metsim_inputs_ERA5_mswep_separatecatchments_loopappend.py`
#

#############################################################################
# Input paths
fcatchment_list = '/export/anthropocene/array-01/pu17449/FUSE_inputs/catchments_p1deg/GBM-p1deg_donorsbest3.txt'
timeperiod = '19790101-20171231'
varnames = ['pr','tasmax','tasmin']
outname = 'metsim-forcing'
datadir = '/export/anthropocene/array-01/pu17449/FUSE_inputs/catchments_p1deg_v3'
tmpdir = '/export/anthropocene/array-01/pu17449/FUSE_inputs/tmp'
outdir = os.path.join(datadir,outname)
if not os.path.exists(outdir):
	os.mkdir(outdir)

if not os.path.exists(tmpdir):
	os.mkdir(tmpdir)



#############################################################################
# Loop over catchmentsn and calculate forcings
with open(fcatchment_list,'r') as f:
	for line in f:
		grdcstr = line.strip()
		grdcid = int(grdcstr)
		outfile = '_'.join([grdcstr,outname,timeperiod+'.nc'])
		fpath_out = os.path.join(outdir,outfile)
		if not os.path.exists(fpath_out):
			catchment_files = []
			for var in varnames:
				varfile = '_'.join([grdcstr,var,timeperiod+'.nc'])
				fpath_var = os.path.join(datadir,var,varfile)
				if var[:3]=='tas':
					#Convert tas  to degC
					ftmp = os.path.join(tmpdir,varfile)
					cdo_cmd1 = ['cdo','subc,273.15',fpath_var,ftmp]
					print(' '.join(cdo_cmd1))
					subprocess.call(cdo_cmd1)
					catchment_files.append(ftmp)
				else:
					catchment_files.append(fpath_var)

			cdo_cmd = ['cdo','merge']+catchment_files+[fpath_out]
			print(' '.join(cdo_cmd))
			ret = subprocess.call(cdo_cmd)
			if not ret==0:
				print('Error, cdo command failed')
				#raise Exception('Error, cdo command failed')

			# NCO command to update tas units to degC
			nco_cmd = ['ncatted','-a','units,tasmax,m,c,degC','-a','units,tasmin,m,c,degC',fpath_out]
			print(' '.join(nco_cmd))
			ret = subprocess.call(nco_cmd)
			if not ret ==0:
				print('Error, nco command failed')
