''' Script to combine tas,pr,pet files for FUSE input
Requires first:
1. Run metsim_merged_MSWEP-p1deg.py
2. running metsim to produce pet

Peter Uhe April 2020
This script uses input files produced on Jasmin, then '''

import os,sys,subprocess
import socket


host = socket.gethostname()

if host=='anthropocene.ggy.bris.ac.uk':
	datadir = '/export/anthropocene/array-01/pu17449/FUSE_inputs/GBM-p1deg-JM'

ftas  = os.path.join(datadir,'ERA5_tas_day_19790101-20171031_GBM-p1deg.nc')
fpet  = os.path.join(datadir,'pet_day_MSWEP-ERA5_p1deg_19790401-20171031.nc')
fpr   = os.path.join(datadir,'MSWEP2-2_010deg_19790101-20171031_GBM-p1deg_v2.nc')
f_out = os.path.join(datadir,'GBM-p1deg_MSWEP-2-2-ERA5_input.nc')
override = True

if not os.path.exists(f_out):

	###########################################################################################
	# NOW set up FUSE inputs:
	tmpdir = os.path.join(datadir,'tmp')
	if not os.path.exists(tmpdir):
		os.mkdir(tmpdir)
	elif override:
		for f in glob.glob(tmpdir+'/*'):
			os.remove(f)

	################################################################
	# tas (input to FUSE): days 91 onwards to match metsim output
	tas_tmp1 = os.path.join(tmpdir,os.path.basename(ftas))
	tas_tmp2 = os.path.join(tmpdir,os.path.basename(ftas)[:-3]+'2.nc')
	if not os.path.exists(tas_tmp2) or override:
		# HACK: ERA5 bias corrected data is in degC already, so skip conversion
		#cdo_cmd = ['cdo','-L','subc,273.15','-seltimestep,91/99999999',tas_path,tas_tmp]
		cdo_cmd = ['cdo','-L','invertlat','-seltimestep,91/99999999',ftas,tas_tmp1]
		print(' '.join(cdo_cmd))
		ret = subprocess.call(cdo_cmd)
		if not ret==0:
			raise Exception('Error with cdo command')
		# Convert lat/lon to float rather than double precision
		ncap_cmd = ['ncap2','-s','lon=float(lon);lat=float(lat)',tas_tmp1,tas_tmp2]
		print(' '.join(ncap_cmd))
		ret = subprocess.call(ncap_cmd)
		if not ret==0:
			raise Exception('Error with ncap command')

	################################################################
	# pr (input to FUSE): days 91 onwards to match metsim output
	pr_tmp  = os.path.join(tmpdir,os.path.basename(fpr))
	if not os.path.exists(pr_tmp) or override:
		cdo_cmd = cdo_cmd = ['cdo','-seltimestep,91/99999999',fpr,pr_tmp]
		print(' '.join(cdo_cmd))
		subprocess.call(cdo_cmd)

	################################################################
	# merge filenames
	cdo_cmd = ['cdo','-L','merge',pr_tmp,tas_tmp2,fpet,f_out]
	print(' '.join(cdo_cmd))
	ret = subprocess.call(cdo_cmd)
	if not ret==0:
		raise Exception('Error with cdo command')

	# Finally rename lat and lon to latitude and longitude to match FUSE domain file
	nco_cmd = ['ncrename','-v','lat,latitude','-d','lat,latitude','-v','lon,longitude','-d','lon,longitude',f_out]
	print(' '.join(nco_cmd))
	ret = subprocess.call(nco_cmd)
	if not ret==0:
			raise Exception('Error with ncrename command')

	################################################################
	# Clean up temp files
	tempfiles = [tas_tmp1,tas_tmp2,pr_tmp]
	for fpath in tempfiles:
		if fpath is not None and os.path.exists(fpath):
			os.remove(fpath)
else:
	print('File already exists, skipping:',f_out)
