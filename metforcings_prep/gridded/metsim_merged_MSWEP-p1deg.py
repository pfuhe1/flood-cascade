# This script produces merged files containing pr,tasmax,tasmin
# These are input files for Metsim to calculate PET
#
# Peter Uhe
# 2020-01-14
#
# Input files are 0.1 degree resolution
# MSWEP v2-2 precipitation and ERA5 temperature (bias corrected to worldclim)


import os,subprocess,sys,glob

adata = '/export/anthropocene/array-01/pu17449/'
outdir    = adata+'FUSE_inputs/GBM-p1deg'
processdir = adata+'FUSE_inputs/tmp'

if not os.path.exists(outdir):
	os.mkdir(outdir)

if not os.path.exists(processdir):
	os.mkdir(processdir)

# Region extent
#grid_x = [73,98]
#grid_y = [22,31.5]
lonlatbox = '73,98,22,31.5'
boxname = 'GBM-p1deg'

# Variables to process
# Note, 'tas' is not used in the merged output file here, but is produced for use in FUSE
varnames = ['pr','tasmax','tasmin','tas']
varfiles = [adata+'isimip_bc/obs/MSWEP2-2/global_daily_010deg_reformat/MSWEP2-2_010deg_*.nc',
adata+'ERA5_regrid/day_corrected/ERA5_tasmax_day_*.nc',
adata+'ERA5_regrid/day_corrected/ERA5_tasmin_day_*.nc',
adata+'ERA5_regrid/day_corrected/ERA5_tas_day_*.nc']
# CDO command to change units (e.g. K to degC)
unitmods = [[],['subc,273.15'],['subc,273.15'],['subc,273.15']]
# ncatted string to change units metadata
unitattmods = ['','units,tasmax,m,c,degC','units,tasmin,m,c,degC','units,tas,m,c,degC']
# name of input variable in netcdf files (to be renamed)
var_orig = ['','mx2t_NON_CDM','mn2t_NON_CDM','']
cat_files = []

###########################################################################################
# Loop over variables and process
for i,var in enumerate(varnames):
	print(var)
	print(varfiles[i])
	flist = []
	var_prefix = os.path.basename(varfiles[i])[:-4]
	fcat = os.path.join(outdir,var_prefix+'19790101-20171031_'+boxname+'.nc')


	if not os.path.exists(fcat):
		###################################################################################
		# Loop over files (for different timeslices)
		for fin in sorted(glob.glob(varfiles[i])):
			fname = os.path.basename(fin)
			fout = os.path.join(processdir,fname) # Temporary file
			if not os.path.exists(fout):

				# Select region and convert units
				cdo_cmd = ['cdo','-L']+unitmods[i]+['-sellonlatbox,'+lonlatbox,fin,fout]
				print(cdo_cmd)
				ret = subprocess.call(cdo_cmd)
				if not ret==0:
					raise Exception('Error with cdo command')

				# Rename variable if needed
				if var_orig[i] != '':
					nco_cmd = ['ncrename','-v',var_orig[i]+','+var,fout]
					print(nco_cmd)
					ret = subprocess.call(nco_cmd)
				if not ret==0:
					raise Exception('Error with ncrename command')

				# rename units to degC or mm/day from K or kg/m2/s
				if unitattmods[i]!='':
					nco_cmd = ['ncatted','-a',unitattmods[i],fout]
					print(nco_cmd)
					ret = subprocess.call(nco_cmd)
				if not ret==0:
					raise Exception('Error with ncatted command')

			flist.append(fout)

		###################################################################################
		# Concatenate files for full timeseries
		cdo_cat = ['cdo','cat']+flist+[fcat]
		print(' '.join(cdo_cat))
		ret = subprocess.call(cdo_cat)
		if not ret==0:
			raise Exception('Error with cdo cat command')
		else: # clean up temporary files
			for f in flist:
				os.remove(f)

	# Append concatenated file to list
	cat_files.append(fcat)

###########
# Manual step to convert MSWEP lat/lon to single precision so the merging works:
# ncap2 -s 'lon=float(lon);lat=float(lat)' f f_v2

###########################################################################################
# Merge variables and select times for metsim state file
f_state = os.path.join(outdir,'metsim-state_day_MSWEP-ERA5_19790101-19790331_'+boxname+'.nc')
if not os.path.exists(f_state):
	cdo_cmd2= ['cdo','-L','merge','-seldate,1979-01-01,1979-03-31','-invertlat',cat_files[1],'-seldate,1979-01-01,1979-03-31','-invertlat',cat_files[2],'-seldate,1979-01-01,1979-03-31',cat_files[0][:-3]+'_v2.nc',f_state]
	print(cdo_cmd2)
	ret = subprocess.call(cdo_cmd2)
	if not ret==0:
		raise Exception('Error with cdo merge command')

###########################################################################################
# Merge variables and select times for metsim forcing file
f_forcing = os.path.join(outdir,'metsim-forcing_day_MSWEP-ERA5_19790401-20171031_'+boxname+'.nc')
if not os.path.exists(f_forcing):
	cdo_cmd3 = ['cdo','-L','merge','-seldate,1979-04-01,2017-10-31','-invertlat',cat_files[1],'-seldate,1979-04-01,2017-10-31','-invertlat',cat_files[2],'-seldate,1979-04-01,2017-10-31',cat_files[0][:-3]+'_v2.nc',f_forcing]
	print(cdo_cmd3)
	ret = subprocess.call(cdo_cmd3)
	if not ret==0:
		raise Exception('Error with cdo merge command')
