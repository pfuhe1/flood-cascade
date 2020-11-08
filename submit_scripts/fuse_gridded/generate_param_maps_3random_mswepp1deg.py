# Script to produce maps of FUSE parameter sets, based on parameter transfer / regionalisation
# Peter Uhe
# August 2019
#
# script '3random' takes the fields from the 3 choices of parameter sets and randomly picks one of them for each grid point.
# This needs to be run AFTER generate_param_maps_3choices_mswepp1deg.py
#
# Script to use anaconda and python3:
# On bluecrystal:
# module load languages/python-anaconda3-2019.03
# source activate petepy

import numpy as np
import pickle,glob,os,sys
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature

np.random.seed(33)

#########################################################################################
# Input files

dec = 900
setup_name = 'GBM-p1deg'
obsversion = 'MSWEP2-2-ERA5'
# for bc3:
#griddir = '/newhome/pu17449/data/fuse/fuse_GBM_v2-2/'
# for bc4:
#griddir = '/mnt/storage/scratch/pu17449/fuse/fuse_GBM_v2-2'
# for bp1:
griddir = '/work/pu17449/fuse/'+setup_name

#########################################################################################

param_maps = {}
param_longname = {}
param_units = {}
lats = None
lons = None
for choice in range(3):
	f_paramset = os.path.join(griddir,'output/'+setup_name+'_'+str(dec)+'_'+obsversion+'-calibrated'+str(choice+1)+'.nc')

	with Dataset(f_paramset,'r') as f:
		if lats is None:
			lats = f.variables['latitude'][:]
		if lons is None:
			lons = f.variables['longitude'][:]
		for varname,var in f.variables.items():
			shp = var.shape
			if len(shp)==2:
				try:
					param_maps[varname][choice,:] = var[:]
				except: # initialise param_map
					param_maps[varname] = np.ones([3,shp[0],shp[1]])*-9999
					param_maps[varname][choice,:] = var[:]
				if not varname in param_longname:
					param_longname[varname] = f.variables[varname].long_name
					param_units[varname] = f.variables[varname].units

			#	if varname == 'MAXFREE_1':
			#		print('maxfree',choice,param_maps[varname][choice,:])

for r in range(20):
	print('donor catchment random choice:',r)
	f_out = os.path.join(griddir,'output/'+setup_name+'_'+str(dec)+'_'+obsversion+'-calibrateRand'+str(r+1).zfill(4)+'.nc')
	choicearr = np.random.choice([0,1,2],size=shp)
	print('writing',f_out,'...')


	# Write grid to output
	with Dataset(f_out,'w') as f_out:

		f_out.createDimension('latitude',len(lats))
		f_out.createVariable('latitude',np.float,('latitude'))
		f_out.variables['latitude'].standard_name = "latitude"
		f_out.variables['latitude'].long_name = "latitude"
		f_out.variables['latitude'].units = "degrees_north"
		f_out.variables['latitude'].axis = "Y"
		f_out.variables['latitude'][:] = lats

		f_out.createDimension('longitude',len(lons))
		f_out.createVariable('longitude',np.float,('longitude'))
		f_out.variables['longitude'].standard_name = "longitude"
		f_out.variables['longitude'].long_name = "longitude"
		f_out.variables['longitude'].units = "degrees_east"
		f_out.variables['longitude'].axis = "X"
		f_out.variables['longitude'][:] = lons

		for param,field in param_maps.items():
			# Simple selection of output params
			#if param[0].isupper()
			var = f_out.createVariable(param,np.float,('latitude','longitude'),fill_value=-9999)
			var.long_name = param_longname[param]
			var.units = param_units[param]

			for j in range(shp[0]):
				for i in range(shp[1]):
					var[j,i] = field[choicearr[j,i],j,i]
