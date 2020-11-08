import pickle,os,glob
from netCDF4 import Dataset
import numpy as np

#########################################################################
# Function to write out catchment forcing to ncfile
#
def write_nc_catchment(fpath_out,f_in,lat,lon,catchment_vars):
	with Dataset(fpath_out,'w') as f_out:

		f_out.createDimension('lat',1)
		f_out.createVariable('lat',np.float,('lat'))
		f_out.variables['lat'].standard_name = "latitude"
		f_out.variables['lat'].long_name = "latitude"
		f_out.variables['lat'].units = "degrees_north"
		f_out.variables['lat'].axis = "Y"
		f_out.variables['lat'][:] = lat

		f_out.createDimension('lon',1)
		f_out.createVariable('lon',np.float,('lon'))
		f_out.variables['lon'].standard_name = "longitude"
		f_out.variables['lon'].long_name = "longitude"
		f_out.variables['lon'].units = "degrees_east"
		f_out.variables['lon'].axis = "X"
		f_out.variables['lon'][:] = lon

		for var in ['elev','mask','frac','area']:
			invar = f_in.variables[var]
			f_out.createVariable(var,invar.dtype,('lat','lon'))
			for att in invar.ncattrs():
				f_out.variables[var].setncattr(att,invar.getncattr(att))
			f_out.variables[var][:] = catchment_vars[var]

#############################################################################
# Input paths

# Pickle file calculated by catchment_forcings_part2 script
f_pkl2 = '/home/pu17449/data2/Discharge Data/GRDC_daily_global/catchments_p1deggrid_indices.pkl'
# Gridded elevations on grid from 60S-85N
fpath_in = '/home/pu17449/data2/metsim_data/domain_global-tiled_p1deg.nc'
# Output directory
outdir = '/home/pu17449/data2/metsim_data/catchment_domains_p1deg/'

if not os.path.exists(outdir):
	os.mkdir(outdir)

#############################################################################
# Open pickle file describing catchment indices
#
indices_set = set()
with open(f_pkl2,'rb') as fdata:
	points_dict = pickle.load(fdata)
	areas_dict = pickle.load(fdata)
	# Indices dict is a dictionary of catchments, containing a list of points used
	# NOTE: indices reference global grid
	indices_dict = pickle.load(fdata)

#############################################################################
# Read in input data
with Dataset(fpath_in,'r') as f_in:
	print(f_in.variables.keys())
	elev = f_in.variables['elev'][:]
	gridlat = f_in.variables['lat'][:]
	gridlon = f_in.variables['lon'][:]

	#############################################################################
	# Loop over catchmentsn and calculate forcings
	print('dictlen',len(indices_dict))
	for grdcid,indices in indices_dict.items():
		print('id',grdcid)
		catchment_elev = 0.
		catchment_area = 0.
		areas = areas_dict[grdcid]
		if len(indices)>0:
			for i,pt in enumerate(indices):
			# Add offset for 60s-85n grid used as input
				glat = pt[0] - 50
				glon = pt[1]
				catchment_elev += elev[glat,glon]*areas[i]*1e6
				catchment_area += areas[i]*1e6
			# normalise weighted sum
			catchment_elev = catchment_elev/catchment_area

			print('elev',catchment_elev)
			# Rough lat/lon
			lat,lon = points_dict[grdcid][0]
			lat2 = gridlat[glat]
			lon2 = gridlon[glon]
			print('check lat,lon:',lat,lon,lat2,lon2)
			catchment_vars = {'elev':catchment_elev,'mask':1,'frac':1,'area':catchment_area}

			fpath_out = os.path.join(outdir,str(grdcid)+'_domain.nc')
			write_nc_catchment(fpath_out,f_in,lat,lon,catchment_vars)
		else:
			print('ERROR, no points for catchment: '+str(grdcid))
