import pickle,os,glob,sys
from netCDF4 import Dataset,MFDataset,num2date
import numpy as np
from multiprocessing import Pool
import datetime

# Peter Uhe Nov 2019:
#
# Read in global datasets of metsim/fuse input variables
# Select data at indices for each catchment
# Write out data for each catchment to netcdf
#

def get_days(timedelta):
	return timedelta.days

#########################################################################
# Function to process catchment
#

def process_catchment(grdcid,f_in,fpath_out,indices,points,areas,var,vardata):
	if not os.path.exists(os.path.dirname(fpath_out)):
		os.mkdir(os.path.dirname(fpath_out))
	print(grdcid)

	catchment_area = 0.
	timevals = np.array([])

	# Now get the data for this catchment
	if len(indices)>0:
		# Rough lat/lon
		lat,lon = points[0]

		# Sum up area of catchment
		for i,pt in enumerate(indices):
			catchment_area += areas[i]

		# Load time and initialise data
		timevar = f_in.variables['time']
		units = timevar.units
		outunits = 'days since 1979-01-01 00:00:00'
		outbase = datetime.datetime(1979,1,1)
		timevals0 = timevar[:]
		timedates = num2date(timevals0,units)
		timevals = list(map(get_days,timedates-outbase))
		nt = len(timevals)
		catchment_data = np.ma.zeros([nt])

		latvals = f_in.variables['lat'][:]
		lonvals = f_in.variables['lon'][:]
		# Check location (first time only)
		if not os.path.exists(fpath_out):
			print('locationcheck',lat,lon,latvals[indices[0][0]],lonvals[indices[0][1]])


		for i,pt in enumerate(indices):
			#print(i,end=':')
			#sys.stdout.flush()
			# hack to reverse ordering of latitudes for era5 data
			if var[:3]=='tas':
				pt0 = len(latvals) - pt[0] -1
				print(latvals[pt0],end=',')
			else:
				pt0 = pt[0]
			catchment_data += vardata[:,pt0,pt[1]]*areas[i]
#			print('data',i,vardata[0,pt0,pt[1]])
			if np.ma.is_masked(vardata[0,pt0,pt[1]]):
				raise Exception('Error, invalid data')
		# normalise weighted sum

		catchment_data = catchment_data/catchment_area
		print('data_weighted',catchment_data[0],catchment_area)
		ncret = write_nc_catchment(fpath_out,f_in,lat,lon,var_map,catchment_data,timevals)
		print('ncret',ncret)
		return fpath_out
	else:
		print('ERROR, no points for catchment: '+str(grdcid))




#########################################################################
# Function to write out catchment forcing to ncfile
#
def write_nc_catchment(fpath_out,f_in,lat,lon,var_map,vals,timevals):

	if os.path.exists(fpath_out):
		mode = 'a'
	else:
		mode = 'w'
	with Dataset(fpath_out,mode) as f_out:
		if mode == 'w':
			print('Creating variables')
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

			timevar = f_in.variables['time']
			#print('intime',timevar)
			f_out.createDimension('time',0)
			f_out.createVariable('time',np.float,('time'))
			for att in timevar.ncattrs():
				f_out.variables['time'].setncattr(att,timevar.__getattr__(att))
			# Override units:
			outunits = 'days since 1979-01-01 00:00:00'
			f_out.variables['time'].setncattr('units',outunits)

			invar = f_in.variables[var_map[var]]
			f_out.createVariable(var,np.float32,('time','lat','lon'),fill_value=invar.__getattr__('_FillValue'))
			for att in invar.ncattrs():
				if att!='_FillValue':
					f_out.variables[var].setncattr(att,invar.__getattr__(att))

		print('Writing data')
		nt = len(f_out.variables['time'])
		f_out.variables['time'][nt:]=timevals
		print('written time')

		f_out.variables[var][nt:] = vals
		print('written',var)
	return True

#############################################################################
# Input paths
f_pkl2 = '/export/anthropocene/array-01/pu17449/FUSE_inputs/catchments_p1deg/catchments_p1deggrid_indices.pkl'
fcatchment_list = '/export/anthropocene/array-01/pu17449/FUSE_inputs/catchments_p1deg/GBM-p1deg_donorsbest3.txt'
timeperiod = '19790101-20171231'

# Set path of input files for each variable
var_files = {}
var_files['pr'] = '/export/anthropocene/array-01/pu17449/isimip_bc/obs/MSWEP2-2/global_daily_010deg_reformat/MSWEP2-2_010deg_??????.nc'
var_files['tas'] = '/export/anthropocene/array-01/pu17449/ERA5_regrid/day_corrected/ERA5_tas_day_*.nc'
var_files['tasmax'] = '/export/anthropocene/array-01/pu17449/ERA5_regrid/day_corrected/ERA5_tasmax_day_*.nc'
var_files['tasmin'] = '/export/anthropocene/array-01/pu17449/ERA5_regrid/day_corrected/ERA5_tasmin_day_*.nc'

# Map of standard names to names in input files (now they do not change)
#var_map= {'pr':'pr','tasmax':'tasmax','tasmin':'tasmin','tas':'tas'}
var_map = {'pr':'pr','tasmax':'mx2t_NON_CDM','tasmin':'mn2t_NON_CDM','tas':'tas'}

outdir = '/export/anthropocene/array-01/pu17449/FUSE_inputs/catchments_p1deg_v3'

#############################################################################
# Open pickle file describing catchment indices
#
indices_set = set()
with open(f_pkl2,'rb') as fdata:
	points_dict = pickle.load(fdata)
	areas_dict = pickle.load(fdata)
	# Indices dict is a dictionary of catchments, containing a list of points used
	# NOTE: indices reference global grid, need to double check the points against the lon/lat of input files
	indices_dict = pickle.load(fdata)


#############################################################################
# Read in input data
for var,fpath_in in var_files.items():
	firsttime = True

	#############################################################################
	# Loop over catchmentsn and calculate forcings
	print('dictlen',len(indices_dict))

	# First read all the data for each input file (hopefully faster although takes more memory)
	for infile in sorted(glob.glob(fpath_in)):
		print(os.path.basename(infile))
		f_in = Dataset(infile,'r')
		vardata = f_in.variables[var_map[var]][:]
		print('read data')
		# Then Loop over catchments and write each timeslice
		with open(fcatchment_list,'r') as f:
			for line in f:
				grdcid = int(line.strip())
				indices = indices_dict[grdcid]
				points  = points_dict[grdcid]
				areas   = areas_dict[grdcid]
				fpath_out = os.path.join(outdir,var,str(grdcid)+'_'+var+'_'+timeperiod+'.nc')

				process_catchment(grdcid,f_in,fpath_out,indices,points,areas,var,vardata)
		firsttime = False
		f_in.close()
		del(vardata)
