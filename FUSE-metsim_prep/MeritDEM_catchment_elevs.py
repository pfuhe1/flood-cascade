import gdal
import numpy as np
import os,glob,shutil
import netCDF4
# Use osgeo (gdal library) to load shapefile
from osgeo import ogr

################################################################################
# Set up paths etc.

# Path for DEM files (cropped for each catchment)
dem_dir = '/export/anthropocene/array-01/pu17449/parameter_transfer_datasets/MeritDEM_catchments'

# output elev_band file for FUSE:
elevbands_dir = '/export/anthropocene/array-01/pu17449/parameter_transfer_datasets/fuse_catchment_elevs'

if not os.path.exists(elevbands_dir):
	os.mkdir(elevbands_dir)

# Shapefile representing catchments used. Requires attributes:
# ID2: format 'GRDC XXXXXXX' specifying the GRDC catchment ID
# StatLat, StatLon: approximate lat-lon for the catchment.
catchments_metafile = '/export/anthropocene/array-01/pu17449/parameter_transfer_datasets/Beck2015_supplementary_data/Catchments.shp'

# Load catchments metadata (for lat/lon of catchment)
driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = driver.Open(catchments_metafile, 0)
layer = dataSource.GetLayer()

# Catchments have a single lat/lon point
n_outx = 1
n_outy = 1
# x,y are first index in lat/lon arrays
y=0
x=0

# Number of elevation bands to use for FUSE (user choice)
nbands = 5
# Note Error occured for some catchments, for those regions nbands = 5 was used

overwrite = False

# Change numpy setting to raise errors instead of printing warnings
np.seterr(all='raise')

################################################################################
# initialise output arrays:
elev = np.zeros([n_outy,n_outx],dtype=np.int)
mask = np.ones([n_outy,n_outx],dtype=np.int64)
frac = np.ones([n_outy,n_outx])

band_elev = np.zeros([nbands,n_outy,n_outx])
band_frac = np.zeros([nbands,n_outy,n_outx])

for dem_file in glob.glob(os.path.join(dem_dir,'dem_*.tif')):
	grdcid = os.path.basename(dem_file).split('_')[1][:-4]

	# Output file name
	elevbands_file = os.path.join(elevbands_dir,'catchment_'+grdcid+'_elev_bands.nc')

	if overwrite or not os.path.exists(elevbands_file):
		print('GRDCID:',grdcid)
		try:

			# find lat/lon from metadata shapefile
			latvals = None
			lonvals = None
			layer.ResetReading()
			for feature in layer:
				fid = feature.GetField('ID2')
				if fid[:4]=='GRDC' and grdcid == fid[5:]:
					latvals = [feature.GetField('StatLat')]
					lonvals = [feature.GetField('StatLon')]
			if latvals is None:
				raise Exception('Error, cant find lat/lon values. GRDCID: '+grdcid)

			################################################################################
			# Read in data to calculate grid box elevations
			f=gdal.Open(dem_file)
			rasterband = f.GetRasterBand(1)
			dem_arr = rasterband.ReadAsArray()
			missingval = -9999. # Todo, could read this from the file
			box = np.ma.masked_values(dem_arr,missingval).compressed()
			npoints = len(box)*1.0
			# Mean elevation
			elev[y,x] = box.mean()
			# Now split elevation into bands
			min_elev = box.min()
			max_elev = box.max()
			bands = np.linspace(min_elev,max_elev,num=nbands+1)
			#print(bands)
			for i in range(nbands):
				subset = box[np.logical_and(box> bands[i],box<bands[i+1])]
				band_elev[i,y,x] = subset.mean()
				if i == nbands-1: # Handle last band separately
					# Hack to make sure fracs sum to 1
					# (otherwise there is a small rounding error)
					band_frac[i,y,x] = 1- band_frac[:-1,y,x].sum()
				else:
					band_frac[i,y,x] = len(subset)/npoints
			print('debug area_sum',band_frac[:,y,x].sum())

		except Exception as e:
			print('Error',e)
			print('min/max elev',min_elev,max_elev)
			print('bands',bands)
			print('subset',i,subset)
			continue
			#raise

		###############################################################################
		# Write output file for FUSE elev bands
		# Follows example of file from fuse_grid example: cesm1-cam5_elev_bands.nc

		with netCDF4.Dataset(elevbands_file,'w') as f_out:
			f_out.createDimension('latitude',n_outy)
			f_out.createVariable('latitude',np.float,('latitude'))
			f_out.variables['latitude'].standard_name = "latitude"
			f_out.variables['latitude'].long_name = "latitude"
			f_out.variables['latitude'].units = "degrees_north"
			f_out.variables['latitude'].axis = "Y"
			f_out.variables['latitude'][:] = latvals

			f_out.createDimension('longitude',n_outx)
			f_out.createVariable('longitude',np.float,('longitude'))
			f_out.variables['longitude'].standard_name = "longitude"
			f_out.variables['longitude'].long_name = "longitude"
			f_out.variables['longitude'].units = "degrees_east"
			f_out.variables['longitude'].axis = "X"
			f_out.variables['longitude'][:] = lonvals

			f_out.createDimension('elevation_band',nbands)
			f_out.createVariable('elevation_band',np.int32,('elevation_band'))
			f_out.variables['elevation_band'].long_name = 'elevation_band'
			f_out.variables['elevation_band'].units = '-'
			f_out.variables['elevation_band'][:] = range(1,nbands+1)

			f_out.createVariable('area_frac',np.float32,('elevation_band', 'latitude', 'longitude'), fill_value=-9999)
			f_out.variables['area_frac'].units = "-" ;
			f_out.variables['area_frac'].long_name = "Fraction of grid cell covered by each elevation band."
			f_out.variables['area_frac'][:] = band_frac

			f_out.createVariable('mean_elev',np.float32,('elevation_band', 'latitude', 'longitude'),fill_value=-9999)
			f_out.variables['mean_elev'].units = "m asl" ;
			f_out.variables['mean_elev'].long_name = "Mean elevation of elevation band." ;
			f_out.variables['mean_elev'][:] = band_elev

			f_out.createVariable('prec_frac',np.float32,('elevation_band', 'latitude', 'longitude'), fill_value=-9999)
			f_out.variables['prec_frac'].units = "-" ;
			f_out.variables['prec_frac'].long_name = "Fraction of cell precipitation that falls on each elevation band." ;
			f_out.variables['prec_frac'].comment = 'PFU note: prec frac is just area_frac for now, need to implement scaling/lapse rate with elevation'
			f_out.variables['prec_frac'][:]=band_frac
