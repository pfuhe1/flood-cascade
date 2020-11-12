# mizuRoute_rivernet.py
#
# Peter Uhe
# May 27 2019
#
# This script uses hydrography data from the MERIT-Hydro dataset.
# Uses stream network shapefiles calculated using TauDEM, with additional GIS processing
# and produces the required river network inputs for the mizuRoute river routing model
# E.g. https://mizuroute.readthedocs.io/en/develop/Input_data.html#river-network-data
#
#
# This script needs requires python libraries for gdal and netCDF4
#
#
###############################################################################################
# Input files:
#
# Some pre-processing was done in QGIS:
#
#Requires outputs from TauDEM: stren_w3d8.tif, stren_net3d8.shp and catchment mask (basins3.tif).
#
#Catchment Mask (processing in QGIS):

#1. Polygonise tif: e.g. gdal_polygonize.py <datadir>/basins3.tif <datadir>/basins3_poly.gpkg" -8 -b 1 -f "GPKG" None DN
#2. Select feature for desired catchment, 'save selected feature as'
#3. Fix geometries


#Catchments:
#1. Clip stren_w3d8.tif by catchment mask (calculated above)
#2. GDAL Polygonize: 'Name of the field to create='LINENO', use 8-connectedness

#Streams:
#1. Select 'extract specific vertices' (vertex 0 for upstream points)
#2. Clip result by catchment mask (calculated able).


###############################################################################################
#
# Load modules
import os,sys,argparse
import numpy as np
from netCDF4 import Dataset
from osgeo import ogr


if __name__ == '__main__':

	###############################################################################################
	# Parse command line arguments
	parser = argparse.ArgumentParser(description='Calculate mizuRoute river network ancillary file')
	parser.add_argument('-s','--streams',help = 'Input shapefile representing stream segments (clipped to desired region with points locations of upstream end of each segment)')
	parser.add_argument('-c','--catchs',help = 'Input shapefile representing catchments corresponding to streams. LINENO attribute corresponds to LINKNO of stream')
	parser.add_argument('-o','--outfile',help = 'Output mizuroute ancillary file representing river network (netcdf file)')
	parser.add_argment('--overwrite',default=False,action='store_true',help = 'Flag to specify overwriting existing output file')

	args = parser.parse_args()
	streams_path = args.streams
	catchs_path  = args.catchs
	fnetcdf      = args.outfile
	overwrite    = args.overwrite

	# Set up Shapefile reading driver
	driver = ogr.GetDriverByName("ESRI Shapefile")


	###############################################################################################

	if not os.path.exists(fnetcdf) or overwrite:

		if os.path.exists(fnetcdf):
			print('Warning, this script will overwrite output file',fnetcdf)

		###############################################################################################
		# Initialise variables
		#
		seg_ids=[]
		lengths = []
		slopes = []
		down_ids = []
		lat_up = []
		lon_up = []


		###############################################################################################
		# Process streams data
		#

		# First load streams
		dataSource = driver.Open(streams_path, 0)
		lines = dataSource.GetLayer()
		print('Opened Lines',streams_path)
		for feature in lines:
			print('id',feature.GetFID())
			seg_ids.append(feature.GetField('LINKNO'))
			lengths.append(feature.GetField('Length'))
			slopes.append(feature.GetField('Slope'))
			down_ids.append(feature.GetField('DSLINKNO'))
			xx,yy,zz = feature.geometry().GetPoint()
			lat_up.append(yy)
			lon_up.append(xx)

		seg_ids = np.array(seg_ids)
		lengths = np.array(lengths)
		slopes = np.array(slopes)
		down_ids = np.array(down_ids)
		lat_up = np.array(lat_up)
		lon_up = np.array(lon_up)

		###############################################################################################
		# Process catchments data
		#

		# Load Catchments
		dataSource = driver.Open(catchs_path, 0)
		lines = dataSource.GetLayer()
		seg_ids_hru=[]
		areas = []
		print('Opened Catchments',catchs_path)
		for feature in lines:
			seg_ids_hru.append(feature.GetField('LINENO'))
			a = feature.GetField('area')
			areas.append(a)
		print('HRUs',len(seg_ids_hru),len(areas))

		seg_ids_hru = np.array(seg_ids_hru,dtype=int)
		hru_ids = np.arange(len(seg_ids_hru),dtype=int)
		areas = np.array(areas)

		#########################################################

		with Dataset(fnetcdf,'w') as f_out:
			f_out.createDimension('segId',len(seg_ids))
			f_out.createVariable('segId',np.int32,('segId'))
			f_out.variables['segId'][:] = seg_ids

			f_out.createVariable('downsegId',np.int32,('segId'),fill_value=-999)
			f_out.variables['downsegId'][:] = down_ids

			f_out.createVariable('slope',np.float32,('segId'))
			f_out.variables['slope'][:] = slopes
			f_out.variables['slope'].units = '-' # m/m

			f_out.createVariable('length',np.float32,('segId'))
			f_out.variables['length'][:] = lengths
			f_out.variables['length'].units = 'm'

			f_out.createVariable('lat_up',np.float32,('segId'))
			f_out.variables['lat_up'][:] = lat_up
			f_out.variables['lat_up'].description = 'Latitude of upstream point of river segment'
			f_out.variables['lat_up'].axis = 'Y'

			f_out.createVariable('lon_up',np.float32,('segId'))
			f_out.variables['lon_up'][:] = lon_up
			f_out.variables['lon_up'].description = 'Longitude of upstream point of river segment'
			f_out.variables['lon_up'].axis = 'X'

			######################

			f_out.createDimension('hruId',len(hru_ids))
			f_out.createVariable('hruId',np.int32,('hruId'))
			f_out.variables['hruId'][:] = hru_ids

			f_out.createVariable('hrusegId',np.int32,('hruId'),fill_value=-999)
			f_out.variables['hrusegId'][:] = seg_ids_hru

			f_out.createVariable('area',np.float32,('hruId'))
			f_out.variables['area'][:] = areas
			f_out.variables['area'].units = 'm2'
	else:
		print('mizuroute ancillary file already exists, skipping',fnetcdf)
