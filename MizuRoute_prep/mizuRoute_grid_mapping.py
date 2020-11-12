# mizuRoute_grid_mapping.py
#
# Peter Uhe
# May 2019
#
# This script takes mizuRoute network, and shapefile representing parts of each catchment split by the runoff model grid (centroids file)
# Produces mapping file that calculates weights for each grid point in the runoff model, that are input for each hydrological response unit/ catchment in the river network.
# E.g. https://mizuroute.readthedocs.io/en/develop/Input_data.html#runoff-mapping-data

# This script requires gdal and netCDF4 libraries

# This script assumes that the features in the 'centroids' shapefile are in order of HRU-ID


########################################################################################
# Inputs
#
# MizuRoute river network file 'network_path': calculated by mizuRoute_rivernet.py

# Netcdf file containing lon-lat co-ordinates of grid for runoff model 'runoff_grid_path'

# Points shapefile, containing centre points 'centroids' of catchments split by grid lines of runoff model
# Pre-processing steps in QGIS:
# 1. 'create grid' (matching the grid of the runoff model above)  <datadir>/'MSWEP2-2_010deg
# 2. 'split by lines <datadir>/lfp-tools/merit_test1/stren_w3d8_GBM_v2area.shp' with grid created above -> <datadir>/lfp-tools/merit_test1/stren_w3d8_GBM_v2_splitp1deg.shp
# 3. Add geometry attributes (ellipsoidal method) -> <datadir>/lfp-tools/merit_test1/stren_w3d8_GBM_v2_splitp1deg-area.shp
# 4. Calculate centroids # Note, some parts of catchments (e.g. single grid boxes have been split off when they probably shouldn't. This version of the script just adds weights/area for any parts of HRUs that fall in each gridcell so this is fine.


###############################################################################################
#
# Load modules
import os,sys,argparse
import numpy as np
from netCDF4 import Dataset
from osgeo import ogr

########################################################################################
# Function to calculate index of closest grid box,
# Inputs are point latitude and longitude coordinates and grid latitude,longitude coordinates
#
def find_closest_1d_v2(p_lat,p_lon,lat_coord,lon_coord):
	ny=lat_coord.shape[0]
	nx=lon_coord.shape[0]
	min_dist=100000000
	minx=-1
	miny=-1
	for j in range(ny):
		dist=(lat_coord[j]-p_lat)**2
		if dist<min_dist:
			miny=j
			min_dist=dist
	min_dist=100000000
	#print 'plon',p_lon
	for i in range(nx):
		dist= (np.abs(lon_coord[i]-p_lon)%360)**2
		if dist<min_dist:
			minx=i
			min_dist=dist
	return miny,minx

###############################################################################################

if __name__=='__main__':

	# Parse command line arguments
	parser = argparse.ArgumentParser(description='Calculate mizuRoute river network ancillary file')
	parser.add_argument('-m','--mizunet',help = 'Mizuroute river network file calculated by mizuRoute_rivernet.py')
	parser.add_argument('-g','--runoffgrid',help = 'Netcdf file containing lon-lat co-ordinates of grid for runoff model')
	parser.add_argument('-c','--centroids',help = 'Shapefile with centroid points of intersections between catchments and runoff grid (see readme)')
	parser.add_argument('-o','--outfile',help = 'Output mizuroute runoff mapping file (netcdf)')
	parser.add_argment('--overwrite',default=False,action='store_true',help = 'Flag to specify overwriting existing output file')

	args = parser.parse_args()
	runoff_grid_path = args.runoffgrid
	network_path     = args.mizunet
	f_out            = args.outfile
	centroids_path   = args.centroids
	overwrite        = args.overwrite

	# Set up Shapefile reading driver
	driver = ogr.GetDriverByName("ESRI Shapefile")

	if not os.path.exists(f_out) or overwrite:

		if os.path.exists(f_out):
			print('Warning, this script will overwrite output file',f_out)

		########################################################################################
		# Open input files

		# open Runoff model grid
		with Dataset(runoff_grid_path,'r') as f:
			lat = f.variables['lat'][:]
			lon = f.variables['lon'][:]

		# Load seg_ids corresponding to HRUID's as a check
		with Dataset(network_path,'r') as f:
			seg_id_hru = f.variables['hrusegId'][:]

		# Open catchments shapefile:
		driver = ogr.GetDriverByName("ESRI Shapefile")
		dataSource = driver.Open(centroids_path, 0)
		layer = dataSource.GetLayer()

		########################################################################################
		# Initialise variables

		nhrus = 3527
		#nparts = 11824
		nparts = 46790 # not including zero size segments

		hruids = np.arange(nhrus,dtype=np.int32)
		noverlaps = np.zeros([nhrus],dtype=np.int32)

		weights = np.zeros([nparts])
		hruids_ofpart = np.zeros([nparts])
		ivals = np.zeros([nparts],dtype=np.int32)
		jvals = np.zeros([nparts],dtype=np.int32)

		# list of ij vals and partid used for each hruId
		# Used to prevent adding multiple entries for the same hruid-gridcell
		ijvals = [] # list ijvals for current HRU
		partids = [] # list partid for current HRU

		########################################################################################
		# Loop over (parts of) catchments/ HRUs
		#
		# The order of this is very important as mizuRoute assumes the hruids are in order
		# (it just happens that the 'centroids' shapefile maintains ordering of hruids)
		# If hruids are not in order, this script would need to be modified.
		i = 0
		hruid=0
		for catchment in layer:
			try:
				oldhruid = hruid
				segid = catchment.GetField('LINENO')
				hruid = int(catchment.GetField('fid')) - 1
				# If we have moved to a different HRU clear list of parts
				if oldhruid !=hruid:
					ijvals = []
					partids = []
					# Do so debug checking on the oldhruid:
					nparts = noverlaps[oldhruid]
					hweights = weights[i-nparts:i]
					print('DEBUG: HRU weights',hweights,hweights.sum())
					if not np.isclose(hweights.sum(),1.0):
						raise Exception('Error, weights for HRU must sum to 1')
				#if hruid == 3: # Debug to check data for first few HRUs
				#	print('Debug exit')
				#	break
				#print('Catchment, seg',i,segid)
				hru_area = catchment.GetField('area') # Area calculated over whole HRU
				part_area = catchment.GetField('area_2') # Area calculated over part of
				if part_area < 7000.: # This is a dodgey segment (less than one grid box). SKIP!
					continue
				weight = part_area/hru_area
				hrucheck = np.where(seg_id_hru==segid)[0][0] # Check HRU-ID
				if not hruid == hrucheck:
					raise Exception('Error, cant match hruid as expected')
				print('HRU-ID',hruid,'number',noverlaps[hruid])

				# Calculate index in runoff grid
				catch_lon,catch_lat,zz = catchment.geometry().GetPoint()
				jlat,ilon = find_closest_1d_v2(catch_lat,catch_lon,lat,lon)
				if (jlat,ilon) in ijvals:
					# We already have matched this grid box
					zz = ijvals.index((jlat,ilon))
					partid = partids[zz]
					print('debug, matched to part',partid)
					print('debug, adding to weight',weights[partid],weight)
					weights[partid] += weight
					if weights[partid]>1.0:
						raise Exception('Error, weight greater than 1',weights[partid])
				else:
					print('Part',i)
					weights[i] = weight
					hruids_ofpart[i] = hruid
					noverlaps[hruid]+=1
					jvals[i] = jlat+1 # plus one for indexing starting at 1
					ivals[i] = ilon+1
					print('lat,lon',jlat,ilon)
					partids.append(i)
					ijvals.append((jlat,ilon))
					i+=1


			except Exception as e:
				print('Error processing catchment:')
				print(e)
				raise

		# cut off end of lists
		weights = weights[:i]
		hruids_ofpart = hruids_ofpart[:i]
		jvals = jvals[:i]
		ivals = ivals[:i]
		nparts = i
		print('nparts used is',i)


		########################################################################################
		# Write out file
		#
		with Dataset(f_out,'w') as f_out:

			# Data on overlapping/ river network hrus split by the runoff model grid
			#
			f_out.createDimension('partID',nparts)
			f_out.createVariable('partID',np.int32,('partID'))
			f_out.variables['partID'][:] = range(nparts)

			f_out.createVariable('weight',np.float32,('partID'))
			f_out.variables['weight'][:] = weights

			f_out.createVariable('i_index',np.int32,('partID'))
			f_out.variables['i_index'][:] = ivals

			f_out.createVariable('j_index',np.int32,('partID'))
			f_out.variables['j_index'][:] = jvals

			f_out.createVariable('hruid_ofpart',np.int32,('partID'))
			f_out.variables['hruid_ofpart'][:] = hruids_ofpart
			f_out.variables['hruid_ofpart'].note = 'Note, this variable is extra metadata, not used by mizuRoute.\nThe hruid for each part is assumed by the order of the hruIds and the nOverlaps'

			# Data on HRU-IDs
			#
			f_out.createDimension('hruId',nhrus)
			f_out.createVariable('hruId',np.int32,('hruId'))
			f_out.variables['hruId'][:] = hruids

			f_out.createVariable('nOverlaps',np.int32,('hruId'))
			f_out.variables['nOverlaps'][:] = noverlaps
	else:
		print('Output file already exists, skipping',f_out)
