# Peter Uhe 7/10/2019
# This script is part of a series of scripts and needs to be run after 'catchment_forcings_part1_p1deg.py'
# This script takes the pickle file containing points in each catchment, and matches them to grid indices in the grid_path file.

# THis script needs to be run from a python console
# (the qgis console doesn't work with netCDF4)


import os,sys,glob,pickle
import numpy as np
from netCDF4 import Dataset

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


# Set paths
grid_path = '/home/pu17449/data2/MSWEP2-2_010deg/MSWEP2-2_010deg_197901.nc'
f_pkl = '/home/pu17449/data2/Discharge Data/GRDC_daily_global/catchments_p1deggrid_indices.pkl'

# open ewembi grid
with Dataset(grid_path,'r') as f:
	lat = f.variables['lat'][:]
	lon = f.variables['lon'][:]

# load pickle file:
with open(f_pkl,'rb') as fdata:
	points_dict = pickle.load(fdata)
	areas_dict = pickle.load(fdata)
indices_dict = {}

# loop over points
for grdcid,points in points_dict.items():
	print('catchment:',grdcid)
	indices = []
	for p_lat,p_lon in points:
		p_index = find_closest_1d_v2(p_lat,p_lon,lat,lon)
		indices.append(p_index)
		print('check latlon- point:',p_lat,p_lon,',grid:',lat[p_index[0]],lon[p_index[1]])
	indices_dict[grdcid]=indices

# write out data
with open(f_pkl,'wb')as fdata:
	pickle.dump(points_dict,fdata,-1)
	pickle.dump(areas_dict,fdata,-1)
	pickle.dump(indices_dict,fdata,-1)
