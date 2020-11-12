# lisflood_setup_bankfullQ.py
# Takes input from mizuRoute maximum discharge at a particular percentile/return period (from 'calc_q_returnperiod.py')
# Also requires running first lisflood_discharge_inputs_qgis.py to produce shapefiles (streamnet)
#
# Peter Uhe
# 2020/01/28
#

# Load modules
import os,sys,glob,pickle,shutil,socket
import numpy as np
from netCDF4 import Dataset,num2date
from osgeo import ogr
import csv
import datetime
from write_bdy_bci import write_bdy, write_bci_v2
import subprocess

import warnings
warnings.filterwarnings("error")

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

########################################################################################
# Function to calculate index of point,
# Inputs are point latitude and longitude coordinates and list of latitude,longitude pairs
#
def find_closest_list(p_lat,p_lon,lat_coord,lon_coord):
	ny=lat_coord.shape[0]
	nx=lon_coord.shape[0]
	min_dist=100000000
	minx=-1
	miny=-1
	for j in range(ny):
		dist=(lat_coord[j]-p_lat)**2 + (lon_coord[j]-p_lon)**2
		if dist<min_dist:
			miny=j
			min_dist=dist
	print(p_lat,p_lon,'min dist',min_dist)
	return miny,min_dist

########################################################################################
# Function to return true if point (xx,yy) is within bounds xbounds,ybounds
def in_bounds(xx,yy,xbound,ybound):
	return xbound[0]<xx and xx<xbound[1] and ybound[0]<yy and yy<ybound[1]


########################################################################################
# Function to check if point is not in list (also ignores -1)
def point_notinlist(point,points_list):
	return (point!=-1 and point not in points_list)

########################################################################################
# Little function to format dates into a string nicely
def format_date(date):
	return date.strftime('%Y-%m-%d')


###############################################################################################
# Start of main script
###############################################################################################
# Set up Shapefile reading driver
driver = ogr.GetDriverByName("ESRI Shapefile")


# Parameters for this domain
###############################################################################################
# extent (for gdal) xmin, ymin, xmax, ymax
# This is for the entire GBM basin
# 3 second dataset
extent = [73.477916667,22.456250000,97.562083333,31.276250000]
res = 0.000833
resname = '3sd8'

# 9 second dataset
#extent = [73.477916667,22.456250000,97.562916667,31.276250000]
#res = 0.0025 # in degrees
#resname = '9sd8'

extentstr = str(extent)[1:-1]

# Specify directories
###############################################################################################

runname = 'GBM-p1deg_MSWEP2-2-ERA5'
#resname = '9sd8'
#res = 0.0025

host = socket.gethostname()
if host[:3]=='bp1':
	# File, storing network attributes
	f_network        = '/home/pu17449/src/mizuRoute/route/ancillary_data/MERIT_mizuRoute_network_meta.nc'
	f_discharge      = '/work/pu17449/mizuRoute/output/q_bankfull_'+runname+'_50ile.nc'
	# Folder with full GBM stream network (output bankfullq files here)
	streamdir        = '/work/pu17449/lisflood/bankfullq_'+resname


# Specify file paths
###############################################################################################
# Shapefiles containing up/downstream points in each river segment/link
vertices_file = os.path.join(streamdir,'strn_network_'+resname+'_vertices.shp')
f_ntdstream    = os.path.join(streamdir,'strn_network_'+resname+'_acc_next-to-downstream.shp')
f_upstream  = os.path.join(streamdir,'strn_network_'+resname+'_acc_upstream.shp')

# Output files
bankfull_shp = os.path.join(streamdir,'strn_network_'+resname+'_'+runname+'_bankfullq.shp')
bankfull_tif = os.path.join(streamdir,'strn_network_'+resname+'_'+runname+'_bankfullq.tif')

###############################################################################################
# Main script: Get points for discharge in river network

# Dict of attributes for each reach/link
# x,y coordinates
points_up  = {}
#points_ntd = {}
#points_ds  = {}
# Upstream and downstream link numbers
dslinks     = {}
uslinks1    = {}
uslinks2    = {}
# Accumulation values
acc_up = {} # upstream accumulation
acc_ntd    = {} # accumulation at next-to-downstream point in reach
# Number of vertices in reach will be equal to the index of the upstream point
index_up  = {}
# Length of stream segment (m)
streamlen = {}
# Slope of reach
#slope      = {}
q_downstream = {}
q_upstream   = {}

nomap = []

# First load upstream points shapefile
dataSource = driver.Open(f_upstream, 0)
points = dataSource.GetLayer()
print('Opened upstream points in domain',f_upstream)
for feature in points:
	#print(feature.keys())
	# Get attributes:
	link = feature.GetField('LINKNO')
	dslinks[link] = feature.GetField('DSLINKNO')
	uslinks1[link] = feature.GetField('USLINKNO1')
	uslinks2[link] = feature.GetField('USLINKNO2')
	acc_up[link] = feature.GetField('ACC_1')
	index_up[link] = feature.GetField('vertex_ind') # index of vertex within link
	streamlen[link] = feature.GetField('Length')
	#slope[link] = feature.GetField('Slope')
	# Get location of point
	xx,yy,zz = feature.geometry().GetPoint()
	points_up[link] = [xx,yy]

links = list(acc_up.keys())

# Load (next-to) downstream points shapefile
dataSource = driver.Open(f_ntdstream, 0)
points = dataSource.GetLayer()
print('Opened next-to-downstream points',f_ntdstream)
for feature in points:
	# Get attributes
	link = feature.GetField('LINKNO')
	acc_ntd[link] = feature.GetField('ACC_1')
	# Get location of point
	#xx,yy,zz = feature.geometry().GetPoint()
	#points_ntd[link] = [xx,yy]



###############################################################################################
# Load Mizuroute river network (latitude and longitude of upstream points of each segment)
# Also load segId (which is not necessarily the same as the index of the segment)
seg_index_map = {}
with Dataset(f_network,'r') as f:
	segids_mizu = f.variables['segId'][:] # mizuroute ID of each river segment
	lats   = f.variables['lat_up'][:] # upstream lat of each river segment
	lons   = f.variables['lon_up'][:] # upstream lon of each river segment
	# Work out the mizuroute segment index corresponding to lisflood river network 'linkno'
	# Match upstream coordinates for each river segment
	#
	seg_index_map = {}
	index_seg_map = {}
	for link,point in points_up.items():
		xx,yy = point
		index,mindist = find_closest_list(yy,xx,lats,lons)
		if mindist < 0.0002:
			seg_index_map[link] = index
			if not index in index_seg_map:
				index_seg_map[index] = link
				print('mizu,taudem link:',segids_mizu[seg_index_map[link]],link)
			else:
			#	raise Exception('Error, multiple links mapping to same mizuroute segment')
				print('Warning, multiple links mapping to same mizuroute segment',index)
				otherlink = index_seg_map[index]
				print('Linknos,len:',otherlink,streamlen[otherlink],link,streamlen[link])
				if streamlen[otherlink]>streamlen[link]:
					del(seg_index_map[link])
					nomap.append(link)
				else:
					del(seg_index_map[otherlink])
					nomap.append(otherlink)
		else:
			print('Error finding link in mizuroute:',link,xx,yy)
		sys.stdout.flush()

###############################################################################################
# Main script: read in bankful discharge for each link
#
fname = os.path.basename(f_discharge)
print(fname)
discharge_dict = {}
# Read mizuroute output file for discharge values
with Dataset(f_discharge,'r') as f_in:
	# Assume the segment ids match those in the mizuroute network
	routed_vals = f_in.variables['IRFroutedRunoff'][0,:]
	print('shape',routed_vals.shape)

# Loop over reaches: set downstream q
for link in links:
	dslink  = dslinks[link]
	uslink1 = uslinks1[link]
	uslink2 = uslinks2[link]
	length = streamlen[link]
	# Set downstream q for this link
	#if length == 0.0: # this link occurs where three tributaries join
	#	q_downstream[link] = 0.0
	if link in nomap:
		# We dont have a mizuroute value for this link (this includes zero length links)
		# Use value from upstream reaches added together
		q_downstream[link] = 0.0
		if uslink1 in seg_index_map:
			q_downstream[link]+=routed_vals[seg_index_map[uslink1]]
		if uslink2 in seg_index_map:
			q_downstream[link]+=routed_vals[seg_index_map[uslink2]]
	else: # Normal link
		index = seg_index_map[link]
		q_downstream[link] = routed_vals[index]


# Loop over reaches: work out upstream q
for link in links:
	# If headwater catchment, also set upstream q for this link
	# Work this out based on fraction of accumulation at upstream point
	uslink1 = uslinks1[link]
	uslink2 = uslinks2[link]
	if uslink1 == -1 and uslink2 == -1:
		upacc = acc_up[link]
		downacc = acc_ntd[link]
		q_upstream[link] = q_downstream[link]*(upacc/downacc)
	else:
		# set q_upstream equal to q_downstream (to have constant q across whole reach)
		# Unless this is greater than the sum of the upstream links, then set to this sum
		q_upstream[link] = min(q_downstream[link], q_downstream[uslink1]+q_downstream[uslink2])

# First copy vertices data to bankfullq shapefile
for ext in ['.dbf','.prj','.shp','.shx','.qpj']:
	shutil.copy(vertices_file[:-4]+ext,bankfull_shp[:-4]+ext)

# Load All Vertices points shapefile
dataSource = driver.Open(bankfull_shp, 1)
points = dataSource.GetLayer()
fields = ['LINKNO', 'DSLINKNO', 'USLINKNO1', 'USLINKNO2', 'DSNODEID', 'strmOrder', 'Length', 'Magnitude', 'DSContArea', 'strmDrop', 'Slope', 'StraightL', 'USContArea', 'WSNO', 'DOUTEND', 'DOUTSTART', 'DOUTMID', 'vertex_pos', 'vertex_ind', 'vertex_par', 'vertex_p_1', 'distance', 'angle']

print('Opened All Vertices',bankfull_shp)
newfield = ogr.FieldDefn('bankfullq', ogr.OFTReal)
points.CreateField(newfield)

# Loop over all vertices, interpolate q between upstream and downstream
for feature in points:
	# Get attributes
	link = feature.GetField('LINKNO')
	dslink = feature.GetField('DSLINKNO')
	# Vertex indexing:
	# Upstream point has index = index_up
	# Downstream point has index 0
	vertex_index = feature.GetField('vertex_ind')
	reachlen = index_up[link]-1 # Length/dist in number of cells (Dont include the downstream point (index 0))
	reachdist = streamlen[link] # length/distance in metres

	# Get upstream and downstream q for this reach
	q_up = q_upstream[link]
	q_down = q_downstream[link]

	# Set bankfull q
	try:
		if vertex_index==0 :
			# Dont assign q for downstream vertices (except for river outlet)
			if dslink == -1:
				thisq = q_down
			else:
				points.DeleteFeature(feature.GetFID())
				continue

		else:
			if reachlen == 0:
				if reachdist == 0.0: # Only a single point skip this
					points.DeleteFeature(feature.GetFID())
					continue
				else: # Set q to q_down
					thisq = q_down
			else:
				# Interpolate between upstream and downstream q
				thisq = q_up + (q_down-q_up) * (reachlen-vertex_index+1)/reachlen
		# Now Write it into the shapefile field
		#print('Link',link,vertex_index,reachlen)
		#print('Q',q_up,q_down,thisq)
		feature.SetField('bankfullq',float(thisq))
		points.SetFeature(feature)
		feature = None # Clear feature so it writes out the data
	except Exception as e:
		print('Error setting q:',e)
		print('Link,Q',link,vertex_index,reachlen,q_up,q_down)

# Finally clean up this shapefile by removing all other fields
#p1 = points.GetNextFeature()
#for field in feature.keys():
#for i in range(feature.GetFieldCount()-1,-1,-1): # Loop from highest to lowest index
#	field = feature.GetField(i)
#	if not field=='xcoord' or field=='ycoord' or field=='bankfullq':
#		points.DeleteField(i)

# Try repacking points layer to ensure deleted features are removed
#dataSource.ExecuteSQL('REPACK ' + points.GetName())
# To write out file, delete dataSource
points = None
dataSource = None


# Finally rasterize the bankfullq
gdal_cmd = ['gdal_rasterize','-a','bankfullq','-tr',str(res),str(res),'-a_nodata','0.0','-te']+extentstr.split(',')+['-ot','Float32','-of','GTiff','-co','COMPRESS=DEFLATE',bankfull_shp,bankfull_tif]
print(' '.join(gdal_cmd))
subprocess.call(gdal_cmd)
