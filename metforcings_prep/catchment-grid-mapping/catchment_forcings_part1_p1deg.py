# Peter Uhe 7/10/2019
# Takes sections of catchment shapefiles intersecting a regular grid, calculated by 'split_catchments_p1deg.py'
# Determine location and areas of each section, and saves to python pickle file

# THis script needs to be run from the qgis console otherwise doesn't seem to work...
# Run by using this command:
# exec(open('/home/pu17449/src/flood-cascade/metforcings_prep/catchment_forcings_part1_p1deg.py').read())

import os,sys,glob,pickle
from qgis.core import (QgsVectorLayer,)
import processing # QGIS processing library

# Set paths
centroids_path = '/home/pu17449/data2/Discharge Data/GRDC_daily_global/shapefiles_processing_p1deg'
f_pkl = '/home/pu17449/data2/Discharge Data/GRDC_daily_global/catchments_p1deggrid_indices.pkl'


# Initialise dictionaries
points_dict = {}
areas_dict = {}

inpath = os.path.join(centroids_path,'*_p1degsplit_centroids.shp')
print(inpath)
# Loop over catchments and split them up
for catchment_split in glob.glob(inpath):
	try:
		cname = os.path.basename(catchment_split[:-4])
		grdcid = int(cname.split('_')[5])
		print(grdcid)

		# Load file
		catchment3 = QgsVectorLayer(catchment_split, "Catchment areas added", "ogr")
		features = catchment3.getFeatures()
		points = []
		areas = []
		indices = []
		for feature in features:
			#print("Feature ID: ", feature.id())
			geom = feature.geometry()
			geomSingleType = QgsWkbTypes.isSingleType(geom.wkbType())
			if geom.type() == QgsWkbTypes.PointGeometry:
				# the geometry type can be of single or multi type
				if geomSingleType:
					x = geom.asPoint()
					p_lon,p_lat = x[0],x[1] # x has lon,lat
					points.append((p_lat,p_lon)) #points has lat,lon (swapped)
			areas.append(feature['area_2']*1e-6)


		# put list from subregions into dictionary by GRDC catchment ID
		points_dict[grdcid]=points
		areas_dict[grdcid] = areas

	except Exception as e:
		print('Error processing file:',catchment)
		print(e)

with open(f_pkl,'wb')as fdata:
	pickle.dump(points_dict,fdata,-1)
	pickle.dump(areas_dict,fdata,-1)
