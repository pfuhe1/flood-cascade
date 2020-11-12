# THis script needs to be run from the qgis console otherwise can't access the QGIS processing library
#
# Run by using this command:
# exec(open('/home/pu17449/src/flood-cascade/metforcings_prep/split_catchments_p1deg.py').read())
import os,sys,glob
from qgis.core import (QgsVectorLayer,)
import processing # QGIS processing library

# Set paths
indir = '/home/pu17449/data2/Discharge Data/GRDC_daily_global/shapefiles'
outdir = '/home/pu17449/data2/Discharge Data/GRDC_daily_global/shapefiles_processing_p1deg'

override=True

if not os.path.exists(outdir):
	os.mkdir(outdir)

# Open grid lines vector layer
lines_path = os.path.join(QgsProject.instance().homePath(), "..", "data2", "MSWEP2-2_010deg/lines_p1deg.shp")
#lines_path = "/home/pu17449/data2/EWEMBI-grid/grid_lines.shp"
lines = QgsVectorLayer(lines_path, "Lines", "ogr")
if not lines.isValid():
    print("Lines Layer failed to load!")
else:
	print('Opened Lines',lines_path)


# Loop over catchments and split them up
for catchment in glob.glob(os.path.join(indir,'*_smoothed*.shp')):
	try:
		clayer = QgsVectorLayer(catchment, "catchment", "ogr")
		cname = os.path.basename(catchment[:-4])
		if not clayer.isValid():
			print("Catchment layer failed to load!",catchment)
		else:
			# First we split catchment using lines
			#algorithmHelp('native:splitwithlines')
			outfile = os.path.join(outdir,cname+'_p1degsplit.shp')
			if override or not os.path.exists(outfile):
				print('Splitting catchment',cname)
				cmd_dict = {'INPUT':clayer,'LINES':lines,'OUTPUT':outfile}
				split_catchment = processing.run('native:splitwithlines',cmd_dict)
				print('Done Split!')

			# Next we want to calculate area fraction of each subset:
			#processing.algorithmHelp('qgis:exportaddgeometrycolumns')
			outfile2 = os.path.join(outdir,cname+'_p1degsplit_areas.shp')
			if override or not os.path.exists(outfile2):
				print('Calculating geometry areas')
				cmd2_dict = { 'INPUT' : outfile, 'CALC_METHOD' : 2, 'OUTPUT' : outfile2 }
				processing.run('qgis:exportaddgeometrycolumns',cmd2_dict)
				print('Done Area Calc!')

			# Calculate centroids (to know which part is in which cell)
			outfile3 = os.path.join(outdir,cname+'_p1degsplit_centroids.shp')
			if override or not os.path.exists(outfile3):
				print('calculating approx location within grids')
				cmd3_dict = { 'INPUT' : outfile2, 'ALL_PARTS' : False, 'OUTPUT' : outfile3 }
				processing.run('native:centroids',cmd3_dict)
				print('Done centroids calculation')

			continue

			# For debugging errors

			# Load outfile3
			catchment3 = QgsVectorLayer(outfile3, "Catchment areas added", "ogr")
			features = catchment3.getFeatures()
			points = []
			areas = []
			for feature in features:
				#print("Feature ID: ", feature.id())
				geom = feature.geometry()
				geomSingleType = QgsWkbTypes.isSingleType(geom.wkbType())
				if geom.type() == QgsWkbTypes.PointGeometry:
					# the geometry type can be of single or multi type
					if geomSingleType:
						x = geom.asPoint()
						#print("Point: ", x)
						points.append(x)
				areas.append(feature['area_2']*1e-6)

				#attrs = feature.attributes()
				# attrs is a list. It contains all the attribute values of this feature
				#print(attrs)

			print('catchment gridcells:',cname)
			print('points (lon,lat)',points)
			print('areas (km2)',areas)
	except Exception as e:
		print('Error processing file:',catchment)
		print(e)
