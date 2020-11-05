import gdal
import numpy as np
import os,glob,shutil,sys
import netCDF4
import matplotlib.pyplot as plt
import subprocess
import multiprocessing

# On BRIDGE servers, gdal is a bit broken, set 'GDAL_DATA='/opt/bridge/CentOS6-64/python/anaconda-5.0.1-2.7/pkgs/libgdal-2.1.0-0/share/gdal'

################################################################################
# Input paths

# Input directory containing DEM tiles (5degree increments)
dem_dir = '/export/anthropocene/array-01/pu17449/MeritDEM/' # files have format e.g. n25e090_dem.tif

out_dir1 = '/export/anthropocene/array-01/pu17449/parameter_transfer_datasets/MeritDEMSlope'
out_dir2 = '/export/anthropocene/array-01/pu17449/parameter_transfer_datasets/MeritDEMSlope_p1deg'
processdir = '/export/anthropocene/array-01/pu17449/parameter_transfer_datasets/tmp'

if not os.path.exists(out_dir2):
	os.mkdir(out_dir2)

numthreads = 6

def processtile(ftile,out_dir1,out_dir2,processdir):
	fname = os.path.basename(ftile)
	print(fname)

	slopefile = os.path.join(out_dir1,fname[:-7]+'slope.tif')
	if not os.path.exists(slopefile):
		# Convert to EPSG:3857 projection (units in of metres not degrees)
		tmpfile = os.path.join(processdir,fname[:-4]+'_epsg3857.tif')
		cmd = ['gdalwarp','-t_srs','EPSG:3857','-r','near','-of','GTiff',ftile,tmpfile]
		#print(cmd)
		ret1 = subprocess.call(cmd)

		# Calculate slope for tile
		cmd = ['gdaldem','slope','-of','GTiff','-b','1','-s','1.0','-co','COMPRESS=DEFLATE',tmpfile,slopefile]
		#print(cmd)
		ret2 = subprocess.call(cmd)
		# Clean up tmpfile
		os.remove(tmpfile)
	else:
		ret1=ret2=0

	# Regrid to 0.5deg
	# For 0.5deg regridding	
	slopefile2 = os.path.join(out_dir2,fname[:-7]+'slope_p1deg.tif')
	if not os.path.exists(slopefile2):	
		cmd = ['gdalwarp', '-t_srs', 'EPSG:4326', '-tr', '0.1', '0.1', '-r', 'average' ,'-of' ,'GTiff','-co','COMPRESS=DEFLATE', slopefile, slopefile2]
		#print(cmd)
		ret3 = subprocess.call(cmd)
	else:
		ret3=0

	if ret1!=0 or ret2!=0 or ret3!=0:
		print('Error processing',ftile)
		return -5
	else:
		return 0

pool = multiprocessing.Pool(processes=numthreads)

for ftile in glob.glob(os.path.join(dem_dir,'*_dem.tif')):
	pool.apply_async(processtile,(ftile,out_dir1,out_dir2,processdir))

# close the pool and make sure the processing has finished
pool.close()
pool.join()
