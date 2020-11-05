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
dem_dir = '/export/anthropocene/array-01/pu17449/parameter_transfer_datasets/global_forest' # files have format e.g. n25e090_dem.tif
out_dir2 = '/export/anthropocene/array-01/pu17449/parameter_transfer_datasets/global_forest_p1deg'
processdir = '/export/anthropocene/array-01/pu17449/parameter_transfer_datasets/tmp'

numthreads = 6

def processtile(ftile,out_dir2,processdir):
	fname = os.path.basename(ftile)
	print(fname)


	# Regrid to 0.5deg
	# For 0.5deg regridding	
	slopefile2 = os.path.join(out_dir2,fname[:-4]+'_p1deg.tif')
	if not os.path.exists(slopefile2):	
		cmd = ['gdalwarp', '-t_srs', 'EPSG:4326', '-tr', '0.1', '0.1', '-r', 'average' ,'-of' ,'GTiff','-co','COMPRESS=DEFLATE', ftile, slopefile2]
		#print(cmd)
		ret3 = subprocess.call(cmd)
		if ret3!=0:
			print('Error processing',ftile)
			return None
		
	return slopefile2

pool = multiprocessing.Pool(processes=numthreads)

for ftile in glob.glob(os.path.join(dem_dir,'*.tif')):
	pool.apply_async(processtile,(ftile,out_dir2,processdir))
	#processtile(ftile,out_dir2,processdir)

# close the pool and make sure the processing has finished
pool.close()
pool.join()
