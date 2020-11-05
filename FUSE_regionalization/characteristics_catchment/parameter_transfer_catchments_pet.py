import numpy as np
import os,glob,shutil,sys
import subprocess
import multiprocessing
import gdal
import pickle

# On BRIDGE servers, gdal is a bit broken, set 'GDAL_DATA='/opt/bridge/CentOS6-64/python/anaconda-5.0.1-2.7/pkgs/libgdal-2.1.0-0/share/gdal'

################################################################################
# Input paths

var = 'pet'

# Input gridded dataset
tile = '/home/pu17449/data2/worldclim_precip/pet_yearmean_v3.tif'
# Shapefiles of catchment boundaries
catchdir = '/home/pu17449/data2/Discharge Data/GRDC_daily_global/shapefiles'

# Output directory
outdir = '/home/pu17449/data2/parameter_transfer_catchments/'+var
# Pickle file containing mean values for slope over each catchment
mean_out = '/home/pu17449/data2/parameter_transfer_catchments/'+var+'.pkl'

# Number of threads for processing in parallel
numthreads = 4

# Create outdir if needed
if not os.path.exists(outdir):
	os.mkdir(outdir)

# Function to do the processing for a single tile/catchment pair
def processtile(ftile,catchment,catchid,tmpfile):
	fname = os.path.basename(ftile)

	#print(fname)
	# Clip tile to catchment region 
	if not os.path.exists(tmpfile):
		cmd = ['gdalwarp','-of','GTiff','-co','COMPRESS=DEFLATE','-cutline',catchment,'-crop_to_cutline',ftile,tmpfile]
		ret = subprocess.call(cmd)
		# If file created successfully, check for valid values
		if ret!=0:
			if os.path.exists(tmpfile):
				os.remove(tmpfile)
			raise Exception('Error with processing files: '+fname+','+os.path.basename(catchment))			
			#return None

	# finally calculate average for catchment
	f = gdal.Open(tmpfile)
	rasterband = f.GetRasterBand(1)
	arr = rasterband.ReadAsArray()
	#missingval = 255 # Todo, could read this from the file
	#missingval = -3.4e+38
	#missingval = -3.3999997493202682e+38
	missingval = -32768
	marr = np.ma.masked_values(arr,missingval)
	return catchid,marr.mean()

# Get lists of tiles and catchments from directories
catchments = glob.glob(os.path.join(catchdir,'*.shp'))

# Initialise variable for mean catchment values
meanvals = {}

# Loop over catchments
pool = multiprocessing.Pool(processes=numthreads)
results = []
for catchment in catchments:

	try:
		catchid = os.path.basename(catchment).split('_')[5][:-4]
		print('Catchment',catchid)
		outfile = os.path.join(outdir,var+'_'+catchid+'.tif')
		#if not os.path.exists(outfile):
		results.append(pool.apply_async(processtile,(tile,catchment,catchid,outfile)))
			
	except Exception as e:
		print('Error processing catchment:',catchid)
		print(e)

# close the pool and make sure the processing has finished
pool.close()
pool.join()

for result in results:
	try:
		catchid,val = result.get()
		meanvals[int(catchid)] = val
	except Exception as e:
		print('Error',e)

	
# Write out means to pickle file
with open(mean_out,'wb') as f:
	pickle.dump(meanvals,f,-1)
	

