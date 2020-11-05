import gdal
import numpy as np
import os,glob,shutil,sys
import netCDF4
import matplotlib.pyplot as plt
import subprocess
import multiprocessing
import gdal
import pickle

# On BRIDGE servers, gdal is a bit broken, set 'GDAL_DATA='/opt/bridge/CentOS6-64/python/anaconda-5.0.1-2.7/pkgs/libgdal-2.1.0-0/share/gdal'

################################################################################
# Input paths


# Relys on first running 'parameter_transfer_p5dataset_MERITslope.py' 

# Input directory containing forestcover tiles (10 degree increments)
tiledir = '/export/anthropocene/array-01/pu17449/parameter_transfer_datasets/global_forest'
# Shapefiles of catchment boundaries
catchdir = '/export/anthropocene/array-01/pu17449/parameter_transfer_datasets/shapefiles'
# Temporary directory for processing data
processdir = '/export/anthropocene/array-01/pu17449/parameter_transfer_datasets/tmp'
# Output directory
outdir = '/export/anthropocene/array-01/pu17449/parameter_transfer_datasets/globalforest_catchments'
# Pickle file containing mean values for slope over each catchment
mean_out = '/export/anthropocene/array-01/pu17449/parameter_transfer_datasets/globalforest_catchments.pkl'
# Number of threads for processing in parallel
numthreads = 4

# Create outdir if needed
if not os.path.exists(outdir):
	os.mkdir(outdir)

# Function to do the processing for a single tile/catchment pair
def processtile(ftile,catchment,catchid,processdir):
	fname = os.path.basename(ftile)

	print(fname)
	# Clip tile to catchment region 
	# Output file will be all missing values if catchment is not in tile
	tmpfile = os.path.join(processdir,fname[:-4]+'_'+catchid+'.tif')
	if not os.path.exists(tmpfile):
		cmd = ['gdalwarp','-of','GTiff','-cutline',catchment,'-crop_to_cutline',ftile,tmpfile]
		ret = subprocess.call(cmd)
		# If file created successfully, check for valid values
		if os.path.exists(tmpfile):
			f = gdal.Open(tmpfile)
			rasterband = f.GetRasterBand(1)
			arr = rasterband.ReadAsArray()
			missingval = 0 # Todo, could read this from the file
			validvals = (arr!=missingval).sum()
			#print('validvals',validvals)
			if validvals == 0:# all values are missing
				os.remove(tmpfile)
				return None
		else:
			raise Exception('Error with processing files: '+fname+','+os.path.basename(catchment))

	return tmpfile

# Get lists of tiles and catchments from directories
tiles = glob.glob(os.path.join(tiledir,'*.tif'))
catchments = glob.glob(os.path.join(catchdir,'*.shp'))

# Initialise variable for mean catchment values
meanvals = {}

# Loop over catchments
for catchment in catchments:

	try:
		catchid = os.path.basename(catchment).split('_')[5][:-4]
		print('Catchment',catchid)
		outfile = os.path.join(outdir,'forestcover_'+catchid+'.tif')

		if not os.path.exists(outfile):
			pool = multiprocessing.Pool(processes=numthreads)

			# Loop over results and find any parts within each tile
			results = []
			outfiles = []
			for ftile in tiles:
				# Multiprocessing
				results.append(pool.apply_async(processtile,(ftile,catchment,catchid,processdir)))
				# Single process:
				#retval = processtile(ftile,catchment,catchid,processdir)
				#if retval is not None:
				#	outfiles.append(retval)
			# close the pool and make sure the processing has finished
			pool.close()
			pool.join()

			# Get output filenames from results (multiprocessing)
			for result in results:
				retval = result.get()
				if retval is not None:
					outfiles.append(retval)
			print(outfiles)

			# Create final result file:
			if len(outfiles)==0:
				raise Exception('No valid results for catchment: '+os.path.basename(catchment))
			elif len(outfiles)==1:
				# simply rename the tmpfile
				os.rename(outfiles[0],outfile)
			else:
				# Merge files:
				cmd = ['gdal_merge.py','-ot','Byte','-of','GTiff','-o',outfile]+outfiles
				ret = subprocess.call(cmd)
				if ret!=0:
					raise Exception('Merging files did not work')
				# Clean up temporary files
				for f in outfiles:
					os.remove(f)

		# finally calculate average for catchment
		f = gdal.Open(outfile)
		rasterband = f.GetRasterBand(1)
		arr = rasterband.ReadAsArray()
		rasterband = f.GetRasterBand(1)
		arr = rasterband.ReadAsArray()
		# PFU, don't mask the landcover array over catchments 
		# (missingval is 0% which we should include)
		#missingval = -9999. # Todo, could read this from the file
		#marr = np.ma.masked_values(arr,missingval)
		meanvals[int(catchid)] = arr.mean()

	except Exception,e:
		print('Error processing catchment:',catchid)
		print(e)
	
# Write out means to pickle file
with open(mean_out,'wb') as f:
	pickle.dump(mean_out,f,-1)
	

