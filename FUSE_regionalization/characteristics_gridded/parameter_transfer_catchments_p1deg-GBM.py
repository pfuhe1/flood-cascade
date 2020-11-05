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

# Temporary directory for processing data
processdir = '/export/anthropocene/array-01/pu17449/parameter_transfer_datasets/tmp'
# Output directory
outdir = '/export/anthropocene/array-01/pu17449/parameter_transfer_datasets/p1deg-GBM/'


# Input directory containing forestcover tiles (10 degree increments)
#tiledir = '/export/anthropocene/array-01/pu17449/parameter_transfer_datasets/global_forest_p1deg'
#varname = 'forestcover'

tiledir = '/export/anthropocene/array-01/pu17449/parameter_transfer_datasets/MeritDEMSlope_p1deg'
varname = 'slope'


# Number of threads for processing in parallel
numthreads = 4

# Create outdir if needed
if not os.path.exists(outdir):
	os.mkdir(outdir)

# Domain extent
xmin = '73'
xmax = '98'
ymin = '22'
ymax = '31.5'

# Function to do the processing for a single tile/catchment pair
def processtile(ftile,processdir):
	fname = os.path.basename(ftile)

	#print(fname)
	# Clip tile to GBM region 
	# Output file will be all missing values if region does not intersect tile
	tmpfile = os.path.join(processdir,fname[:-4]+'.tif')
	if not os.path.exists(tmpfile):
		cmd = ['gdal_translate','-of','GTiff','-projwin',xmin,ymax,xmax,ymin,'-a_nodata','255',ftile,tmpfile]
		ret = subprocess.call(cmd)
		# If file created successfully, check for valid values
		if ret==0:
			f = gdal.Open(tmpfile)
			rasterband = f.GetRasterBand(1)
			arr = rasterband.ReadAsArray()
			missingval = 255 
			validvals = (arr!=missingval).sum()
			#print('validvals',validvals)
			if validvals == 0:# all values are missing
				os.remove(tmpfile)
				return None
		else:
			if os.path.exists(tmpfile):
				os.remove(tmpfile)
			#raise Exception('Error with processing files: '+fname+','+os.path.basename(catchment))			
			return None

	return tmpfile

# Get lists of tiles and catchments from directories
tiles = glob.glob(os.path.join(tiledir,'*.tif'))

# Initialise variable for mean catchment values
meanvals = {}

try:
	outfile = os.path.join(outdir,varname+'_p1deg-GBM.tif')

	if not os.path.exists(outfile):
		pool = multiprocessing.Pool(processes=numthreads)

		# Loop over results and find any parts within each tile
		results = []
		outfiles = []
		for ftile in tiles:
			# Multiprocessing
			results.append(pool.apply_async(processtile,(ftile,processdir)))
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
			raise Exception('No valid results')
		elif len(outfiles)==1:
			# simply rename the tmpfile
			os.rename(outfiles[0],outfile)
		else:
			# Merge files:
			cmd = ['gdal_merge.py','-ot','Byte','-of','GTiff','-co','COMPRESS=DEFLATE','-o',outfile]+outfiles
			ret = subprocess.call(cmd)
			if ret!=0:
				raise Exception('Merging files did not work')
			# Clean up temporary files
			for f in outfiles:
				os.remove(f)

except Exception as e:
	print(e)
	

