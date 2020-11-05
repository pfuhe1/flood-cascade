import gdal
import numpy as np
import os,glob,shutil,sys
import subprocess

# On BRIDGE servers, gdal is a bit broken, set 'GDAL_DATA='/opt/bridge/CentOS6-64/python/anaconda-5.0.1-2.7/pkgs/libgdal-2.1.0-0/share/gdal'

################################################################################

# Input paths
infiles = {}
infiles['aridity']='/home/pu17449/data2/worldclim_precip/wc2.0_30s_aridity/wc2.0_30s_aridity.tif'
infiles['clyppt'] = '/home/pu17449/data2/soilgrids/CLYPPT_M_ave_250m_ll.tif'
infiles['pr'] = '/home/pu17449/data2/worldclim_precip/wc2.0_30s_prec/wc2.0_30s_prec_year.tif'
infiles['tas'] = '/home/pu17449/data2/worldclim_precip/wc2.0_30s_tavg/wc2.0_30s_tavg_yearmean_v2.tif'
infiles['pet'] = '/home/pu17449/data2/worldclim_precip/pet_yearmean_v3.tif'
infiles['snowcover'] = '/home/pu17449/data2/modis_snowcover/modis_snowcover_mean_2001-2014.tif'


# Output directory
outdir = '/home/pu17449/data2/parameter_transfer_catchments/p1deg-GBM/'
# Create outdir if needed
if not os.path.exists(outdir):
	os.mkdir(outdir)

# Domain extent
xmin = '73'
xmax = '98'
ymin = '22'
ymax = '31.5'


# Function to do the processing for a single tile/catchment pair
def processtile(ftile,tmpfile):
	fname = os.path.basename(ftile)

	#print(fname)
	# Clip tile to GBM region 
	# Output file will be all missing values if region does not intersect tile
	if not os.path.exists(tmpfile):
		cmd = ['gdalwarp', '-t_srs', 'EPSG:4326','-te',xmin,ymin,xmax,ymax, '-tr', '0.1', '0.1', '-r', 'average' ,'-of' ,'GTiff','-co','COMPRESS=DEFLATE',ftile,tmpfile]
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
				raise Exception('Error with processing files: '+fname+','+os.path.basename(catchment))			

	return tmpfile

# Loop over files and processs
for varname,infile in infiles.items():
	outfile = os.path.join(outdir,varname+'_p1deg-GBM.tif')
	processtile(infile,outfile)

