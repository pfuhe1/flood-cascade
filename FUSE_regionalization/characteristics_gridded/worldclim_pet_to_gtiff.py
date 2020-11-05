# Process PET from worldclim temperature
# Peter Uhe
# 25/3/2019
#
# This script uses python3, set up using conda environment gdal_env2

import os,sys,glob
import gdal
import netCDF4

inpath = '/home/pu17449/data2/worldclim_precip/'
ncfile = '/home/pu17449/data2/worldclim_precip/pet_yearmean_v3.nc'
f_tiffout = '/home/pu17449/data2/worldclim_precip/pet_yearmean_v3.tif'
	
###############################################################################
# Write output file for modsim: 

# Follows format of variables needed in 'domain' file for metsim e.g. '/home/bridge/pu17449/src/MetSim/metsim/data/domain.nc'
with netCDF4.Dataset(ncfile,'r') as fin:
	data_arr = fin.variables['pet'][:]

# Try writing out to geotiff (currently not working)
gtiffDriver = gdal.GetDriverByName('GTiff')
f_template = os.path.join(inpath,'wc2.0_30s_prec','wc2.0_30s_prec_01.tif')
f1 = gdal.Open(f_template)
f2 = gtiffDriver.CreateCopy(f_tiffout,f1,strict=0,options=["COMPRESS=DEFLATE"])
rasterband = f2.GetRasterBand(1)
rasterband.WriteArray(data_arr)

# Close tif datasets
#f1 = None
#f2 = None
print('done')
