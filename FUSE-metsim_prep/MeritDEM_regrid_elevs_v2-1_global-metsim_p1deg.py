import gdal
import numpy as np
import os,glob,shutil,sys
import netCDF4
import matplotlib.pyplot as plt
import multiprocessing,subprocess

# On BRIDGE servers, gdal is a bit broken, set 'GDAL_DATA='/opt/bridge/CentOS6-64/python/anaconda-5.0.1-2.7/pkgs/libgdal-2.1.0-0/share/gdal'

################################################################################
# Input paths

# Input directory containing DEM tiles (5degree increments)
dem_dir = '/export/anthropocene/array-01/pu17449/MeritDEM/' # files have format e.g. n25e090_dem.tif
processdir = '/export/anthropocene/array-01/pu17449/parameter_transfer_datasets/tmp'
mergelist_path = os.path.join(processdir,'mergelist.txt')
numthreads = 4

# File containing output domain. Values are 0 or masked (-9999.)
grid_template = '/export/anthropocene/array-01/pu17449/isimip_bc/obs/MSWEP2-2/global_daily_010deg_reformat/MSWEP2-2_010deg_199008.nc'
# Netcdf file, created by running command 'cdo gridarea' on grid template file
grid_area = '/export/anthropocene/array-01/pu17449/isimip_bc/obs/MSWEP2-2/global_010deg_gridarea.nc'

#####################################
# Ouptut paths

# output 'domain' file for metsim
outdomain_file = '/home/bridge/pu17449/src/MetSim/global/domain_global-tiled_p1deg.nc'

# output elev_band file for FUSE:
#elevbands_file = '/home/bridge/pu17449/src/MetSim/global/global-tiled_elev_bands_p1deg.nc'

##################################

# open grid info
f_grid = netCDF4.Dataset(grid_template)
latvals = f_grid.variables['lat'][:]
lonvals = f_grid.variables['lon'][:]
#gridmask =  f_grid.variables['BCmask'][0,:]
n_outx = len(lonvals)


# Grid extent for EWEMBI 0.5 deg grid
#grid_x = [73,98]
#grid_y = [22,31.5]

# Merit data only spans 85N - 60S
latvals_out = latvals[50:1500]
n_outy = len(latvals_out)

# Use grid area calculated by cdo 
with netCDF4.Dataset(grid_area,'r') as f_area:
	area = f_area.variables['cell_area'][50:1500,:]
print('shape',area.shape)

# Regrid to 0.1x0.1deg grid
def processtile(ftile,processdir):
	fname = os.path.basename(ftile)
	outfile = os.path.join(processdir,fname[:-4]+'_p1deg.tif')
	if not os.path.exists(outfile):
		cmd = ['gdalwarp', '-t_srs', 'EPSG:4326', '-tr', '0.1', '0.1', '-r', 'average' ,'-of' ,'GTiff', ftile,outfile]
		subprocess.call(cmd)
	return outfile

##################################

outfile = os.path.join(processdir,'meritdem_p5deg_merged.tif')
if not os.path.exists(outfile):
	pool = multiprocessing.Pool(processes=numthreads)
	results = []
	for ftile in glob.glob(os.path.join(dem_dir,'*_dem.tif')):
		results.append(pool.apply_async(processtile,(ftile,processdir)))

	# close the pool and make sure the processing has finished
	pool.close()
	pool.join()	

	# Get output filenames from returnvals
	outfiles = []
	with open(mergelist_path,'w') as outlist:
		for result in results:
			retval = result.get()
			if retval is not None:
				outfiles.append(retval)
				outlist.write(retval+'\n')


	cmd = ['gdal_merge.py','-ot','Float32','-of','GTiff','-o',outfile,'-a_nodata','-9999.','--optfile',mergelist_path]
	ret = subprocess.call(cmd)
	if ret!=0:
		raise Exception('Merging files did not work')

#print('debug, exiting without writing netcdf')
#sys.exit(1)

# Open merged file
f = gdal.Open(outfile)
rasterband = f.GetRasterBand(1)
elev = rasterband.ReadAsArray()

print('shapes',elev.shape,area.shape)


################################################################################
# initialise output arrays:
mask = elev!=-9999. # true is used, false is outside region
frac = mask # TODO, could consider making this into a fraction instead of true/false
	
###############################################################################
# Write output file for modsim: 

# Follows format of variables needed in 'domain' file for metsim e.g. '/home/bridge/pu17449/src/MetSim/metsim/data/domain.nc'
with netCDF4.Dataset(outdomain_file,'w') as f_out_modsim:

	f_out_modsim.createDimension('lat',n_outy)
	f_out_modsim.createVariable('lat',np.float,('lat'))
	f_out_modsim.variables['lat'].standard_name = "latitude"
	f_out_modsim.variables['lat'].long_name = "latitude"
	f_out_modsim.variables['lat'].units = "degrees_north"
	f_out_modsim.variables['lat'].axis = "Y"
	f_out_modsim.variables['lat'][:] = latvals_out

	f_out_modsim.createDimension('lon',n_outx)
	f_out_modsim.createVariable('lon',np.float,('lon'))
	f_out_modsim.variables['lon'].standard_name = "longitude"
	f_out_modsim.variables['lon'].long_name = "longitude"
	f_out_modsim.variables['lon'].units = "degrees_east"
	f_out_modsim.variables['lon'].axis = "X"
	f_out_modsim.variables['lon'][:] = lonvals

	f_out_modsim.createVariable('elev',np.int32,('lat','lon'),fill_value=-9999)
	#f_out_modsim.variables['elev']._FillValue = -9999
	f_out_modsim.variables['elev'].long_name = 'gridcell_elevation'
	f_out_modsim.variables['elev'].units = 'm'
	f_out_modsim.variables['elev'][:] = elev

	f_out_modsim.createVariable('mask',np.int64,('lat','lon'))
	f_out_modsim.variables['mask'].long_name = 'domain_mask'
	f_out_modsim.variables['mask'].comment = '0 indicates cell is not active'
	f_out_modsim.variables['mask'][:] = mask

	f_out_modsim.createVariable('frac',np.float,('lat','lon'))
	f_out_modsim.variables['frac'].long_name = 'fraction of grid cell that is active'
	f_out_modsim.variables['frac'].units = '1'
	f_out_modsim.variables['frac'][:] = frac

	f_out_modsim.createVariable('area',np.float,('lat','lon'))
	f_out_modsim.variables['area'].standard_name = 'area'
	f_out_modsim.variables['area'].long_name = 'area of grid cell'
	f_out_modsim.variables['area'].units = 'm2'
	f_out_modsim.variables['area'][:] = area


