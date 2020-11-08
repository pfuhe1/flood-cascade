import gdal
import numpy as np
import os,glob,shutil,sys
import netCDF4
import matplotlib.pyplot as plt

# Create FUSE elev_bands grid files for the GBM.
# Also creates a 'domain' file for metsim (https://github.com/UW-Hydro/MetSim)
# Assumes elevation data is provided by MeritDEM and FUSE grid is 0.1 degree resolution.
# Original file name: MeritDEM_regrid_elevs_GBM-p1deg.py

griddir = '/export/anthropocene/array-01/pu17449/FUSE_inputs/grids'

################################################################################
# Input paths

# Input directory containing DEM tiles (5degree increments), downloaded from http://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/
# files have format e.g. n25e090_dem.tif
dem_dir = '/export/anthropocene/array-01/pu17449/MeritDEM/'

# File containing output domain. Values are 0 or masked (-9999.)
grid_template = os.path.join(griddir,'GBM-Grid_p1deg_v2.nc')
# Netcdf file, created by running command 'cdo gridarea' on grid template file
grid_area = os.path.join(griddir,'gridarea_GBM-p1deg.nc')

#####################################
# Ouptut paths

# output 'domain' file for metsim
outdomain_file = os.path.join(griddir,'domain_GBM-p1deg.nc')

# output elev_band file for FUSE:
elevbands_file = os.path.join(griddir,'elev-bands_GBM-p1deg.nc')

##################################

# open grid info
f_grid = netCDF4.Dataset(grid_template)
latvals = f_grid.variables['lat'][:]
lonvals = f_grid.variables['lon'][:]
gridmask = f_grid.variables['Band1'][:]
n_outx = len(lonvals)
n_outy = len(latvals)

# Grid extent for EWEMBI 0.5 deg grid
#grid_x = [73,98]
#grid_y = [22,31.5]


# number of DEM points per grid box
stepx = 120
stepy = 120

# Number of elevation bands to use for FUSE (user choice)
nbands = 16

################################################################################
# initialise output arrays:
elev = np.zeros([n_outy,n_outx],dtype=np.int)
mask = np.ma.logical_not(gridmask).filled(0) # true is used, false is outside region
frac = mask # TODO, could consider making this into a fraction instead of true/false

band_elev = np.zeros([nbands,n_outy,n_outx],dtype=np.float32)
band_frac = np.zeros([nbands,n_outy,n_outx],dtype=np.float32)

# calculate (approx) grid box area
#area_p5deg = (111000/2.)**2.
# calculate weights (convert degrees to radians)
#area1D = area_p5deg * np.cos(latvals/180.*np.pi)
# repeat the 1D array across all longitudes
#area = np.repeat(area1D[:,np.newaxis],n_outx,axis=1)

# Use grid area calculated by cdo
with netCDF4.Dataset(grid_area,'r') as f_area:
	area = f_area.variables['cell_area'][:]
print('shape',area.shape)

tile_last = None

for x,lon in enumerate(lonvals):
	for y,lat in enumerate(latvals):
		print('lon,lat',lon,lat,mask[y,x])

		if mask[y,x]:
			tilex = str(int(5 * (lon // 5))).zfill(3)
			tiley =str(int(5 * (lat // 5))).zfill(2)
			f_tile = os.path.join(dem_dir,'n'+tiley+'e'+tilex+'_dem.tif')
			print(os.path.basename(f_tile))

			################################################################################
			# Read in data to calculate grid box elevations
			if not tile_last == f_tile: # check if file is already open
				f=gdal.Open(f_tile)
				rasterband = f.GetRasterBand(1)
				dem_arr = rasterband.ReadAsArray()
			else:
				tile_last = f_tile

			# intial indices for gridbox in tile
			# (number of 0.1 degrees past the nearest 5 degrees) (Should0-49)
#			yinit = 49-int((lat%5 - 0.05)*2) # Test flipping direction
			yinit = 49-int((lat%5 - 0.05)*10) # Test flipping direction
			xinit = int((lon%5 - 0.05)*10)
			print('x,y',xinit,yinit)
			#box = dem_arr[y*stepy:(y+1)*stepy,x*stepx:(x+1)*stepx]
			box = dem_arr[yinit*stepy:(yinit+1)*stepy,xinit*stepx:(xinit+1)*stepx]
			if lon == 73.75 and lat == 24.75:
				plt.figure()
				plt.pcolor(box)
				plt.colorbar()
				plt.show()
			# Try with masked array
			box = np.ma.masked_values(box,-9999.)
			box = box.compressed()
			if len(box)==0: # mask out watery boxes
				elev[y,x] = -9999.
				mask[y,x] = 0
				frac[y,x] = 0
				band_frac[0,y,x] = -9999.
				band_elev[0,y,x] = -9999.
				continue
			# Mean elevation
			elev[y,x] = box.mean()
			# Now split elevation into bands
			min_elev = box.min()
			#min_elev = max(0,min_elev)
			max_elev = box.max()
			bands = np.linspace(min_elev,max_elev,num=nbands+1)

			for i in range(nbands)[::-1]:
				subset = box[np.logical_and(box> bands[i],box<bands[i+1])]
				#band_elev[i,y,x] = subset.mean()
				band_elev[i,y,x] = (bands[i]+bands[i+1])/2. # halfway point of band
				band_frac[i,y,x] = len(subset)
			# Normalise band_frac
			band_frac[:,y,x] = band_frac[:,y,x]/band_frac[:,y,x].sum()
			print('debug area_sum',band_frac[:,y,x].sum())
		else:
			# Region not within domain
			elev[y,x] = -9999.
			band_frac[:,y,x] = -9999.
			band_elev[:,y,x] = -9999.


#print('debug, exiting without writing files')
#sys.exit(1)

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
	f_out_modsim.variables['lat'][:] = latvals

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



###############################################################################
# Write output file for FUSE elev bands
# Follows example of file from fuse_grid example: cesm1-cam5_elev_bands.nc

with netCDF4.Dataset(elevbands_file,'w') as f_out:
	f_out.createDimension('latitude',n_outy)
	f_out.createVariable('latitude',np.float,('latitude'))
	f_out.variables['latitude'].standard_name = "latitude"
	f_out.variables['latitude'].long_name = "latitude"
	f_out.variables['latitude'].units = "degrees_north"
	f_out.variables['latitude'].axis = "Y"
	f_out.variables['latitude'][:] = latvals

	f_out.createDimension('longitude',n_outx)
	f_out.createVariable('longitude',np.float,('longitude'))
	f_out.variables['longitude'].standard_name = "longitude"
	f_out.variables['longitude'].long_name = "longitude"
	f_out.variables['longitude'].units = "degrees_east"
	f_out.variables['longitude'].axis = "X"
	f_out.variables['longitude'][:] = lonvals

	f_out.createDimension('elevation_band',nbands)
	f_out.createVariable('elevation_band',np.int32,('elevation_band'))
	f_out.variables['elevation_band'].long_name = 'elevation_band'
	f_out.variables['elevation_band'].units = '-'
	f_out.variables['elevation_band'][:] = range(1,nbands+1)

	f_out.createVariable('area_frac',np.float32,('elevation_band', 'latitude', 'longitude'), fill_value=-9999)
	f_out.variables['area_frac'].units = "-" ;
	f_out.variables['area_frac'].long_name = "Fraction of grid cell covered by each elevation band."
	f_out.variables['area_frac'][:] = band_frac

	f_out.createVariable('mean_elev',np.float32,('elevation_band', 'latitude', 'longitude'),fill_value=-9999)
	f_out.variables['mean_elev'].units = "m asl" ;
	f_out.variables['mean_elev'].long_name = "Mean elevation of elevation band." ;
	f_out.variables['mean_elev'][:] = band_elev

	f_out.createVariable('prec_frac',np.float32,('elevation_band', 'latitude', 'longitude'), fill_value=-9999)
	f_out.variables['prec_frac'].units = "-" ;
	f_out.variables['prec_frac'].long_name = "Fraction of cell precipitation that falls on each elevation band." ;
	f_out.variables['prec_frac'].comment = 'PFU note: prec frac is just area_frac for now, need to implement scaling/lapse rate with elevation'
	f_out.variables['prec_frac'][:]=band_frac
