# Process PET from worldclim temperature
# Peter Uhe
# 25/3/2019
#
# This script uses python3, set up using conda environment gdal_env2

import os,sys,glob
import gdal
import numpy as np
import datetime,calendar
import netCDF4

sys.path.append('/home/pu17449/gitsrc/PyETo')
import pyeto

inpath = '/home/pu17449/data2/worldclim_precip/'
outfile = '/home/pu17449/data2/worldclim_precip/pet_yearmean_v3.nc'
#f_tiffout = '/home/pu17449/data2/worldclim_precip/pet_yearmean_v2.tif'

# Construct lat and lon arrays (centrepoints)
step = 0.25/30.
lons = np.arange(-180,180,step)+step/2.
lats = np.arange(90,-90,-1*step)-step/2.
nlats = len(lats)
nlons = len(lons)
lats_rad = pyeto.deg2rad(lats)
#lat2D = np.repeat(lats_rad[:,np.newaxis],len(lons),axis=1)

# Data array to store time/lon slices of PET
pet_timeslice = np.zeros([12,nlons],dtype=np.float32)
pet = np.zeros([nlons],dtype=np.float32) #yearly mean
	
###############################################################################
# Write output file for modsim: 

# Follows format of variables needed in 'domain' file for metsim e.g. '/home/bridge/pu17449/src/MetSim/metsim/data/domain.nc'
with netCDF4.Dataset(outfile,'w') as f_out:

	f_out.createDimension('lat',nlats)
	f_out.createVariable('lat',np.float,('lat'))
	f_out.variables['lat'].standard_name = "latitude"
	f_out.variables['lat'].long_name = "latitude"
	f_out.variables['lat'].units = "degrees_north"
	f_out.variables['lat'].axis = "Y"
	f_out.variables['lat'][:] = lats

	f_out.createDimension('lon',nlons)
	f_out.createVariable('lon',np.float,('lon'))
	f_out.variables['lon'].standard_name = "longitude"
	f_out.variables['lon'].long_name = "longitude"
	f_out.variables['lon'].units = "degrees_east"
	f_out.variables['lon'].axis = "X"
	f_out.variables['lon'][:] = lons

	#f_out.createDimension('time',None)
	#f_out.createVariable('time',np.float,('time'))
	#f_out.variables['time'].standard_name = "time"
	#f_out.variables['time'].long_name = "time"
	#f_out.variables['time'].units = "days since 2000-01-01"
	#f_out.variables['time'].axis = "T"
	#f_out.variables['time'][:] = [15]

	f_out.createVariable('pet',np.float32,('lat','lon'))
	f_out.variables['pet'].long_name = 'Yearly average Potential Evapotranspiration'
	f_out.variables['pet'].comment = 'Calculated by Hargraeves Method, using WorldClim2 monthly data and pyeto (https://pyeto.readthedocs.io)'
	f_out.variables['pet'].units = 'mm/day'

	# Try writing out to geotiff (currently not working)
	#gtiffDriver = gdal.GetDriverByName('GTiff')
	#f_template = os.path.join(inpath,'wc2.0_30s_tavg','wc2.0_30s_tavg_01.tif')
	#f1 = gdal.Open(f_template)
	#f2 = gtiffDriver.CreateCopy(f_tiffout,f1,strict=0,options=["COMPRESS=DEFLATE"])
	#rasterband = f2.GetRasterBand(1)


	# Monthdays: Jan Feb    Mar Apr May Jun Jul Aug Sep Oct Nov Dec
	monthdays = [31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

	for i,lat in enumerate(lats_rad):
		print(i)
		#pet = np.zeros([12,nlons],dtype=np.float32)
		for month in range(1,13):
			

			# Test for Jan (this is a climatology project, just use date of year 1999)
			# choose (approx) middle of month
			tmp,day_in_month = calendar.monthrange(1999,month)
			day_of_year = datetime.date(1999,month,(day_in_month+1)//2).timetuple().tm_yday
			sol_dec = pyeto.sol_dec(day_of_year)            # Solar declination
			sha = pyeto.sunset_hour_angle(lat, sol_dec)
			ird = pyeto.inv_rel_dist_earth_sun(day_of_year)
			et_rad = pyeto.et_rad(lat, sol_dec, sha, ird)   # Extraterrestrial radiation

			# Load data files
			f_tavg = os.path.join(inpath,'wc2.0_30s_tavg','wc2.0_30s_tavg_'+str(month).zfill(2)+'.tif')
			#print('Loading tavg:',f_tavg)
			f=gdal.Open(f_tavg)
			rasterband = f.GetRasterBand(1)
			tavg = rasterband.ReadAsArray(yoff=i,win_ysize=1)
			#print(tavg.shape)
			#f.close()

			f_tmax = os.path.join(inpath,'wc2.0_30s_tmax','wc2.0_30s_tmax_'+str(month).zfill(2)+'.tif')
			#print('Loading tmax:',f_tmax)
			f=gdal.Open(f_tmax)
			rasterband = f.GetRasterBand(1)
			tmax = rasterband.ReadAsArray(yoff=i,win_ysize=1)
			#f.close()

			f_tmin = os.path.join(inpath,'wc2.0_30s_tmin','wc2.0_30s_tmin_'+str(month).zfill(2)+'.tif')
			#print('Loading tmin:',f_tmax)
			f=gdal.Open(f_tmin)
			rasterband = f.GetRasterBand(1)
			tmin = rasterband.ReadAsArray(yoff=i,win_ysize=1)
			#f.close()
			
			# pyeto gives mm/day, multiply by monthdays to get mm/month
			pet_timeslice[month-1,:] = pyeto.hargreaves(tmin,tmax,tavg, et_rad)*monthdays[month-1]

		# Sum over months to get mm/year, and write out
		pet = pet_timeslice.sum(0,keepdims=True)
		f_out.variables['pet'][i,:] = pet[0,:]

		#print(rasterband)
		# Also write data out to tiff format
		#rasterband.WriteArray(pet,yoff=i)

# Close tif datasets
#f1 = None
#f2 = None
print('done')
