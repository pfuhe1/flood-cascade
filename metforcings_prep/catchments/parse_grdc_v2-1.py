import os,sys,glob,socket
import datetime
import numpy as np
from netCDF4 import Dataset

def parse_grdc_file(fname,netcdf_out=None,override=False):
	# Missing values: convert from -999 to -9999
	in_missing_val = -999.
	out_missing_val = -9999.

	# Data arrays
	dates = [] # format is days since firstDate
	vals = [] # discharge data (mm/day)

	#
	firstDate = None
	notes = ''

	lineno = 1
	# Open file and loop over lines
	with open(fname,'br') as f:
		for l in f:
			# Convert textfile encoding to ASCII
			line = l.decode('ISO-8859-15')
			if line[0]=='#':
				# Header line
				notes+=line
				if lineno == 13:
					lat = float(line.split(':')[-1])
				if lineno == 14:
					lon = float(line.split(':')[-1])
				if lineno == 15:
					area = float(line.split(':')[-1]) * 1e6 # convert from km2 to m2

			else:
				# Parse line for date and value
				datetime_str,val = line.strip().split()
				date,hrmn,tmp = datetime_str.split(';')
				yr,mon,day = date.split('-')
				try:
					d = datetime.datetime(int(yr),int(mon),int(day))
					if firstDate is None:
						firstDate = d
					delta = (d - firstDate).days
					dates.append(delta)
					vals.append(float(val))
				except:
					print('error parsing line',line)
					continue
			lineno+=1

	print('area,lat,lon:',area,lat,lon)
	# hack for catchment with no area
	if os.path.basename(fname)=='1434820_Q_Day.Cmd.txt':
		area = 4758398216.865 # m2, calculated from shapefile area
	#if area == -999.*1e6:
	if area < 0.0:
			raise Exception('Error, negative area in discharge file',fname)

	# Convert data into arrays
	dates = np.array(dates)
	vals = np.ma.masked_values(vals,in_missing_val)
	# Convert vals from m3/s to mm/day
	vals = vals / area * 86400. *1e3


	# Output netcdf file
	if netcdf_out is not None and (not os.path.exists(netcdf_out) or override):
		print('Writing out to netcdf:',netcdf_out)
		with Dataset(netcdf_out,'w') as f:
			f.createDimension('time',len(dates))
			f.createVariable('time',np.float32,('time'))
			f.variables['time'].standard_name = 'time'
			f.variables['time'].long_name = 'time'
			f.variables['time'].units = 'days since '+str(firstDate).split()[0]
			f.variables['time'][:] = dates

#			f.createDimension('lat',1)
#			f.createVariable('lat',np.float32,('lat'))
#			f.variables['lat'][:]=[lat]
#			f.createDimension('lon',1)
#			f.createVariable('lon',np.float32,('lon'))
#			f.variables['lon'][:]=[lon]

			f.createVariable('q_obs',np.float32,('time',),fill_value = out_missing_val)
			f.variables['q_obs'].long_name = 'Mean observed daily discharge'
			f.variables['q_obs'].units = 'mm/day'
			f.variables['q_obs'][:] = vals.filled(out_missing_val)
			f.GRDC_header = notes


	return firstDate,dates,vals

if __name__ == '__main__':

	override = False

	# Set output path
	fpath_in = '/home/pu17449/data2/Discharge Data/GRDC_daily_global/textfiles/*_Q_Day.Cmd.txt'
	fpath_out = '/home/pu17449/data2/Discharge Data/GRDC_daily_global/netcdf_v2-1'

	if not os.path.exists(fpath_out):
		os.mkdir(fpath_out)

	# Loop over GRDC text files and process
	for fname in glob.glob(fpath_in):
		try:
			bname = os.path.basename(fname)
			print(fname)
			fname_netcdf = os.path.join(fpath_out,bname.split('.')[0]+'.nc')
			if not os.path.exists(fname_netcdf) or override:
				firstDate,dates,vals = parse_grdc_file(fname,fname_netcdf,override=override)
			# Use for debugging without writing out files:
			#firstDate,dates,vals = parse_grdc_file(fname,netcdf_out=None,override=False)

			print('Missing values',vals.mask.sum(),'of',len(vals))
		except Exception as e:
			print(e)
