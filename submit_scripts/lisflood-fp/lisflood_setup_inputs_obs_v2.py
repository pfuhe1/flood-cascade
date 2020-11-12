# lisflood_setup_inputs_obs_v2.py
# Takes input from mizuRoute discharge, river network information and produces lisflood boundary conditions
# Requires running first lisflood_discharge_inputs_qgis_clip.py to produce shapefiles (streamnet)
# Requires running first main LISFLOOD-FP_prep script calling LFPTools scripts to generate lisflood-fp model files
#
# Peter Uhe
# 2020/01/28
#

# Load modules
import os,sys,glob,pickle,shutil,socket
import numpy as np
from netCDF4 import Dataset,num2date
from osgeo import ogr
import csv
import datetime, calendar
from write_bdy_bci import write_bdy, write_bci_v2
import subprocess

import warnings
warnings.filterwarnings("error")

########################################################################################
# Function to calculate index of closest grid box,
# Inputs are point latitude and longitude coordinates and grid latitude,longitude coordinates
#
def find_closest_1d_v2(p_lat,p_lon,lat_coord,lon_coord):
	ny=lat_coord.shape[0]
	nx=lon_coord.shape[0]
	min_dist=100000000
	minx=-1
	miny=-1
	for j in range(ny):
		dist=(lat_coord[j]-p_lat)**2
		if dist<min_dist:
			miny=j
			min_dist=dist
	min_dist=100000000
	#print 'plon',p_lon
	for i in range(nx):
		dist= (np.abs(lon_coord[i]-p_lon)%360)**2
		if dist<min_dist:
			minx=i
			min_dist=dist
	return miny,minx

########################################################################################
# Function to calculate index of point,
# Inputs are point latitude and longitude coordinates and list of latitude,longitude pairs
#
def find_closest_list(p_lat,p_lon,lat_coord,lon_coord):
	ny=lat_coord.shape[0]
	nx=lon_coord.shape[0]
	min_dist=100000000
	minx=-1
	miny=-1
	for j in range(ny):
		dist=(lat_coord[j]-p_lat)**2 + (lon_coord[j]-p_lon)**2
		if dist<min_dist:
			miny=j
			min_dist=dist
	print(p_lat,p_lon,'min dist',min_dist)
	return miny,min_dist

########################################################################################
# Function to return true if point (xx,yy) is within bounds xbounds,ybounds
def in_bounds(xx,yy,xbound,ybound):
	return xbound[0]<xx and xx<xbound[1] and ybound[0]<yy and yy<ybound[1]


########################################################################################
# Function to check if point is not in list (also ignores -1)
def point_notinlist(point,points_list):
	return (point!=-1 and point not in points_list)

########################################################################################
# Little function to format dates into a string nicely
def format_date(date):
	return date.strftime('%Y-%m-%d')


###############################################################################################
# Start of main script
###############################################################################################
# Set up Shapefile reading driver
driver = ogr.GetDriverByName("ESRI Shapefile")


# Parameters for this domain
###############################################################################################
# extent xmin, xmax, ymin, ymax
extent = [87.628,92.7379,21.13,26.681]
regname = 'bangladesh-chansolverMSWEP_2020-11-02'

ds_buffer = 0.3 # buffer around each reach to set the downstream boundary

res = 0.0025 # in degrees
resname = '9sd8'

#res = 0.000833
#resname = '3sd8'

clipname = regname+'_'+resname

# Pattern of name of mizuroute simulations to process (can contain wildcards '*','?' etc)
mizuRuns = 'GBM-p1deg_900_MSWEP2-2-ERA5-calibrated1_MSWEP2-2-ERA5/q_*-1-1.nc'

startmon = 4 # april (1)
endmon = 10 # october (31)
ntimes = None

sublist = []

dryrun = False
skip_existing = 'started' # 'no', 'completed'

# Setup for file paths
host = socket.gethostname()
if host[:3]=='bp1':
	mizuroute_outdir = '/work/pu17449/mizuRoute/output/'
	# File, storing network attributes
	f_network        = '/home/pu17449/src/mizuRoute/route/ancillary_data/MERIT_mizuRoute_network_meta.nc'
	# Lisflood simulation directories
	lfdir = '/work/pu17449/lisflood/lisfloodfp_'+clipname
	logdir = '/work/pu17449/lisflood/logs'

	# Qsub file for HPC queue
	qsub_script = '/home/pu17449/src/setup_scripts/lisflood_discharge/call_pythonscript_bp.sh'
	control_script = '/home/pu17449/src/setup_scripts/lisflood_discharge/lisflood-fp_wrapper.py'
	exe_file = '/home/pu17449/bin/lisflood_double_rect_r688'
	ncpus = 24 # (node size is 24) # number of processors per job
	jobsize = 24 # number of processors per simulation
	simsperjob = 1
	qsub_cmd = ['qsub','-l','select=1:ncpus='+str(ncpus)+':ompthreads='+str(jobsize)+':mem='+str(4*ncpus)+'gb','-v','CONTROL_FLIST,EXE_FILE,LOGDIR,CONTROL_SCRIPT,LISFLOOD_DIR,RESULTDIR',qsub_script]

# Specify directories (relative paths)
###############################################################################################
lf_dischargedir = os.path.join(lfdir,'dischargeobs')
lf_pardir = os.path.join(lfdir,'parfilesobs')
resultdir = os.path.join(lfdir,'results')
streamdir = os.path.join(lfdir,'streamnet')

if not os.path.exists(lf_dischargedir):
	os.mkdir(lf_dischargedir)
if not os.path.exists(lf_pardir):
	os.mkdir(lf_pardir)

# Specify file paths (relative paths)
###############################################################################################
# Shapefiles containing up/downstream points in each river segment/link
# Use d8 discharge locations (needed if running model converted to d4)
resname2 = resname[:2]+'d8'
clipname2 = regname.split('-')[0]+'_'+resname2
f_downstream    = os.path.join(streamdir,clipname2+'_acc_downstream.shp')
f_ntdstream    = os.path.join(streamdir,clipname2+'_acc_next-to-downstream.shp')
f_upstream_clip = os.path.join(streamdir,clipname2+'_acc_upstream.shp')
f_upstream_all  = os.path.join(streamdir,'strn_network_'+resname2+'_acc_upstream.shp')

# parameter file template for lisflood
f_par_template = os.path.join(lfdir,'077_template.par')
f_bci = os.path.join(lf_dischargedir,clipname+'.bci')
print('bcifile:',f_bci)

# Export common environment variables
os.environ['LISFLOOD_DIR'] = lfdir
os.environ['EXE_FILE'] = exe_file
os.environ['LOGDIR'] = logdir
os.environ['CONTROL_SCRIPT'] = control_script
os.environ['NCPUS']=str(ncpus)
os.environ['OMP_NUM_THREADS'] = str(jobsize)
os.environ['RESULTDIR'] = resultdir


###############################################################################################
# Main script: Get points for discharge in river network

# Dict of attributes for each reach/link
# x,y coordinates
points_upall  = {}
points_upclip  = {}
points_ntd = {}
points_ds  = {}
# Upstream and downstream link numbers
dslinks     = {}
uslinks1    = {}
uslinks2    = {}
# Accumulation values
acc_upclip = {} # upstream accumulation (clipped to domain)
acc_upall  = {} # upstream accumulation (some points outside domain)
acc_ntd    = {} # accumulation at next-to-downstream point in reach
# Slope of reach
slope      = {}

# First load upstream points shapefile
dataSource = driver.Open(f_upstream_clip, 0)
try:
	points = dataSource.GetLayer()
except:
	print('Error, failed to open streamnet file: '+f_upstream_clip)
	raise
print('Opened upstream points in domain',f_upstream_clip)
for feature in points:
	# Get attributes:
	link = feature.GetField('LINKNO')
	dslinks[link] = feature.GetField('DSLINKNO')
	uslinks1[link] = feature.GetField('USLINKNO1')
	uslinks2[link] = feature.GetField('USLINKNO2')
	acc_upclip[link] = feature.GetField('ACC_1')
	slope[link] = feature.GetField('Slope')
	# Get location of point
	xx,yy,zz = feature.geometry().GetPoint()
	points_upclip[link] = [xx,yy]

links = list(points_upclip.keys())

# Load (next-to) downstream points shapefile
dataSource = driver.Open(f_ntdstream, 0)
points = dataSource.GetLayer()
print('Opened next-to-downstream points',f_ntdstream)
for feature in points:
	# Get attributes
	link = feature.GetField('LINKNO')
	acc_ntd[link] = feature.GetField('ACC_1')
	# Get location of point
	xx,yy,zz = feature.geometry().GetPoint()
	points_ntd[link] = [xx,yy]

# Load (next-to) downstream points shapefile
dataSource = driver.Open(f_ntdstream, 0)
points = dataSource.GetLayer()
print('Opened next-to-downstream points',f_ntdstream)
for feature in points:
	# Get attributes
	link = feature.GetField('LINKNO')
	# Get location of point
	xx,yy,zz = feature.geometry().GetPoint()
	points_ds[link] = [xx,yy]

# Load upstream points shapefile (all points)
dataSource = driver.Open(f_upstream_all, 0)
points = dataSource.GetLayer()
print('Opened upstream points (all)',f_upstream_all)
for feature in points:
	# Get attributes
	link = feature.GetField('LINKNO')
	if link in links: # Only get point for reaches in domain
		acc_upall[link] = feature.GetField('ACC_1')
		xx,yy,zz = feature.geometry().GetPoint()
		points_upall[link] = [xx,yy]

###############################################################################################
# Work out downstream boundary conditions,
# Boundary should be closed except at downstream boundaries
# create buffer around each downstream point for
# use slope of reach to prescribe water surface slope at boundary
# TODO:	possible issue if buffers from different reaches overlap.
# 		check implications of this
ds_bcs = []
for link in links:
	#####################################################################
	# Check if ds point is on the boundary
	xd,yd = points_ds[link]
	# Check if the downstream link is exiting into the ocean
	# TODO: generalise to allow different fixed height or variable height (for tides/ storm surge)
	dslink = dslinks[link]
	if dslink<0: # ocean boundary
		bc_str = 'P '+str(xd)+' '+str(yd)+' HFIX 0.0'
		ds_bcs.append(bc_str)
	elif xd - res*2 < extent[0]: # West
		# create string setting west boundary condition
		bc_str = 'W '+str(yd-ds_buffer)+' '+str(yd+ds_buffer)+' FREE '+str(slope[link])
		ds_bcs.append(bc_str)
	elif xd + res*2 > extent[1]: # East
		# create string setting east boundary condition
		bc_str = 'E '+str(yd-ds_buffer)+' '+str(yd+ds_buffer)+' FREE '+str(slope[link])
		ds_bcs.append(bc_str)
	elif yd-res*2 < extent[2]: # South
		# create string setting south boundary condition
		bc_str = 'S '+str(xd-ds_buffer)+' '+str(xd+ds_buffer)+' FREE '+str(slope[link])
		ds_bcs.append(bc_str)
	elif yd+res*2 > extent[3]: # North:
		# create string setting north boundary condition
		bc_str = 'N '+str(xd-ds_buffer)+' '+str(xd+ds_buffer)+' FREE '+str(slope[link])
		ds_bcs.append(bc_str)

###############################################################################################
# Load Mizuroute river network (latitude and longitude of upstream points of each segment)
# Also load segId (which is not necessarily the same as the index of the segment)
seg_index_map = {}
with Dataset(f_network,'r') as f:
	segids_mizu = f.variables['segId'][:] # mizuroute ID of each river segment
	lats   = f.variables['lat_up'][:] # upstream lat of each river segment
	lons   = f.variables['lon_up'][:] # upstream lon of each river segment
	# Work out the mizuroute segment index corresponding to lisflood river network 'linkno'
	# Match upstream coordinates for each river segment
	#
	seg_index_map = {}
	for link,point in points_upall.items():
		xx,yy = point
		index,mindist = find_closest_list(yy,xx,lats,lons)
		if mindist < 0.0001:
			seg_index_map[link] = index
			print('mizu,taudem link:',segids_mizu[seg_index_map[link]],link)
		else:
			print('Error finding link in mizuroute:',link,xx,yy)

###############################################################################################
# Main script: Loop over discharge files and write out discharge values
#
fpattern = os.path.join(mizuroute_outdir,mizuRuns)
for f_discharge in glob.glob(fpattern):
	#f_discharge = '/home/pu17449/data2/mizuRoute/merithydro/q_GBM_MERIT-Hydro_1988-1-1.nc'
	fname = os.path.basename(f_discharge)
	print(fname)
	runname = fname[2:-7]
	year = int(runname[-4:])
	print(runname,year)
	f_csv = os.path.join(lf_dischargedir,runname+'.csv')
	f_bdy = os.path.join(lf_dischargedir,runname+'.bdy')
	f_par = os.path.join(lf_pardir,runname+'.par')

	outdir = os.path.join(resultdir,runname)
	maxtif = os.path.join(outdir,runname+'-max.tif') # Only generated after the run is finished
	if os.path.exists(maxtif) and not skip_existing == 'no':
		print('Result already exists, skipping')
		continue
	if os.path.exists(outdir) and skip_existing == 'started':
		print('Result already started, skipping')
		continue

	if not os.path.exists(f_bdy):

		discharge_dict = {}
		# Read mizuroute output file for discharge values
		with Dataset(f_discharge,'r') as f_in:
			linkids = f_in.variables['reachID'][:]
			runoff_vals = f_in.variables['dlayRunoff']
			routed_vals = f_in.variables['IRFroutedRunoff']
			times = f_in.variables['time']
			dates = num2date(times[:],times.units)
			print('debug',)

			# Subset data to be within startdate-endate
			startdate = datetime.datetime(year,startmon,1)
			d0,d1 = calendar.monthrange(year,endmon)
			enddate   = datetime.datetime(year,endmon,d1)
			startindex = np.where(dates==startdate)[0][0]
			endindex = np.where(dates==enddate)[0][0]
			routed_vals = routed_vals[startindex:endindex+1,:]
			runoff_vals = runoff_vals[startindex:endindex+1,:]
			dates = dates[startindex:endindex+1]
			datestrings = map(format_date,dates)
			ntimes = len(dates)

			print('shapes',dates.shape,routed_vals.shape,runoff_vals.shape)


			print('writing to csv file')
			with open(f_csv,'w') as f:
				fwriter = csv.writer(f)
				fwriter.writerow(['linkno(mizu)','x','y','linkno(taudem)']+list(datestrings))
				# Loop over reaches:
				for link in links:
					dslink = dslinks[link]
					uslink = uslinks1[link]
					try:
						seg_index = seg_index_map[link]
					except:
						print('Skipping link not found in mizuroute',link)
						continue

					if uslink == -1: # headwater catchment (no upstream link)
						# split up discharge depending on up/downstream accumulation
						upacc = acc_upclip[link]
						downacc = acc_ntd[link]
						runoff = runoff_vals[:,seg_index]
						top_flow = runoff*(upacc/downacc)
						bottom_flow = runoff*(1-upacc/downacc)

						# Write flow at top of reach
						x,y = points_upclip[link]
						row = [str(segids_mizu[seg_index])+'_top',x,y,link]+list(top_flow)
						fwriter.writerow(row)

						# Write flow at downstream of reach
						x,y = points_ntd[link]
						row = [segids_mizu[seg_index],x,y,link]+list(bottom_flow)
						fwriter.writerow(row)

					elif uslink not in links:
						# Only part of this catchment is within the domain
						# Set upstream flow to the routed discharge/runoff
						# The routed flow is put in at the upstream boundary
						# A fraction of the runoff is put in at the upstream/downstream boundaries
						runoff = runoff_vals[:,seg_index]
						routed = routed_vals[:,seg_index]
						# use the upstream (outside and clipped) accumulation to determine the fraction of catchment outside of the domain
						upacc1 = acc_upclip[link]
						upacc2 = acc_upall[link]
						downacc = acc_ntd[link]
						top_flow = routed - runoff * (downacc - upacc1) / (downacc - upacc2)
						bottom_flow = runoff * (downacc - upacc1) / (downacc - upacc2)

						# Write flow at top of reach
						x,y = points_upclip[link]
						row = [str(segids_mizu[seg_index])+'_top',x,y,link]+list(top_flow)
						fwriter.writerow(row)

						# Write flow at next-to-downstream point of reach
						x,y = points_ntd[link]
						row = [segids_mizu[seg_index],x,y,link]+list(bottom_flow)
						fwriter.writerow(row)

					else: # This is a normal reach within the domain
						# Add the runoff from this reach to the downstream point
						runoff = runoff_vals[:,seg_index]
						# Write flow at next-to-downstream point of reach
						x,y = points_ntd[link]
						row = [segids_mizu[seg_index],x,y,link]+list(runoff)
						fwriter.writerow(row)


		print('writing to bdy file')
		write_bdy(f_bdy,f_csv,ntimes)
	else:
		print('file exists, skipping:',f_bdy)

	# Write bci file (just needs the list of points, so can use any csv file, only create once)
	if not os.path.exists(f_bci):
		print('writing to bci file')
		write_bci_v2(f_bci,f_csv,ds_bcs)

	bci_relative = f_bci[len(lfdir)+1:]
	bdy_relative = f_bdy[len(lfdir)+1:]

	if ntimes is None:
		with open(f_csv,'r') as f: tmp = f.readline()
		ntimes = len(tmp.split(','))-4
	sim_time = str(ntimes*86400)


	if not os.path.exists(f_par):
		# copy template and modify
		shutil.copy(f_par_template,f_par)
		sed_expr = 's#<RUNNAME>#'+runname+'#g; '
		sed_expr += 's#<SIMTIME>#'+sim_time+'#g; '
		sed_expr += 's#<BCI>#'+bci_relative+'#g; '
		sed_expr += 's#<BDY>#'+bdy_relative+'#g'
		print(sed_expr)
		cmd = ['sed','-i','-e',sed_expr ,f_par]
		subprocess.call(cmd)


	# add par to sublist
	sublist.append(f_par)

	if len(sublist)>=simsperjob:
		print('Submitting jobs',len(sublist))
		os.environ['CONTROL_FLIST']=':'.join(sublist)
		print(os.environ['CONTROL_FLIST'])
		print(' '.join(qsub_cmd))
		if not dryrun:
			subprocess.call(qsub_cmd)
		else:
			print('Dry run, not submitting')
		# reset sublist
		sublist = []

# submit any that are left
if len(sublist)>0:
	print('Submitting jobs',len(sublist))
	os.environ['CONTROL_FLIST']=':'.join(sublist)
	print(os.environ['CONTROL_FLIST'])
	print(' '.join(qsub_cmd))
	if not dryrun:
		subprocess.call(qsub_cmd)
	else:
		print('Dry run, not submitting')
else:
	print('No more simulations to submit!')
