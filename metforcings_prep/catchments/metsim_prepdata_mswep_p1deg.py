import os,sys,subprocess,glob,datetime
from netCDF4 import Dataset,num2date
import multiprocessing
from catchment_obsforcings_func import merge_obs

# Loops over files for individual catchments
# Run Metsim using the inputs from script `metsim_splitstate_ERA5_mswep_separatecatchments.py`
# Then creates FUSE input files

def write_metsim_conf(state,forcing,domain,f_conf,outdir,outprefix):

	# First get start and end of forcing data
	with Dataset(forcing,'r') as f:
		timevar = f.variables['time']
		s = num2date(timevar[0],timevar.units)
		start = s.strftime('%Y-%m-%d')
		start2 = s.strftime('%Y%m%d')
		e   = num2date(timevar[-1],timevar.units)
		end = e.strftime('%Y-%m-%d')
		end2 = e.strftime('%Y%m%d')

	with open(f_conf,'w') as fout:
		fout.write('[MetSim]\n\n')

		fout.write('# Time step in minutes\n')
		fout.write('time_step = 1440\n')
		fout.write('#Forcings begin here (year/month/day:hour\n')
		fout.write('start = '+start+'\n')
		fout.write('# Forcings end at this date (year/month/day)\n')
		fout.write('stop = '+end+'\n\n')

		fout.write('# Input specification\n')
		fout.write('forcing = '+forcing+'\n')
		fout.write('state  = '+state+'\n')
		fout.write('domain = '+domain+'\n')
		fout.write('forcing_fmt = netcdf\ndomain_fmt = netcdf\nstate_fmt = netcdf\n\n')

		fout.write('#Output specification\n')
		fout.write('out_fmt = netcdf\n')
		fout.write('out_dir = '+outdir+'\n')
		fout.write('out_prefix = '+outprefix+'\n')
		fout.write('out_precision = f8\n')
		fout.write('out_vars = [\'pet\']\n')
		fout.write('utc_offset = True\n\n')

		fout.write('# How to disaggregate\n')
		fout.write('method = mtclim\n\n')

		fout.write('[forcing_vars]\n')
		fout.write('prec = pr\n')
		fout.write('t_max = tasmax\n')
		fout.write('t_min = tasmin\n\n')

		fout.write('[state_vars]\n')
		fout.write('prec = pr\n')
		fout.write('t_max = tasmax\n')
		fout.write('t_min = tasmin\n\n')

		fout.write('[domain_vars]\n')
		fout.write('lat = lat\n')
		fout.write('lon = lon\n')
		fout.write('mask = mask\n')
		fout.write('elev = elev\n\n')

		fout.write('[chunks]\n')
		fout.write('lat = 1\n')
		fout.write('lon = 1\n')

	return os.path.join(outdir,outprefix+'_'+start2+'-'+end2+'.nc')

################################################################
#
# Function to process a single run
def process_run(grdcid,datadir):

	# Directories
	indir = os.path.join(datadir,'metsim-forcing')
	domaindir = os.path.join(datadir,'metsim-domain')
	confdir = os.path.join(datadir,'metsim-conf')
	outdir = os.path.join(datadir,'pet')
	fusedir = os.path.join(datadir,'fuse-forcing')
	qdir = '/export/anthropocene/array-01/pu17449/FUSE_inputs/catchment_obsdischarge'
	obsname = 'mswep-p1deg'

	clim_start = datetime.datetime(1979,4,1)
	clim_end = datetime.datetime(2017,10,31)

	# Final Output file
	#f_fuse = os.path.join(fusedir,'grdc_'+grdcid+'_mswep-p1deg_input.nc')
	#if not os.path.exists(f_fuse):
	fpattern_out = os.path.join(fusedir,grdcid+'_catchment_'+obsname+'_*.nc')
	print(fpattern_out,len(glob.glob(fpattern_out)))
	return
	if len(glob.glob(fpattern_out))<1:

		################################################################
		# Run metsim to produce PET

		# Input files
		fforcing = os.path.join(indir,grdcid+'_'+'forcing_19790401-20171031.nc')
		fstate = os.path.join(indir,grdcid+'_state_19790101-19790331.nc')
		fdomain = os.path.join(domaindir,grdcid+'_domain.nc')

		# Output files
		fconf = os.path.join(confdir,'metsim_catchment-'+grdcid+'.conf') #File to write metsim conf to
		if not os.path.exists(confdir):
			os.mkdir(confdir)

		pet_prefix = 'pet_'+grdcid+'_MSWEP-ERA5-p1deg' # Prefix for pet output file
		fpet = write_metsim_conf(fstate,fforcing,fdomain,fconf,outdir,pet_prefix)

		if not os.path.exists(fpet):
			if not os.path.exists(outdir):
				os.mkdir(outdir)
			metsim_cmd = ['ms',fconf,'-v','-n','1']
			print(' '.join(metsim_cmd))
			subprocess.call(metsim_cmd)

		###########################################################################################
		# NOW set up FUSE inputs:
		tmpdir = os.path.join(datadir,'tmp')
		if not os.path.exists(tmpdir):
			os.mkdir(tmpdir)

		################################################################
		# tas (input to FUSE): days 91 onwards to match metsim output
		tasdir = os.path.join(datadir,'tas')
		ftas = grdcid+'_tas_19790101-20171231.nc'
		tas_path = os.path.join(tasdir,ftas)
		tas_tmp = os.path.join(tmpdir,ftas)
		cdo_cmd = ['cdo','-L','subc,273.15','-seltimestep,91/99999999',tas_path,tas_tmp]
		print(cdo_cmd)
		subprocess.call(cdo_cmd)
		nco_cmd = ['ncatted','-a','units,tas,m,c,degC',tas_tmp]
		print(nco_cmd)
		subprocess.call(nco_cmd)

		################################################################
		# pr (input to FUSE): days 91 onwards to match metsim output
		prdir  = os.path.join(datadir,'pr')
		fpr  = grdcid+'_pr_19790101-20171231.nc'
		pr_path = os.path.join(prdir,fpr)
		pr_tmp  = os.path.join(tmpdir,fpr)
		cdo_cmd = cdo_cmd = ['cdo','-seltimestep,91/99999999',pr_path,pr_tmp]
		print(cdo_cmd)
		subprocess.call(cdo_cmd)

		################################################################
		# merge files for input to FUSE
		if not os.path.exists(fusedir):
			os.mkdir(fusedir)

		qobs_path = os.path.join(qdir,grdcid+'_Q_Day.nc')
		f_fuse = merge_obs(grdcid,obsname,qobs_path,tas_tmp,pr_tmp,fpet,fpattern_out,fusedir,tmpdir,clim_start,clim_end)

		################################################################
		# Clean up temp files
		tempfiles = [tas_tmp,pr_tmp]
		for fpath in tempfiles:
			if fpath is not None and os.path.exists(fpath):
				os.remove(fpath)
	else:
		print('File already exists, skipping:',fpattern_out)
		return

################################################################
#
# Main part of the script
if __name__ == '__main__':

	#############################################################################
	# Input paths
	fcatchment_list = '/export/anthropocene/array-01/pu17449/FUSE_inputs/catchments_p1deg/GBM-p1deg_donorsbest3.txt'
	datadir = '/export/anthropocene/array-01/pu17449/FUSE_inputs/catchments_p1deg_v3'

	#############################################################################
	# Loop over catchmentsn and calculate forcings
	pool = multiprocessing.Pool(processes=3)
	with open(fcatchment_list,'r') as f:
		for line in f:
			grdcid = line.strip()
			print(grdcid)
			ret = pool.apply_async(process_run,[grdcid,datadir])
			#break # for debugging a single run

	pool.close()
	pool.join()
	print(ret.get())
