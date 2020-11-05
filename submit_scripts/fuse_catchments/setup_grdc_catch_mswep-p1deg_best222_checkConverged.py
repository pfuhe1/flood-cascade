# script to setup FUSE calibration of GRDC catchments
# Peter Uhe May 15 2019
#
# Assumes catchment_obsforcings_toBC3.py script has been run to copy the elev_bands and catchment_obsforcing files for each catchment
# Also requires list of catchment ids to process (f_catchlist)

import os,glob,subprocess,sys,shutil,socket
import datetime
import queue

# Check see output file:
# Relies that the second last line of a successful calibration starts with 'CRITERION'
def success_sce_output(f_sce,nfromend=4):
	lprev=queue.Queue(maxsize=nfromend-1)
	with open(f_sce,'r') as f:
		for line in f:
			if lprev.full():
				checkstr = lprev.get()
			lprev.put(line.strip())
	#print(checkstr)
	if checkstr == '*** PRINT THE FINAL PARAMETER ESTIMATE AND ITS CRITERION VALUE ***':
		return True
	else:
		return False

def check_converged(f_sce):
	# Search for 'EXCEEDED' in sce output
	cmd = ['grep','EXCEEDED',f_sce]
	ret = subprocess.call(cmd) # Returns 0 if found exceeed
	return  ret

# User choice
fuse_decision_id = 902


host=socket.gethostname()
# Define paths (so scripts can be set up on different clusters)
if host[:3]=='bp1': # blue pebble
	templatedir = '/home/pu17449/src/fuse_templates'
	settingsdir = os.path.join(templatedir,'settings_dir')
	data_dir    = '/work/pu17449/fuse/p1deg_catchments/'
	# List of catchments to calibrate (produced by fuse_regionalization script)
	f_catchlist = os.path.join(templatedir,'GBM-p1deg_donorsbest3.txt')
	qsub_script = '/home/pu17449/src/setup_scripts/fuse_catchments/call_pythonscript_bp1.sh'
	# NOTE: I compiled an optimized version of fuse.exe, but that seems to cause problems with the sce parameter estimation
	#fuse_exe    = '/home/pu17449/src/fuse_dev/bin/fuse.exe'
	fuse_exe    = '/home/pu17449/src/fuse_dev/bin/fuse_update_noconsistency.exe'
	ncpus       = 1
	qsub_command    = ['qsub','-l','select=1:ncpus='+str(ncpus)+':mem='+str(ncpus)+'gb','-v','FM_FLIST,NCPUS,DATADIR,FUSE_EXE']

sublist = []

with open(f_catchlist,'r') as f:
	for line in f:
		grdcid=line.strip()
		sim_name0 = 'grdc_'+grdcid

		print(grdcid)
		sim_dir = os.path.join(data_dir,'fuse_grdc_'+grdcid)
		print(sim_dir)
		fm_template = os.path.join(templatedir,'fm_catch_grdc_template_25000loops_mswep-p1deg.txt')
		fm_file = os.path.join(sim_dir,'fm_grdc_'+grdcid + '_dec_'+str(fuse_decision_id)+'_mswep-p1deg.txt')
		# Paths for catchment
		idir = os.path.join(sim_dir,'input')
		if not os.path.exists(idir): os.mkdir(idir)
		odir = os.path.join(sim_dir,'output')
		if not os.path.exists(odir):
			os.mkdir(odir)
		else: # Check if the calibration has already been run scucessfully
			fsce_out = os.path.join(odir,sim_name0+'_'+str(fuse_decision_id)+'_mswep-p1deg_sce_output.txt')
			print(fsce_out)
			if os.path.exists(fsce_out) and success_sce_output(fsce_out):
				if check_converged(fsce_out):
					print('Calibration already completed, skipping')
					continue
				else: # Try more iterations
					print('Calibration previously did not converge')
					continue # go to next catchment
		# Overwrite existing FM_FILE (needed to do a clean run)
		if os.path.exists(fm_file):
			os.remove(fm_file)
		sdir = os.path.join(sim_dir,'settings')
		if not os.path.exists(sdir): shutil.copytree(settingsdir,sdir)
		ldir = os.path.join(sim_dir,'logs')
		if not os.path.exists(ldir): os.mkdir(ldir)

		# Input files
		elevs_file = os.path.join(sim_dir,'catchment_'+grdcid+'_elev_bands.nc')
		elevs_file_ln = os.path.join(idir,sim_name0+'_elev_bands.nc')
		force_file = glob.glob(os.path.join(sim_dir,grdcid+'_catchment_mswep-p1deg_*.nc'))[0]
		force_file_ln = os.path.join(idir,sim_name0+'_mswep-p1deg_input.nc')

		sdate,edate = os.path.basename(force_file).split('_')[-1][:-3].split('-')
		sdatetime = datetime.datetime(int(sdate[:4]),int(sdate[4:6]),int(sdate[6:8]))
		sdatestring = sdatetime.strftime("%Y-%m-%d")
		edatetime = datetime.datetime(int(edate[:4]),int(edate[4:6]),int(edate[6:8]))
		edatestring = edatetime.strftime("%Y-%m-%d")
		#print('Full period',sdatetime,edatetime)
		# For calibration/evaluation period, just split in two
		days = (edatetime-sdatetime).days
		s_eval = sdatetime+datetime.timedelta(days/2)
		s_evalstring = s_eval.strftime("%Y-%m-%d")
		e_eval = edatetime
		e_evalstring = e_eval.strftime("%Y-%m-%d")
		#print('Eval period',s_eval,e_eval)

		# copy fuse filemanager template and modify
		if not os.path.exists(fm_file):
			shutil.copy(fm_template,fm_file)
			sed_expr = 's/<grdcid>/'+grdcid+'/g; '
			sed_expr += 's/<decid>/'+str(fuse_decision_id)+'/g; '
			sed_expr += 's/<s_sim>/'+sdatestring+'/g; '
			sed_expr += 's/<e_sim>/'+edatestring+'/g; '
			sed_expr += 's/<s_eval>/'+s_evalstring+'/g; '
			sed_expr += 's/<e_eval>/'+e_evalstring+'/g'
			cmd = ['sed','-i','-e',sed_expr ,fm_file]
			subprocess.call(cmd)
			#cmd = ['sed','-i','s/<decid>/'+str(fuse_decision_id)+'/g',fm_file]
			#subprocess.call(cmd)
			#cmd = ['sed','-i',,fm_file]
			#subprocess.call(cmd)
			#cmd = ['sed','-i','s/<e_sim>/'+edatestring+'/g',fm_file]
			#subprocess.call(cmd)
			#cmd = ['sed','-i','s/<s_eval>/'+s_evalstring+'/g',fm_file]
			#subprocess.call(cmd)
			#cmd = ['sed','-i','s/<e_eval>/'+e_evalstring+'/g',fm_file]
			#subprocess.call(cmd)

		# Rename input files to what fuse expects
		if not os.path.exists(elevs_file_ln):
			os.symlink(elevs_file,elevs_file_ln)
		if not os.path.exists(force_file_ln):
			os.symlink(force_file,force_file_ln)
		# Extra command to rename lat/lon if necessary
		#cmd = ['ncrename','-v','lat,latitude','-v','lon,longitude','-d','lat,latitude','-d','lon,longitude',force_file]
		#subprocess.call(cmd)

		# add fm_file to sublist
		sublist.append(fm_file)

		# submit ncpu files at a time
		if len(sublist)==ncpus:
			print('Submitting jobs',len(sublist))
			# Set list of fm files as environment variable
			os.environ['FM_FLIST']=':'.join(sublist)
			os.environ['NCPUS'] = str(ncpus)
			os.environ['FUSE_EXE'] = fuse_exe
			os.environ['DATADIR'] = data_dir
			print(os.environ['FM_FLIST'])
			subprocess.call(qsub_command+[qsub_script])
			# Reset submit list
			sublist = []


# submit any that are left
if len(sublist)>0:
	print('Submitting jobs',len(sublist))
	# First export environment variables used in the job
	os.environ['FM_FLIST']=':'.join(sublist)
	os.environ['NCPUS'] = str(ncpus)
	os.environ['FUSE_EXE'] = fuse_exe
	os.environ['DATADIR'] = data_dir
	print(os.environ['FM_FLIST'])
	#print(' '.join(qsub_command+[qsub_script]))
	subprocess.call(qsub_command+[qsub_script])
else:
	print('No more simulations to submit!')
