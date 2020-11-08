# script to setup and submit FUSE gridded runs
# Peter Uhe Aug 21 2019
#

# First call script to produce different calibration files (scripts in fuse_param_maps folder)
# Unless running with default parameter sets, then set calib_choices = ['rundef']

import os,glob,subprocess,sys,shutil
import datetime
import queue
import socket


# IDs of FUSE structures/decisions to use
fuse_decision_ids = [900,902,904]

# Name given to FUSE regional model domain
setup_name = 'GBM-p1deg'

# Name given to observational dataset used as input
obsversion = 'MSWEP2-2-ERA5'

# Name given to observational dataset used for calibration (usually the same as obsversion)
calibversion = 'MSWEP2-2-ERA5'

# Start and end date of simulations
ssim = '1979-04-01'
esim = '2017-10-31'


calib_choices = ['rundef',calibversion+'-calibrated1',calibversion+'-calibrated2',calibversion+'-calibrated3']

#calib_choices = [obsversion+'-calibrateRand0001',obsversion+'-calibrateRand0002',obsversion+'-calibrateRand0003',obsversion+'-calibrateRand0004',obsversion+'-calibrateRand0005',obsversion+'-calibrateRand0006',obsversion+'-calibrateRand0007',obsversion+'-calibrateRand0008',obsversion+'-calibrateRand0009',obsversion+'-calibrateRand0010',obsversion+'-calibrateRand0011',obsversion+'-calibrateRand0012',obsversion+'-calibrateRand0013',obsversion+'-calibrateRand0014',obsversion+'-calibrateRand0015',obsversion+'-calibrateRand0016',obsversion+'-calibrateRand0017',obsversion+'-calibrateRand0018',obsversion+'-calibrateRand0019',obsversion+'-calibrateRand0020']]

override = True

###################################################################################
# Define computer specific paths/variables
if host[:3] == 'bp1':# blue pebble cluster
	fm_template = '/home/pu17449/src/fuse_templates/fm_genforcing_template.txt'
	basedir = '/work/pu17449/fuse/'+setup_name
	ncpus = 1
	nperjob = 1
	fuse_exe =  '/home/pu17449/src/fuse_dev/bin/fuse_update_noconsistency.exe'
	qsub_script = '/home/pu17449/src/setup_scripts/fuse_GBM/call_pythonscript_bp1.sh'
	qsub_command    = ['qsub','-l','select=1:ncpus='+str(ncpus)+':mem='+str(ncpus*8)+'gb','-v','FM_FLIST,LOGDIR,FUSE_EXE',qsub_script] # NCPUs environment variable is set by PBS

###################################################################################
# Fuse folder paths

settingsdir = os.path.join(basedir,'settings')
inputdir    = os.path.join(basedir,'input')
outputdir    = os.path.join(basedir,'output')
logdir    = os.path.join(basedir,'logs')
if not os.path.exists(logdir):
	os.mkdir(logdir)

###################################################################################
# Export environment variables used by all QSUB commands
os.environ['NCPUS'] = str(ncpus)
os.environ['FUSE_EXE'] = fuse_exe
os.environ['LOGDIR'] = logdir


###################################################################################
sublist = []
# Loop over runs to submit and add fm_file to sublist for submission
for dec in fuse_decision_ids:
	for calib in calib_choices:
		sim_name = setup_name+'_'+str(dec)+'_'+calib+'_'+obsversion
		print(sim_name)

		if calib == 'rundef':
			out_file = os.path.join(outputdir,sim_name+'runs_def.nc')
		else:
			out_file = os.path.join(outputdir,sim_name+'runs_pre_dist.nc')

		# Check if this run has already been computed
		if not os.path.exists(out_file) or override:

			fm_file = os.path.join(settingsdir,'fm_'+sim_name+'.txt')

			if os.path.exists(fm_file) and override:
				os.remove(fm_file)

			# copy fuse filemanager template and modify
			if not os.path.exists(fm_file):
				shutil.copy(fm_template,fm_file)
				sed_expr =  's#<settings_dir>#'+settingsdir+'#g; '
				sed_expr += 's#<input_dir>#'+inputdir+'#g; '
				sed_expr += 's#<output_dir>#'+outputdir+'#g; '
				sed_expr += 's#<decid>#'+str(dec)+'#g; '
				sed_expr += 's#<calibid>#'+calib+'#g; '
				sed_expr += 's#<forcing>#'+obsversion+'#g; '
				sed_expr += 's#<s_sim>#'+ssim+'#g; '
				sed_expr += 's#<e_sim>#'+esim+'#g; '
				sed_expr += 's#<s_eval>#'+ssim+'#g; '
				sed_expr += 's#<e_eval>#'+esim+'#g'
				cmd = ['sed','-i','-e',sed_expr ,fm_file]
				print(cmd)
				subprocess.call(cmd)

			# add fm_file to sublist
			sublist.append(fm_file)

			# submit nperjob simulations at a time
			if len(sublist)==nperjob:
				print('Submitting jobs',len(sublist))
				# First export environment variables used in the job
				os.environ['FM_FLIST']=':'.join(sublist)
				print(os.environ['FM_FLIST'])
				subprocess.call(qsub_command)
				# Reset submit list
				sublist = []

		else:
			print('Output file already exists, skipping:',out_file)


# submit any that are left
if len(sublist)>0:
	print('Submitting jobs',len(sublist))
	os.environ['FM_FLIST']=':'.join(sublist)
	print(os.environ['FM_FLIST'])
	# Submit to compute queue
	subprocess.call(qsub_command)
	#subprocess.call([qsub_script]) # call on login node for debugging
else:
	print('No more simulations to submit!')
