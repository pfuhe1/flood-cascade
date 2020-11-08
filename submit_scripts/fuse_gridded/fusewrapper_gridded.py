#!/cm/shared/languages/python-3.3.2/bin/python
# submit script for FUSE GBM runs
# Peter Uhe 29 Oct 2019
# 

# call this script from 'setup_GBM_runs_mswep_p5degcalib' which creates a qsub job to submit to the HPC queue
# This script is actually called from 'call_pythonscript.sh' (which is needed to load the modules before calling the script)

import os,glob,subprocess,sys,shutil,multiprocessing
import datetime

def call_subproc(cmd,logfile):
	with open(logfile,'w') as flog:
		# Write command into the log file before running
		flog.write('Running command: '+' '.join(cmd)+'\n\n')
		subprocess.call(cmd,stdout=flog,stderr=subprocess.STDOUT)

# Print start time 
print('Starting',datetime.datetime.now())

#fm_files = sys.argv[1:]
fm_files = os.environ['FM_FLIST'].split(':')
logdir = os.environ['LOGDIR']
ncpus = int(os.environ['NCPUS'])
fuseexe = os.environ['FUSE_EXE']
print('Environment: NCPUS',ncpus,'LOGDIR',logdir,'fuseexe',fuseexe)
print('running simulations',len(fm_files))
print(os.environ['FM_FLIST'])
pool = multiprocessing.Pool(processes=ncpus)

for fm_file in fm_files:
	# Todo, could add check if this simulation has already been run
	fname = os.path.basename(fm_file)
	settingsdir = os.path.dirname(fm_file)
	sim_name = fname[3:-4]
	# Sim name has format FUSEsetup_FUSEdecison_FUSEcalib_<GCM-runstring>
#	setup,dec,calib = sim_name.split('_')
	tmp = sim_name.split('_')
	setup = tmp[0]
	calib = tmp[2]
	calib_file = '_'.join(tmp[:3])+'.nc'
#	calib_file = sim_name+'.nc' # Calib file is in output directory
	logfile = os.path.join(logdir,sim_name+'.log')
	if calib == 'rundef':
		cmd = ['time',fuseexe,fm_file,setup,'run_def']
	else:
		cmd = ['time',fuseexe,fm_file,setup,'run_pre_dist',calib_file]
	print('command',cmd)
	print('log',logfile)
	#ret = pool.apply_async(subprocess.call,cmd,{'stdout':open(logfile,'w') ,'stderr':subprocess.STDOUT})
	#subprocess.call(cmd,stdout=open(logfile,'w'),stderr=subprocess.STDOUT)
	ret = pool.apply_async(call_subproc,[cmd,logfile])
pool.close()
pool.join()

# Print end time 
print('Finished',datetime.datetime.now())
