#!/cm/shared/languages/python-3.3.2/bin/python
# submit script for submission of mizuRoute simualtions
# Peter Uhe Oct 29 2019
#

# call this script from 'run_mizuRoute_templated_mswep050calib.py which creates a qsub job to submit to the HPC queue
# This script is actually called from 'call_pythonscript.sh' (which is needed to load modules before calling the script)

import os,glob,subprocess,sys,shutil,multiprocessing
import datetime

def call_subproc(cmd,logfile):
	subprocess.call(cmd,stdout=open(logfile,'w'),stderr=subprocess.STDOUT)

# Print start time
print('Starting:',datetime.datetime.now())

# Get environment variables
control_files = os.environ['CONTROL_FLIST'].split(':')
logdir = os.environ['LOGDIR']
ncpus = int(os.environ['NCPUS'])
mizuexe = os.environ['MIZU_EXE']

print('running simulations',len(control_files))
print(os.environ['CONTROL_FLIST'])
pool = multiprocessing.Pool(processes=ncpus)

for control_file in control_files:
	# Todo, could add check if this simulation has already been run
	fname = os.path.basename(control_file)
	sim_name =fname[8:-4]
	logfile = os.path.join(logdir,sim_name+'.log')
	cmd = ['time',mizuexe,control_file]
	print('command',cmd)
	print('log',logfile)
	#ret = pool.apply_async(subprocess.call,cmd,{'stdout':open(logfile,'w') ,'stderr':subprocess.STDOUT})
	#subprocess.call(cmd,stdout=open(logfile,'w'),stderr=subprocess.STDOUT)
	ret = pool.apply_async(call_subproc,[cmd,logfile])
pool.close()
pool.join()

print('Finished:',datetime.datetime.now())
