#!/cm/shared/languages/python-3.3.2/bin/python
# submit script for FUSE calibration of GRDC catchments
# Peter Uhe May 15 2019
#

# call this script from 'setup_grdc_catch_best73.py which creates a qsub job to submit to the HPC queue
# This script is actually called from 'call_pythonscript.sh' (which is needed to load the python module before calling the script)

import os,glob,subprocess,sys,shutil,multiprocessing
import datetime

# Function to call subprocess and send output to log file
def call_subproc(cmd,logfile):
	subprocess.call(cmd,stdout=open(logfile,'w'),stderr=subprocess.STDOUT)

print('Start',datetime.datetime.now())

#fm_files = sys.argv[1:]
fm_files = os.environ['FM_FLIST'].split(':')
ncpus = int(os.environ['NCPUS'])
datadir = os.environ['DATADIR']
fuseexe = os.environ['FUSE_EXE']
print('running simulations',len(fm_files))
print(os.environ['FM_FLIST'])
pool = multiprocessing.Pool(processes=ncpus)

for fm_file in fm_files:
	# Todo, could add check if this simulation has already been run
	fname = os.path.basename(fm_file)
	grdcid = fname.split('_')[2]
	sim_name = 'grdc_'+grdcid
	logfile = os.path.join(datadir,'fuse_'+sim_name,'logs',fname[3:-4]+'.log')
	cmd = ['time',fuseexe,fm_file,sim_name,'calib_sce']
	print('command',cmd)
	print('log',logfile)
	#ret = pool.apply_async(subprocess.call,cmd,{'stdout':open(logfile,'w') ,'stderr':subprocess.STDOUT})
	#subprocess.call(cmd,stdout=open(logfile,'w'),stderr=subprocess.STDOUT)
	ret = pool.apply_async(call_subproc,[cmd,logfile])
pool.close()
pool.join()

print('Finished',datetime.datetime.now())
