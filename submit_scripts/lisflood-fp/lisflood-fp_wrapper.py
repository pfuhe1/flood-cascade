#!/cm/shared/languages/python-3.3.2/bin/python
# submit script for running lisflood runs in parallel
# Peter Uhe May 15 2019
# V3 2020/02/04: using cuncurrent.futures.ProcessPoolExecutor

import os,subprocess
# trying concurrent.futures to handle lisflood segfaults etc. better
from concurrent.futures import ProcessPoolExecutor
import datetime
from convert_output_totiff import convert_to_tif_v4
from contextlib import redirect_stdout

def call_subproc(cmd,sim_name,logfile):
	print('command',cmd)
	print('log',logfile)
	ret = subprocess.call(cmd,stdout=open(logfile,'w'),stderr=subprocess.STDOUT)
	# converts to tiff, relies on output being in folder relative to working directory
	convert_to_tif_v2('results/'+sim_name,sim_name)
	return ret

# v2: specify jobsize for multiprocess conversion to tiff
#     specify outdir (where the output files should be)
#	  redirects stdout to logfile
def call_subproc_v2(cmd,sim_name,outdir,logfile,jobsize):
	print('command',cmd)
	print('log',logfile)
	ret = subprocess.call(cmd,stdout=open(logfile,'w'),stderr=subprocess.STDOUT)
	# Append convert_to_tiff output to logfile
	with open(logfile,'a') as f:
		with redirect_stdout(f):
			convert_to_tif_v4(outdir,sim_name,jobsize=jobsize)


# Print start time
print('Start:',datetime.datetime.now())

control_files = os.environ['CONTROL_FLIST'].split(':')
exefile = os.environ['EXE_FILE']
# Work out number of concurrent processes to run, based on node size and number of processes per job
nodesize = int(os.environ['NCPUS'])
jobsize = int(os.environ['OMP_NUM_THREADS'])
logdir = os.environ['LOGDIR']
resultdir = os.environ['RESULTDIR'] # The results for each simulations will be a subdirectory
numprocesses = int(nodesize/jobsize)

print('Running',numprocesses,'simultaneous processes using',jobsize,'cores each')
print('running simulations',len(control_files))
print(os.environ['CONTROL_FLIST'])
with ProcessPoolExecutor(max_workers=numprocesses) as pool:# since node has 16 cores, and OMP_NUM_THREADS=2

	for control_file in control_files:
		# Todo, could add check if this simulation has already been run
		fname = os.path.basename(control_file)
		sim_name =fname[:-4]
		outdir = os.path.join(resultdir,sim_name)
		logfile = os.path.join(logdir,sim_name+'.log')
		cmd = ['time',exefile,'-v',control_file]
		#ret = pool.apply_async(subprocess.call,cmd,{'stdout':open(logfile,'w') ,'stderr':subprocess.STDOUT})
		#subprocess.call(cmd,stdout=open(logfile,'w'),stderr=subprocess.STDOUT)
		ret = pool.submit(call_subproc_v2,cmd,sim_name,outdir,logfile,jobsize)

print(ret.result())
print('Finished:',datetime.datetime.now())
