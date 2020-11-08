# script to setup mizuroute simulations and submit to queue
# This version for MSWEP2-2 p1deg simulations
# Peter Uhe April 23 2020
#

import os,glob,subprocess,sys,shutil,socket
#exec(open('/newhome/pu17449/src/setup_scripts/module_python3_init').read())
#module('load','apps/cdo-1.9.3')

host = socket.gethostname()
##################################################################################
# Inputs from FUSE simulations:

runlength = 38 # years to run (for checking complete output)


# List of decision IDS and calibration choices to loop over
# These need to correspond to FUSE simulations
# FUSE simulations have the naming format: <setup-name>_<decision-id>_<calib-choice>
setup_name = 'GBM-p1deg'
calibversion = 'MSWEP2-2-ERA5' # observations used in calibration
obsversion    = 'MSWEP2-2-ERA5'      # obsversion data used in FUSE run (can be obs or modeled)
fuse_decision_ids = [900,902,9042]
#calib_choices = ['']
#calib_choices = ['rundef','calibrated1','calibrated2','calibrated3']
#calib_choices  = ['longcalib1','longcalib2','longcalib3']
calib_choices = ['calibrateRand0001','calibrateRand0002','calibrateRand0003','calibrateRand0004','calibrateRand0005', 'calibrateRand0006','calibrateRand0007','calibrateRand0008','calibrateRand0009','calibrateRand0010', 'calibrateRand0011','calibrateRand0012','calibrateRand0013','calibrateRand0014','calibrateRand0015', 'calibrateRand0016','calibrateRand0017','calibrateRand0018','calibrateRand0019']

# initialise dictionary of input runs (name and input file)
input_runs = {}

##################################################################################
# Mizuroute configuration and paths

if host[:3]=='bp1': #blue pebble
	mizudir = '/work/pu17449/mizuRoute'
	inbase  = '/work/pu17449/fuse/GBM-p1deg/output/'
	control_template = '/home/pu17449/src/mizuRoute/route/settings/GBM-MERIT_MSWEPp1deg.control_template'
	qsub_script = '/home/pu17449/src/setup_scripts/mizuRoute/call_pythonscript_bp1.sh' # use this one so we can load modules
	ncpus       = 1
	nperjob  = 1
	qsub_command    = ['qsub','-l','select=1:ncpus='+str(ncpus)+':mem='+str(ncpus)+'gb','-v','CONTROL_FLIST,LOGDIR,MIZU_EXE',qsub_script]
	mizuexe = '/home/pu17449/src/mizuRoute/route/bin/mizuRoute'

outbase    = os.path.join(mizudir,'output')
controldir = os.path.join(mizudir,'control_files')
logdir     = os.path.join(mizudir,'logs')


if not os.path.exists(outbase):
	os.mkdir(outbase)
if not os.path.exists(controldir):
	os.mkdir(controldir)
if not os.path.exists(logdir):
	os.mkdir(logdir)

# User choice
override = False

# Initialize list of control files
sublist = []


##################################################################################
#
# Loop over decision ids and calib choices
for dec in fuse_decision_ids:
	for calib in calib_choices:
		# First do some preprocessing to remove extra dimension and select q_instnt variable
		# TODO: for large sets of runs, could move this command to the qsub script rather than here
		#
		#if calib == '':# runs def (first version)
			#sim_name = setup_name+'_'+str(dec)+'_'+obsversion
			#sim_name = sim_name[:-1]# remove last '_' from name
		if calib == 'rundef': # FUSE default parameters
			sim_name = setup_name+'_'+str(dec)+'_rundef_'+obsversion
			infile = os.path.join(inbase,sim_name+'_runs_def.nc')
			infile2 = os.path.join(inbase,sim_name+'_runs_def_qinst.nc')
		else:
			sim_name = setup_name+'_'+str(dec)+'_'+calibversion+'-'+calib+'_'+obsversion
			infile = os.path.join(inbase,sim_name+'_runs_pre_dist.nc')
			infile2 = os.path.join(inbase,sim_name+'_runs_pre_dist_qinst.nc')

		if os.path.exists(infile2):
			ret = 0
		elif os.path.exists(infile) and not os.path.exists(infile2):
			cdo_cmd = ['cdo','--reduce_dim','selvar,q_instnt',infile,infile2]
			print(cdo_cmd)
			ret = subprocess.call(cdo_cmd)
		else:
			print('Error, infile doesnt exist',infile)
			ret = -1
		if ret == 0:
			input_runs[sim_name] = infile2
		else:
			print('Error running cdo command on input file!')


##################################################################################
#
# Loop over input files and set up mizuRoute simulation
for runname,inpath in input_runs.items():
		print(runname)

		if not os.path.exists(inpath):
			print('Error, input file doesnt exist',runname,inpath)
			continue

		outdir = os.path.join(outbase,runname)+'/'
		if not os.path.exists(outdir):
			os.mkdir(outdir)
		indir = os.path.dirname(inpath)+'/'
		fin   = os.path.basename(inpath)
		outname = 'q_'+runname

		control_file = os.path.join(controldir,'control_'+runname+'.txt')

		# Check if the calibration has already been run scucessfully
		outfiles = glob.glob(os.path.join(outdir,outname+'*.nc'))
		if len(outfiles) >= runlength and not override:
			print('Run',runname,'already completed, skipping')
			continue

		# copy template and modify
		if override and os.path.exists(control_file):
			os.remove(control_file)
		if not os.path.exists(control_file):
			shutil.copy(control_template,control_file)
		sed_expr = 's#<OUTDIR>#'+outdir+'#g; '
		sed_expr += 's#<INDIR>#'+indir+'#g; '
		sed_expr += 's#<FIN>#'+fin+'#g; '
		sed_expr += 's#<OUTNAME>#'+outname+'#g'
		print(sed_expr)
		cmd = ['sed','-i','-e',sed_expr ,control_file]
		subprocess.call(cmd)

		# add fm_file to sublist
		sublist.append(control_file)

		if len(sublist)==nperjob:
			print('Submitting jobs',len(sublist))
			os.environ['CONTROL_FLIST']=':'.join(sublist)
			os.environ['LOGDIR'] = logdir
			os.environ['MIZU_EXE'] = mizuexe
			print(os.environ['CONTROL_FLIST'])
			print(' '.join(qsub_command))
			subprocess.call(qsub_command)
			# clear sublist
			sublist = []

# submit any that are left
if len(sublist)>0:
	print('Submitting jobs',len(sublist))
	os.environ['CONTROL_FLIST']=':'.join(sublist)
	os.environ['LOGDIR'] = logdir
	os.environ['MIZU_EXE'] = mizuexe
	print(os.environ['CONTROL_FLIST'])
	print(' '.join(qsub_command))
	subprocess.call(qsub_command)
else:
	print('No more simulations to submit!')
