#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -j oe

# NOTE: some resources are passed through by the command line
# (-l select=1,ncpus=X,ompthreads=Y:mem=12gb)
export
echo "Starting script"

# Load modules for blue pebble
module load lib/gdal/2.4.2 
module load lang/python/anaconda/3.7-2019.03
source activate pyenv
which python
# ENV vars used by this script (exported by submit script):
# LISFLOOD_DIR and CONTROL_SCRIPT

# LISFLOOD_DIR is the working directory (env var exported by submit script)
echo "CDing to lisflood_dir $LISFLOOD_DIR"
cd $LISFLOOD_DIR

# Call python CONTROL_SCRIPT
# ENV VARS needed (exported by submit script):  CONTROL_FLIST, EXE_FILE, LOGDIR
# ENV VARS needed (set by QSUB):                OMP_NUM_THREADS, NCPUS
echo "Calling python control script $CONTROL_SCRIPT"
python $CONTROL_SCRIPT 
echo "Done"
