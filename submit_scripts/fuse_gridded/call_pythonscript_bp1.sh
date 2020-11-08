#!/bin/sh
#PBS -l walltime=11:55:00
#PBS -j oe

# NOTE: some resources are passed through by the command line
# e.g. (-l select=1,ncpus=X,ompthreads=Y:mem=12gb)
module load lib/netcdf/4.70
module load lang/python/anaconda/3.7-2019.03
module list
python /home/pu17449/src/flood-cascade/setup_scripts/fuse_gridded/fusewrapper_gridded.py
