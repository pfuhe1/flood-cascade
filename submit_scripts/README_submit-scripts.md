### README file for submit scripts.

Code to be run on a cluster environment to prepare and submit simulations for each of the model components:

### fuse_catchments
 - `setup_grdc_catch_mswep-p5deg_best73_checkConverged.py`: Example submit script to calibrate catchments based on 0.5 degree grid regionalization.
 - `setup_grdc_catch_mswep-p1deg_best222_checkConverged.py`: Example submit script to calibrate catchments based on 0.1 degree grid regionalization
 - `qsub_grdc_catch.py`: python wrapper script, using a multiprocessing queue to loop over multiple simulations submitted as a single batch job.
 - `call_pythonscript_bp1.sh`: bash script specifying resources to a batch job scheduler and set up environment before calling python wrapper script. Different scripts will be used for different computing environments.

### fuse_param_maps
Calculate parameter set maps for the fuse gridded model, based on the calibrated catchments and the mappings determined by the fuse_regionalization scripts

### fuse_gridded


### mizuRoute


### lisflood-fp

## Additional configuration
There are also additional template files and settings for FUSE in the folder:
 - fuse_templates
