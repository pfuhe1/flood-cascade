## README file for submit scripts.

Code to be run on a cluster environment to prepare and submit simulations for each of the model components:

### fuse_catchments
 - `setup_grdc_catch_mswep-p5deg_best73_checkConverged.py`: Example submit script to calibrate catchments based on 0.5 degree grid regionalization.
 - `setup_grdc_catch_mswep-p1deg_best222_checkConverged.py`: Example submit script to calibrate catchments based on 0.1 degree grid regionalization

 Other scripts to note:
 - `fusewrapper_grdc_catch.py`: python wrapper script, using a multiprocessing queue to loop over multiple simulations submitted as a single batch job.
 - `call_pythonscript_bp1.sh`: bash script specifying resources to a batch job scheduler and set up environment before calling python wrapper script. Different scripts will be used for different computing environments.

### fuse_gridded
1. First calculate parameter set maps for the fuse gridded model, based on the calibrated catchments and the mappings determined by the fuse_regionalization scripts. Example for 0.1 degree resolution using MSWEP2.2 calibration simulations.
 - `generate_param_maps_3choicesmswepp1deg_v2.py`: Produce three parameter set maps corresponding to the closest donor parameter sets for each grid-cell.
 - `generate_param_maps_3random_mswepp1deg.py`: Takes the three parameter set maps calculated above and produces new maps by randomly sampling different parameters sets at different grid-cells.
2. Then submit simulations.
 - `setup_GBM_runs.py`: Example script to setup and submit gridded simulations to cluster. Loops over FUSE structures and parameter-sets to run batches of simulations. `fuse_templates/fm_template_genforcing.txt` is used as a template file for the above script

 Other scripts to note:
 - `fusewrapper_gridded.py`:  python wrapper script, using a multiprocessing queue to loop over multiple simulations submitted as a single batch job.
 - `call_pythonscript_bp1.sh`: bash script specifying resources to a batch job scheduler and set up environment before calling python wrapper script. Different scripts will be used for different computing environments.

### mizuRoute
`run_mizuRoute_templated_mswep010.py`: Example script to setup and submit mizuRoute simulations to a cluster. Loops over FUSE structures and parameter-sets to run batches of simulations.

Other scripts to note:
- `GBM-MERIT_p1deg_v2.control_template`: A template control file used in the above script to generate control files for individual simulations.
- `mizuRoute_wrapper.py`:  python wrapper script, using a multiprocessing queue to loop over multiple simulations submitted as a single batch job.
- `call_pythonscript_bp1.sh`: bash script specifying resources to a batch job scheduler and set up environment before calling python wrapper script. Different scripts will be used for different computing environments.

### lisflood-fp
**Determine bankfullQ from mizuRoute simulations**

These are used in the `LISFLOOD-FP_prep` scripts to calibrate the 1D channel parameters in the LISFLOOD-FP model to be consistent with the discharge simulated by mizuRoute.
- `calc_q_returnperiod.py`: Calculate a particular percentile flow of yearly maximum discharge correesponding to bankfull (e.g. 50 percentile). Reads in mizuRoute simulations, outputs netcdf file matching mizuRoute output with variable name 'IRFroutedRunoff'
- `lisflood_setup_bankfullQ.py`: Calculates bankfull discharge at every point along the lisflood river network channels. Takes output file calculated above, and shapefiles representing the river network (calculated by `rivernet_prep/lisflood_discharge_inputs_qgis.py`)
The output files from these scripts are then used to finalize the LISFLOOD-FP model before running the simulations.

**Set-up simulations based on mizuRoute discharge**

 `lisflood_setup_inputs_obs_v2.py`: Sets up and submits LISFLOOD-FP simulations to a cluster. Requires inputs:
 - Discharge from mizuRoute simulations
 - mizuRoute ancillary files
 - Model set-up produced in `LISFLOOD-FP-prep` folder

Other scripts to note:
 - `lisflood-fp_wrapper.py`: python wrapper to submit lisflood-fp simulations and post-process output into geotiff format.
 - `convert_output_totiff.py`: python functions for post-processing output. NOTE: function 'convert_to_tif_v4' converts timeseries into a single file with each time-step as a single BAND. Previous versions produce one file per time-step.
 - `077_template.par`: template parameter file used by `lisflood_setup_inputs_obs_v2.py` to generate simulations with different inputs.


### Additional configuration
There are also additional template files and settings for FUSE in the folder:
 - fuse_templates
