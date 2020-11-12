## Readme for metforcings_prep
Example scripts to prepare meteorological forcing data for MetSIM and FUSE models

### Gridded model inputs
1. `metsim_merged_MSWEP-p1deg.py`: Takes daily precipitation, minimum temperature and maximum temperature and produces input files for MetSIM. Example MetSim configuration file is given (EXAMPLE_metsim.conf)
2. `fuse_prepdata_mswep_GBM-p1deg.py`: Takes daily precipitation, mean temperature and potential evapotranspiration (calculated by MetSIM), and produces an input file for FUSE.

### Catchment to grid mapping
Requires shapefiles for catchments to use for calibration. Shapefiles provided by GRDC: https://www.bafg.de/GRDC.
Steps 1. and 2. use QGIS processing, so are run using the QGIS python console.
1. `split_catchments_p1deg.py`: Splits up catchment shapefiles into segments corresponding to each FUSE grid-box. Adds area information to output shapefiles.
2. `catchment_forcings_part1_p1deg.py`: Reads shapefiles output by step 1. Outputs python pickle file containing lists of grid points and grid areas for each catchment
3. `catchment_forcings_part2_p1deg.py`: Reads pickle file output by step 2, and netcdf file containing the global grid that is used for the meteorological forcings. Add lists of grid indices for each catchment to the pickle file output by step 2.

### Catchment based model inputs:
Calculate FUSE inputs for calibration of catchments from the GRDC archive.

Uses global gridded data matching the inputs to the gridded model. Additionally processes observed daily river discharge timeseries from GRDC.
1. `parse_grdc_v2-1.py`: Converts GRDC discharge textfiles into netcdf format with variable 'q_obs'.
2. `metsim_inputs_ERA5_mswep_separatecatchments_loopappend.py`: Loops over list of catchments and produces averages of daily pr,tas,tasmin,tasmax variables for each catchment. Requires global gridded data, list of catchment IDs to process and pickle file containing catchment indices (pickle file calculated above).
3. `metsim_merge_ERA5_mswep_separatecatchments.py`: merges pr,tasmin,tasmax variables into a single file for each catchment
4. `metsim_splitstate_ERA5_mswep_separatecatchments.py`: Splits up file from step 3 into two files required by MetSIM (first 90 days go into 'state' file, and rest of time series into 'forcing' file.)
5. `metsim_prepdata_mswep_p1deg.py`:  Runs MetSim, using the inputs from step 4. Constructs FUSE input files for each catchment, combining pet, tas and q_obs
