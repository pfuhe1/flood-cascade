## flood-cascade
### A model cascade to translate meteorological drivers to flood hazard.
This modelling chain consists of distinct models to be run in sequence. The code is split up into a number of components relating to preparing or running a specific model. Each of these components has its own README file.

### Components:
**metforcings_prep:**
Scripts to prepare meteorological forcings to run through the MetSIM and FUSE software.

**rivernet_prep:**
Example configuration file for creating a river network in LFP-Tools.
Scripts to downsample a river network to lower resolution (e.g. 3arc-seconds to 9arc-seconds)

**FUSE-metsim_prep:**
Scripts to generate domain ancillary files for the MetSIM and FUSE software.

**FUSE_regionalization:**
Scripts to process datasets for parameter regionalization
Code to determine similarity of characteristics between model grid-cells and GRDC catchments

**MizuRoute_prep:**
Scripts to generate mizuroute ancillary files
Creates river network file based on output of rivernet_prep code
Creates grid mapping file based on output of rivernet_prep and FUSE model grid

**LISFLOOD_prep:**
Example scripts to generate a LISFLOOD-FP model using LFP-Tools
Script to determine location of lisflood discharge inputs from river network
TODO: Add code for calculating bankfull discharge

**submit_scripts**
Code to be run on a cluster environment to prepare and submit simulations for each of the model components:
1. fuse_catchments
2. fuse_gridded
3. mizuRoute
4. lisflood-fp

**submodels:**
Code from other projects used in this model chain
TODO:
 - setup/submit scripts (for running on cluster) - OR move up?
 - LFP-tools
 - FUSE
 - MizuRoute
