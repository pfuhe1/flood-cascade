## flood-cascade
### A model cascade to translate meteorological drivers to flood hazard.
This modelling chain consists of distinct models to be run in sequence. The code is split up into a number of components relating to preparing or running a specific model. Each of these components has its own README file.
***
### Components:
**[metforcings_prep](metforcings_prep/README_metforcings-prep.md):**
Scripts to prepare meteorological forcings to run through the MetSIM and FUSE software.

**[rivernet_prep](rivernet_prep/README_rivernet-prep.md):**
Example configuration file for creating a river network in LFP-Tools.
Scripts to downsample a river network to lower resolution (e.g. 3arc-seconds to 9arc-seconds)

**[FUSE-metsim_prep](FUSE-metsim_prep/README_fuse-metsim_prep.md):**
Scripts to generate domain ancillary files for the MetSIM and FUSE software.

**[FUSE_regionalization](FUSE_regionalization/README_fuse-regionalization.md):**
Scripts to process datasets for parameter regionalization
Code to determine similarity of characteristics between model grid-cells and GRDC catchments

**[MizuRoute_prep](MizuRoute_prep/Readme_mizuroute-prep.md):**
Scripts to generate mizuroute ancillary files
Creates river network file based on output of rivernet_prep code
Creates grid mapping file based on output of rivernet_prep and FUSE model grid

**[LISFLOOD-FP_prep](LISFLOOD-FP_prep/README_lisflood-prep.md):**
Example scripts to generate a LISFLOOD-FP model using LFP-Tools
Script to determine location of lisflood discharge inputs from river network
TODO: Add code for calculating bankfull discharge

**[submit_scripts](submit_scripts/README_submit-scripts.md)**
Code to be run on a cluster environment to prepare and submit simulations for each of the model components:
1. fuse_catchments
2. fuse_gridded
3. mizuRoute
4. lisflood-fp

***
### Submodels

Code from other projects used in this model chain, included using git submodules. Tagged at versions tested and working with this code.
 - LFP-tools
 - fuse
 - mizuRoute

These submodules can be downloaded by:
1. Cloning this repository using this option: `git clone --recurse-submodules https://github.com/pfuhe1/flood-cascade.git` OR
2. Running `git submodule update --init` after cloning

**LISFLOOD-FP**:

In addition to the above submodels, the LISFLOOD-FP and channel solver code are licensed for non-commercial use only. The model code used in this model are available through zenodo [DOI:10.5281/zenodo.4268656](https://doi.org/10.5281/zenodo.4268656).

***
### Dependencies

- python
- [QGIS](https://qgis.org): used for some of the processing of spatial data
- [cdo](https://code.mpimet.mpg.de/projects/cdo/): Climate data operators (installable through anaconda).
- [nco](http://nco.sourceforge.net): Netcdf operators (installable through anaconda).
- [MetSim](https://github.com/UW-Hydro/MetSim): Python code to calculate daily potential evapotranspiration (FUSE input).
- [PyETo](https://pyeto.readthedocs.io/en/latest/): Python code to calculate monthly potential evapotranpiration (used for parameter regionalization datasets)


Other python libraries (e.g. installable through anaconda):
 - numpy
 - gdal
 - netCDF4
 - matplotlib and cartopy (for plotting routines only)

Dependencies for LFPTools:
- scipy
- [Pandas](https://pandas.pydata.org/)
- [Geopandas](http://geopandas.org/)
- [xarray](http://xarray.pydata.org/en/stable/)
- [Cython](https://cython.org/)
- [TauDEM](http://hydrology.usu.edu/taudem/taudem5/index.html)
- [gdalutils](https://github.com/jsosa/gdalutils.git)
- [statsmodels](https://www.statsmodels.org) (for fixelevs script only, not used when running channel solver)
