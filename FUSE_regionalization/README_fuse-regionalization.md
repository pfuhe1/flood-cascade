## README file for FUSE parameter regionalization.

### Catchment characteristics

Eight catchment characteristics were used to determine the (dis)similarity between cells in the gridded model domain and catchments from the GRDC archive:

- pr: yearly mean precipitation from worldclim v2 (http://worldclim.org/version2)
- tas: yearly mean temperature from worldclim v2
- PET: calculated from worldclim v2 data using Hargreaves method and PyETO code (https://pyeto.readthedocs.io)
- aridity index: calculated as pr/PET (from data listed above)
- Forest cover: global forest dataset (http://earthenginepartners.appspot.com/science-2013-global-forest/download_v1.6.html)
- Snow cover: Modis MOD10CM v6 dataset (https://search.earthdata.nasa.gov/search?q=MOD10CM&ok=MOD10CM)
- Soil Clay %: Soilgrids250m dataset (https://files.isric.org/soilgrids/)
- Terrain slope: Calculated from MERIT-DEM (http://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/)

Characteristics have to be calculated for both the gridded datasets to match the FUSE gridded model domain, and for individual catchments.

Data processing scripts for these characteristics are provided as examples:
- **characteristics_gridded**: Characteristics calculated over the GBM 0.1 degree grid, output in geotiff format
- **characteristics_catchment**: Characteristics calculated over catchments. Output in python pickle files, where values are stored in a python dictionary with catchment IDs as the keys.

### Calculating similarity
Scripts calculate an array of GRDC catchment IDs  matched to each grid-cell (shape [len(catchments),ny,nx]).

For each grid point the catchments IDs are sorted (from more to less similar), up to 10 catchments. Note: for the modelling, only parameter sets for the first 3 catchments were used.

- 0.5 degree resolution grid: `parameter_transfer_p5deg-GBM_distances_v1-1.py`
- 0.1 degree resolution grid: `parameter_transfer_p1deg-GBM_distances_v1-1.py`

Output files are in python pickle format, to be used later by the script: `submit_scripts/fuse_gridded/generate_param_maps_3choices.py`.
