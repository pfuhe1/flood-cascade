# Readme file for creating elev_bands and domain files for FUSE and metsim:

FUSE and MetSim require grid/catchment specification files containing elevation information.
They are calculated for gridded domains (at 0.5 and 0.1 degree resolution), and catchment domains (used for parameter regionalization).

### Inputs:
- MERIT DEM elevations,
- netcdf mask file specifying domain for gridded model (masked regions are not included in simulations).
- Shapefiles for extent of catchments used for parameter regionalization
- Metsim domain files for catchments additionally require pickle files (mapping catchments to gridcells) calculated in `metforcings_prep/catchment-grid-mapping/catchment_forcings_part1_p1deg.py`

### Gridded files:

Scripts produce both FUSE 'elev_bands' file and MetSIM 'domain' file.
 - 0.5 degrees: `FUSE-elevbands_MeritDEM_p5deg.py`
 - 0.1 degrees: `FUSE-elevbands_MeritDEM_p1deg.py`

### Catchment files

1. FUSE elev_bands files:

  1. First clip MERIT DEM files to the extent of each catchment
	`python parameter_transfer_catchments_MERITDEM-v1-1.py`
  2. Then create elev_bands files for each catchment:
	`python MeritDEM_catchment_elevs.py`

2. MetSIM domain files:

  1. Create global (60S-85N) elevation file at specified resolution:
	- 0.5 degrees: `MeritDEM_regrid_elevs_v2-1_global-metsim_v2.py`
	- 0.1 degrees: `MeritDEM_regrid_elevs_v2-1_global-metsim_p1deg.py`

 2. Create domain files for each catchment:
	- 0.5 degrees: `metsim_inputs-domain_separatecatchments_v2.py`
	- 0.1 degrees: `metsim_inputs-domain_separatecatchments_p1deg.py`
