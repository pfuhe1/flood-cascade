# Readme for Mizuroute river network preparation:

########################################################################################################
## Mizuroute river network ancillary file

Requires outputs from TauDEM: stren_w3d8.tif, stren_net3d8.shp and catchment mask (basins3.tif).

**Catchment Mask (processing in QGIS):**
1. Polygonise tif:
`	gdal_polygonize.py <datadir>/basins3.tif <datadir>/basins3_poly.gpkg" -8 -b 1 -f "GPKG" None DN`
2. Select feature for desired catchment, 'save selected feature as'
3. Fix geometries


**Catchments:**
1. Clip stren_w3d8.tif by catchment mask (calculated above)
2. GDAL Polygonize: 'Name of the field to create='LINENO', use 8-connectedness

**Streams:**
1. Select 'extract specific vertices' (vertex 0 for upstream points)
2. Clip result by catchment mask (calculated able).

**Calculate mizuRoute ancillary file:**
`python mizuRoute_rivernet.py -c <catchments> -s <streams> -o <output_file>`


########################################################################################################
## Runoff mapping data calculation

**Reqires:**
1. netcdf file with grid of latitude and longitude points (for grid cell centres) used by runoff model
2. Shapefile representing grid boundaries matching above grid (greate in qgis by 'create grid' tool)
3. Polygonised catchments file (created above)
4. mizuroute river network file (created above)

**Preprocessing steps:**
1. Add geometry attributes (ellipsoidal method) to catchments file (adds area attribute)
2. Split by lines: file from step 1 split using shapefile grid.
3. Add geometry attributes (ellipsoidal method), to file from step 2 (adds area_2 attribute)
4. Calculate centroids from file in step 3. 

**Generate grid mapping file:**
`python mizuRoute_grid_mapping.py`

########################################################################################################

# Running Mizuroute:
A template control file 'GBM-MERIT_p1deg_v2.control_template is included with specifications matching the input files produced above.
Setup scripts to submit simulations to a cluster are included separately.

Code tested from mizuRoute develop branch:
commit 8992396f52dabed55699b03fb8b52c103e76da34
Merge: a72ef46 7b9a9e7
Author: Naoki Mizukami <mizukami@ucar.edu>
Date:   Wed May 8 13:04:29 2019 -0600
