### River network generation.

1. Use LFPTools (and TauDEM) to generate river network.
 - 1a. `lfp-prepdata -i GBM_splitd8.cfg`
 - 1b. `lfp-split -i GBM_splitd8.cfg`
 - EXAMPLE configuration file to define Ganges-Brahmaputra-Meghna (GBM) stream network using MERIT-hydro as input: `GBM_splitd8.cfg`. NOTE basin 077 was manually determined to correspond to the GBM basin by checking values in 'basins3.tif'.
 - See installation notes in LFPtools repository

2. (Optional): Downsample river network to lower horizontal resolution for LISFLOOD-FP simualtions. (scripts originally from https://github.com/pfuhe1/downsample_hydro).
  - 2a. `python downsample_hydro.py --nwindow 3 --count_thresh 2 --max_start_acc 510 --datadir <DATADIR>`
  - 2b. `python call_streamnet_downsample.py --res '9s' --datadir <DATADIR>`
  - DATADIR is the location of the river basin extracted by the 'lfp-split' script e.g. outdir in `GBM_splitd8.cfg`


 3. `lisflood_discharge_inputs_qgis.py`:

 Run in qgis python console to create shapefiles representing upstream accumulation for vertices in the river network. Used by `calc_q_returnperiod.py` and `lisflood_setup_bankfullQ.py` in `submit_scripts/lisflood-fp` folder.
