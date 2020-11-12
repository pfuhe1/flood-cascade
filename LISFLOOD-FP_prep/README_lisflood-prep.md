### Required steps:
- First produce river network using rivernet_prep scripts
- Once mizuRoute baseline simulations have been run: from setup_scripts/lisflood-fp, run ` calc_q_returnperiod.py` and `lisflood_setup_bankfullQ.py`

### Option using Mannings equation for channel depths
- Run main script calling LFPtools routines to generate lisflood-fp model files for a particular domain e.g. `LISFLOOD-FP_prep/077_main_allres_manningdepths.py`

### Option using channel solver:
- Run main script calling LFPtools routines to generate lisflood-fp model files for a particular domain e.g. `LISFLOOD-FP_prep/077_main_allres_rectlarger_solverchannels.py`. set `chansolv_out = None` to run first time. This outputs csv file `joinf` which is the input for the channel solver.
- Run channel solver script `run_solver.m`, specifying joinf file as input. The channel solver outputs a csv file with the optimized channel parameters.
- Run the main script again, but now set the variable `chansolv_out` to the file output by the channel solver.

### Calculate discharge points
- `lisflood_discharge_inputs_qgis_clip.py`: Clips streamnet to domain used for LISFLOOD-FP simulations, then creates points shapefiles representing locations where discharge will be input, (and also downstream boundaries).
- NOTE: due to differences in clipping rasters and vector datasets, the domain 'extent' specified here needs to be adjusted slightly from the domain 'regbound' used in the main scripts above. This was done manually by checking the actual extent of one of the clipped rasters calculated above.
- Output used as 'streamnet' in `setup_scripts/lisflood-fp/lisflood_setup_inputs_obs_v2.py`
