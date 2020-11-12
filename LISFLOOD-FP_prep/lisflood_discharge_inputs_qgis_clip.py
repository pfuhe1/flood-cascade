# Peter Uhe 2020/01/27
# Extract upstream and downstream points in river network for lisflood inputs.
# THis script needs to be run from the qgis console otherwise doesn't seem to work.
# Run by using this command e.g:
# exec(open('/home/pu17449/src/flood-cascade/LISFLOOD-FP_prep/lisflood_discharge_inputs_qgis_clip.py').read())
#
# These files then need to be copied to the location of the mizuRoute discharge files
# Then the lisflood input files are calculated by the script: e.g. lisflood_setup_inputs_obs_v2.py

import os,socket
import processing # (qgis processing)

####################################################
# BASE PATH:
datadir = '/home/pu17449/data2/lfp-tools/splitd8_v2'

####################################################
# INPUT FILES:
# Files for different resolutions:
res = '3s'
resname = '3sd8'
streamnetpath = os.path.join(datadir,'077/strn_network_3s_d8/strn_network_3sd8_077.shp')
accpath = os.path.join(datadir,'077/077_acc.tif')

#res = '9s'
#resname = res+'d8'
#streamnetpath = os.path.join(datadir,'077/strn_network_'+resname+'.out/strn_network_'+resname+'.shp')
#accpath = os.path.join(datadir,'077/acc_downsample_'+res+'.tif')


####################################################
# Parameters for this domain

# Give this domain a name
regname = 'rectlarger'
# Extent should (roughly) match the extent given in the LFPtools_scripts (e.g. 077_main_allres_rectlarger_3s_2020-07-13.py)
regbound = [89.081,24.,90.3,26.5]
# Because of rounding error, we need to slightly adjust the extent at the East and south boundaries.
# Also order of co-ordinates are switched (to [xmin,xmax,ymin,ymax])
if res == '9s':
	extent = [89.081,90.297,24.002,26.5]
elif res == '3s':
	extent = [89.081,90.3,24.0,26.5]
clipname = regname+'_'+resname

extentstr = str(extent)[1:-1]
overwrite = False

####################################################
# OUTPUT PATH
outdir = os.path.join(datadir,clipname)
if not os.path.exists(outdir):
	os.mkdir(outdir)

####################################################
# Main Script
###########################################################################
# First clip to extent
# Algorithm clipvectorbyextent (see help):
# processing.algorithmHelp('gdal:clipvectorbyextent')
clipfile = os.path.join(outdir,clipname+'_streamnet.shp')
if not os.path.exists(clipfile) or overwrite:
	print('Clipping streamnet')
	cmd_dict = {'INPUT':streamnetpath,'EXTENT':extentstr,'OUTPUT':clipfile}
	processing.run('gdal:clipvectorbyextent',cmd_dict)


#NOTE after this step I manually made a small change to the clipped streamnet
# Deleted a small section of one of the reaches which was split off
# This was because it went out of the domain then back in
# After this change, I deleted the subsequent files and reran this script.

###########################################################################
# next extract specific vertices: upstream,downstream, next-to-downstream
# Use algorithm: processing.algorithmHelp('qgis:extractspecificvertices')
upall = os.path.join(outdir,'strn_network_'+resname+'_upstream.shp')
if not os.path.exists(upall) or overwrite:
	print('Extracting vertices: upstream all')
	cmd_dict = {'INPUT':streamnetpath,'VERTICES':'-1','OUTPUT':upall}
	processing.run('qgis:extractspecificvertices',cmd_dict)


upfile = os.path.join(outdir,clipname+'_upstream.shp')
if not os.path.exists(upfile) or overwrite:
	print('Extracting vertices: upstream domain')
	cmd_dict = {'INPUT':clipfile,'VERTICES':'-1','OUTPUT':upfile}
	processing.run('qgis:extractspecificvertices',cmd_dict)

downfile = os.path.join(outdir,clipname+'_downstream.shp')
if not os.path.exists(downfile) or overwrite:
	print('Extracting vertices: downstream')
	cmd_dict = {'INPUT':clipfile,'VERTICES':'0','OUTPUT':downfile}
	processing.run('qgis:extractspecificvertices',cmd_dict)

nextdownfile = os.path.join(outdir,clipname+'_next-to-downstream.shp')
if not os.path.exists(nextdownfile) or overwrite:
	print('Extracting vertices: next-to-downstream')
	cmd_dict = {'INPUT':clipfile,'VERTICES':'1','OUTPUT':nextdownfile}
	processing.run('qgis:extractspecificvertices',cmd_dict)

###########################################################################
# next sample accumulation at these points
# algorithm: processing.algorithmHelp('qgis:rastersampling')
acc_up_file = os.path.join(outdir,clipname+'_acc_upstream.shp')
if not os.path.exists(acc_up_file) or overwrite:
	print('Sampling accumulation values: upstream domain')
	cmd_dict = {'INPUT':upfile,'RASTERCOPY':accpath,'COLUMN_PREFIX':'ACC','OUTPUT':acc_up_file}
	processing.run('qgis:rastersampling',cmd_dict)

acc_upall_file = os.path.join(outdir,'strn_network_'+resname+'_acc_upstream.shp')
if not os.path.exists(acc_upall_file) or overwrite:
	print('Sampling accumulation values: upstream all')
	cmd_dict = {'INPUT':upall,'RASTERCOPY':accpath,'COLUMN_PREFIX':'ACC','OUTPUT':acc_upall_file}
	processing.run('qgis:rastersampling',cmd_dict)

acc_ntd_file = os.path.join(outdir,clipname+'_acc_next-to-downstream.shp')
if not os.path.exists(acc_ntd_file) or overwrite:
	print('Sampling accumulation values: next-to-downstream')
	cmd_dict = {'INPUT':nextdownfile,'RASTERCOPY':accpath,'COLUMN_PREFIX':'ACC','OUTPUT':acc_ntd_file}
	processing.run('qgis:rastersampling',cmd_dict)

acc_ds_file = os.path.join(outdir,clipname+'_acc_downstream.shp')
if not os.path.exists(acc_ds_file) or overwrite:
	print('Sampling accumulation values: downstream')
	cmd_dict = {'INPUT':downfile,'RASTERCOPY':accpath,'COLUMN_PREFIX':'ACC','OUTPUT':acc_ds_file}
	processing.run('qgis:rastersampling',cmd_dict)
