# Peter Uhe 2020/01/27
# Extract upstream and downstream points in river network for lisflood inputs.
# THis script needs to be run from the qgis console otherwise doesn't seem to work...
# Run by using this command e.g:
# exec(open('/home/pu17449/src/flood-cascade/rivernet_prep/lisflood_discharge_inputs_qgis.py').read())
#
# The output of this script is used along with the mizuRoute output to calculate bankfull Q at each point along river network.

import os,socket
import processing # (qgis processing)

####################################################
# BASE PATH:
datadir = '/home/pu17449/data2/lfp-tools/splitd8_v2'

####################################################
# INPUT FILES:
# Files for different resolutions:
#resname = '3sd8'
#streamnetpath = os.path.join(datadir,'077/strn_network_3s_d8/strn_network_3sd8_077.shp')
#accpath = os.path.join(datadir,'077/077_acc.tif')

res = '9s'
resname = res+'d8'
streamnetdir = os.path.join(datadir,'077/strn_network_'+resname+'.out')
streamnetpath = os.path.join(streamnetdir,'strn_network_'+resname+'.shp')
accpath = os.path.join(datadir,'077/acc_downsample_'+res+'.tif')

overwrite = False

####################################################
# Main Script
###########################################################################

###########################################################################
# extract specific vertices: upstream,downstream, next-to-downstream
# Use algorithm: processing.algorithmHelp('qgis:extractspecificvertices')
upall = os.path.join(streamnetdir,'strn_network_'+resname+'_upstream.shp')
if not os.path.exists(upall) or overwrite:
	print('Extracting vertices: upstream all')
	cmd_dict = {'INPUT':streamnetpath,'VERTICES':'-1','OUTPUT':upall}
	processing.run('qgis:extractspecificvertices',cmd_dict)

nextdownfile = os.path.join(streamnetdir,'strn_network_'+resname+'_next-to-downstream.shp')
if not os.path.exists(nextdownfile) or overwrite:
	print('Extracting vertices: next-to-downstream')
	cmd_dict = {'INPUT':streamnetpath,'VERTICES':'1','OUTPUT':nextdownfile}
	processing.run('qgis:extractspecificvertices',cmd_dict)

###########################################################################
# next sample accumulation at these points
# algorithm: processing.algorithmHelp('qgis:rastersampling')
acc_upall_file = os.path.join(streamnetdir,'strn_network_'+resname+'_acc_upstream.shp')
if not os.path.exists(acc_upall_file) or overwrite:
	print('Sampling accumulation values: upstream all')
	cmd_dict = {'INPUT':upall,'RASTERCOPY':accpath,'COLUMN_PREFIX':'ACC','OUTPUT':acc_upall_file}
	processing.run('qgis:rastersampling',cmd_dict)

acc_ntd_file = os.path.join(streamnetdir,'strn_network_'+resname+'_acc_next-to-downstream.shp')
if not os.path.exists(acc_ntd_file) or overwrite:
	print('Sampling accumulation values: next-to-downstream')
	cmd_dict = {'INPUT':nextdownfile,'RASTERCOPY':accpath,'COLUMN_PREFIX':'ACC','OUTPUT':acc_ntd_file}
	processing.run('qgis:rastersampling',cmd_dict)

acc_ds_file = os.path.join(streamnetdir,clipname+'_acc_downstream.shp')
if not os.path.exists(acc_ds_file) or overwrite:
	print('Sampling accumulation values: downstream')
	cmd_dict = {'INPUT':downfile,'RASTERCOPY':accpath,'COLUMN_PREFIX':'ACC','OUTPUT':acc_ds_file}
	processing.run('qgis:rastersampling',cmd_dict)

###########################################################################
# Create points shapefile from streamnet
# algorithm: processing.algorithmHelp('native:extractvertices')
vertices_file = os.path.join(streamnetdir,'strn_network_'+resname+'_vertices.shp')
if not os.path.exists(vertices_file) or overwrite:
	print('Creating points (vertices) file')
	cmd_dict = {'INPUT':streamnetpath,'OUTPUT':vertices_file}
	processing.run('native:extractvertices',cmd_dict)
