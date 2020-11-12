#!/usr/bin/env python

# inst: university of bristol
# auth: Peter Uhe
# mail: Peter.Uhe@bristol.ac.uk / pete.uhe@gmail.com

import os
import lfptools as lfp
from glob import glob
from shutil import copyfile
from subprocess import call
import gdalutils
import numpy as np

# Define domain
res = '9s'
regname = 'rectlarger'
buildname = 'chansolverMSWEP'
regbound = [89.081,24.,90.3,26.5]

indir = '/home/pu17449/data2/lfp-tools/splitd8_v2/077'
outdir1 = os.path.join(indir,'lfptools_'+buildname+'_'+res)
outdir2 = os.path.join(indir,'lisfloodfp_'+buildname+'_'+res+'_'+regname)

basinno = '077'
overwrite=False

####### Inputs ######

# Inputs at full res
dem_masked = indir + '/'+basinno+'_MERITDEM-masked.tif'
wthf  = '/home/pu17449/data2/grwl_widths/GBM-tiles_grwl_width.vrt'
#wthf = indir+'/'+basinno+'_wth.tif' # Merit hydro widths

# File produced by lisflood_setup_bankfullQ.py script
bfqf = indir+'/../bankfullq_'+res+'d8/strn_network_'+res+'d8_GBM-p1deg_MSWEP2-2-ERA5_bankfullq.tif'

# NOTE: this is needed for model builder but overridden for actual simulations.
discharge_csv = 'dummy_discharge.csv'

if res == '3s':
	#dem = indir+'/'+basinno+'_dem.tif'
	dem = indir+'/'+basinno+'_MERITDEM.tif' # Merit dem
	recf = indir+'/'+basinno+ '_rec.csv'
	netf = indir+'/'+basinno+'_net.tif'
	dirf = indir+'/'+basinno+'_dir.tif'
else:
	####### Inputs at low  res ######
	# Files produced by downsample_hydro script
	netf = indir+'/net_downsample_'+res+'.tif'
	dirf = indir+'/dir_d8_downsample_'+res+'.tif'
	dem = indir+'/dem_downsample_'+res+'.tif' # Downsample by simple mean

	# produced by call_streamnet_downsample script
	recf = indir+'/rec_downsample_'+res+'.csv'

# After running this once and getting channel solver output, then run again with chansolv_out set
#chansolv_out = None
chansolv_out = os.path.join(indir,'../rectlarger_channelsolver/csvfiles','chansolverMSWEP_brahm_solv_thin5_smooth500.csv')

##### Outputs #####

# Outputs from first set of scripts (without file extensions)
# shapefiles and tif files are both produced
wdt_shp    = os.path.join(outdir1,basinno+'_wdt_'+res)
slp_shp    = os.path.join(outdir1,basinno+'_slp_'+res)
bnk_shp    = os.path.join(outdir1,basinno+'_bnk_'+res)
bfq_shp    = os.path.join(outdir1,basinno+'_bfq_'+res)
dpt_shp    = os.path.join(outdir1,basinno+'_dpt_'+res)
bed_shp    = os.path.join(outdir1,basinno+'_bed_'+res)

# Outputs (buildmodel):
# '.asc. files are also produced by the buildmodel script
parf    = os.path.join(outdir2,basinno+'.par')
bcif    = os.path.join(outdir2,basinno+'.bci')
bdyf    = os.path.join(outdir2,basinno+'.bdy')
evapf   = os.path.join(outdir2,basinno+'.evap')
gaugef  = os.path.join(outdir2,basinno+'.gauge') # Turned off
stagef  = os.path.join(outdir2,basinno+'.stage') # Turned off
dembnkf = os.path.join(outdir2,basinno+'_dembnk.tif')
dembnk1D = os.path.join(outdir2,basinno+'_dembnk_1D.tif')

# CSV file with combined channel, used for running Neal channel solver
joinf = os.path.join(outdir2,'recjoin_'+buildname+'_'+res+'_'+regname+'.csv')

wdtf_2     = os.path.join(outdir1,basinno+'_wdt-solver_'+res+'_'+regname)
bnkf_2     = os.path.join(outdir1,basinno+'_bnk-solver_'+res+'_'+regname)
bedf_2     = os.path.join(outdir1,basinno+'_bed-solver_'+res+'_'+regname)

#############################################################################################
# Running LFPtool Scripts

# Creating a folder to save lfptools files
try:
    os.makedirs(outdir1)
except FileExistsError:
    pass

# Calling lfp-getbankfullq
if not os.path.exists(bfq_shp+'.shp') or overwrite:
	lfp.getbankfullq(thresh=0.002, output=bfq_shp, recf=recf, netf=netf,
		proj='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',
		fbankfullq=bfqf)

# Calling lfp-getwidths
if not os.path.exists(wdt_shp+'.shp') or overwrite:
	lfp.getwidths(output=wdt_shp,
		recf=recf,
		netf=netf,
		method = 'var_intreach',
		fbankfullq = bfq_shp+'.shp',
		proj='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',
		fwidth=wthf)

# Calling bankelevs
if not os.path.exists(bnk_shp+'.shp') or overwrite:
	lfp.getbankelevs(outlier='yes',
		method='near',
		hrnodata=-9999,
		thresh=0.05,
		output=bnk_shp,
		recf = recf,
		netf = netf,
		hrdemf=dem,
		proj='+proj=longlat + ellps=WGS84 + datum=WGS84 + no_defs')

# Calling lfp-getslopes
if not os.path.exists(slp_shp+'.shp') or overwrite:
	lfp.getslopes(output=slp_shp,
		recf=recf,
		netf=netf,
		proj='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',
		source=bnk_shp+'.shp',
		step = 50)

if not os.path.exists(dpt_shp+'.shp') or overwrite:
	lfp.getdepths(proj='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',
		netf=netf,
		method='depth_manning',
		output=dpt_shp,
		wdtf=wdt_shp+'.shp',
		slpf=slp_shp+'.shp',
		qbnkf = bfq_shp+'.shp',
		n=0.035) # Choose q to be consistent with SGCn in lisflood config

#############################################################################################
# Crop domain for region

# First clip floats:
output_tifs = [wdt_shp+'.tif',dem,slp_shp+'.tif',bfq_shp+'.tif',dpt_shp+'.tif']
for tif in output_tifs:
	# Snap extent to match input tif grid cells
	geo = gdalutils.get_geo(tif)
	xmin0,ymin0,xmax0,ymax0 = regbound
	# Geo has format [xmin, ymin, xmax, ymax, xn, yn, xres, yres, ....]
	xmin = geo[0] + np.floor((xmin0 - geo[0])/geo[6])*geo[6]
	ymin = geo[1] + np.floor((ymin0 - geo[1])/geo[7])*geo[7]
	xmax = geo[2] + np.floor((xmax0 - geo[2])/geo[6])*geo[6]
	ymax = geo[3] + np.floor((ymax0 - geo[3])/geo[7])*geo[7]
	print('Clipextent:',xmin,xmax,ymin,ymax)
	tifclip = outdir1+'/'+os.path.basename(tif[:-4])+'_'+regname+'.tif'
	if not os.path.exists(tifclip) or overwrite:
		call(["gdalwarp", '-ot', 'Float32', "-te", str(xmin), str(ymin), str(xmax), str(ymax), "-overwrite", "-co", "BIGTIFF=YES", "-co", 'COMPRESS=DEFLATE', tif,tifclip])

	# Clip shp files too!
	shp0    = tif[:-4]+'.shp'
	shpclip = outdir1+'/'+os.path.basename(tif[:-4])+'_'+regname+'.shp'
	fname0 = os.path.basename(tif)[:-4]
	if not os.path.exists(shpclip) or overwrite:
		cmd = ['ogr2ogr','-spat',str(xmin),str(ymax),str(xmax),str(ymin),'-clipsrc','spat_extent',shpclip,shp0,fname0,'-f','ESRI Shapefile']
		print(' '.join(cmd))
		call(cmd)

# Then clip ints (dirn)
output_tifs = [dirf]
for tif in output_tifs:
	# Snap extent to match input tif grid cells
	geo = gdalutils.get_geo(tif)
	xmin0,ymin0,xmax0,ymax0 = regbound
	# Geo has format [xmin, ymin, xmax, ymax, xn, yn, xres, yres, ....]
	xmin = geo[0] + np.floor((xmin0 - geo[0])/geo[6])*geo[6]
	ymin = geo[1] + np.floor((ymin0 - geo[1])/geo[7])*geo[7]
	xmax = geo[2] + np.floor((xmax0 - geo[2])/geo[6])*geo[6]
	ymax = geo[3] + np.floor((ymax0 - geo[3])/geo[7])*geo[7]
	tifclip = outdir1+'/'+os.path.basename(tif[:-4])+'_'+regname+'.tif'
	if not os.path.exists(tifclip) or overwrite:
		call(["gdalwarp", '-ot', 'Int16', "-te", str(xmin), str(ymin), str(xmax), str(ymax), "-overwrite", "-co", "BIGTIFF=YES", "-co", 'COMPRESS=DEFLATE', tif,tifclip])

# Creating a folder to save LISFLOOD-FP files
try:
    os.makedirs(outdir2)
except FileExistsError:
    pass

# join river network files into single csv for input to channel solver
if not os.path.exists(joinf) or overwrite:
	lfp.join_chan(output=joinf,recf=recf,bankelevf=bnk_shp,bankfixf=bnk_shp,widthf=wdt_shp,slopef=slp_shp,bfqf=bfq_shp,depthf=dpt_shp,extent=[xmin,xmax,ymin,ymax])

# Now try getting channel output:
if os.path.exists(chansolv_out):
	lfp.getsolverchan(chansolv_out,tifclip,'+proj=longlat + ellps=WGS84 + datum=WGS84 + no_defs',bnkf_2,bedf_2,wdtf_2)
else:
	raise Exception('Error, run channel solver using file: '+joinf+'to produce: '+chansolv_out)

#############################################################################################
# Running LFPtool Buildmodel
# Above scripts output tif files, which then need to be converted into 'asc' files

# Calling buildmodel
lfp.buildmodel(parlfp=parf,
               bcilfp=bcif,
               bdylfp=bdyf,
               evaplfp=evapf,
               gaugelfp=gaugef,
               stagelfp=stagef,
               dembnktif=dembnkf,
               dembnktif_1D=dembnk1D,
               bedtif=bedf_2+'.tif',
               wdttif=wdtf_2+'.tif',
               demtif=outdir1+'/'+os.path.basename(dem[:-4])+'_'+regname+'.tif',
               fixbnktif=bnkf_2+'.tif',
               dirtif=outdir1+'/'+os.path.basename(dirf[:-4])+'_'+regname+'.tif',
               runcsv=discharge_csv,
               reccsv=recf,
               date1='1988-01-01',
               date2='1988-12-31',
               d8dirn=True)

# Copy asc files to final output folder
tmp = glob(outdir1+'/*.asc')
print(tmp)
for i in tmp:
	os.rename(i, os.path.join(outdir2,os.path.basename(i)))
