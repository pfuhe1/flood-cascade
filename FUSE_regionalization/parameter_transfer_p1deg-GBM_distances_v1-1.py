import gdal
import numpy as np
import os,glob,shutil,sys,pickle
import subprocess
import matplotlib.pyplot as plt
import matplotlib.cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from netCDF4 import Dataset
import copy


# Peter Uhe Nov 2019
# Run on pc using conda environment: gdal_env2

################################################################################
# Set up stuff
varnames = ['aridity','clyppt','pr','tas','pet','snowcover','forestcover','slope']

# Input paths
gridpath = '/home/pu17449/data2/parameter_transfer_catchments/p1deg-GBM'
catchpath = '/home/pu17449/data2/parameter_transfer_catchments/catchment_vals'
maskfile = '/home/pu17449/data2/MSWEP2-2_010deg/GBM-Grid_p1deg_v2.nc'

# Output
fpickle = '/home/pu17449/data2/parameter_transfer_catchments/GBM-p1deg_distances_GBM-hisnowIQR.pkl'
fpickle_reduced = '/home/pu17449/data2/parameter_transfer_catchments/GBM-p1deg_distances_GBM-reduced.pkl'
ftxt_reduced3 = '/home/pu17449/data2/parameter_transfer_catchments/GBM-p1deg_donorsbest3.txt'
override = False

# Inter quartile ranges for each characteristic
iqr = {}
# Beck 2016 IQR values:
#iqr['aridity']=0.88 # ratio pet to pr
#iqr['clyppt'] = 13.77 # %
#iqr['pr'] = 743. # mm/year
#iqr['tas'] = 26.49 # deg C
#iqr['pet'] = 1054. # mm/year
#iqr['snowcover'] = 57. # %
#iqr['forestcover'] = 45. # %
#iqr['slope'] = 1.08 # degrees (from Beck)


# IQR values calculated over the GBM region:
iqr['aridity']=1.3 # ratio pet to pr
iqr['clyppt'] = 16 # %
iqr['pr'] = 943. # mm/year
iqr['tas'] = 17.6 # deg C
iqr['pet'] = 882. # mm/year
iqr['snowcover'] = 4. # %
iqr['forestcover'] = 20. # %
iqr['slope'] = 12 # degrees

# Try larger iqr value for snowcover:
iqr['snowcover'] = 57.

#########################################################
# Function to do the processing for a single tile/catchment pair
def process_gridpoint(ingrids,i,j,incatchs,iqr):
	print(i,j)
	distances = []
	for grdcid in catchments:
		vardists = np.zeros([len(varnames)])
		for k,var in enumerate(varnames):
			vardists[k] = np.abs(ingrids[var][j,i]-incatchs[var][grdcid])/iqr[var]
		dist = vardists.mean()
		distances.append((grdcid,dist,vardists))

	return sorted(distances,key = lambda x : x[1])

if not os.path.exists(fpickle) or override:

	#########################################################
	# Load data

	# Gridded files with data for each characteristic
	ingrids = {}
	for var in varnames:
		fpath = os.path.join(gridpath,var+'_p1deg-GBM.tif')
		f = gdal.Open(fpath)
		rasterband = f.GetRasterBand(1)
		ingrids[var] = rasterband.ReadAsArray()

	ny,nx = ingrids['aridity'].shape

	# Catchment files are pickle files with dictionary indexed by GRDC catchment ID
	incatchs = {}
	for var in varnames:
		fpath = os.path.join(catchpath,var+'.pkl')
		with open(fpath,'rb') as fdata:
			incatchs[var] = pickle.load(fdata,encoding='latin1')
			print('var:',var)
			print('ncatchments',len(incatchs[var].keys()))
			print('dtype',type(next(iter(incatchs[var].keys()))))

	catchments = incatchs['aridity'].keys()

	#########################################################
	# Calculate distances

	bestdist = np.zeros([ny,nx])
	catcharray = np.zeros([len(catchments),ny,nx],dtype=int)
	distarray = np.zeros([len(catchments),ny,nx])
	vdistarray = np.zeros([len(catchments),len(varnames),ny,nx])
	# Loop over files and processs
	for i in range(nx):
		for j in range(ny):
			distances = process_gridpoint(ingrids,i,j,incatchs,iqr)
			catchlist,distlist,vardists = zip(*distances)
			catcharray[:,j,i] = catchlist
			distarray[:,j,i] = distlist
			vdistarray[:,:,j,i] = np.array(vardists)
	#		bestdist[j,i] = distances[0[-1]

	#########################################################
	# Write out results to pickle file

	with open(fpickle,'wb') as f:
			pickle.dump(catcharray,f,-1)
			pickle.dump(distarray,f,-1)
			pickle.dump(vdistarray,f,-1)
else:
	#########################################################
	# Load existing  results from pickle file

	with open(fpickle,'rb') as f:
		catcharray = pickle.load(f)
		distarray = pickle.load(f)
		vdistarray = pickle.load(f)
	nc,ny,nx = catcharray.shape

#########################################################
# Plot stuff

with Dataset(maskfile,'r') as f:
	#mask = f.variables['Band1'][:]==-9999.
	mask = f.variables['Band1'][:].mask
mask10 = np.repeat(mask[np.newaxis,:],10,axis=0)

# Print number of non-masked cells
print('active cells',np.logical_not(mask).sum())

#lats = [31.25, 30.75, 30.25, 29.75, 29.25, 28.75, 28.25, 27.75, 27.25, 26.75,
#    26.25, 25.75, 25.25, 24.75, 24.25, 23.75, 23.25, 22.75, 22.25]
#lons = [73.25, 73.75, 74.25, 74.75, 75.25, 75.75, 76.25, 76.75, 77.25, 77.75,
#    78.25, 78.75, 79.25, 79.75, 80.25, 80.75, 81.25, 81.75, 82.25, 82.75,
#    83.25, 83.75, 84.25, 84.75, 85.25, 85.75, 86.25, 86.75, 87.25, 87.75,
#    88.25, 88.75, 89.25, 89.75, 90.25, 90.75, 91.25, 91.75, 92.25, 92.75,
#    93.25, 93.75, 94.25, 94.75, 95.25, 95.75, 96.25, 96.75, 97.25, 97.75]

# p1deg grid:
lats = np.arange(31.45,22,-.1)
lons = np.arange(73.05,98,.1)

fig=plt.figure()
ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
ax.set_extent([73.25, 97.75, 22.25, 31.25], crs=ccrs.PlateCarree())
#ax.pcolormesh(lons,lats,distarray[0,::-1,:],vmin=0,vmax=1,cmap='jet')
marr = np.ma.masked_where(mask,distarray[:10,:,:].mean(0))
cm = ax.pcolormesh(lons,lats,marr,cmap='jet',vmin=0,vmax=1,transform=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cfeature.BORDERS)
plt.colorbar(cm,ax=ax,shrink=0.3)
plt.title('Dissimilarity between best 10 donor catchments and grid cells')
#plt.show()
plt.savefig('figs/p1deg_dissimilarity_map.png')

for i,var in enumerate(varnames):
	fig=plt.figure()
	ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
	ax.set_extent([73.25, 97.75, 22.25, 31.25], crs=ccrs.PlateCarree())
	#ax.pcolormesh(lons,lats,distarray[0,::-1,:],vmin=0,vmax=1,cmap='jet')
	marr = np.ma.masked_where(mask,vdistarray[:10,i,:,:].mean(0))
	cm = ax.pcolormesh(lons,lats,marr,cmap='jet',vmin=0,vmax=2,transform=ccrs.PlateCarree())
	ax.coastlines()
	ax.add_feature(cfeature.BORDERS)
	plt.colorbar(cm,ax=ax,shrink=0.3)
	plt.title('Dissimilarity: '+var+', best 10 GRDC catchments')
	#plt.show()
	plt.savefig('figs/p1deg_dissimilarity_map_'+var+'.png')

fig=plt.figure()
ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
ax.set_extent([73.25, 97.75, 22.25, 31.25], crs=ccrs.PlateCarree())
#ax.pcolormesh(lons,lats,distarray[0,::-1,:],vmin=0,vmax=1,cmap='jet')
marr = np.ma.masked_where(mask,catcharray[0,:,:])
cm = ax.pcolormesh(lons,lats,marr,vmin=catcharray.min(),vmax = catcharray.max(),cmap='jet',transform=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cfeature.BORDERS)
plt.colorbar(cm,ax=ax,shrink=0.3)
plt.title('Best donor catchments (GRDC ID)')
#plt.show()
plt.savefig('figs/p1deg_bestdonor_map.png')

#########################################################
# Evaluate number of uses by donor catchments
donor_uses = {}
bestlot = np.ma.masked_where(mask10[:10],catcharray[:10,:])
#for use in catcharray[:10,:].flatten():
for use in bestlot.compressed():
	if use in donor_uses:
		donor_uses[use]+=1
	else:
		donor_uses[use]=1

print('catchments used',len(donor_uses))
plt.figure()
#plt.hist(donor_uses.values(),bins=len(donor_uses),density=True)
plt.hist(donor_uses.values(),bins=25)
plt.title('Best 10 donor catchments for each gridcell in GBM region (total: '+str(len(donor_uses))+')')
plt.xlabel('Number of grid cells in GBM region')
plt.ylabel('Number of donor catchments')
#plt.show()
plt.savefig('figs/p1deg_donor_histogram.png')


####################################################################################
# To reduce the amount of computation,
# we might want to exlude donor catchments that only have a couple of recipients.
# The following code tests the effect of that
exclude_thresh = 50
removelist = []
for catch,uses in donor_uses.items():
	if uses<=exclude_thresh:
		removelist.append(catch)
carray_reduced = copy.deepcopy(catcharray[:10,:]) # Best 10 donor catchments
print('removed catchments with only '+str(exclude_thresh)+' recipient gridcell',len(removelist))
for catch in removelist:
	catchmatch = carray_reduced==catch
	carray_reduced[catchmatch] = -1

fig=plt.figure()
ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
ax.set_extent([73.25, 97.75, 22.25, 31.25], crs=ccrs.PlateCarree())
#ax.pcolormesh(lons,lats,distarray[0,::-1,:],vmin=0,vmax=1,cmap='jet')
cmap = matplotlib.cm.jet
cmap.set_under('k')
marr = np.ma.masked_where(mask,(carray_reduced!=-1).sum(0)) # number of non nans
cm = ax.pcolormesh(lons,lats,marr,vmin=1,vmax=10,cmap=cmap,transform=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cfeature.BORDERS)
plt.colorbar(cm,ax=ax,shrink=0.3,extend='min')
plt.title('Number of donor catchments, after removing '+str(exclude_thresh)+' uses')
#plt.show()
plt.savefig('figs/p1deg_donor_number_exclude'+str(exclude_thresh)+'.png')

####################################################################################
# Second option (more aggressive, but enforce at least 3 donors for each grid point)
exclude_thresh = 50
min_donors = 3
num_removed = 0
carray_reduced = np.ma.masked_where(mask10,copy.deepcopy(catcharray[:10,:])) # Best 10 donor catchments
catchs_used = list(donor_uses.keys())
for thresh in range(1,exclude_thresh+1):
	removelist = []
	for catch,uses in donor_uses.items():
		if uses==thresh:
			removelist.append(catch)
	print('catchments with only '+str(thresh)+' recipient gridcell',len(removelist))
	for catch in removelist:
		donormin = (carray_reduced!=-1).sum(0)==3
		catchmatch = carray_reduced == catch
		if not np.any(np.logical_and(catchmatch.sum(0)>0,donormin)): # we would reduce donors to less than 3
			carray_reduced[catchmatch] = -1
			num_removed += 1
			catchs_used.remove(catch)
	print('removed',str(num_removed))

print('total donors removed',str(num_removed))
print('left',len(catchs_used))

####################################################################################
# Write out reduced list:
with open(fpickle_reduced,'wb') as f:
	pickle.dump(carray_reduced,f,-1)


####################################################################################
# Make even more reduced list
carray_3best = np.ones([3,ny,nx],dtype=int)
donor_set = set()
for j in range(ny):
	for i in range(nx):
		best3 = np.ma.masked_values(carray_reduced[:,j,i],-1).compressed()[:3]
		for donor in best3:
			donor_set.add(donor)
		if len(best3)==3:
			carray_3best[:,j,i] = donor
		else:
			if not mask[j,i]:
				print('Error, there arent three donors for grid',j,i)

print(donor_set)
print(len(donor_set))
with open(ftxt_reduced3,'w') as f:
	for donor in donor_set:
		f.write(str(donor)+'\n')

####################################################################################
fig=plt.figure()
ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
ax.set_extent([73.25, 97.75, 22.25, 31.25], crs=ccrs.PlateCarree())
#ax.pcolormesh(lons,lats,distarray[0,::-1,:],vmin=0,vmax=1,cmap='jet')
cmap = matplotlib.cm.jet
cmap.set_under('k')
marr = np.ma.masked_where(mask,(carray_reduced!=-1).sum(0)) # number of non nans
cm = ax.pcolormesh(lons,lats,marr,vmin=1,vmax=10,cmap=cmap,transform=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cfeature.BORDERS)
plt.colorbar(cm,ax=ax,shrink=0.3,extend='min')
plt.title('Number of donor catchments, after removing '+str(exclude_thresh)+' uses')
#plt.show()
plt.savefig('figs/p1deg_donor_number_exclude'+str(exclude_thresh)+'_keep3.png')

#############################################################
# plot dissimilarity map after reducing donors

fig=plt.figure()
ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
ax.set_extent([73.25, 97.75, 22.25, 31.25], crs=ccrs.PlateCarree())
#ax.pcolormesh(lons,lats,distarray[0,::-1,:],vmin=0,vmax=1,cmap='jet')
marr = np.ma.masked_where(carray_reduced==-1,distarray[:10,:,:])
print('masked',marr.mask.sum())
cm = ax.pcolormesh(lons,lats,marr.mean(0),cmap='jet',vmin=0,vmax=1,transform=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cfeature.BORDERS)
plt.colorbar(cm,ax=ax,shrink=0.3)
plt.title('Dissimilarity between donor catchments and grid cells, after excluding some')
#plt.show()
plt.savefig('figs/p1deg_dissimilarity_map_excluded.png')

#############################################################
# plot rank of best donor after reducing donors used
ranks,ys,xs = np.where(carray_reduced!=-1)
ny,nx = carray_reduced[0,:].shape
rankarray = np.ones([ny,nx],dtype=int)*-1
for i,r in enumerate(ranks):
	yy = ys[i]
	xx = xs[i]
	if rankarray[yy,xx] == -1:
		#print(i,r,yy,xx)
		rankarray[yy,xx] = r + 1

####################################################################################
fig=plt.figure()
ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
ax.set_extent([73.25, 97.75, 22.25, 31.25], crs=ccrs.PlateCarree())
#ax.pcolormesh(lons,lats,distarray[0,::-1,:],vmin=0,vmax=1,cmap='jet')
marr = np.ma.masked_where(rankarray==-1,rankarray)
marr = np.ma.masked_where(mask,marr)
print('masked',marr.mask.sum())
cm = ax.pcolormesh(lons,lats,marr,cmap='jet',vmin=0,vmax=5,transform=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cfeature.BORDERS)
plt.colorbar(cm,ax=ax,shrink=0.3)
plt.title('Rank of best donor, after excluding some')
#plt.show()
plt.savefig('figs/p1deg_rank_map_excluded.png')
