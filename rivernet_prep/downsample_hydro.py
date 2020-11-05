#!/usr/bin/env python

# inst: university of bristol
# auth: Peter Uhe
# mail: peter.uhe@bristol.ac.uk

import os,argparse
import gdalutils # From https://github.com/jsosa/gdalutils
import numpy as np
from skimage.measure import block_reduce

################################################################################################
# Function to trace down from the upstream of the river network
# To the bottom of the network, or to part of the network that has already been processed
#
def remove_trib(dirn,j,i):
	"""
	Remove upstream section of river until getting to a point which no further upstream points
	point j,i is the first point to remove, and needs to have its upstream link removed
	otherwise this function will not remove any points
	"""
	dirindex = np.array([[0,0],[0,1],[-1,1],[-1,0],[-1,-1],[0,-1],[1,-1],[1,0],[1,1]],dtype=np.int16)
	reverse_dirn = np.array([-1,5,6,7,8,1,2,3,4],dtype=np.int16)
	counter = 0
	while True:
		#print('removing point j,i,dirn',j,i,dirn)
			# Remove this cell
		dirn[j,i] = 0
		counter += 1
		# Check if this cell has any upstream points (if not, exit)
		foundupstream=False
		# Loop over d8 directions
		for d in range(1,9):
			# If the neighbor points to this cell
			if dirn[j+dirindex[d,0],i+dirindex[d,1]] == reverse_dirn[d]:
				# Go to next cell
				nextj = j + dirindex[d,0]
				nexti = i + dirindex[d,1]
				foundupstream = True
		# Exit if no upstream point
		if not foundupstream:
			print('Removed trib with size',counter)
			return dirn
		else:
			j = nextj
			i = nexti

################################################################################################
# Function to trace down from the upstream of the river network
# To the bottom of the network, or to part of the network that has already been processed
#
def trace_down(j,i,acc_arr,ord_arr,process_mask,dir_arr,max_step_acc):

	# Arrays to reference the 8 possible directions
	# d8 dirs gives the steps in j,i directions
	d8_dirs = np.array([ [0,1] , [-1, 1] , [-1, 0] , [-1, -1], [0,-1], [ 1, -1], [1, 0], [1, 1] ])
	# dir_vals gives the direction value (1-8) corresponding with the directions
	# 1 -East, 2 - Northeast, 3 - North, 4 - Northwest, 5 - West, 6 - Southwest, 7 - South, 8 - Southeast.
	dir_vals = np.arange(1,9)

	# Counter for the number of cells in this tributary
	counter = 1

	# Loop over cells from upstream to downstream
	while True:
		process_mask[j,i] = True # mask this cell so not processed again
		try:
			# Work out the values for the neighboring cells
			nextij = tuple(np.transpose(d8_dirs + [j,i]))
			nextacc = acc_arr[nextij]
			nextord = ord_arr[nextij]
			nextmask = process_mask[nextij]
			# Possible items must have greater accumulation (and same or greater strahler number)
			possible_items = np.logical_and(nextord >= ord_arr[j,i] , nextacc > acc_arr[j,i] )
#			print('Debug, possible_items',possible_items,nextord[possible_items])
			if sum(possible_items)==0:
				# Try special case of the next cell having the same acc: (due to bug in merit-hydro)
				nextacc.mask=False
				possible_items = np.logical_and(nextord>=ord_arr[j,i] , np.abs(nextacc - acc_arr[j,i]) < 1e-6)
				if sum(possible_items) >= 1:
					zz = np.where(possible_items)[0][0]
					print('Equal accumulations, continuing...',acc_arr[j,i],zz)
				else:
					dir_arr[j,i] = -1
					print('reached bottom of network. len=',counter)
					return
			elif min(nextord[possible_items]) == ord_arr[j,i] and min(nextacc[possible_items])< acc_arr[j,i]+max_step_acc:
				# Normal case (there is just a small increase in accumulation to the next cell)
				# mask nextacc by possible items then take the lowest value
				nextacc = np.ma.masked_where(np.logical_not(possible_items),acc_arr[nextij])
				zz = np.ma.argmin(nextacc,fill_value=1.e10)
#				print(nextacc,'dir',zz+1)
			else:
				nextacc.mask = False
				possible_items = np.logical_and(nextord > ord_arr[j,i] , nextacc > acc_arr[j,i] ) # Try merging with larger river stem
				if sum(possible_items)==0:
					# Continue on same stream (assume smaller tributary has merged)
					possible_items = np.logical_and(nextord==ord_arr[j,i] ,nextacc>acc_arr[j,i])
					# mask nextacc by possible items then take the lowest value
					nextacc = np.ma.masked_where(np.logical_not(possible_items),acc_arr[nextij])
					zz = np.ma.argmin(nextacc,fill_value=1.e10)
				elif sum(possible_items)==1:
					# Merge into the main tributary, only one possible cell
					zz = np.where(possible_items)[0][0]
				else:
					# Merge into the main tributary,
					# Need to determine which point in the main stem the tributary should join
					zvals = np.where(possible_items)[0]
					nextacc = acc_arr[nextij][zvals]
					# If the larger of the next acc values is greater by more than acc0, then chose the greater one
					# First sort from highest to lowest accumulation
					pairs =  sorted(zip(zvals,nextacc), key= lambda x: x[1],reverse = True)
					match = pairs[-1] # default choice is lowest accumulation
					for y in range(len(zvals)-1):
						if pairs[y][1] > pairs[y+1][1] + acc_arr[j,i]:
							match = pairs[y]
							break
					zz = match[0]

			# Set value of directions array, and index of next downstream cell
			dir_arr[j,i] = dir_vals[zz]
			nextj = d8_dirs[zz][0]+j
			nexti = d8_dirs[zz][1]+i

			if process_mask[nextj,nexti]:
				#The downstream cell has already been processed
				#if counter == 1:
				#	# Dont include single cells as tributaries
				#	dir_arr[j,i] = 0
				#	print('single cell tributary, removing')
				#if ord_arr[nextj,nexti] == ord_arr[j,i] and acc_arr[j,i] < max_start_acc:
				#if counter < 10 and acc_arr[j,i] < max_start_acc:
				# Look 2/3 ahead and check if ord is still the sample
				dirnext = dir_arr[nextj,nexti]
				next2j = d8_dirs[dirnext-1][0]+nextj
				next2i = d8_dirs[dirnext-1][1]+nexti
				dirnext2 = dir_arr[next2j,next2i]
				next3j = d8_dirs[dirnext2-1][0]+next2j
				next3i = d8_dirs[dirnext2-1][1]+next2i
				if ord_arr[next3j,next3i] == ord_arr[j,i] and acc_arr[j,i] < max_start_acc:
				# small tributary created by cutting off a loop. Remove this
					print('Removing tributary which shouldnt be in network')
					print('Debug: ORD012,ACC',ord_arr[j,i],ord_arr[nextj,nexti],ord_arr[next2j,next2i],acc_arr[j,i])
					dir_arr = remove_trib(dir_arr,j,i)
				else:
					print('reached end of trib: len=',counter)
				return
			else:
				j = nextj
				i = nexti
				counter+=1
		except Exception as e: # We've probably gone outside the boundary
			print ('error: probably at boundary',nextij)
			print(e)
			dir_arr[j,i] = -1
			return


#########################################################################################
# Function to loop over tributaries and trace out each one, creating direction maps as we go
# Inputs are ac_arr (accumulation) and ord_arr (strahler order)
#
def create_directions(acc_arr,ord_arr,max_start_acc=500,max_step_acc = 40):
	# Use process mask to keep track of river cells we have already traced.
	process_mask = ord_arr<=0
	dir_arr = np.zeros(ord_arr.shape,dtype=np.int16)
	# Loop over tributaries, start with smallest accumulations
	while True:
		# Mask out processed accumulations to work out lowest yet to be processed
		acc_arr = np.ma.masked_where(process_mask,acc_arr)
		# Work out index of smallest accumulation
		start_ij = np.unravel_index(np.ma.argmin(acc_arr,fill_value = 1.e10),acc_arr.shape)
		acc_arr.mask = False # Remove mask again...
		# Only trace tributaries with accumulation less than max_start_acc
		if acc_arr[start_ij[0],start_ij[1]] < max_start_acc:
			print('starting tracedown',start_ij[0],start_ij[1],acc_arr[start_ij[0],start_ij[1]])
			trace_down(start_ij[0],start_ij[1],acc_arr,ord_arr,process_mask,dir_arr,max_step_acc)
		else:
			return dir_arr


#########################################################################################
# Simple function to return number of values > 0 for use in block_reduce
#
def count(data,axis=None):
	return (data>0).sum(axis=axis)


#########################################################################################
# Main code
if __name__ == '__main__':

	#########################################################################################
	# Input parameters for algorithm

	parser = argparse.ArgumentParser(description='Downsample river network to lower horizontal resolution')
	parser.add_argument('-n','--nwindow',default = '3',type=int,help = 'Number of cells to accumulate into a single pixel')
	parser.add_argument('-c','--count_thresh',default = '2',type=int,help = 'Minimum river cells needed within a (nwindow x nwindow) block to include')
	parser.add_argument('-m','--max_start_acc',default = '510',type=int,help = 'Maximum accumulation that will be considered as an upstream point of a tributary.')
	parser.add_argument('-d','--datadir',help = 'Directory containing river network data (input and output)')
	args = parser.parse_args()

	nwindow = args.nwindow
	count_thresh = args.count_thresh
	max_start_acc = args.max_start_acc
	datadir = args.datadir


	# nwindow is block of ( nwindow x nwindow ) cells to aggregate over

	# count thresh is the minimum number of river cells within each window (otherwise this window is not included in the downsampled river network

	# max_start_acc is the maximum accumulation that will be considered for the upstream point in a tributary. If the most upstream point of the tributary is greater than this it wont be included.
	# This is to prevent including cells from the original network being included as separate tributaries, if they are not traced by the algorithm
	# min acc is 250km^2, so will include tributaries with lowest accumulation from 250 - max_start_acc

	# max_step_acc is the maximum accumulation between cells before we assume that a tributary has merged

	max_step_acc = 40 # Could potentially increase this to 250?


	# Recommended for 9s resolution
	# nwindow = 3
	# count_thresh = 2
	# max_start_acc = 510

	# Recommended for 15s resolution
	# nwindow = 5
	# count_thresh = 3
	# max_start_acc = 550


	# Recommended for 30s resolution
	# nwindow = 10
	# count_thresh = 5
	# max_start_acc = 580

	#########################################################################################
	# Set up paths

	# Input files for 3s river network
	#datadir = '/home/pu17449/data2/lfp-tools/splitd8_v2/077/'

	# Accumulation at each point at the river network, masked for accumulation > 250 km^2
	acctif = os.path.join(datadir,'077_acc.tif')
	# Digital elevation model to downsample (not used to calculate directions)
	demtif = os.path.join(datadir,'077_dem.tif')
	# Strahler order for each point on the river network (calculated by TauDEM)
	ordtif = os.path.join(datadir,'077_ord.tif')

	# Output name assumes that input files are 3s resolution
	dem_downsample = os.path.join(datadir,'dem_downsample_'+str(nwindow*3)+'s.tif')
	acc_downsample = os.path.join(datadir,'acc_downsample_'+str(nwindow*3)+'s.tif')
	ord_downsample = os.path.join(datadir,'ord_downsample_'+str(nwindow*3)+'s.tif')
	f_dir          = os.path.join(datadir,'dir_d8_downsample_'+str(nwindow*3)+'s.tif')
	f_net          = os.path.join(datadir,'net_downsample_'+str(nwindow*3)+'s.tif')
	f_outlets      = os.path.join(datadir,'outlets_downsample_'+str(nwindow*3)+'s.tif')

	#########################################################################################
	# Get geometry information from acctif, assume all arrays have the same geometry
	#
	geo = gdalutils.get_geo(acctif)
	print(geo)
	# modify number of cells: divide by nwindow and round up (block_reduce pads the arrays)
	geo[4] = int(np.ceil(geo[4]/nwindow))
	geo[5] = int(np.ceil(geo[5]/nwindow))
	# modify resolution: multiply by nwindow
	geo[6] = geo[6]*nwindow
	geo[7] = geo[7]*nwindow

	#########################################################################################
	# Downsample dem array
	if not os.path.exists(dem_downsample):
		data = gdalutils.get_data(demtif)
		print('inshape',data.shape)
		downsample_dem = block_reduce(data,block_size=(nwindow,nwindow),func = np.mean,cval=-9999)
		print('downsampled dem',downsample_dem.shape)
		gdalutils.write_raster(downsample_dem, dem_downsample, geo, 'Float32', -9999)

	#########################################################################################
	# Downsample ord and acc arrays (for calculation of directions)
	#
	if not os.path.exists(ord_downsample) or not os.path.exists(acc_downsample):
		data = gdalutils.get_data(ordtif)
		downsample = block_reduce(data,block_size=(nwindow,nwindow),func = np.max,cval=-32767)
		downsample_count = block_reduce(data,block_size=(nwindow,nwindow),func = count,cval=-32767)
		downsample_ord = np.where(downsample_count>=count_thresh,downsample,0)
		print('downsampled ord',downsample.shape)
		gdalutils.write_raster(downsample_ord, ord_downsample, geo, 'Int16', -32767)

		data = gdalutils.get_data(acctif)
		downsample = block_reduce(data,block_size=(nwindow,nwindow),func = np.max,cval=-9999)
		downsample_acc = np.where(downsample_count>=count_thresh,downsample,0)
		print('downsampled acc',downsample.shape)
		gdalutils.write_raster(downsample_acc, acc_downsample, geo, 'Float32', -9999)
	else:
		downsample_ord = gdalutils.get_data(ord_downsample)
		downsample_acc = gdalutils.get_data(acc_downsample)


	#########################################################################################
	# Create  directions array along stream network
	#
	if not os.path.exists(f_dir):
		dir_arr = create_directions(downsample_acc,downsample_ord,max_start_acc=max_start_acc,max_step_acc=max_step_acc)

		# Create outlets array based on directions (use special value of -1 at bottom of river or domain boundary). Setting missing value to 0 allows processing in qgis by 'raster pixels to points'
		#
		outlets_arr = dir_arr == -1
		gdalutils.write_raster(outlets_arr, f_outlets, geo, 'Int16', 0)

		# Remove outlets from final outlets_arr, write out to file
		dir_arr[outlets_arr] = 0
		gdalutils.write_raster(dir_arr, f_dir, geo, 'Int16', -32767)

		# Create network array of points in river network
		#
		net_arr = dir_arr > 0
		gdalutils.write_raster(net_arr, f_net, geo, 'Int16', -32767)

	print('Finished!')
