import os,sys,subprocess,glob
from check_forcing_dates import check_dates
import datetime
import numpy as np

def merge_obs(catch,obsname,f_qobs,f_tas,f_pr,f_pet,fpattern_out,outdir,tmpdir,clim_start,clim_end):

	clim_start_str = clim_start.strftime('%Y-%m-%d')
	clim_end_str = clim_end.strftime('%Y-%m-%d')


	#if not os.path.exists(f_out):
	infile = glob.glob(fpattern_out)
	if len(infile)==0: # if outfile doesn't exist
		if os.path.exists(f_qobs) and os.path.exists(f_tas) and os.path.exists(f_pr) and os.path.exists(f_pet):

			# First check dates of f_obs:
			qstart,qend = check_dates(f_qobs)
			qstart_str = qstart.strftime('%Y-%m-%d')
			qend_str = qend.strftime('%Y-%m-%d')
			# Work out the start and end of combined time series
			if (qstart - clim_start).days > 0:
				start_all = qstart
			else:
				start_all = clim_start
			if (qend - clim_end).days < 0:
				end_all = qend
			else:
				end_all = clim_end
			start_str = start_all.strftime('%Y%m%d')
			end_str = end_all.strftime('%Y%m%d')
			f_out = os.path.join(outdir,catch+'_catchment_'+obsname+'_'+start_str+'-'+end_str+'.nc')
			# debug statements for time periods
			#print('q',qstart_str,qend_str)
			#print('clim',clim_start_str,clim_end_str)
			#print('combined',start_str,end_str)

			# Firstly merge pet file and qobs then use expr to extract q
			# This tricks cdo into adding the lat/lon coordinates to the qobs variable
			tmp1 = os.path.join(tmpdir,'merge_'+catch+'_tmp1.nc')
			# Set missing value to -9999.0 for pet (to match q)
			cmd1 = ['cdo','-L','merge','-seldate,'+clim_start_str+','+clim_end_str,f_qobs,'-seldate,'+qstart_str+','+qend_str,'-setmissval,-9999.',f_pet,tmp1]
			print(' '.join(cmd1))
			subprocess.call(cmd1)
			tmp2 = os.path.join(tmpdir,'qobs_'+catch+'_tmp2.nc')
			cmd2 = ['cdo','expr,q_obs3d=pet*0+q_obs',tmp1,tmp2]
			print(' '.join(cmd2))
			subprocess.call(cmd2)
			os.remove(tmp1)
			# re-add attributes
			cmd3 = ['ncatted','-a','longname,q_obs3d,c,c,"Mean observed daily discharge"','-a','units,q_obs3d,c,c,mm/day',tmp2]
			subprocess.call(cmd3)
			# Now merge the 4 variables we want:
			cmd4 = ['cdo','-L','merge','-seldate,'+qstart_str+','+qend_str,'-selvar,pr',f_pr,'-seldate,'+qstart_str+','+qend_str,f_tas,'-seldate,'+qstart_str+','+qend_str,f_pet,tmp2,f_out]
			print(' '.join(cmd4))
			subprocess.call(cmd4)
			os.remove(tmp2)
		else:
			print('ERROR files dont all exist:',f_qobs,f_tas,f_pr,f_pet,'\n')
			return -1
		return f_out
	else:
#			print('exists,skipping')
		print('exists, renaming lat/lon')
		cmd = ['ncrename','-v','lat,latitude','-d','lat,latitude','-v','lon,longitude','-d','lon,longitude', infile[0] ]
		print(cmd)
		subprocess.call(cmd)
		return infile[0]
