import netCDF4 as nc
import glob
import geospatialtools.netcdf_tools as netcdf_tools
import numpy as np
import os
import geospatialtools.gdal_tools as gdal_tools
import pickle
import datetime
import time
import geospatialtools.upscaling_tools_fortran as upscaling_fortran
import random
import dateutil.relativedelta as relativedelta
import rasterio

def Create_Output_Files(metadata,rank,size,vars,startdate,enddate):

 #startdate = metadata['idate']
 #enddate = metadata['fdate']
 rdir = metadata['rdir']
 workspace = '%s/experiments/simulations/%s/postprocess/workspace' % (rdir,metadata['experiment'])
 output_dir = '%s/experiments/simulations/%s/postprocess/output_dir' % (rdir,metadata['experiment'])
 #nstripes = info['nstripes']
 #ncores = info['ncores']
 #output_dir = info['output_dir']
 date = startdate
 dt = 1
 #nt = (enddate - startdate).days + 1
 #nt_out = 24
 nt = 24*((enddate - startdate).days + 1)

 #Create dates array
 dates = []
 date = startdate
 while date <= enddate:
  dates.append(date)
  date = date + relativedelta.relativedelta(days=1)
 dates = np.array(dates)
 
 #Define the dimensions
 file = '%s/experiments/simulations/%s/postprocess/mapping.tif' % (rdir,metadata['experiment'])
 metadata = gdal_tools.retrieve_metadata(file)
 #nlon_chunk = int(metadata['nx']/size**0.5)/2
 #nlat_chunk = int(metadata['ny']/size**0.5)/2
 nlon_chunk = int(metadata['nx']/10)
 nlat_chunk = int(metadata['ny']/10)
 #ntime_chunk = nt_out
 dims = {'nlat':metadata['ny'],
         'nlon':metadata['nx'],
         'res':metadata['resx'],
         'minlon':metadata['minx'] + metadata['resx']/2,
         'minlat':metadata['miny'] + metadata['resy']/2,
         'undef':-9999.0,
         'chunksize':(-9999,nlat_chunk,nlon_chunk)}

 for date in dates[rank::size]:

  print("Creating the file for ",date,flush=True)
  #Update the number of days
  dtt = relativedelta.relativedelta(days=1)
  nt_out = 24#((date + dtt) - date).days
  dims['chunksize'] = nt_out

  #Define the file
  ncfile = '%s/%04d%02d%02d.nc' % (output_dir,date.year,date.month,date.day)

  #Create the netcdf file
  #nlat = md['nlat']
  md = {}
  md['file'] = ncfile
  md['vars'] = vars
  md['vars_info'] = vars
  md['tstep'] = '%dhr' % dt
  md['nt'] = nt_out
  md['nlat'] = dims['nlat']
  md['nlon'] = dims['nlon']
  md['res'] = dims['res']
  md['minlon'] = dims['minlon']
  md['minlat'] = dims['minlat']
  md['undef'] = dims['undef']
  md['chunksize'] = dims['chunksize']
  md['tinitial'] = date
  md['tinitial_all'] = startdate
  #fp = netcdf_tools.Create_NETCDF_File(dims,ncfile,vars,vars,date,'%dday' % dt,nt_out)
  fp = netcdf_tools.Create_NETCDF_File(md)

  #Iterate through all the files
  files = glob.glob('%s/workspace/%04d%02d%02d/*.pck' % (output_dir,date.year,date.month,date.day))
  for file in files:
   output = pickle.load(open(file,'rb'))
   ilats = output['coords']['ilats']
   ilons = output['coords']['ilons']
   vars_keys = output['vars'].keys()
   for var in vars_keys:
    fp.variables[var][:,ilats[0]:ilats[-1]+1,ilons[0]:ilons[-1]+1] = output['vars'][var]
   del output

  #Close the file
  fp.close()

  #Compress the file
  #print("Compressing the file for ",date,flush=True)
  #tmp = '%s/workspace/%04d%02d%02d/tmp.nc' % (output_dir,date.year,date.month,date.day)
  #os.system('nccopy -k 3 -d 4 %s %s' % (ncfile,tmp))
  #os.system('mv %s %s' % (tmp,ncfile))

  #Remove the workspace
  os.system('rm -rf %s/workspace/%04d%02d%02d' % (output_dir,date.year,date.month,date.day))

  #Update the time step
  #date = date + datetime.timedelta(hours=nt_out)
  #date = date + relativedelta.relativedelta(days=nt_out)

 #Create the control file
 #if rank == 0:
 # tstep = '1dy'
 # file_template = '^%s.nc' % ('%y4%m2',)
 # ctl_file = '%s/hydrobloks_output.ctl' % output_dir
 # netcdf_tools.Update_Control_File('nc',startdate,dims,nt,tstep,file_template,ctl_file)

 return

def Create_Virtual_Rasters(metadata):

 dir = metadata['dir']
 output_dir = metadata['output_dir']
 ncatch = metadata['ncatch']

 #Create the virtual raster of the hsus
 os.system('rm -f %s/workspace/files.txt' % output_dir)
 files = []
 for icatch in xrange(ncatch):
  #file = '%s/catch_%d/workspace/hsu_mapping_latlon.tif' % (dir,icatch)
  file = '%s/hru/%d.tif' % (dir,icatch)
  os.system('echo %s >> %s/files.txt' % (file,dir))
 os.system('gdalbuildvrt -srcnodata -9999 -input_file_list %s/files.txt %s/hsu_map.vrt' % (dir,dir))

 #Create the virtual raster of the catchment ids
 os.system('rm -f %s/workspace/files.txt' % output_dir)
 files = []
 for icatch in xrange(ncatch):
  #file = '%s/catch_%d/workspace/icatch_latlon.tif' % (dir,icatch)
  file = '%s/cid/%d.tif' % (dir,icatch)
  os.system('echo %s >> %s/files.txt' % (file,output_dir))
 os.system('gdalbuildvrt -srcnodata -9999 -input_file_list %s/files.txt %s/icatch_map.vrt' % (dir,dir))

 return

def Create_Finescale_Maps(metadata):

 dir = metadata['dir']
 bbox = metadata['coverage']['bbox']
 minlat = bbox['minlat']
 maxlat = bbox['maxlat']
 minlon = bbox['minlon']
 maxlon = bbox['maxlon']
 output_dir = metadata['output_dir']
 #hsu map
 vrt_in = '%s/hru/hsu_map.vrt' % dir
 file_out = '%s/workspace/hsu_map_region.tif' % output_dir
 os.system('rm -f %s' % file_out)
 os.system('gdalwarp -te %f %f %f %f %s %s' % (minlon,minlat,maxlon,maxlat,vrt_in,file_out))
 metadata['hru_map'] = file_out
 #cid map
 vrt_in = '%s/cid/icatch_map.vrt' % dir
 file_out = '%s/workspace/cid_map_region.tif' % output_dir
 os.system('rm -f %s' % file_out)
 os.system('gdalwarp -te %f %f %f %f %s %s' % (minlon,minlat,maxlon,maxlat,vrt_in,file_out))
 metadata['cid_map'] = file_out

 return metadata

def Create_Upscale_Template(metadata):
 
 res = metadata['upscaling']['res']
 rdir = metadata['rdir']
 file_cid = '%s/experiments/simulations/%s/postprocess/cids.vrt' % (rdir,metadata['experiment'])
 file = '%s/experiments/simulations/%s/postprocess/mapping.tif' % (rdir,metadata['experiment'])
 os.system('rm -rf %s' % file)
 #Create the upscaled grid -> summary per cell info
 os.system('gdalwarp %s -srcnodata -9999 -dstnodata -9999 -tr %.16f %.16f %s' % (file_cid,res,res,file))

 return

def Create_Upscale_Mapping(metadata,rank,bbox):
 
 rdir = metadata['rdir']
 file_cid = '%s/experiments/simulations/%s/postprocess/cids.vrt' % (rdir,metadata['experiment'])
 file_hrus = '%s/experiments/simulations/%s/postprocess/hrus.vrt' % (rdir,metadata['experiment'])
 workspace = '%s/experiments/simulations/%s/postprocess/workspace' % (rdir,metadata['experiment'])
 os.system('mkdir -p %s' % workspace)
 

 #Read in the upscaled mapping
 lats_upscale = bbox['lats_upscale']
 lons_upscale = bbox['lons_upscale']
 ilats_upscale = bbox['ilats_upscale']
 ilons_upscale = bbox['ilons_upscale']
 res_upscale = bbox['res_upscale']

 #Read in the fine scale catchment map
 lats_finescale = bbox['lats_finescale']
 lons_finescale = bbox['lons_finescale']
 ilats_finescale = bbox['ilats_finescale']
 ilons_finescale = bbox['ilons_finescale']
 res_finescale = bbox['res_finescale']
 
 #Read in the fine scale data
 dims = {'nx':lons_finescale.size,'ny':lats_finescale.size,
         'ixmin':ilons_finescale[0],'iymin':ilats_finescale[0]}
 #print(dims,ilons_finescale[0].dtype,ilats_finescale[0].dtype,lons_finescale.size,lats_finescale.size)
 #icatch_finescale = gdal_tools.read_raster_subarea(file_cid,dims)
 #hrus_finescale = gdal_tools.read_raster_subarea(file_hrus,dims)
 #window = rasterio.windows.Window(dims['iymin'],dims['ixmin'],dims['nx'],dims['ny'])
 window = rasterio.windows.Window(dims['ixmin'],dims['iymin'],dims['nx'],dims['ny'])
 icatch_finescale = rasterio.open(file_cid).read(1,window=window)
 hrus_finescale = rasterio.open(file_hrus).read(1,window=window)
 icatch_finescale[icatch_finescale < 0] = -9999
 hrus_finescale[icatch_finescale < 0] = -9999

 #Iterate through all upscaled cells and find their mappings
 lats_upscale = list(lats_upscale)
 lons_upscale = list(lons_upscale)
 print(rank,"Determining the mapping",flush=True)
 Output = {}
 icell = 0
 for lat in lats_upscale:
  ilat = lats_upscale.index(lat)
  #print "cell %d" % icell
  mask_lats = np.where((lats_finescale >= lats_upscale[ilat]-res_upscale/2) & 
                      (lats_finescale < lats_upscale[ilat]+res_upscale/2))[0]
  #for ilon in ilons_upscale[0:-1]:
  for lon in lons_upscale:
   ilon = lons_upscale.index(lon)
   #Determine the fine scale cells that fall into each coarse cell
   mask_lons = np.where((lons_finescale >= lons_upscale[ilon]-res_upscale/2) & 
                       (lons_finescale < lons_upscale[ilon]+res_upscale/2))[0]
   #Memorize the HRUs and their catchments for this coarse cell
   icatchs = icatch_finescale[mask_lats[0]:mask_lats[-1]+1,mask_lons[0]:mask_lons[-1]+1]
   hrus = hrus_finescale[mask_lats[0]:mask_lats[-1]+1,mask_lons[0]:mask_lons[-1]+1]
   #If the are no suitable candidates then skip the current cell
   mask_map = (icatchs >= 0) & (hrus >= 0)
   icatchs = icatchs[mask_map]
   hrus = hrus[mask_map]
   #icatchs = np.unique(icatchs[mask_map])
   if icatchs.size != 0:
    #Iterate through the catchments
    ncells = icatchs.size
    info = {}
    for icatch in np.unique(icatchs).astype(np.int):
     mask_icatch = icatchs == icatch
     info[icatch] = {'hru':[],'pct':[]}
     for hru in np.unique(hrus[mask_icatch]).astype(np.int):
       mask_hru = hrus[mask_icatch] == hru
       pct = float(np.sum(mask_hru))/float(ncells)
       info[icatch]['hru'].append(hru)
       info[icatch]['pct'].append(pct)
     info[icatch]['hru'] = np.array(info[icatch]['hru'])
     info[icatch]['pct'] = np.array(info[icatch]['pct'])
    Output[icell] = {'coords':{'ilat_local':ilat,'ilon_local':ilon},'info':info}
       
   #Update the cell id
   icell += 1

 #Output the mapping
 print(rank,"Writing out the mapping info",flush=True)
 pickle.dump(Output,open('%s/mapping_%d.pck' % (workspace,rank),'wb'),pickle.HIGHEST_PROTOCOL)
 
 return

def Map_Model_Output_Full(info,var,MPI):

 ncatch = info['ncatch']
 nt_in = info['nt_in']
 startdate = info['startdate']
 enddate = info['enddate']
 dt = info['dt']
 output_dir = info['output_dir']
 dir = info['dir']
 nt_out = info['nt_out']
 rank = MPI.COMM_WORLD.Get_rank()
 size = MPI.COMM_WORLD.Get_size()

 #Define the file
 file = '%s/%s_%d_%d.nc' % (output_dir,var,startdate.year,enddate.year)

 #Define the variables
 vars = [var,]

 #Create the netcdf file
 fp = h5py.File(file,'a',driver='mpio',comm=MPI.COMM_WORLD)

 #Open access to all the catchments
 fp_in = nc.Dataset('%s/catch_%d/of.nc' % (dir,0))

 #Read in the mapping
 hmll = fp_in.groups['latlon_mapping'].variables['hmll'][:].astype(np.int32)
 tmp = np.zeros(hmll.shape)
 mask = hmll >= 0
 hmll = hmll[mask]

 #Iterate through time steps
 tmp[:] = -9999.0
 for itime in np.arange(info['nt_out'])[rank::size]:
  print(rank,itime)
  #Compute the map
  data = fp_in.groups['catchment'].variables[var][itime,:]
  tmp[mask] = data[hmll]
  #Write the data
  fp[var][itime,:,:] = np.flipud(tmp[:])

 #Close the file
 fp.close()

 return

def Map_Model_Output(metadata,vars,rank,bbox,startdate,enddate):

 '''ncatch = info['ncatch']
 nt_in = info['nt_in']
 startdate = info['startdate']
 enddate = info['enddate']
 dt = info['dt']
 output_dir = info['output_dir']
 dir = info['dir']
 nt_out = info['nt_out']'''
 #startdate = metadata['idate']
 #enddate = metadata['fdate']
 #nt_in = metadata['nt_in']
 #nt_out = metadata['nt_out']
 rdir = metadata['rdir']
 workspace = '%s/experiments/simulations/%s/postprocess/workspace' % (rdir,metadata['experiment'])
 output_dir = '%s/experiments/simulations/%s/postprocess/output_dir' % (rdir,metadata['experiment'])
 os.system('mkdir -p %s' % output_dir)
 file_cid = '%s/experiments/simulations/%s/postprocess/cids.vrt' % (rdir,metadata['experiment'])

 #Read in the mapping info
 print("Reading in the mapping info",flush=True)
 mapping = pickle.load(open('%s/mapping_%d.pck' % (workspace,rank),'rb'))

 #Determine the unique catchments in the box
 #Read in the fine scale data
 dims = {'nx':bbox['lons_finescale'].size,'ny':bbox['lats_finescale'].size,
         'ixmin':bbox['ilons_finescale'][0],'iymin':bbox['ilats_finescale'][0]}
 #icatch_finescale = gdal_tools.read_raster_subarea(info['cid_map'],dims)
 #window = rasterio.windows.Window(dims['iymin'],dims['ixmin'],dims['ny'],dims['nx'])
 window = rasterio.windows.Window(dims['ixmin'],dims['iymin'],dims['nx'],dims['ny'])
 icatch_finescale = rasterio.open(file_cid).read(1,window=window)
 icatchs = np.unique(icatch_finescale[icatch_finescale >= 0]).astype(np.int)

 nlat = bbox['lats_upscale'].size#-1#metadata_upscale['ny']
 nlon = bbox['lons_upscale'].size#-1#metadata_upscale['nx']
 mask = np.zeros((nlat,nlon))
  
 #Extract true location info
 ilats_upscale = bbox['ilats_upscale']
 ilons_upscale = bbox['ilons_upscale']
 ilats_upscale_flipped = bbox['ilats_upscale_flipped']

 #Open access to all the catchments
 fps = {}
 for cid in icatchs:
  #file_output = '%s/catch_%d/output.nc' % (dir,icatch)
  file_output = '%s/experiments/simulations/%s/output_data/%d/%04d-%02d-%02d.nc' % (rdir,metadata['experiment'],cid,startdate.year,startdate.month,startdate.day)
  fps[cid] = nc.Dataset(file_output)

 #Determine nt_out
 nt_out = fps[cid]['data'].variables['trad'].shape[0]

 #Initialize the output
 output = {}
 for var in vars:
  output[var] = np.zeros((nt_out,nlat,nlon))

 #Iterate through all the cachments
 print(rank,"Begin: Reading and preparing the output",time.ctime(),flush=True)
 for icatch in icatchs:
  #flag_catchment = True
  data_catchment = {}
  for var in vars:
   #try:
   #data_catchment[var] = fps[icatch].groups['catchment'].variables[var][:,:]
   data_catchment[var] = fps[icatch]['data'].variables['%s' % var][:,:]
   #except:
   # flag_catchment = False
  #if flag_catchment == False:continue

  #Iterate through all the cell ids
  #cid = 0
  cids = mapping.keys()
  for cid in cids:
   info = mapping[cid]['info']
   ilat = mapping[cid]['coords']['ilat_local']
   ilon = mapping[cid]['coords']['ilon_local']
   if icatch in info:
    hrus = info[icatch]['hru']
    pcts = info[icatch]['pct']
    for var in vars:
     #data = np.sum(pcts*fps[icatch].groups['catchment'].variables[var][0:nt_in,hrus],axis=1)
     data = np.sum(pcts*data_catchment[var][:,hrus],axis=1)
     #Save the output
     output[var][:,ilat,ilon] += upscaling_fortran.time_average(data,data.size)
    #Add to the mask
    mask[ilat,ilon] += 1
 print(rank,"End: Reading and preparing the output",time.ctime(),flush=True)

 #Clear up the undefined
 for var in vars:
  output[var][:,mask == 0] = -9999.0
 
 #Create the dates array (monthly)
 '''dates_monthly = []
 date = datetime.datetime(startdate.year,startdate.month,startdate.day,0,0,0)
 fdate = datetime.datetime(enddate.year,enddate.month,enddate.day,23,0,0)
 timedelta = relativedelta.relativedelta(months=1)
 while date <= fdate:
  dates_monthly.append(date)
  date = date + timedelta
 dates_monthly = np.array(dates_monthly)'''
 #Create the dates array (daily)
 dates_daily = []
 timedelta = datetime.timedelta(hours=24)
 date = startdate
 while date <= enddate:
  dates_daily.append(date)
  date = date + timedelta
 dates_daily = np.array(dates_daily)
 #Create the dates array (hourly)
 dates_hourly = []
 #timedelta = datetime.timedelta(hours=24)
 timedelta = datetime.timedelta(hours=1)
 date = startdate
 while date <= enddate:
  dates_hourly.append(date)
  date = date + timedelta
 dates_hourly = np.array(dates_hourly)
 #Write the data
 for var in output:
  output[var] = np.fliplr(output[var][:])
 print(rank,"Begin: Writing the data",time.ctime())
 #Define the file
 random.shuffle(dates_daily)
 timedelta = relativedelta.relativedelta(days=1)
 for date in dates_daily:
  daily_dir = '%s/workspace/%04d%02d%02d' % (output_dir,date.year,date.month,date.day)
  os.system('mkdir -p %s' % daily_dir)
  file = '%s/%d.pck' % (daily_dir,rank)
  #print date,date+timedelta
  #mask = (dates_daily >= date) & (dates_daily < date + timedelta)
  mask = (dates_hourly >= date) & (dates_hourly < date + timedelta)
  #output_daily = {'coords':{'ilats':ilats_upscale_flipped,'ilons':ilons_upscale},
  output_hourly = {'coords':{'ilats':ilats_upscale_flipped,'ilons':ilons_upscale},
                  'vars':{}}
  for var in vars:
   output_hourly['vars'][var] = output[var][mask,:,:]

  #Save the data
  pickle.dump(output_hourly,open(file,'wb'),pickle.HIGHEST_PROTOCOL)
 print(rank,"End: Writing the data",time.ctime())

 return
