import netCDF4 as nc
import h5py
import netcdf_tools
import numpy as np
import os
import gdal_tools
import cPickle as pickle
import datetime
import time
import upscaling_fortran

def Initialize_Output_Files(info):

 startdate = info['startdate']
 enddate = info['enddate']
 nt_out = info['nt_out']
 dt = info['dt']
 vars = info['vars']
 nstripes = info['nstripes']
 ncores = info['ncores']
 output_dir = info['output_dir']
 date = startdate
 dt = 1
 nt_out = 24
 
 #Define the dimensions
 metadata = gdal_tools.retrieve_metadata('%s/workspace/mapping.tif' % info['output_dir'])
 nlon_chunk = int(metadata['nx']/ncores**0.5)/2
 nlat_chunk = int(metadata['ny']/ncores**0.5)/2
 ntime_chunk = nt_out
 dims = {'nlat':metadata['ny'],
         'nlon':metadata['nx'],
         'res':metadata['resx'],
         'minlon':metadata['minx'] + metadata['resx']/2,
         'minlat':metadata['miny'] + metadata['resy']/2,
         'undef':-9999.0,
         'chunksize':(nt_out,nlat_chunk,nlon_chunk)}

 #Define the directory
 #vardir = '%s/%s' % (info['output_dir'],var)
 #os.system('mkdir -p %s' % vardir)
 #os.system('lfs setstripe -c %d %s' % (nstripes,vardir))

 while date <= enddate:

  print "Creating the file for ",date
  #Define the file
  file = '%s/%04d%02d%02d.nc' % (output_dir,date.year,date.month,date.day)

  #Create the netcdf file
  fp = netcdf_tools.Create_NETCDF_File(dims,file,vars,vars,startdate,'%dhour' % dt,nt_out)

  #fill in the data
  for var in vars:
   fp.variables[var][:] = 0.1
  fp.close()

  #Update the time step
  date = date + datetime.timedelta(hours=nt_out)

 return

def Create_Virtual_Rasters(metadata):

 dir = metadata['dir']
 output_dir = metadata['output_dir']
 ncatch = metadata['ncatch']

 #Create the virtual raster of the hsus
 os.system('rm -rf %s/workspace/files.txt' % output_dir)
 files = []
 for icatch in xrange(ncatch):
  file = '%s/catch_%d/workspace/hsu_mapping_latlon.tif' % (dir,icatch)
  os.system('echo %s >> %s/workspace/files.txt' % (file,output_dir))
 os.system('gdalbuildvrt -srcnodata -9999 -input_file_list %s/workspace/files.txt %s/workspace/hsu_map.vrt' % (output_dir,output_dir))

 #Create the virtual raster of the catchment ids
 os.system('rm -rf %s/workspace/files.txt' % output_dir)
 files = []
 for icatch in xrange(ncatch):
  file = '%s/catch_%d/workspace/icatch_latlon.tif' % (dir,icatch)
  os.system('echo %s >> %s/workspace/files.txt' % (file,output_dir))
 os.system('gdalbuildvrt -srcnodata -9999 -input_file_list %s/workspace/files.txt %s/workspace/icatch_map.vrt' % (output_dir,output_dir))

 return

def Create_Upscale_Template(metadata):
 
 res = metadata['res']
 res = res/3600.0
 output_dir = metadata['output_dir']
 file = '%s/workspace/mapping.tif' % output_dir
 os.system('rm -rf %s' % file)
 #Create the upscaled grid -> summary per cell info
 os.system('gdalwarp %s/workspace/icatch_map.vrt -srcnodata -9999 -dstnodata -9999 -tr %.16f %.16f %s' % (output_dir,res,res,file))

 return

def Create_Upscale_Mapping(info,rank,bbox):
 
 output_dir = info['output_dir']

 #Read in the upscaled mapping
 #metadata_upscale = gdal_tools.retrieve_metadata('%s/workspace/mapping.tif' % output_dir)
 lats_upscale = bbox['lats_upscale']#np.linspace(metadata_upscale['miny'],metadata_upscale['maxy'],metadata_upscale['ny']+1)
 lons_upscale = bbox['lons_upscale']#np.linspace(metadata_upscale['minx'],metadata_upscale['maxx'],metadata_upscale['nx']+1)
 ilats_upscale = bbox['ilats_upscale']
 ilons_upscale = bbox['ilons_upscale']
 res_upscale = bbox['res_upscale']

 #Read in the fine scale catchment map
 #metadata_finescale = gdal_tools.retrieve_metadata('%s/workspace/icatch_map.vrt' % output_dir)
 lats_finescale = bbox['lats_finescale']
 lons_finescale = bbox['lons_finescale']
 ilats_finescale = bbox['ilats_finescale']
 ilons_finescale = bbox['ilons_finescale']
 res_finescale = bbox['res_finescale']
 #lats_finescale = np.linspace(metadata_finescale['miny'],metadata_finescale['maxy'],metadata_finescale['ny']+1)
 #lons_finescale = np.linspace(metadata_finescale['minx'],metadata_finescale['maxx'],metadata_finescale['nx']+1)
 
 #Read in the fine scale data
 dims = {'nx':lons_finescale.size,'ny':lats_finescale.size,
         'ixmin':ilons_finescale[0],'iymin':ilats_finescale[0]}
 icatch_finescale = gdal_tools.read_raster_subarea('%s/workspace/icatch_map.vrt' % output_dir,dims)
 hrus_finescale = gdal_tools.read_raster_subarea('%s/workspace/hsu_map.vrt' % output_dir,dims)
 icatch_finescale[icatch_finescale < 0] = -9999
 hrus_finescale[icatch_finescale < 0] = -9999

 #Iterate through all upscaled cells and find their mappings
 lats_upscale = list(lats_upscale)
 lons_upscale = list(lons_upscale)
 print rank,"Determining the mapping"
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
   #Save the cell id number
   #data_upscale[ilat,ilon] = icell
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
 print rank,"Writing out the mapping info"
 pickle.dump(Output,open('%s/workspace/mapping_%d.pck' % (output_dir,rank),'wb'),pickle.HIGHEST_PROTOCOL)
 
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
 fp_in = nc.Dataset('%s/catch_%d/output_data.nc' % (dir,0))

 #Read in the mapping
 hmll = fp_in.groups['latlon_mapping'].variables['hmll'][:].astype(np.int32)
 tmp = np.zeros(hmll.shape)
 mask = hmll >= 0
 hmll = hmll[mask]

 #Iterate through time steps
 tmp[:] = -9999.0
 for itime in np.arange(info['nt_out'])[rank::size]:
  print rank,itime
  #Compute the map
  data = fp_in.groups['catchment'].variables[var][itime,:]
  tmp[mask] = data[hmll]
  #Write the data
  fp[var][itime,:,:] = np.flipud(tmp[:])

 #Close the file
 fp.close()

 return

def Map_Model_Output(info,vars,MPI,bbox):

 ncatch = info['ncatch']
 nt_in = info['nt_in']
 startdate = info['startdate']
 enddate = info['enddate']
 dt = info['dt']
 output_dir = info['output_dir']
 dir = info['dir']
 nt_out = info['nt_out']
 rank = MPI.COMM_WORLD.Get_rank()

 #Define the file
 #vardir = '%s/%s' % (info['output_dir'],var)
 #file = '%s/%s_%d_%d.nc' % (vardir,var,startdate.year,enddate.year)

 #Define the variables
 #vars = [var,]

 #Create the netcdf file
 #fp = h5py.File(file,'a',driver='mpio',comm=MPI.COMM_WORLD)
 #fp.atomic = True

 #Read in the mapping info
 print "Reading in the mapping info"
 mapping = pickle.load(open('%s/workspace/mapping_%d.pck' % (output_dir,rank)))

 #Determine the unique catchments in the box
 #Read in the fine scale data
 dims = {'nx':bbox['lons_finescale'].size,'ny':bbox['lats_finescale'].size,
         'ixmin':bbox['ilons_finescale'][0],'iymin':bbox['ilats_finescale'][0]}
 icatch_finescale = gdal_tools.read_raster_subarea('%s/workspace/icatch_map.vrt' % output_dir,dims)
 icatchs = np.unique(icatch_finescale[icatch_finescale >= 0]).astype(np.int)

 #metadata_upscale = gdal_tools.retrieve_metadata('%s/workspace/mapping.tif' % output_dir)
 nlat = bbox['lats_upscale'].size#-1#metadata_upscale['ny']
 nlon = bbox['lons_upscale'].size#-1#metadata_upscale['nx']
 mask = np.zeros((nlat,nlon))
  
 #Extract true location info
 ilats_upscale = bbox['ilats_upscale']
 ilons_upscale = bbox['ilons_upscale']
 ilats_upscale_flipped = bbox['ilats_upscale_flipped']

 #Open access to all the catchments
 fps = {}
 for icatch in icatchs:#xrange(ncatch):
  fps[icatch] = nc.Dataset('%s/catch_%d/output_data.nc' % (dir,icatch))

 #Initialize the output
 output = {}
 for var in vars:
  output[var] = np.zeros((nt_out,nlat,nlon))

 #Iterate through all the cachments
 print rank,"Reading and preparing the output"
 for icatch in icatchs:#xrange(ncatch):
  #print 'catchment: %d' % icatch
  data_catchment = {}
  for var in vars:
   data_catchment[var] = fps[icatch].groups['catchment'].variables[var][:,:]

  #Iterate through all the cell ids
  #cid = 0
  for cid in mapping:
   info = mapping[cid]['info']
   ilat = mapping[cid]['coords']['ilat_local']
   ilon = mapping[cid]['coords']['ilon_local']
   if icatch in info:
    hrus = info[icatch]['hru']
    pcts = info[icatch]['pct']
    for var in vars:
     data = np.sum(pcts*data_catchment[var][0:nt_in,hrus],axis=1)
     #Save the output
     output[var][:,ilat,ilon] += upscaling_fortran.time_average(data,nt_out)
    #Add to the mask
    mask[ilat,ilon] += 1

    #Update the cell id
    cid += 1

 #Clear up the undefined
 for var in vars:
  output[var][:,mask == 0] = -9999.0
 
 dates = []
 #Create the dates array
 date = datetime.datetime(startdate.year,startdate.month,startdate.day,0,0,0)
 fdate = datetime.datetime(enddate.year,enddate.month,enddate.day,0,0,0)
 timedelta = datetime.timedelta(hours=1)
 while date <= fdate:
  dates.append(date)
  date = date + timedelta
 dates = np.array(dates)
 #Write the data
 #fp[var][:,ilats_upscale[0]:ilats_upscale[-1]+1,ilons_upscale[0]:ilons_upscale[-1]+1] = np.fliplr(output[:])
 for var in output:
  output[var] = np.fliplr(output[var][:])
 print rank,"Writing the data"
 #Define the file
 #vardir = '%s/%s' % (info['output_dir'],var)
 timedelta = datetime.timedelta(hours=24)
 date = startdate
 while date <= enddate:
  print rank,date
  file = '%s/%04d%02d%02d.nc' % (output_dir,date.year,date.month,date.day)
  mask = (dates >= date) & (dates < date + timedelta)
  #Open the netcdf file
  #fp = h5py.File(file,'a',driver='mpio',comm=MPI.COMM_WORLD)
  #fp.atomic = True
  fp = h5py.File(file,'a')
  for var in vars:
   #dset = fp[var]
   fp[var][:,ilats_upscale_flipped[0]:ilats_upscale_flipped[-1]+1,ilons_upscale[0]:ilons_upscale[-1]+1] = output[var][mask,:,:]
  #Close the file
  fp.close()
  #Update the time step
  date = date + timedelta

 return

'''startdate = datetime.datetime(2004,1,1,0)
enddate = datetime.datetime(2004,12,31,23)
#enddate = datetime.datetime(2004,1,31,0)
dt = 24 #hours
ncatch = 9
nt_in = 24*((enddate - startdate).days+1)
nt_out = nt_in/dt
dir = '/home/freeze/nchaney/LittleWashitaRegion'
res = 90 #arcsec
output_dir = '%s/%darsec' % (dir,res)
os.system('rm -rf %s' % output_dir) 
os.system('mkdir %s' % output_dir)

#Create the virtual rasters
#Create_Virtual_Rasters()

#Create the upscaling template
#Create_Upscale_Template(res)

#Create upscaling mapping
#Create_Upscale_Mapping()

#Map the data
vars = ['smc1',]#'prcp','smc1','lh','sh','qout_surface','qout_subsurface']
for var in vars:
 print var
 Map_Model_Output(ncatch,nt_in,startdate,enddate,dt,var,output_dir,nt_out)

#Output the upscaled data
#Write_Upscaled_Data(nt,startdate,dt)

#Split area into n processes
#Have each process work on its own segment -> use parallel netcdf to write the data'''
