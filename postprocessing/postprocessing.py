import warnings
warnings.filterwarnings('ignore')
import postprocessing.upscaling_python as upscaling_python
import geospatialtools.gdal_tools as gdal_tools
import numpy as np
import datetime
import sys

def Read_Metadata_File(file):

 import json
 metadata = json.load(open(file))

 return metadata

def driver(comm,metadata_file):

 size = comm.Get_size()
 rank = comm.Get_rank()
 #Read in the metadata
 #metadata_file = '%s/metadata.json' % edir
 metadata = Read_Metadata_File(metadata_file)
 metadata['idate'] = datetime.datetime(metadata['startdate']['year'],
                           metadata['startdate']['month'],
                           metadata['startdate']['day'],0)
 metadata['fdate'] = datetime.datetime(metadata['enddate']['year'],
                           metadata['enddate']['month'],
                           metadata['enddate']['day'],0) + datetime.timedelta(days=1) - datetime.timedelta(seconds=metadata['dt'])
 startdate = metadata['idate']
 enddate = metadata['fdate']
 rdir = metadata['rdir']
 edir = '%s/experiments/simulations/%s' % (rdir,metadata['experiment'])

 #Create the upscale template
 if rank == 0:
  upscaling_python.Create_Upscale_Template(metadata) #laura, uncommented
 #comm.Barrier()

 #Determine the sites box
 print(rank,"Determing the bounding box",flush=True)
 bbox_metadata = Determine_Bounding_Box(metadata,rank,size)

 #Create upscaling mapping
 print(rank,"Creating the upscaling mapping",flush=True)
 upscaling_python.Create_Upscale_Mapping(metadata,rank,bbox_metadata) #laura, uncommented
 
 #Map the data
 vars = metadata['upscaling']['vars']
 #metadata['nt_in'] = 365*24
 #metadata['nt_out'] = 365*24
 for year in range(metadata['idate'].year,metadata['fdate'].year+1):
  if year == metadata['idate'].year:
   startdate = metadata['idate']
   enddate = datetime.datetime(year,12,31,23)
  elif year == metadata['fdate'].year:
   startdate = datetime.datetime(year,1,1,0)
   enddate = metadata['fdate']
  else:
   startdate = datetime.datetime(year,1,1,0)
   enddate = datetime.datetime(year,12,31,23)

  upscaling_python.Map_Model_Output(metadata,vars,rank,bbox_metadata,startdate,enddate)
  #print('upscale_python_map',flush=True)
  #Pause until all files have been processed
  comm.Barrier()

  #Create files
  print(rank,"Creating the output files (%d)" % year,flush=True)
  upscaling_python.Create_Output_Files(metadata,rank,size,vars,startdate,enddate)
  #Pause until all files have been processed
  comm.Barrier()

 return

def Determine_Bounding_Box(metadata,rank,size):
 print(rank,size,flush=True)
 sys.exit('para')
 rdir = metadata['rdir']
 file_cid = '%s/experiments/simulations/%s/postprocess/cids.vrt' % (rdir,metadata['experiment'])
 file_mapping = '%s/experiments/simulations/%s/postprocess/mapping.tif' % (rdir,metadata['experiment'])
 #Retrieve metadata for entire region
 metadata_upscale = gdal_tools.retrieve_metadata(file_mapping)
 res_upscale = metadata_upscale['resx']
 lats_upscale = np.linspace(metadata_upscale['miny']+res_upscale/2,metadata_upscale['maxy']-res_upscale/2,metadata_upscale['ny'])
 lons_upscale = np.linspace(metadata_upscale['minx']+res_upscale/2,metadata_upscale['maxx']-res_upscale/2,metadata_upscale['nx'])
 lats_upscale_flipped = np.flipud(lats_upscale)
 metadata_finescale = gdal_tools.retrieve_metadata(file_cid)
 res_finescale = metadata_finescale['resx']
 lats_finescale = np.linspace(metadata_finescale['miny']+res_finescale/2,
                  metadata_finescale['maxy']-res_finescale/2,metadata_finescale['ny'])
 lons_finescale = np.linspace(metadata_finescale['minx']+res_finescale/2,
                  metadata_finescale['maxx']-res_finescale/2,metadata_finescale['nx'])

 #Determine its bounding region
 #comm = MPI.COMM_WORLD
 #size = comm.size
 #rank = comm.rank
 #ibox = int(rank) / int(size**0.5)
 #jbox = int(rank % size**0.5)
 split_size = int(np.ceil(lons_upscale.size/size))
 ilons_upscale_min = rank*split_size
 ilons_upscale_max = (rank+1)*split_size
 if rank == size-2:
  ilons_upscale_max = ilons_upscale_max - 1
 if rank == size-1:
  ilons_upscale_min = ilons_upscale_min - 1
 if ilons_upscale_max > lons_upscale.size:ilons_upscale_max = lons_upscale.size

 #Determine the lats/lons and ilats/ilons for the bounding box (coarsescale)
 lats_upscale_box = lats_upscale#[ibox*nlat_upscale:(ibox+1)*nlat_upscale]#+1]
 #lons_upscale_box = lons_upscale[jbox*nlon_upscale:(jbox+1)*nlon_upscale]#+1]
 lons_upscale_box = lons_upscale[ilons_upscale_min:ilons_upscale_max+1]
 ilats_upscale_box = np.where(np.in1d(lats_upscale,lats_upscale_box))[0]
 ilons_upscale_box = np.where(np.in1d(lons_upscale,lons_upscale_box))[0]
 ilats_upscale_flipped_box = np.where(np.in1d(lats_upscale_flipped,lats_upscale_box))[0]
 
 print(rank,np.where((lons_finescale+res_finescale/2 <= lons_upscale_box[-1]+res_upscale/2) & (lons_finescale-res_finescale/2 >= lons_upscale_box[0]-res_upscale/2))[0],flush=True)
 sys.exit('para') 

 #Determine the lats/lons and ilats/ilons for the bounding box (finescale)
 ilats_finescale_box = np.where((lats_finescale+res_finescale/2 <= lats_upscale_box[-1]+res_upscale/2) &
                                (lats_finescale-res_finescale/2 >= lats_upscale_box[0]-res_upscale/2))[0]
 ilons_finescale_box = np.where((lons_finescale+res_finescale/2 <= lons_upscale_box[-1]+res_upscale/2) &
                                (lons_finescale-res_finescale/2 >= lons_upscale_box[0]-res_upscale/2))[0]
 #Skips one?
 lats_finescale_box = lats_finescale[ilats_finescale_box]
 lons_finescale_box = lons_finescale[ilons_finescale_box]
 bbox_metadata = {'lats_upscale':lats_upscale_box,'ilats_upscale':ilats_upscale_box,
                  'lons_upscale':lons_upscale_box,'ilons_upscale':ilons_upscale_box,
                  'res_upscale':res_upscale,'ilats_upscale_flipped':ilats_upscale_flipped_box,
                  'lats_finescale':lats_finescale_box,'ilats_finescale':ilats_finescale_box,
                  'lons_finescale':lons_finescale_box,'ilons_finescale':ilons_finescale_box,
                  'res_finescale':res_finescale}

 return bbox_metadata
