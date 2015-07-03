from mpi4py import MPI
import sys
import pickle
import gdal_tools
import upscaling_python
import numpy as np

def Determine_Bounding_Box(info,MPI):

 output_dir = info['output_dir']

 #Retrieve metadata for entire region
 metadata_upscale = gdal_tools.retrieve_metadata('%s/workspace/mapping.tif' % output_dir)
 res_upscale = metadata_upscale['resx']
 lats_upscale = np.linspace(metadata_upscale['miny']+res_upscale/2,metadata_upscale['maxy']-res_upscale/2,metadata_upscale['ny'])
 lons_upscale = np.linspace(metadata_upscale['minx']+res_upscale/2,metadata_upscale['maxx']-res_upscale/2,metadata_upscale['nx'])
 lats_upscale_flipped = np.flipud(lats_upscale)
 metadata_finescale = gdal_tools.retrieve_metadata('%s/workspace/icatch_map.vrt' % output_dir)
 res_finescale = metadata_finescale['resx']
 lats_finescale = np.linspace(metadata_finescale['miny']+res_finescale/2,
                  metadata_finescale['maxy']-res_finescale/2,metadata_finescale['ny'])
 lons_finescale = np.linspace(metadata_finescale['minx']+res_finescale/2,
                  metadata_finescale['maxx']-res_finescale/2,metadata_finescale['nx'])

 #Determine its bounding region
 comm = MPI.COMM_WORLD
 size = comm.size
 rank = comm.rank
 ibox = int(rank) / int(size**0.5)
 jbox = int(rank % size**0.5)

 #Determine the lats/lons and ilats/ilons for the bounding box (coarsescale)
 nlat_upscale = lats_upscale.size/size**0.5
 nlon_upscale = lons_upscale.size/size**0.5
 lats_upscale_box = lats_upscale[ibox*nlat_upscale:(ibox+1)*nlat_upscale]#+1]
 lons_upscale_box = lons_upscale[jbox*nlon_upscale:(jbox+1)*nlon_upscale]#+1]
 ilats_upscale_box = np.where(np.in1d(lats_upscale,lats_upscale_box))[0]
 ilons_upscale_box = np.where(np.in1d(lons_upscale,lons_upscale_box))[0]
 ilats_upscale_flipped_box = np.where(np.in1d(lats_upscale_flipped,lats_upscale_box))[0]

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

#Determine the rank and size
comm = MPI.COMM_WORLD
rank = comm.rank
size = comm.size

#Read in the directory
dir = sys.argv[1]

#Read in the metadata
metadata = pickle.load(open('%s/workspace/metadata.pck' % dir))

if metadata['type'] == 'regional':

 #Determine the sites box
 bbox_metadata = Determine_Bounding_Box(metadata,MPI)

 #Create upscaling mapping
 upscaling_python.Create_Upscale_Mapping(metadata,rank,bbox_metadata)

 #Map the data
 #vars = metadata['vars']#['smc1',]#'prcp','smc1','lh','sh','qout_surface','qout_subsurface']
 vars = metadata['vars']
 #for var in vars:
 upscaling_python.Map_Model_Output(metadata,vars,MPI,bbox_metadata)

elif metadata['type'] == 'catchment':

 #Map the data
 vars = metadata['vars']#['smc1',]#'prcp','smc1','lh','sh','qout_surface','qout_subsurface']
 for var in vars:
  print rank,var
  upscaling_python.Map_Model_Output_Full(metadata,var,MPI) 
