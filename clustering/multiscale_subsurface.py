import warnings
warnings.filterwarnings('ignore')
import numpy as np
import netCDF4 as nc
import os
import sys
import pickle
from geospatialtools import gdal_tools
from geospatialtools import terrain_tools
from mpi4py import MPI
import time

def build_adjacency_from_labels(labels, nodata=-9999):
  """Vectorized adjacency builder that counts shared edges between labeled regions.
  Returns a square matrix of shape (max_label, max_label) with counts of shared
  boundary pixels between labeled regions. Assumes labels are 1-based integers.
  """
  # Work on a plain ndarray
  if isinstance(labels, np.ma.MaskedArray):
    arr = labels.filled(nodata)
  else:
    arr = labels

  unique = np.unique(arr)
  unique = unique[unique != nodata]
  if unique.size == 0:
    return np.zeros((0, 0), dtype=np.int64)

  max_label = int(unique.max())
  adj = np.zeros((max_label, max_label), dtype=np.int64)

  # Define 4-neighbour shifts (up, down, left, right)
  # For each shift, slice the overlapping window and accumulate pairs
  # Up
  c = arr[1:, :]
  n = arr[:-1, :]
  mask = (c != nodata) & (n != nodata) & (c != n)
  if mask.any():
    centers = c[mask].astype(int) - 1
    neighs = n[mask].astype(int) - 1
    np.add.at(adj, (centers, neighs), 1)

  # Down
  c = arr[:-1, :]
  n = arr[1:, :]
  mask = (c != nodata) & (n != nodata) & (c != n)
  if mask.any():
    centers = c[mask].astype(int) - 1
    neighs = n[mask].astype(int) - 1
    np.add.at(adj, (centers, neighs), 1)

  # Left
  c = arr[:, 1:]
  n = arr[:, :-1]
  mask = (c != nodata) & (n != nodata) & (c != n)
  if mask.any():
    centers = c[mask].astype(int) - 1
    neighs = n[mask].astype(int) - 1
    np.add.at(adj, (centers, neighs), 1)

  # Right
  c = arr[:, :-1]
  n = arr[:, 1:]
  mask = (c != nodata) & (n != nodata) & (c != n)
  if mask.any():
    centers = c[mask].astype(int) - 1
    neighs = n[mask].astype(int) - 1
    np.add.at(adj, (centers, neighs), 1)

  return adj

def get_neighboring_area_ids(raster, specific_area_id):
  """Vectorized neighbor-ID finder for a specific area id.
  Only checks the four direct neighbors for cells that belong to the area.
  """
  neighboring_ids = set()
  # operate on filled array if masked
  if isinstance(raster, np.ma.MaskedArray):
    arr = raster.filled(-9999)
  else:
    arr = raster

  coords = np.argwhere(arr == specific_area_id)
  directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
  nrows, ncols = arr.shape
  for i, j in coords:
    for dx, dy in directions:
      x, y = i + dx, j + dy
      if 0 <= x < nrows and 0 <= y < ncols:
        neighbor_id = int(arr[x, y])
        if neighbor_id != specific_area_id and neighbor_id != -9999:
          neighboring_ids.add(neighbor_id)
  return neighboring_ids
   
def generate_risfu_maps(regunits,cid):
 '''Function to copy and recategorize the type of regional units {'clusters','basins'}'''
 #Create units with consecutive names
 risfu_ids   = np.unique(regunits);
 risfu_ids   = risfu_ids[risfu_ids!=-9999]
 recat_risfu = np.copy(regunits)
 recat = 1; #idx
 for unit_id in risfu_ids:
  if unit_id != -9999:
   recat_risfu[regunits == unit_id] = recat
   recat += 1;
 #Create regional units with a larger id based on cid.
 recat_risfu[recat_risfu != -9999] = recat_risfu[recat_risfu != -9999] + cid*1000
 return recat_risfu

#-----------------------------------------------------------------------------------------------------------------------------------
# Define Main Functions
def generate_regional_and_intermediate_units_maps(cid,edir,rdir,hydroblocks_info):
 '''Function to create maps of regional units and vrt'''
 #(1) Define variables
 input_dir = '%s/%d' % (edir,cid)
 data_dir = "%s/data/cids/%d" % (rdir,cid)
 metadata = gdal_tools.retrieve_metadata('%s/mask_latlon.tif' % data_dir)
 metadata['nodata'] = -9999.0
 
 #(2) Compute the regional units based on cluster of watersheds or basins defined in HMC
 if hydroblocks_info['multiscale_subsurface']['RISFU_type']=='clusters':
  basin_clusters = gdal_tools.read_data('%s/basin_clusters_latlon.tif' % input_dir)
  recat_risfu = generate_risfu_maps(basin_clusters.data,cid)
 elif hydroblocks_info['multiscale_subsurface']['RISFU_type']=='basins':
  basins = gdal_tools.read_data('%s/basins_latlon.tif' % input_dir)
  recat_risfu = generate_risfu_maps(basins.data,cid)
 
 #(3) Write out the RISFU maps for multiscale scheme
 os.makedirs('%s/groundwater' % input_dir, exist_ok=True) #create directory to store subsurface files 
 file_ca = '%s/groundwater/risfu_map.tif' % input_dir
 gdal_tools.write_raster(file_ca,metadata,recat_risfu) #write the file out
 
 return

def create_vrt_risfu(workspace,edir,size):
 os.makedirs(workspace, exist_ok=True)
 odir = '%s/regunits.vrt' % workspace
 list_files = '%s/file_paths.txt' % workspace #Create list of risfu files for each cid
 with open(list_files,'w') as f:
  for sub in range(size):
   file = '%s/%d/groundwater/risfu_map.tif' % (edir,sub+1)
   f.write(file + '\n')
 os.system('gdalbuildvrt -input_file_list %s %s' % (list_files,odir))
 return
  
def create_risfu_mask(cid,edir,rdir,workspace):
 #(1) Define variables
 data_dir = "%s/data/cids/%d" % (rdir,cid)
 metadata = gdal_tools.retrieve_metadata('%s/mask_latlon.tif' % data_dir)
 metadata['nodata'] = -9999.0
 #(2) Create cliped files for subdomain interaction
 minx = metadata['minx']
 miny = metadata['miny']
 maxx = metadata['maxx']
 maxy = metadata['maxy']
 res  = abs(metadata['resx'])
 lproj    = metadata['proj4']
 file_in  = '%s/regunits.vrt' % workspace
 file_out = '%s/risfu_%s.tif' % (workspace,cid)
 #clip the regional units from the entire domain
 os.system('gdalwarp -t_srs \'%s\' -dstnodata -9999 -tr %.16f %.16f -te %.16f %.16f %.16f %.16f -q -overwrite %s %s' % (lproj,res,res,minx,miny,maxx,maxy,file_in,file_out))
 
 #Creating the clipped array
 mask = gdal_tools.read_data('%s/mask_latlon.tif' % data_dir).data
 rows , cols = np.where(mask == int(cid)) #Get the rows and columns where the conditions is true
 
 miny_bf = rows.min() - 1; #Create the mask with the buffer for cid 1. 1 cells will do
 maxy_bf = rows.max() + 1;
 minx_bf = cols.min() - 1;
 maxx_bf = cols.max() + 1;
  
 risfu_clip = gdal_tools.read_data(file_out).data #Get the regional/intermediate units with the neighbouring areas from other subdomains
 risfu_clip = np.ma.masked_array(risfu_clip,risfu_clip==-9999)

 mask_risfu = risfu_clip[miny_bf : maxy_bf + 1, minx_bf : maxx_bf + 1]; #+1 to account for python indexing
 clip_file = '%s/risfu_%s.pck' % (workspace,cid)
 pickle.dump(mask_risfu,  open(clip_file, 'wb'));
 return mask_risfu

def generate_risfu_aggregated_properties(cid,edir,rdir,workspace,hydroblocks_info):
 '''Function to compute hru aggregated properties for the risfu areas'''
 #(1) Define variables
 input_dir = '%s/%d' % (edir,cid)
 data_dir = "%s/data/cids/%d" % (rdir,cid)
 metadata = gdal_tools.retrieve_metadata('%s/mask_latlon.tif' % data_dir)
 metadata['nodata'] = -9999.0
 dz = hydroblocks_info['dz']
 mask_object = gdal_tools.read_data('%s/mask_latlon.tif' % data_dir)
 terrain_tools.calculate_area(mask_object)
 resx = np.mean(mask_object.area**0.5) #all pixels in the subdomain have the same resolution in x and y; still not ideal a but much better than resx = 90...
 hrus = gdal_tools.read_data('%s/hru_mapping_latlon.tif' % input_dir).data
 risfu_org = gdal_tools.read_data('%s/groundwater/risfu_map.tif' % input_dir).data
 
 #(2) Get the number of hrus and regional/intermediate units
 gwunits_id = np.unique(risfu_org) #ids for regional/intermediate units
 gwunits_id = gwunits_id[gwunits_id!=-9999]
 num_gwus   = gwunits_id.size
 hru_ids    = np.unique(hrus) #ids of hrus
 hru_ids    = hru_ids[hru_ids!=-9999]
 num_hrus   = hru_ids.size
 #print('num gw units: ', num_gwus,flush=True)
 #print('num hrus: '    , num_hrus,flush=True)
 
 #(3) Counting how many hrus there are inside each regional/intermediate unit
 gwunits_count = {};
 for i in gwunits_id:
  gwunits_count['gw_units_%d' % i] = {}
  hrus_risfu = hrus[risfu_org == i]; #hrus in risfu
  gwunits_count['gw_units_%d' % i]['data_hrus'] = hrus_risfu[:];
  unique, counts = np.unique(hrus_risfu, return_counts=True)
  gwunits_count['gw_units_%d' % i]['counts_hrus'] = dict(zip(unique, counts)) 
 #print('Counting HRUS inside basins done',flush=True)
 
 #(4) Computing the fraction of area covered by each hrus within each groundwater unit
 areas_gws = np.zeros((num_gwus,num_hrus), dtype=np.float64);
 n = 0 #the array structure follows the gwunits_id array
 for i in gwunits_id:
  area_risfu = sum(gwunits_count['gw_units_%d' % i]['counts_hrus'].values()) #the area of each regional/intermediate unit
  for nhru, hru in enumerate(hru_ids):
   try:
    areas_gws[n,nhru] = gwunits_count['gw_units_%d' % i]['counts_hrus'][hru] / area_risfu; #compute the percentage of area covered by the hru
   except:
    areas_gws[n,nhru] = .0; #no hru present in basin
  n += 1 #move on to the next regional/intermediate unit

 #print('Area of HRUS inside basins done',flush=True)
 pickle.dump(areas_gws,  open('%s/groundwater/area_hrus.pck' % input_dir, 'wb')) #save file

 #(5) Compute the weighted average of the soil properties for the regional/intermediate areas
 total_area = np.sum(areas_gws, axis = 1); #total area per regional/intermediate units
 units      = np.count_nonzero(total_area); #number of regional/intermediate units with hrus 
 area       = np.zeros((units,areas_gws.shape[1]));
 area[:]    = areas_gws[total_area > 0]; # area of regional/intermediate units
 
 #file = '%s/input_file_basins.nc' % experiment
 file = '%s/input_file.nc' % input_dir  #Open input file with the parameters of each hru
 fp = nc.Dataset(file)
    
 bb = fp['parameters']['BB'][:].astype(np.float64)
 sp = fp['parameters']['SATPSI'][:].astype(np.float64)
 ks = fp['parameters']['SATDK'][:].astype(np.float64)
 tr = fp['parameters']['DRYSMC'][:].astype(np.float64)
 ts = fp['parameters']['MAXSMC'][:].astype(np.float64)
 m  = fp['parameters']['m'][:].astype(np.float64)
 
 fp.close()

 names = ['bb', 'sp', 'ks', 'tr', 'ts', 'm']  # name of parameters
 param_gws = {}
 tmp = np.ones(area.shape[1])
 params = {'bb': bb, 'sp': sp, 'ks': ks, 'tr': tr, 'ts': ts, 'm': m}
 for var in names:
   if var == 'm':
     param_gws[var] = np.dot(params[var][:] * area, tmp)
   else:
     parameters = np.zeros((num_gwus, len(dz)), dtype=np.float64)
     for il in range(len(dz)):
       parameters[:, il] = np.dot(params[var][:, il] * area, tmp)
     param_gws[var] = parameters
  #print('Hrus aggregation of subsurface parameters',flush=True)
 pickle.dump(param_gws,  open('%s/groundwater/param_aggregation.pck' % input_dir, 'wb')) #save file
 
 #(7) Computing the mean elevation for each unit
 dem = gdal_tools.read_data('%s/dem_latlon.tif' % data_dir).data #Open the elevation file
 dem = np.ma.masked_array(dem,dem==-9999)

 gw_elv = np.zeros(num_gwus); #Compute the elevation for each unit
 n = 0
 for i in gwunits_id:
  gw_elv[n] = np.mean(dem[risfu_org == i]);
  n += 1

 dem_gws = np.copy(risfu_org); #Assing the new elevation values to the raster of regional units
 for i in gwunits_id:
  dem_gws[risfu_org == i] = np.mean(dem[risfu_org == i]);
  
 file_ca = '%s/groundwater/dem_reg_map.tif' % input_dir
 gdal_tools.write_raster(file_ca,metadata,dem_gws) #write the file out

 pickle.dump(gw_elv,open('%s/groundwater/elv.pck' % input_dir,'wb')); #save file
 
 return risfu_org,gwunits_id,num_gwus,resx,input_dir, gwunits_count
 
def compute_connection_matrix_risfu(risfu_org,gwunits_id,ngwus,resx,input_dir,gwunits_count):
 '''Function to compute the connection matrix (reduced) between the regional units for intermediate subsurface flow'''
 #(1) Defining the shared length between units and computing the distance between centroids
 recat_gw = np.copy(risfu_org) # copy array with the risfu areas to change its names
 nid = 1;
 for i in gwunits_id:
  recat_gw[risfu_org == i] = nid
  nid += 1
  
  recat_cp = np.copy(recat_gw)
  gwus_counts = build_adjacency_from_labels(recat_cp, nodata=-9999)
  w_gws = gwus_counts * resx  # Converting the number of cells to distance in meters
 #print('Defining the shared length between units',flush=True)

 #To compute dx matrix we need to compute the areas of each gw unit and then dived the area by the shared length between them.
 gw_area = np.zeros((ngwus));  #Computing area of regional/intermediate units
 n = 0
 for i in gwunits_id:
  gw_area[n] = sum(gwunits_count['gw_units_%d' % i]['counts_hrus'].values())*resx*resx
  n += 1
 #print('area for',rank + 1,gw_area.shape,flush=True)
 
 #The final dx matrix assumes that the distance between units is equal to the total area divided by the shared length. The results are then further 
 #divided by 2 assuming a centroid and then, the centroid of the coneected unit is added to the distance.
 dx_gws = 0.5*((gw_area/w_gws)+(gw_area/w_gws).T);
 #print('Defining the distance between centroids',flush=True)
 pickle.dump(gw_area,  open('%s/groundwater/area.pck' % input_dir, 'wb')); #save file
 
 #(7) Reshaping matrix of connections to match changes in hb for a reduced matrix
 #print('Reshaping the matrix of connections',flush=True)
 rows, cols = np.where(w_gws > 0) #Obtaining the indices for the connections from the big matrix 
 maxcox     = np.count_nonzero(w_gws, axis = 1).max()

 connections = np.zeros((ngwus,maxcox), dtype = np.int32);
 wreg        = np.zeros((ngwus,maxcox));
 dxreg       = np.zeros((ngwus,maxcox));

 for i in range(ngwus):
  tmp = cols[rows==i]
  if tmp.size != ngwus-1:
   tmp = np.append(tmp, np.zeros(maxcox-tmp.size)+i);
  connections[i,:] = tmp[:];
  wreg[i,:] = w_gws[i,connections[i,:]]; #converting all the filled values with widht of the coneection
  dxreg[i,:] = dx_gws[i,connections[i,:]]; #converting all the filled values with widht of the coneection

 pickle.dump(wreg,  open('%s/groundwater/w_reg.pck' % input_dir, 'wb')); #save file
 pickle.dump(dxreg, open('%s/groundwater/dx_reg.pck' % input_dir, 'wb')); #save file
 pickle.dump(connections,open('%s/groundwater/conx_reg.pck' % input_dir, 'wb')); #save file
 
 return

def create_files_regional_interaction(cid,mask_risfu,resx,input_dir,edir):
 '''Function to compute the connection matrix between the risfu at the boundaries'''
 #(1)Count the neighbours that actually shared a border with the subdomain in consideration
 # Recategorize the groundwater units. Only with the number of the cids
 gw_ids = np.unique(mask_risfu).data.astype(int)
 if -9999 in gw_ids:
  gw_ids = np.delete(gw_ids, np.where(gw_ids == -9999))
 
 #print('gw ids: %s' % str(cid), gw_ids ,flush = True)
 cids = gw_ids // 1000
 recat_cids = np.copy(mask_risfu)

 nid = 0;
 for i in gw_ids:
  recat_cids[mask_risfu == i] = cids[nid]
  nid += 1

 recat_cids = np.ma.masked_array(recat_cids,recat_cids==-9999)
 #print('cid %s: ' %str (cid),np.unique(cids) ,flush = True)
  
 # Specific area ID to check for neighbors
 specific_area_id =  cid
 # Get neighboring area IDs for the specific area
 neighboring_ids = get_neighboring_area_ids(recat_cids, specific_area_id)
 print(f'Neighboring area IDs for area {specific_area_id}: {neighboring_ids}',flush=True)
 
 neighboring_ids.add(int(cid))
 new = np.copy(recat_cids);
 new[:,:] = -9999
 for i in neighboring_ids:
  new[recat_cids == i] = i
  
 copy_risfu = np.zeros(mask_risfu.shape) - 9999
 copy_risfu[new != -9999] = mask_risfu[new != -9999]
 copy_risfu = np.ma.masked_array(copy_risfu,copy_risfu==-9999)
 
 #(2) Recategorize the groundwater units. This allow us to count the gridcells and putting the total count in an array
 gw_ids   = np.unique(copy_risfu).data.astype(int)
 recat_gw = np.copy(copy_risfu)

 if -9999 in gw_ids:
  gw_ids = np.delete(gw_ids, np.where(gw_ids == -9999))
  
 nid = 1;
 for i in gw_ids:
  recat_gw[mask_risfu == i] = nid
  nid += 1

 recat_gw = np.ma.masked_array(recat_gw,recat_gw==-9999)
 
 #(3) Count the neighbours using the basins with the new id recat_gw
 nunits = gw_ids.size;
 # Build adjacency counts using the vectorized helper
 recat_cp = np.copy(recat_gw)
 gw_counts = build_adjacency_from_labels(recat_cp, nodata=-9999)
 #print('counts for cid: %s' % str(rank+1), flush = True)
 
 #(4) Create the matrix of widths and the distance between centroids
 w_gws = gw_counts*resx; #assuming 30 m pixel

 #To compute the distance between centroids we must bring the corresponding area of the groundwater units stored in the corresponding directory 

 #find the domains  involved
 reg_conx = {};
 cids = np.unique(gw_ids // 1000);

 #create a dictionary with the following structure {cid: [array containing the index of the groundwater units from that subdomain]}
 for cid in cids:
  reg_conx[cid] = gw_ids[np.logical_and(gw_ids>=cid*1000, gw_ids<((cid+1)*1000))].data - cid*1000;
  
 pickle.dump(reg_conx,  open('%s/groundwater/reg_ids.pck' % input_dir, 'wb')); #save file
 #print('reg conx for rank: %s' % str(rank+1), reg_conx, flush = True)

 #Compute the distance to the centroid
 gw_area = [] #bring the area from the other subdomains and append it to the list
 for cid in cids:
  file = '%s/%s/groundwater/area.pck' % (edir,str(cid))
  while not os.path.exists(file):
    print(f"Waiting for file: {file}",flush=True)
    time.sleep(5)  # Check every 5 seconds
  #open the area file
  area      = pickle.load(open(file,'rb'))
  #Units of that cid interacting with the subdomain
  units_cid = reg_conx[cid].astype(int)
  #Extract the area for those units
  tmp       = area[units_cid - 1] #-1 due to python indexing
  gw_area   = np.append(gw_area, tmp) #an array with size equal to units involved and values corresponding to areas

 dx_gws  = 0.5*((gw_area/w_gws)+(gw_area/w_gws).T); #assuming a rectangle with same area
 #print('w and dx for cid: %s' % str(rank+1), flush = True)
 
 #(5) Creating the matrix of connections
 rows, cols = np.where(w_gws > 0) #Obtaining the indices for the connections from the big matrix 
 maxcox     = np.count_nonzero(w_gws, axis = 1).max()

 connections = np.zeros((nunits,maxcox), dtype = np.int32);
 w_cids      = np.zeros((nunits,maxcox));
 dx_cids     = np.zeros((nunits,maxcox));

 for i in range(nunits):
  tmp = cols[rows==i]
  if tmp.size != nunits-1:
   tmp = np.append(tmp, np.zeros(maxcox-tmp.size)+i);
  connections[i,:] = tmp[:];
  w_cids[i,:] = w_gws[i,connections[i,:]]; #converting all the filled values with widht of the connection
  dx_cids[i,:] = dx_gws[i,connections[i,:]]; #converting all the filled values with widht of the connection

 #connections = gw_ids[connections]; #use the real ids for the groundwater units

 pickle.dump(w_cids,  open('%s/groundwater/w_cids.pck' % input_dir, 'wb')); #save files
 pickle.dump(dx_cids, open('%s/groundwater/dx_cids.pck' % input_dir, 'wb'));
 pickle.dump(connections,open('%s/groundwater/conx_cids.pck' % input_dir, 'wb'));

 #print('finished for cid: %s' % str(rank+1), flush = True)
 
 return

def multiscale_subsurface_preprocessing(comm,edir,rdir,hydroblocks_info,cids):
 '''Main function to generate the files for multiscale scheme for subsurface flow'''
 rank = comm.Get_rank()
 size = comm.Get_size()
 workspace = '%s/workspace' % edir #name of directory
 for cid in cids[rank::size]:
  print(f'mssubsurface from {cids[rank::size]} risfu:{rank} {cid} {len(cids)}',flush=True)
  #(1) Create the maps of regional/intermediate areas
  generate_regional_and_intermediate_units_maps(cid,edir,rdir,hydroblocks_info)
 comm.Barrier()
 #(2) Create virtual raster 
 file = '%s/regunits.vrt' % workspace
 if rank == 0:
  if os.path.exists(file) == False:
   create_vrt_risfu(workspace,edir,len(cids))
   print(f"Creating file: {file}",flush=True)
 comm.Barrier()
 while not os.path.exists(file):
  print(f"Waiting for file: {file}",flush=True)
  time.sleep(5)  # Check every 5 seconds
 mask_risfu_list = {}
 resx_list = {}
 input_dir_list = {}
 for cid in cids[rank::size]:
  print(f'mssubsurface connections:{cids[rank::size]}',flush=True)
  #(2) Create mask regional and intermediate units
  mask_risfu = create_risfu_mask(cid,edir,rdir,workspace)
  #(2) Generate files of aggregated soil hydraulic properties
  (risfu_org,gwunits_id,ngwus,resx,input_dir,gwunits_count) = generate_risfu_aggregated_properties(cid,edir,rdir,workspace,hydroblocks_info)
  mask_risfu_list[cid] = mask_risfu
  resx_list[cid] = resx
  input_dir_list[cid] = input_dir
  #(3) Generate matrix of connections for intermediate subsurface flow
  compute_connection_matrix_risfu(risfu_org,gwunits_id,ngwus,resx,input_dir,gwunits_count)
 comm.Barrier()
 #for cid in range(start_idx, end_idx):
 for cid in cids[rank::size]:
  #(4) Generate matrix of connections for regional subsurface flow
  create_files_regional_interaction(cid,mask_risfu_list[cid],resx_list[cid],input_dir_list[cid],edir)
 print(f'mssubsurface completed:{rank}',flush=True)
 comm.Barrier()
 return