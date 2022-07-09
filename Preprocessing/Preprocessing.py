import warnings
warnings.filterwarnings('ignore')
import sys
sys.path.append('Tools')
import pickle
import datetime
import numpy as np
import scipy.sparse as sparse
import scipy.stats as stats
from cv2 import dilate, getStructuringElement, MORPH_RECT, filter2D
#import model_tools as mt
import os
import netCDF4 as nc
import time
import glob
from geospatialtools import gdal_tools
#from geospatialtools import terrain_tools
import terrain_tools as terrain_tools
import gc
from scipy.interpolate import griddata
#import numba

dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append('%s/../HydroBlocks/pyHWU/' % dir )
#import management_funcs as mgmt_funcs

def plot_data(data, *arg_mask):

 if arg_mask:
   mask = arg_mask[0]
   data = np.ma.masked_array(data,mask==0)

 import matplotlib.pyplot as plt
 data = np.ma.masked_array(data,data==-9999)
 plt.figure(figsize=(10,10))
 plt.imshow(data)
 plt.colorbar()
 #plt.savefig('tmp.png')
 plt.show()

 return


def Prepare_Model_Input_Data(hydroblocks_info):

 #Prepare the info dictionary
 info = {}

 #Define the start/end dates
 info['time_info'] = {}
 info['time_info']['startdate'] = hydroblocks_info['idate']
 info['time_info']['enddate'] = hydroblocks_info['fdate']
 info['time_info']['dt'] = hydroblocks_info['dt']

 #Define the workspace
 workspace = hydroblocks_info['workspace']

 #Define the model input data directory
 input_dir = workspace

 #Create the dictionary to hold all of the data
 output = {}

 #Create the Latin Hypercube (Clustering)
 nhru = hydroblocks_info['nhru']
 icatch = hydroblocks_info['icatch']

 #Get metadata
 md = gdal_tools.retrieve_metadata('%s/mask_latlon.tif' % workspace)
 
 #Prepare the input file
 wbd = {}
 wbd['bbox'] = {'minlat':md['miny'],'maxlat':md['maxy'],
                'minlon':md['minx'],'maxlon':md['maxx'],
                'res':abs(md['resx'])}

 wbd['files'] = {
  'WLTSMC':'%s/theta1500_ea.tif' % workspace,
  'TEXTURE_CLASS':'%s/texture_class_ea.tif' % workspace,
  #'cslope':'%s/cslope_ea.tif' % workspace,
  'MAXSMC':'%s/thetas_ea.tif' % workspace,
  'BB':'%s/bb_ea.tif' % workspace,
  'DRYSMC':'%s/thetar_ea.tif' % workspace,
  #'fdir':'%s/fdir_ea.tif' % workspace,
  'QTZ':'%s/qtz_ea.tif' % workspace,
  'SATDW':'%s/dsat_ea.tif' % workspace,
  'REFSMC':'%s/theta33_ea.tif' % workspace,
  'mask':'%s/mask_ea.tif' % workspace,
  'SATPSI':'%s/psisat_ea.tif' % workspace,
  'lc':'%s/lc_ea.tif' % workspace,
  #'carea':'%s/carea_ea.tif' % workspace,
  #'ti':'%s/ti_ea.tif' % workspace,
  'ndvi':'%s/ndvi_ea.tif' % workspace,
  'F11':'%s/f11_ea.tif' % workspace,
  'SATDK':'%s/ksat_ea.tif' % workspace,
  'dem':'%s/dem_ea.tif' % workspace,
  #'demns':'%s/demns_ea.tif' % workspace,
  'sand':'%s/sand_ea.tif' % workspace,
  'clay':'%s/clay_ea.tif' % workspace,
  'silt':'%s/silt_ea.tif' % workspace,
  'om':'%s/om_ea.tif' % workspace,
  'bare30':'%s/bare30_ea.tif' % workspace,
  'water30':'%s/water30_ea.tif' % workspace,
  'tree30':'%s/tree30_ea.tif' % workspace,
  'irrig_land':'%s/irrig_land_ea.tif' % workspace,
  'dbedrock':'%s/dbedrock_ea.tif' % workspace
  }
 if hydroblocks_info['water_management']['hwu_agric_flag']:
   wbd['files']['irrig_land'] = '%s/irrig_land_ea.tif' % workspace
   wbd['files']['start_growing_season'] = '%s/start_growing_season_ea.tif' % workspace
   wbd['files']['end_growing_season']   = '%s/end_growing_season_ea.tif' % workspace

 wbd['files_meteorology'] = {
  'lwdown':'%s/lwdown.nc' % workspace,
  'swdown':'%s/swdown.nc' % workspace,
  'tair':'%s/tair.nc' % workspace,
  'precip':'%s/precip.nc' % workspace,
  'psurf':'%s/psurf.nc' % workspace,
  'wind':'%s/wind.nc' % workspace,
  'spfh':'%s/spfh.nc' % workspace,
  }

 if hydroblocks_info['water_management']['hwu_flag'] == True:
  wbd['files_water_use'] = {}
  if hydroblocks_info['water_management']['hwu_domest_flag']:
   wbd['files_water_use']['domestic']   = '%s/domestic.nc' % workspace
  if hydroblocks_info['water_management']['hwu_indust_flag']:
   wbd['files_water_use']['industrial'] = '%s/industrial.nc' % workspace
  if hydroblocks_info['water_management']['hwu_lstock_flag']:
   wbd['files_water_use']['livestock']  = '%s/livestock.nc' % workspace

 #Create the clusters and their connections
 output, hydroblocks_info = Create_Clusters_And_Connections(workspace,wbd,output,input_dir,nhru,info,hydroblocks_info)
 #import matplotlib.pyplot as plt
 #plt.imshow(output['hru_map']); plt.show()
 
 #Extract the meteorological forcing
 print("Preparing the meteorology", flush=True)
 Prepare_Meteorology_Semidistributed(workspace,wbd,output,input_dir,info,hydroblocks_info)

 #Extract the water use demands
 print("Preparing the water use", flush=True)
 if hydroblocks_info['water_management']['hwu_flag'] == True:
  Prepare_Water_Use_Semidistributed(workspace,wbd,output,input_dir,info,hydroblocks_info) 

 #Write out the files to the netcdf file
 fp = hydroblocks_info['input_fp']
 data = output

 #Write out the metadata
 grp = fp.createGroup('metadata')
 grp.latitude = (wbd['bbox']['minlat'] + wbd['bbox']['maxlat'])/2
 lon = (wbd['bbox']['minlon'] + wbd['bbox']['maxlon'])/2 
 if lon < 0:lon += 360
 grp.longitude = lon
 metadata = gdal_tools.retrieve_metadata(wbd['files']['mask']) 
 grp.dx = metadata['resx']

 #Write the HRU mapping
 #CONUS conus_albers metadata
 metadata['nodata'] = -9999.0
 #Save the conus_albers metadata
 grp = fp.createGroup('conus_albers_mapping')
 grp.createDimension('nx',metadata['nx'])
 grp.createDimension('ny',metadata['ny'])
 hmca = grp.createVariable('hmca','f4',('ny','nx'))#,zlib=True) 
 hmca.gt = metadata['gt']
 hmca.projection = metadata['projection']
 hmca.description = 'HRU mapping (conus albers)'
 hmca.nodata = metadata['nodata']
 #Save the conus albers mapping
 hru_map = np.copy(output['hru_map'])
 hru_map[np.isnan(hru_map) == 1] = metadata['nodata']
 hmca[:] = hru_map

 #Write out the mapping
 file_ca = '%s/hru_mapping_ea.tif' % workspace
 gdal_tools.write_raster(file_ca,metadata,hru_map)
 nhru_ea = np.unique(hru_map[ hru_map != -9999]).size

 #Map the mapping to regular lat/lon
 file_ll = '%s/hru_mapping_latlon.tif' % workspace
 os.system('rm -f %s' % file_ll)
 res = wbd['bbox']['res']
 minlat = wbd['bbox']['minlat']
 minlon = wbd['bbox']['minlon']
 maxlat = wbd['bbox']['maxlat']
 maxlon = wbd['bbox']['maxlon']
 log = '%s/log.txt' % workspace
 lon_0 = (maxlon+minlon)/2.0
 lat_0 = (maxlat+minlat)/2.0
 #os.system('gdalwarp -tr %f %f -dstnodata %f -r mode -t_srs \'+proj=longlat +lon_0=%f +lat_0=%f +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +no_defs \' -te %.16f %.16f %.16f %.16f %s %s >> %s 2>&1' % (res,res,metadata['nodata'],lon_0,lat_0,minlon,minlat,maxlon,maxlat,file_ca,file_ll,log))
 os.system('gdalwarp -overwrite -tr %.16f %.16f -dstnodata %f -t_srs \'+proj=longlat \' -te %.16f %.16f %.16f %.16f %s %s >> %s 2>&1' % (res,res,metadata['nodata'],lon_0,lat_0,minlon,minlat,maxlon,maxlat,file_ca,file_ll,log)) # Noemi 
 tmp = gdal_tools.read_raster(file_ll)
 nhru_latlon = np.unique(tmp[tmp != -9999]).size
 if nhru_ea != nhru_latlon: exit('Catch: %s - nhru in hru_mapping_ea.tif (%s) and hru_mapping_latlon.tif (%s) do not match -- check resampling' % (str(icatch),nhru_ea,nhru_latlon))
 
 #Write a map for the catchment id
 file_icatch = '%s/icatch_latlon.tif' % workspace
 metadata = gdal_tools.retrieve_metadata(file_ll)
 metadata['nodata'] = -9999.0
 tmp = gdal_tools.read_raster(file_ll)
 tmp[tmp >= 0] = hydroblocks_info['icatch']
 gdal_tools.write_raster(file_icatch,metadata,tmp)

 #Write the connection matrices
 #width
 wmatrix = output['cmatrix']['width']
 nconnections = wmatrix.data.size
 grp = fp.createGroup('wmatrix')
 grp.createDimension('connections_columns',wmatrix.indices.size)
 grp.createDimension('connections_rows',wmatrix.indptr.size)
 grp.createVariable('data','f4',('connections_columns',))
 grp.createVariable('indices','f4',('connections_columns',))
 grp.createVariable('indptr','f4',('connections_rows',))
 grp.variables['data'][:] = wmatrix.data
 grp.variables['indices'][:] = wmatrix.indices
 grp.variables['indptr'][:] = wmatrix.indptr

 
 #Write the model parameters
 grp = fp.createGroup('parameters')
 vars = ['slope','area_pct','land_cover','channel',
        'dem','soil_texture_class','ti','carea','area',
        'BB','F11','SATPSI','SATDW','QTZ','clay',
        'WLTSMC','MAXSMC','DRYSMC','REFSMC','SATDK',
        'm','hand']

 if hydroblocks_info['water_management']['hwu_agric_flag']:
  for var in ['centroid_lats', 'centroid_lons', 'irrig_land', 'start_growing_season', 'end_growing_season']:
    vars.append(var)

 for var in vars:
  grp.createVariable(var,'f4',('hru',))#,zlib=True)
  grp.variables[var][:] = data['hru'][var]

 if hydroblocks_info['water_management']['hwu_flag']:
  grp.createVariable('hru_min_dist','f4',('hru','hru'))#,zlib=True)
  grp.variables['hru_min_dist'][:] = data['hru']['hru_min_dist']

 #Remove info from output
 del output['hru']

 #Add in the catchment info
 output['wbd'] = wbd

 #Close the file
 fp.close()

 return output

def Compute_HRUs_Semidistributed_HMC(covariates,mask,hydroblocks_info,wbd,eares):

 #PARAMETERS (NEED TO GO OUTSIDE)
 #eares = 30 #meters
 
 print("catchment area: %.2f km2" % ((np.sum(mask)*eares*eares)/(10**6)),flush=True)

 #Define the parameters for the hierarchical multivariate clustering
 ncatchments = hydroblocks_info['hmc_parameters']['number_of_characteristic_subbasins']
 dh = hydroblocks_info['hmc_parameters']['average_height_difference_between_bands']
 nclusters = hydroblocks_info['hmc_parameters']['number_of_intraband_clusters']

 #Pre-process DEM
 dem = covariates['dem']
 #Remove pits in dem
 print("Removing pits in dem", flush=True)
 
 # normalize DEM -- Noemi
 invalid_dem = dem == -9999
 min_dem = np.min(dem[~invalid_dem])
 dem = dem - min_dem
 dem[invalid_dem] = -9999
 delta_dem = np.max(dem)

 demns = terrain_tools.ttf.remove_pits_planchon(dem,eares)
 dem = spatial_imputation(dem,-9999.0,'nearest',mask)
 demns = spatial_imputation(demns,-9999.0,'nearest',mask)
 covariates['dem'] = dem 
 covariates['demns'] = demns

 #Calculate slope and aspect
 print("Calculating slope and aspect", flush=True)
 res_array = np.copy(demns)
 res_array[:] = eares
 (slope,aspect) = terrain_tools.ttf.calculate_slope_and_aspect(demns,res_array,res_array)
 
 #Compute accumulated area
 m2 = np.copy(demns)
 m2[:] = 1
 
 print("Calculating accumulated area", flush=True)
 (area,fdir) = terrain_tools.ttf.calculate_d8_acc(demns,m2,eares)

 #Calculate channel initiation points (2 parameters) -- units? 
 C = area/eares*slope**2
 basin_area_threshold = 0.5*(10**6) # 0.1km2,  1km2 = 10**6 m2
 ipoints = ((C > 200) & (area > basin_area_threshold)).astype(np.int32)
 ipoints[ipoints == 0] = -9999

 #Create area for channel delineation
 (ac,fdc) = terrain_tools.ttf.calculate_d8_acc_wipoints(demns,m2,ipoints,eares)
 ac[ac != 0] = area[ac != 0]

 #Compute the channels
 print("Defining channels", flush=True)
 channels = terrain_tools.ttf.calculate_channels_wocean(ac,basin_area_threshold,basin_area_threshold,fdc,m2)
 channels = np.ma.masked_array(channels,channels<=0)
 #plot_data(channels,mask)
 #plot_data(ac)
 
 #If the dem is undefined then set to undefined
 channels[dem == -9999] = -9999
 
 #Compute the basins
 print("Defining basins", flush=True)
 basins = terrain_tools.ttf.delineate_basins(channels,m2,fdir)
 # Remove channel artifacts from basin delineation
 basins[(channels > 0) & (mask == 1)] = -9999
 #plot_data(basins,mask)
 basins = spatial_imputation(basins.astype(float),-9999.0,'nearest',mask).astype(np.int32)
 #plot_data(basins)

 # Create coarser sub-basins -- Noemi. This should be removed or solved with new routing implementation
 coarse_basins_threshold = 25*(10**6) # 25 km2
 c_ipoints = ((C > 200) & (area > coarse_basins_threshold)).astype(np.int32)
 c_ipoints[c_ipoints == 0] = -9999
 (c_ac,c_fdc) = terrain_tools.ttf.calculate_d8_acc_wipoints(demns,m2,c_ipoints,eares)
 c_ac[c_ac != 0] = area[c_ac != 0]
 c_channels = terrain_tools.ttf.calculate_channels_wocean(c_ac,coarse_basins_threshold,coarse_basins_threshold,c_fdc,m2)
 c_channels = np.ma.masked_array(c_channels,c_channels<=0)
 c_channels[dem == -9999] = -9999
 c_basins = terrain_tools.ttf.delineate_basins(c_channels,m2,fdir)
 c_basins[(c_channels > 0) & (mask == 1)] = -9999
 c_basins = spatial_imputation(c_basins.astype(float),-9999.0,'nearest',mask).astype(np.int32)
 #plot_data(c_basins,mask)
 
 '''
 # merge tiny c_basins clusters -- this will be unecessary after using dem to delineate catchments
 tmp = np.copy(c_basins)
 tmp[mask == 0] = -9999
 ubcs = np.unique(tmp)
 ubcs = ubcs[ubcs!=-9999]
 for ubc in ubcs:
  m = tmp == ubc
  ubc_size = float(eares*eares*np.sum(m))/float(10**6) # km2
  print("sub-basins %i area %.2f km2" % (ubc,ubc_size),flush=True)
  if ubc_size < 0.1 * (coarse_basins_threshold/float(10**6)): # if smaller than 10% of the coarse_basins_threshold 
    c_basins[m] = -9999
 c_basins = spatial_imputation(c_basins.astype(float),-9999.0,'nearest',mask).astype(np.int32)
 '''

 # Identify lakes and wetlands
 lakes = cluster_lakes(covariates,mask)
 wetlands = cluster_wetlands(covariates,mask)

 # set channels with lakes and wetlands
 channels_w_lakes = np.copy(channels)
 for lake in np.unique(lakes[lakes!=-9999]):
   m = (lakes == lake)
   channels_w_lakes[m] = np.max(channels)+1
 for wetland in np.unique(wetlands[wetlands!=-9999]):
   m = (wetlands == wetland)
   channels_w_lakes[m] = np.max(channels)+1
 #plot_data(channels_w_lakes,mask)

 #Calculate the height above nearest drainage area
 print("Computing height above nearest drainage area", flush=True)
 #hand = terrain_tools.ttf.calculate_depth2channel(channels,basins,fdir,demns) # dem should be normalized
 hand = terrain_tools.ttf.calculate_depth2channel(channels_w_lakes,basins,fdir,demns) # dem should be normalized
 hand[ hand < 0 ] = -9999.0
 hand[ channels > 0 ] = 0  # Noemi
 hand[ channels_w_lakes > 0 ] = 0  # Noemi
 hand[ lakes > 0 ] = 0.0
 hand[ wetlands > 0 ] = 0.0
 hand[ covariates['lc'] == 17 ] = 0.0
 hand[ covariates['lc'] == 11 ] = 0.0
 hand = spatial_imputation(hand,-9999.0,'nearest',mask)
 #plot_data(hand,mask)

 # Calculate coarse hand
 #c_hand = terrain_tools.ttf.calculate_depth2channel(c_channels,c_basins,fdir,demns) # dem should be normalized
 #c_hand[ c_hand < 0 ] = -9999.0
 #c_hand[ c_channels > 0 ] = 0.0
 #c_hand = spatial_imputation(c_hand,-9999.0,'nearest',mask)
 #plot_data(c_hand,mask)
 
 # calculate smooth coarse hand
 #smooth_c_hand = np.copy(c_hand)
 #filter_size = int(np.around(250./eares,0)) # 1000m
 #kernel = np.ones((filter_size,filter_size), np.float32)/float(filter_size*filter_size)
 #smooth_c_hand = filter2D(smooth_c_hand,-1,kernel)
 #smooth_c_hand[c_channels > 0 ] = 0.0
 #c_hand = np.copy(smooth_c_hand)
 #plot_data(smooth_c_hand,mask)

 
 # Identify large flat lands
 flatlands = cluster_flatlands(slope,mask)  
 # set hand at the flat lands
 for flat in np.unique(flatlands[flatlands!=-9999]):
   m = (flatlands == flat)
   tmp = hand[m][hand[m]!=-9999]
   if len(tmp) > 0: hand[m] = np.min(tmp)

 # Set hand, dem, demns at the lakes
 for lake in np.unique(lakes[lakes!=-9999]):
   tmp = dem[m][dem[m]!=-9999]                   
   if len(tmp) > 0: dem[m] = np.min(tmp)
   tmp = demns[m][demns[m]!=-9999]      
   if len(tmp) > 0: demns[m] = np.min(tmp)
   #tmp = basins[m][basins[m]!=-9999]
   #if len(tmp) > 0: basins[m] = np.max(basins)+1
   #tmp = c_basins[m][c_basins[m]!=-9999]
   #if len(tmp) > 0: c_basins[m] = np.max(c_basins)+1

 #Calculate topographic index
 print("Computing topographic index")
 ti = np.copy(area)
 m = (area != -9999) & (slope != -9999) & (slope != 0.0)
 ti[m] = np.log(area[m]/eares/slope[m])
 ti[slope == 0] = 15.0
 ti[slope < 0.0001] = 15.0

 # cleanup
 slope[mask != 1] = -9999
 aspect[mask != 1] = -9999
 area[mask != 1] = -9999
 channels[mask != 1] = -9999
 channels = np.ma.masked_array(channels,channels<=0)
 basins[mask != 1] = -9999
 c_basins[mask != 1] = -9999
 ti[mask != 1] = -9999
 #plot_data(mask)

 covariates['dem'] = dem
 covariates['demns'] = demns
 covariates['slope'] = slope
 covariates['aspect'] = aspect
 covariates['carea'] = area
 #covariates['fdir'] = fdir
 covariates['ti'] = ti


 '''
 #Calculate the subbasin properties
 print("Assembling the subbasin properties")
 hp_in = terrain_tools.calculate_basin_properties_updated(basins,eares,covariates,hydroblocks_info['hmc_parameters']['subbasin_clustering_covariates'])


 #Clustering the basins
 print("Clustering the basins")
 #Assemble input data
 cvs = {}
 for var in hydroblocks_info['hmc_parameters']['subbasin_clustering_covariates']:
  tmp = np.copy(hp_in[var])
  cvs[var] = {'min':np.min(hp_in[var][hp_in[var]!=-9999]),
              'max':np.max(hp_in[var][hp_in[var]!=-9999]),
              't':-9999,
              'd':tmp}

 #(basin_clusters,) = terrain_tools.cluster_basins_updated(basins,cvs,hp_in,ncatchments)
 '''
 basin_clusters = renumber(c_basins).astype(np.int32)
 print('number of sub-basins clusters: %i' % (np.unique(basin_clusters[basin_clusters != -9999]).size), flush=True)
 #plot_data(basin_clusters,mask)

 '''
 # print each sub-basin area
 ubcs = np.unique(basin_clusters[basin_clusters!=-9999])
 sub_basins_areas =[]
 for ubc in ubcs:
   ubc_size = float(eares*eares*np.sum( (basin_clusters == ubc) & (mask ==1) ))/float(10**6) # km2
   sub_basins_areas.append(ubc_size)
   print("sub-basins %i area %.2f km2" % (ubc,ubc_size),flush=True)
 print("sub-basins mean area %.2f km2" % (np.mean(sub_basins_areas)),flush=True)
 ''' 
 #plot_data(basin_clusters,mask)
 
 # remove tiny basins clusters
 ubcs = np.unique(basin_clusters)
 ubcs = ubcs[ubcs!=-9999]
 for ubc in ubcs:
  m = basin_clusters == ubc
  valid_vals = (hand[m] != -9999) 
  #if the basin cluster has less than 50% of valid hand values
  if np.sum(valid_vals) > 0:
    missing_ratio = 1.0 - (np.sum(valid_vals)/float(np.sum(m)))  
  else:
    missing_ratio = 1.0
  if missing_ratio > 0.95 : 
    basin_clusters[m] = -9999
    basins[m] = -9999
 if len(basin_clusters[basin_clusters != -9999]) < 1 : 
  #plot_data(basin_clusters)
  exit('Error_basin_clustering: hand_full_of_nans %s' % hydroblocks_info['icatch']) 
 

 #Divide each subbasin into height bands
 (tiles,new_hand,tile_position) = terrain_tools.create_basin_tiles(basin_clusters,hand,demns,basins,dh)   
 #plot_data(tiles,mask)
 #plot_data(new_hand,mask) 

 # De-normalize DEM
 dem = dem + min_dem
 demns = demns + min_dem
 covariates['dem'] = dem
 covariates['demns'] = demns
 
 #Disagregate land cover - create binary covariate for each land cover type
 intraband_clust_vars = hydroblocks_info['hmc_parameters']['intraband_clustering_covariates']
 if 'lc' in intraband_clust_vars: 
  intraband_clust_vars.remove('lc')
  disag = [i for i in covariates.keys() if 'lc_' in i]
  intraband_clust_vars = intraband_clust_vars + disag

 # where slope is flat, set aspect -4 to avoid noise at the flat lands
 flats = ( covariates['slope'] < 10./100. ) & ( covariates['slope'] != -9999 )
 covariates['aspect'][flats] = -3.14*5

 #Calculate the hrus (kmeans on each tile of each basin)
 cvs = {}
 for var in intraband_clust_vars:
  cvs[var] = {'min':np.min(covariates[var][covariates[var]!=-9999]),
              'max':np.max(covariates[var][covariates[var]!=-9999]),
              't':-9999,
              'd':covariates[var]}
 
 #print("Clustering the height bands into %d clusters" % nclusters)
 hrus = terrain_tools.create_hrus_hydroblocks(basin_clusters,tiles,channels,cvs,nclusters)
 hrus[hrus!=-9999] = hrus[hrus!=-9999] - 1
 nhru = np.unique(hrus[hrus!=-9999]).size

 # Add lakes as independent HRUs
 count = np.max(hrus)+1
 for lake in np.unique(lakes[lakes!=-9999]):
   m = lakes == lake 
   hrus[m] = count
   count = count+1
   new_hand[m] = 0.0
 #for wetland in np.unique(wetlands[wetlands!=-9999]):
 #  m = wetlands == wetland
 #  hrus[m] = count
 #  count = count+1
 #  new_hand[m] = 0.0
 unic = np.unique(hrus[hrus!=-9999])                        
 new_hrus = np.ones(hrus.shape)*(-9999)                         
 for nclust, oclust in enumerate(unic, 0):                      
   new_hrus[hrus == oclust] = nclust
 hrus = np.copy(new_hrus)
 nhru = np.unique(hrus[hrus!=-9999]).size
 #plot_data(hrus,mask)

 ngrids = np.sum(mask==1)
 print('Statistics - grids:%i  HRUs:%i  grids_to_hrus_ratio:%i  delta-dem:%im ' % (ngrids,nhru,float(ngrids)/float(nhru),delta_dem), flush=True) 

 # add a buffer at the channels and lakes in tile position, this helps channels connectivity at the valleys
 tile_position = add_buffer_to_channels(tile_position,channels,new_hand,mask)
 tile_position = add_buffer_to_water_and_wetlands(tile_position,covariates,mask)
 #plot_data(tile_position)


 #Construct HMC info for creating connections matrix
 HMC_info = {}
 #HMC_info['basins'] = basins
 HMC_info['basins'] = c_basins # Noemi -- use coarse scale sub-basins or set ridge connection to true
 HMC_info['tile_position'] = tile_position


 return (hrus.astype(np.float32),nhru,new_hand,HMC_info,covariates)

def add_buffer_to_water_and_wetlands(tile_position,covariates,mask):
  binary = np.zeros(mask.shape,dtype=np.uint8)
  binary[ ((covariates['lc'] == 17) | (covariates['lc'] == 11)) & (mask ==1) ] = 1

  diskv = getStructuringElement(MORPH_RECT,(3,3))
  buff = dilate(binary.astype(np.uint8),diskv)

  tile_position[(buff == 1) & (mask == 1)] = 0
  return tile_position
 

def add_buffer_to_channels(tile_position,channels,hand,mask,kernel_size=3):
  binary = np.zeros(channels.shape,dtype=np.uint8)
  binary[((channels > 0) | (hand == 0.0) ) & (mask == 1)] = 1
 
  diskv = getStructuringElement(MORPH_RECT,(kernel_size,kernel_size))
  buff = dilate(binary.astype(np.uint8),diskv)
  
  tile_position[(buff == 1) & (mask == 1)] = 0
  return tile_position


def renumber(data):
  new = np.zeros(data.shape,dtype=int,order='F')
  new[:] = -9999
  ub = np.unique(data)
  if -9999 in ub: ub = ub[1:]
  fill = np.arange(len(ub))
  for i,b in zip(fill,ub):
    new[data == b] = i
  return new

def split_subbasins(basins,mask):
  clust_map = np.copy(basins)
  clust_map[mask==1] = -9999

  basins[mask==0] = -9999
  ub = np.unique(basins)
  if -9999 in ub: ub = ub[1:]

  from sklearn.cluster import DBSCAN
  for b in ub:
    pos = np.where(basins == b)
    X = list(zip(*pos))
    clustering = DBSCAN(eps=1.0).fit(X)
    unic = np.unique(clustering.labels_)
    if -1 in unic: unic = unic[1:]
    if len(unic) > 0:
      increment = 1./float(len(unic))
      for i,c in enumerate(unic):
        m = clustering.labels_ == c
        ml,mc = pos[0][m],pos[1][m]
        clust_map[ml,mc] = b + float(i)*increment
  clust_map = renumber(clust_map)
  return clust_map

  
def cluster_lakes(covariates,mask):
 #Identify water bodies with area
 lc = covariates['lc']                           
 lc[mask != 1] = -9999 
 pos = np.where(lc == 17)
 clust_map = np.ones(mask.shape)*(-9999)

 # If there are at least 10 grids
 if len(pos[0]) < 10.: return clust_map 
         
 X = list(zip(*pos))              
 from sklearn.cluster import DBSCAN                                 
 clustering = DBSCAN(eps=1).fit(X)                

 unic = np.unique(clustering.labels_)                               
 unic = unic[unic>=0]            
 for c in unic:
  m = clustering.labels_ == c     
  if np.sum(m) > 20.:                            
   ml,mc = pos[0][m],pos[1][m]    
   clust_map[ml,mc] = c  
 return clust_map

def cluster_wetlands(covariates,mask):
 #Identify water bodies with area
 lc = covariates['lc']
 lc[mask != 1] = -9999
 pos = np.where(lc == 11)
 clust_map = np.ones(mask.shape)*(-9999)

 # If there are at least 10 grids
 if len(pos[0]) < 10.: return clust_map

 X = list(zip(*pos))
 from sklearn.cluster import DBSCAN
 clustering = DBSCAN(eps=1).fit(X)

 unic = np.unique(clustering.labels_)
 unic = unic[unic>=0]
 for c in unic:
  m = clustering.labels_ == c
  if np.sum(m) > 20.:
   ml,mc = pos[0][m],pos[1][m]
   clust_map[ml,mc] = c
 return clust_map

def cluster_flatlands(slope,mask):
 #Identify areas (> 0.1 km2) where the slope is too flat for channels 
 slope[mask != 1] = -9999
 pos = np.where(((slope >= 0) & (slope < 0.0005)))
 clust_map = np.ones(mask.shape)*(-9999)
 if len(pos[0]) < 20: return clust_map

 X = list(zip(*pos))
 from sklearn.cluster import DBSCAN
 clustering = DBSCAN(eps=1,n_jobs=-1).fit(X)

 unic = np.unique(clustering.labels_)
 unic = unic[unic>=0]
 #print(unic)
 for c in unic:
  m = clustering.labels_ == c
  if np.sum(m) > 100: # 1000 grids ~ 1km2
   ml,mc = pos[0][m],pos[1][m]
   clust_map[ml,mc] = c
 #print(np.unique(clust_map[clust_map!=-9999]))
 return clust_map



def Assign_Parameters_Semidistributed(covariates,metadata,hydroblocks_info,OUTPUT,cluster_ids,mask):

 nhru = hydroblocks_info['nhru']
 #Initialize the arrays
 vars = ['area','area_pct','BB','DRYSMC','F11','MAXSMC','REFSMC','SATPSI',
         'SATDK','SATDW','WLTSMC','QTZ','slope','ti','dem','carea','channel',
         'land_cover','soil_texture_class','clay','sand','silt',
         'm','hand']

 if hydroblocks_info['water_management']['hwu_agric_flag']:
  for var in ['centroid_lats', 'centroid_lons', 'irrig_land', 'start_growing_season', 'end_growing_season']:
    vars.append(var)

 OUTPUT['hru'] = {}
 if hydroblocks_info['water_management']['hwu_flag']: OUTPUT['hru']['hru_min_dist'] = np.zeros((nhru,nhru))
 for var in vars:
   OUTPUT['hru'][var] = np.zeros(nhru)


 #Metadata
 for hru in np.arange(nhru):
  #Set indices
  idx = np.where(cluster_ids == hru)

  # Check for spurious data
  for v in covariates.keys():                     
   vtmp = covariates[v][idx]
   if v[:3] == 'lc_': continue     
   if np.min(vtmp) == -9999: exit("Spurious values in %s hru: %i" % (v, hru))   

  #Calculate area per hru
  OUTPUT['hru']['area'][hru] = metadata['resx']**2*idx[0].size # units?
  #Calculate area percentage per hru
  OUTPUT['hru']['area_pct'][hru] = 100*OUTPUT['hru']['area'][hru]/(metadata['resx']**2*mask[mask].size)
  #Soil properties
  for var in ['BB','DRYSMC','F11','MAXSMC','REFSMC','SATPSI','SATDK','SATDW','WLTSMC','QTZ','clay','sand','silt']:
   if var in ['SATDK','SATDW']:
    #print('HMEAN:',var,covariates[var][idx][covariates[var][idx] <= 0 ])
    OUTPUT['hru'][var][hru] = stats.mstats.hmean(covariates[var][idx])/3600.0/1000.0 #mm/hr -> m/s
   else:
    OUTPUT['hru'][var][hru] = np.mean(covariates[var][idx])
  #Average Slope
  OUTPUT['hru']['slope'][hru] = np.nanmean(covariates['slope'][idx])
  #Topographic index
  OUTPUT['hru']['ti'][hru] = np.nanmean(covariates['ti'][idx])
  #HAND
  OUTPUT['hru']['hand'][hru] = np.nanmean(covariates['hand'][idx])
  if OUTPUT['hru']['hand'][hru] < 0.1: OUTPUT['hru']['hand'][hru] = 0.0
  #OUTPUT['hru']['hand'][hru] = np.nanmean(covariates['combined_hand_dem'][idx])
  #DEM
  OUTPUT['hru']['dem'][hru] = np.nanmean(covariates['demns'][idx])
  #if OUTPUT['hru']['hand'][hru] == 0.0: OUTPUT['hru']['dem'][hru] = np.nanmin(covariates['demns'][idx])

  #Average Catchment Area
  OUTPUT['hru']['carea'][hru] = np.nanmean(covariates['carea'][idx])
  #Channel?
  #OUTPUT['hru']['channel'][hru] = stats.mode(covariates['channels'][idx])[0]

  #Land cover type 
  tmp = covariates['lc'][idx]
  tmp = tmp[tmp>=1]
  if len(tmp) >= 1 :
   OUTPUT['hru']['land_cover'][hru] = stats.mode(tmp)[0][0]
  else:
   OUTPUT['hru']['land_cover'][hru] = 17  # if there is no valid value, set to water #Noemi
      
  #Soil texture class
  OUTPUT['hru']['soil_texture_class'][hru] = stats.mode(covariates['TEXTURE_CLASS'][idx])[0][0]

  #Define the estimate for the model parameters
  OUTPUT['hru']['m'][hru] = np.nanmean(covariates['dbedrock'][idx]) #0.1 #Form of the exponential decline in conductivity (0.01-1.0)
  
  #f11, qtz = get_f11_and_qtz_from_soil_type(OUTPUT['hru']['soil_texture_class'][hru])
  #if np.isnan(f11) == False: OUTPUT['hru']['F11'][hru] = f11
  #if np.isnan(qtz) == False: OUTPUT['hru']['QTZ'][hru] = qtz

  # Limit parameters to possible range
  #if OUTPUT['hru']['BB'][hru] > 11.55: OUTPUT['hru']['BB'][hru] = 11.55
  #if OUTPUT['hru']['BB'][hru] < 2.79: OUTPUT['hru']['BB'][hru] = 2.79
  #if OUTPUT['hru']['DRYSMC'][hru] > 0.138: OUTPUT['hru']['DRYSMC'][hru] = 0.138
  #if OUTPUT['hru']['SATPSI'][hru] > 0.617: OUTPUT['hru']['SATPSI'][hru] = 0.617
  #if OUTPUT['hru']['SATPSI'][hru] < 0.036: OUTPUT['hru']['SATPSI'][hru] = 0.036
  #if OUTPUT['hru']['SATDK'][hru] < 9.74e-7: OUTPUT['hru']['SATDK'][hru] = 9.74e-7
  #if OUTPUT['hru']['SATDK'][hru] > 4.66e-5: OUTPUT['hru']['SATDK'][hru] = 4.66e-5
  #if OUTPUT['hru']['SATDW'][hru] < 0.608e-6: OUTPUT['hru']['SATDW'][hru] = 0.608e-6
  #if OUTPUT['hru']['SATDW'][hru] > 0.239e-4: OUTPUT['hru']['SATDW'][hru] = 0.239e-4


  # Recalcualte the soil paramters based on the HRU mean values - Noemi
  #lamda_campbell = (np.log(OUTPUT['hru']['REFSMC'][hru])-np.log(OUTPUT['hru']['WLTSMC'][hru]))/(np.log(1500.)-np.log(33.))
  #psisat_campbell = 33*np.power(OUTPUT['hru']['REFSMC'][hru]/OUTPUT['hru']['MAXSMC'][hru], 1/lamda_campbell) #kPa
  #OUTPUT['hru']['BB'][hru] = 1./lamda_campbell
  #OUTPUT['hru']['SATPSI'][hru] = psisat_campbell*(10./100) # KPa --> cm --> m
  #OUTPUT['hru']['SATDW'][hru] = OUTPUT['hru']['BB'][hru] * OUTPUT['hru']['SATDK'][hru] * OUTPUT['hru']['SATPSI'][hru] / OUTPUT['hru']['MAXSMC'][hru] # m2/hr
 
                           
  # Define soil types over urban areas 
  if OUTPUT['hru']['land_cover'][hru] == 13:
    # if there is urban areas in the river channel, set to grassland
    #if OUTPUT['hru']['hand'][hru] == 0 :
    #  OUTPUT['hru']['land_cover'][hru] == 10
    #  #else, set soil to bedrock
    #else:
      OUTPUT['hru']['BB'][hru] = 2.79    
      OUTPUT['hru']['MAXSMC'][hru] = 0.20
      OUTPUT['hru']['REFSMC'][hru] = 0.10
      OUTPUT['hru']['WLTSMC'][hru] = 0.006 
      OUTPUT['hru']['DRYSMC'][hru] = 0.004
      OUTPUT['hru']['SATPSI'][hru] = 0.069 
      #OUTPUT['hru']['SATDW'][hru] = 0.136e-3  # bedrock
      #OUTPUT['hru']['SATDK'][hru] = 1.41e-4   # bedrock
      #OUTPUT['hru']['QTZ'][hru] = 0.07    
      #OUTPUT['hru']['soil_texture_class'][hru] = 15 
      


  # Water Management Variables
  if hydroblocks_info['water_management']['hwu_agric_flag'] == True:
   # Irrigation: 1 Irrigated, 2 paddy crop
   irrig_vec = np.copy(covariates['irrig_land'][idx].flatten())
   nii = np.nansum( irrig_vec == 1 )
   npi = np.nansum( irrig_vec == 2 )
   ttt = np.nansum( irrig_vec >= 0 )
   iratio = 0
   if ttt > 0: 
    iratio = float((nii+npi))/ttt
    if   iratio >= 0.50 and nii > npi: OUTPUT['hru']['irrig_land'][hru] = 1
    elif iratio >= 0.50 and nii < npi: OUTPUT['hru']['irrig_land'][hru] = 2
    else: OUTPUT['hru']['irrig_land'][hru] = 0
   else: OUTPUT['hru']['irrig_land'][hru] = 0
 
    # Crop Calendar
   OUTPUT['hru']['start_growing_season'][hru] = int(stats.mode(covariates['start_growing_season'][idx])[0][0])
   OUTPUT['hru']['end_growing_season'][hru] = int(stats.mode(covariates['end_growing_season'][idx])[0][0])
  
  if hydroblocks_info['water_management']['hwu_flag'] == True:
   #HRU Centroids for water management
   OUTPUT['hru']['centroid_lons'][hru] = np.nanmean(covariates['lons'][idx])
   OUTPUT['hru']['centroid_lats'][hru] = np.nanmean(covariates['lats'][idx])
   #print OUTPUT['hru']['centroid_lons']  

 #if hydroblocks_info['water_management']['hwu_flag'] == True: 
  #for hru in np.arange(nhru):
   ##HRU distance between the centroids of the hru and all the other hrus  -- need to fix skimage dependences
   #OUTPUT['hru']['hru_min_dist'][hru,:] = mgmt_funcs.calculate_min_distance(hru, nhru, cluster_ids, covariates['lats'], covariates['lons'], OUTPUT['hru']['centroid_lats'], OUTPUT['hru']['centroid_lons'])
   #OUTPUT['hru']['hru_min_dist'][hru,hru] = 0.0
   
 return OUTPUT

def Determine_HMC_Connectivity_old(h1,h2,b1,b2,tp1,tp2,hmc):

 if (h2 == -9999):return False
 if (h1 == h2):return True

 # Noemi
 #if ((tp1 == tp2) & (tp1 == 0) & (b1 != b2) & (hmc['intervalley_connectivity'] == True)):return True
 #if (b1 != b2) & (hmc['interridge_connectivity'] == False):return False
 #if (np.abs(tp1 - tp2) != 1) & (hmc['intraband_connectivity'] == False):return False
 
 if (b1 != b2) & ((tp1 == itp2) | (tp1 == 0)): return True
 if (b1 != b2) : return False

 return True


def Determine_HMC_Connectivity(h1,h2,b1,b2,tp1,tp2,hmc):

 ivc = hmc['intervalley_connectivity']
 irc = hmc['interridge_connectivity']
 ibc = hmc['intraband_connectivity']

 if (h2 == -9999): return False
 if (h1 == h2): return True
 if (b1 == b2): return True

 if (b1 != b2):
  
   if ((tp1 == tp2) & (tp1 == 0) & (ivc)): return True
   if (irc == False): return False
   if (irc == True): return True

 #if (np.abs(tp1 - tp2) == 1) & (ibc == False): return False

 print("Missing connectivity case:", h1,h2,b1,b2,tp1,tp2, flush=True)

 return True



def Calculate_HRU_Connections_Matrix_HMC(covariates,cluster_ids,nhru,dx,HMC_info,hydroblocks_info):

 #Add pointers for simplicity
 tile_position = HMC_info['tile_position']
 basins = HMC_info['basins'] 

 #plot_data(tile_position)
 #plot_data(cluster_ids)
 #plot_data(basins)
  
 #Define spatial resolution
 res = dx
 
 horg = []
 hdst = []
 #Count the connections
 for i in range(cluster_ids.shape[0]):
  for j in range(cluster_ids.shape[1]):
   h1 = cluster_ids[i,j]
   b1 = basins[i,j]
   tp1 = tile_position[i,j]
   if h1 == -9999: continue
   #up
   if (i+1) < cluster_ids.shape[0]:
    h2 = cluster_ids[i+1,j]
    b2 = basins[i+1,j]
    tp2 = tile_position[i+1,j]
    if Determine_HMC_Connectivity(h1,h2,b1,b2,tp1,tp2,hydroblocks_info['hmc_parameters']):
     horg.append(h1)
     hdst.append(h2)
   #down
   if (i-1) > 0:
    h2 = cluster_ids[i-1,j]
    b2 = basins[i-1,j]
    tp2 = tile_position[i-1,j]
    if Determine_HMC_Connectivity(h1,h2,b1,b2,tp1,tp2,hydroblocks_info['hmc_parameters']):
     horg.append(h1)
     hdst.append(h2)
   #left
   if (j-1) > 0:
    h2 = cluster_ids[i,j-1]
    b2 = basins[i,j-1]
    tp2 = tile_position[i,j-1]
    if Determine_HMC_Connectivity(h1,h2,b1,b2,tp1,tp2,hydroblocks_info['hmc_parameters']):
     horg.append(h1)
     hdst.append(h2)
   #right
   if (j+1) < cluster_ids.shape[1]:
    h2 = cluster_ids[i,j+1]
    b2 = basins[i,j+1]
    tp2 = tile_position[i,j+1]
    if Determine_HMC_Connectivity(h1,h2,b1,b2,tp1,tp2,hydroblocks_info['hmc_parameters']):
     horg.append(h1)
     hdst.append(h2)

   
   # Include diagonal connections. This needs to be revised: widths should be account for 8D or channel delineation should be done with 4D only.
   #up right
   if (i+1) < cluster_ids.shape[0] and (j+1) < cluster_ids.shape[1]:
    h2 = cluster_ids[i+1,j+1]
    b2 = basins[i+1,j+1]
    tp2 = tile_position[i+1,j+1]
    if Determine_HMC_Connectivity(h1,h2,b1,b2,tp1,tp2,hydroblocks_info['hmc_parameters']):
     horg.append(h1)
     hdst.append(h2)
   #down right 
   if (i-1) > 0 and (j+1) < cluster_ids.shape[1]:
    h2 = cluster_ids[i-1,j+1]
    b2 = basins[i-1,j+1]
    tp2 = tile_position[i-1,j+1]
    if Determine_HMC_Connectivity(h1,h2,b1,b2,tp1,tp2,hydroblocks_info['hmc_parameters']):
     horg.append(h1)
     hdst.append(h2)
   #up left
   if (i+1) < cluster_ids.shape[0] and (j-1) > 0:
    h2 = cluster_ids[i+1,j-1]
    b2 = basins[i+1,j-1]
    tp2 = tile_position[i+1,j-1]
    if Determine_HMC_Connectivity(h1,h2,b1,b2,tp1,tp2,hydroblocks_info['hmc_parameters']):
     horg.append(h1)
     hdst.append(h2)
   #down left
   if (i-1) > 0 and (j-1) > 0 :
    h2 = cluster_ids[i-1,j-1]
    b2 = basins[i-1,j-1]
    tp2 = tile_position[i-1,j-1]
    if Determine_HMC_Connectivity(h1,h2,b1,b2,tp1,tp2,hydroblocks_info['hmc_parameters']):
     horg.append(h1)
     hdst.append(h2)
  

 print("HRUs connections: %i orgin %i dest" % (len(horg),len(hdst)),flush=True)
 horg = np.array(horg)
 hdst = np.array(hdst)

 #Prepare the sparse matrix
 #cmatrix = sparse.coo_matrix((np.ones(hdst.size),(horg,hdst)),dtype=np.float32)
 cmatrix = sparse.coo_matrix((np.ones(hdst.size),(horg,hdst)),shape=(nhru,nhru),dtype=np.float32)
 
 #import matplotlib.pyplot as plt
 #from matplotlib.colors import LogNorm
 #plt.imshow(cmatrix.toarray()+1, norm=LogNorm())
 #plt.colorbar()
 #plt.show()

 cmatrix = cmatrix.tocsr()
 
 #Prepare length, width, and ksat matrices
 wmatrix = cmatrix.copy()
 wmatrix[:] = res*wmatrix[:]

 #Prepare output dictionary
 cdata = {'width':wmatrix.T,}

 return cdata

def Create_and_Curate_Covariates(wbd,hydroblocks_info):

 covariates = {}
 #Read in and curate all the covariates
 for file in wbd['files']:
  if os.path.isfile(wbd['files'][file]): 
   covariates[file] = gdal_tools.read_raster(wbd['files'][file])

 # check if lc is a covariates, and disagregate it in classes
 if 'lc' in hydroblocks_info['hmc_parameters']['intraband_clustering_covariates']:
  for lc in np.unique(covariates['lc'][covariates['mask'].astype(np.bool)]):
   if lc >= 0 :
    vnam = u'lc_%i' % lc
    masklc = (covariates['lc'] == lc)
    covariates[vnam] = np.zeros(covariates['lc'].shape)
    covariates[vnam][masklc] = 1.0
    hydroblocks_info['covariates'][vnam] = 'n'
    
 #Create lat/lon grids
 lats = np.linspace(wbd['bbox']['minlat']+wbd['bbox']['res']/2,wbd['bbox']['maxlat']-wbd['bbox']['res']/2,covariates['dem'].shape[0])
 lons = np.linspace(wbd['bbox']['minlon']+wbd['bbox']['res']/2,wbd['bbox']['maxlon']-wbd['bbox']['res']/2,covariates['dem'].shape[1])

 #Need to fix so that it doesn't suck up all the clustering:
 lats, lons = np.meshgrid(lats, lons)
 covariates['lats'] = lats.T
 covariates['lons'] = lons.T

 #Define the mask
 mask = np.copy(covariates['mask'])
 mask[mask >= 0] = 1
 mask[mask < 0] = 0
 mask = mask.astype(np.bool)
 
 # Identify missing values
 for var in covariates:
  mask0 = (np.isinf(covariates[var]) == 1) | (np.isnan(covariates[var]) == 1)
  covariates[var][mask0] = -9999.0

 # Spatial imputation on -9999 
 for var in covariates:

  m2 = ( mask == 1 ) & (covariates[var] != -9999.0)
  missing_ratio = 1.0 - np.sum(m2)/float(np.sum(mask))
  missing_ratio_threshold = 0.95
  if missing_ratio > missing_ratio_threshold : 
    print("Warning: Covariate %s in catchment %s has %.2f %% of nan's" % (var,hydroblocks_info['icatch'],100*missing_ratio), flush=True) # Noemi insert
    if var == 'lc': 
      mlc = (covariates[var] == -9999) & mask
      covariates[var][mlc] = 17  # Water
    if var == 'dbedrock':
      # fill no-data dbedrock values at water -- this should be improved
      covariates['dbedrock'][ (covariates['dbedrock'] == -9999) & (covariates['lc'] == 17) ] = 50
      m2 = ( mask == 1 ) & (covariates[var] != -9999.0)
      missing_ratio = 1.0 - np.sum(m2)/float(np.sum(mask))
    if missing_ratio > missing_ratio_threshold and var in ['dem','fdir','sand','clay','silt','dbedrock','TEXTURE_CLASS','SATDK']:
      exit('Error_clustering: %s_full_of_nans %s' % (var,hydroblocks_info['icatch']))

  if var not in ['mask']:
   if var in ['fdir','nlcd','TEXTURE_CLASS','lc','irrig_land','bare30','water30','tree30','start_growing_season','end_growing_season','SATDK','SATDW']:
    #print(var, np.sum( covariates[var] != -9999 ), flush=True)
    covariates[var] = spatial_imputation(covariates[var],-9999.0,'nearest',mask)
   elif var[:3] == 'lc_': pass
   elif var in ['dem']: pass  # do imputation afterwards
   else:
    #print(var, np.sum( covariates[var] != -9999 ), flush=True)
    covariates[var] = spatial_imputation(covariates[var],-9999.0,'nearest',mask) #faster
    #covariates[var] = spatial_imputation(covariates[var],-9999.0,'linear',mask)

 #Set everything outside of the mask to -9999
 for var in covariates:
  if var in hydroblocks_info['hmc_parameters']['subbasin_clustering_covariates']:continue
  covariates[var][mask <= 0] = -9999.0

 return (covariates,mask)

from scipy.stats import find_repeats
def faster_mode1D(a,fill_val):                             
    arr = np.asarray(a) # would be _chk_array
    if len(arr) < 1: return fill_val

    v, c = find_repeats(arr)
    #print(arr,v,c,flush=True) 
    if c.size == 0:                           
        arr.sort() # mimic first value behavior   
        return arr[0]#, 1.                            
    else:       
        pos = c.argmax()                        
        return v[pos]#, c[pos]

def spatial_imputation(array,fill_val,method,mask_catch):

 array = np.ma.filled(array,fill_value=fill_val)
 array[np.isnan(array)] = fill_val

 # If nothing to be filled return
 if np.sum(array==fill_val) == 0: return array  # if not  missing return

 # Fill outside the catchment with mode:
 pos_fill = (array==fill_val) & (~mask_catch)
 nmissing = np.sum(pos_fill)
 if nmissing > 0:
  fill_with = (array!=fill_val) & (~mask_catch)
  array[pos_fill] = faster_mode1D(array[array!=fill_val].ravel(),fill_val)

 # Fill inside the catchment 
 pos_fill = (array==fill_val) & mask_catch
 nmissing = np.sum(pos_fill)
 #print(nmissing,flush=True) 
 if nmissing > 0:

  # If there are way too many missing values inside the catchmentment 
  if nmissing > (30*30)*50: # If missing more than 50 km2, return with the mode
   fill_with = (array!=fill_val) & (mask_catch)
   array[pos_fill] = faster_mode1D(array[fill_with].ravel(),fill_val)
   del fill_with
   return array

  # Else use nearest for filling inside the catchment 
  pos_mask = np.zeros(array.shape,dtype=np.uint8)
  pos_mask[pos_fill] = 1
  if method == 'nearest': diskv = getStructuringElement(MORPH_RECT,(3,3))
  else: diskv = getStructuringElement(MORPH_RECT,(5,5))
  mask_good = dilate(pos_mask.astype(np.uint8),diskv)
  mask_good[pos_fill]=0
  mask_good[~mask_catch]=0
  del pos_mask; gc.collect()

  array[pos_fill]=np.nan
  del pos_fill; gc.collect()
  array2 = np.copy(array)
  array2[(mask_good == 0)] = np.nan
   
  x = np.arange(0, array.shape[1])
  y = np.arange(0, array.shape[0])
  array = np.ma.masked_invalid(array)
  array2 = np.ma.masked_invalid(array2)
  xx, yy = np.meshgrid(x, y)
  x1 = xx[~array2.mask] #train vals
  y1 = yy[~array2.mask] #train vals
  newarr = array[~array2.mask] #train vals
  xf = xx[array.mask] # fill vals
  yf = yy[array.mask] #fill vals
  #print(x1,y1,flush=True)
  #print(xf,yf,flush=True)

  # If there is data in the souroundings to train:
  if len(x1) > 0: 
    filled  = griddata((x1, y1), newarr.ravel(),(xf, yf), method=method)
    array[array.mask] = filled  # replace fill vals
  else: # fill with the mode
    fill_with = (array!=fill_val) & (mask_catch)
    array[array.mask] = faster_mode1D(array[fill_with].ravel(),fill_val)

  if method in ['linear','cubic']:
   # fill the missing boundaries with nearest
   array = np.ma.masked_invalid(array)
   xf = xx[array.mask] # fill vals
   yf = yy[array.mask] #fill vals
   filled  = griddata((x1, y1), newarr.ravel(),(xf, yf), method='nearest')
   array[array.mask] = filled

  array = np.ma.filled(array)
 return array
 
def Create_Clusters_And_Connections(workspace,wbd,output,input_dir,nhru,info,hydroblocks_info):
 
 #Retrieve some metadata
 metadata = gdal_tools.retrieve_metadata(wbd['files']['mask'])
 resx = metadata['resx']

 print("Creating and curating the covariates", flush=True)
 (covariates,mask) = Create_and_Curate_Covariates(wbd,hydroblocks_info)
 
 #Determine the HRUs (clustering if semidistributed; grid cell if fully distributed)
 print("Computing the HRUs", flush=True)
 (cluster_ids,nhru,hand,HMC_info,covariates) = Compute_HRUs_Semidistributed_HMC(covariates,mask,hydroblocks_info,wbd,resx)
 covariates['hand'] = hand
 hydroblocks_info['nhru'] = nhru
  
 #Create the netcdf file
 file_netcdf = hydroblocks_info['input_file']
 hydroblocks_info['input_fp'] = nc.Dataset(file_netcdf, 'w', format='NETCDF4')
 

 #Create the dimensions (netcdf)
 idate = hydroblocks_info['idate']
 fdate = hydroblocks_info['fdate']
 dt = hydroblocks_info['dt']
 ntime = 24*3600*((fdate - idate).days+1)/dt
 nhru = hydroblocks_info['nhru']
 hydroblocks_info['input_fp'].createDimension('hru',nhru)
 hydroblocks_info['input_fp'].createDimension('time',ntime)
 
 #Create the groups (netcdf)
 hydroblocks_info['input_fp'].createGroup('meteorology')
 hydroblocks_info['input_fp'].createGroup('water_use')

 #Prepare the hru connections matrix (darcy clusters)
 print("Calculating the connections between HRUs", flush=True)
 cmatrix = Calculate_HRU_Connections_Matrix_HMC(covariates,cluster_ids,nhru,resx,HMC_info,hydroblocks_info)

 #Define the metadata
 metadata = gdal_tools.retrieve_metadata(wbd['files']['dem'])

 #Make the output dictionary for the basin
 OUTPUT = {'hru':{},'metadata':metadata,'mask':mask,'cmatrix':cmatrix}

 #Remember the map of hrus
 OUTPUT['hru_map'] = cluster_ids

 #Assign the model parameters
 print("Assigning the model parameters", flush=True)
 OUTPUT = Assign_Parameters_Semidistributed(covariates,metadata,hydroblocks_info,OUTPUT,cluster_ids,mask)

 #Add the new number of clusters
 OUTPUT['nhru'] = nhru
 OUTPUT['mask'] = mask

 return OUTPUT, hydroblocks_info

def Prepare_Meteorology_Semidistributed(workspace,wbd,OUTPUT,input_dir,info,hydroblocks_info):

 #Define the mapping directory
 mapping_info = {}
 #Calculate the fine to coarse scale mapping
 for data_var in wbd['files_meteorology']:
  
  #Define the variable name
  var = data_var
  mapping_info[var] = {}

  #Read in the coarse and fine mapping
  file_coarse = '%s/%s_latlon_coarse.tif' % (workspace,data_var)
  file_fine = '%s/%s_ea_fine.tif' % (workspace,data_var)
  mask_coarse = gdal_tools.read_raster(file_coarse)
  mask_fine = gdal_tools.read_raster(file_fine)
  nlat = mask_coarse.shape[0]
  nlon = mask_coarse.shape[1]

  #Compute the mapping for each hru
  for hru in np.arange(hydroblocks_info['nhru']):
   idx = OUTPUT['hru_map'] == hru
   #print "Catch:",hydroblocks_info['icatch'], "HRU: ", mask_fine[idx].astype(np.int)
   icells = np.unique(mask_fine[idx][mask_fine[idx] != -9999.0].astype(np.int))   # Add != -9999 for unique and bicount - Noemi
   counts = np.bincount(mask_fine[idx][mask_fine[idx] != -9999.0].astype(np.int))
   coords,pcts = [],[]
   for icell in icells:
    ilat = int(np.floor(icell/mask_coarse.shape[1]))
    jlat = icell - ilat*mask_coarse.shape[1]
    #ilat = int(mask_coarse.shape[0] - ilat - 1) #CAREFUL
    pct = float(counts[icell])/float(np.sum(counts))
    coords.append([ilat,jlat])
    pcts.append(pct)
   pcts = np.array(pcts)
   coords = list(np.array(coords).T)
   if len(pcts) == 0: 
     print("Error - Preprocessing: no meteorological data at HRU:", hru, flush=True)
   mapping_info[var][hru] = {'pcts':pcts,'coords':coords}

 #Iterate through variable creating forcing product per HRU
 idate = info['time_info']['startdate']
 fdate = info['time_info']['enddate']
 dt = info['time_info']['dt']
 nt = int(3600*24/dt)*((fdate - idate).days+1)
 #Create structured array
 meteorology = {}
 for data_var in wbd['files_meteorology']:
  meteorology[data_var] = np.zeros((nt,hydroblocks_info['nhru']))
 #Load data into structured array
 for data_var in wbd['files_meteorology']:
  var = data_var#data_var.split('_')[1]
  date = idate
  file = wbd['files_meteorology'][data_var]
  fp = nc.Dataset(file)
  #fp = h5py.File(file)
  
  #Determine the time steps to retrieve
  #fidate = ' '.join(fp.variables['t'].units.split(' ')[2::])
  #dates = nc.num2date(fp.variables['t'][:],units='hours since %s' % (fidate))  # num2date doesnt work for step of 3h.
  nc_step = int(fp.variables['t'].units.split(' ')[0].split('h')[0])
  nc_idate = np.array(fp.variables['t'].units.split(' ')[2].split('-'))
  nc_nt = len(fp.variables['t'][:])
  dates = [datetime.datetime(int(nc_idate[0]),int(nc_idate[1]),int(nc_idate[2]))]
  for it in range(1,nc_nt): dates.append(dates[0] + datetime.timedelta(hours=it*nc_step))
  dates=np.array(dates)
  startdate = info['time_info']['startdate']
  enddate = info['time_info']['enddate']
  #print(startdate,enddate)
  mask_dates = (dates >= startdate) & (dates <= enddate)
  #print(mask_dates)
  #print(fp.variables[var][:,:,:].shape)
  data = np.ma.getdata(fp.variables[var][mask_dates,:,:])
  fp.close()

  #Assing to hrus
  for hru in mapping_info[var]:
   #print( hru, var, data.shape, hru,mapping_info[var][hru]['pcts'], mapping_info[var][hru]['coords'], flush=True)
   pcts = mapping_info[var][hru]['pcts']
   coords = mapping_info[var][hru]['coords']
   coords[0][coords[0] >= data.shape[1]] = data.shape[1] - 1
   coords[1][coords[1] >= data.shape[2]] = data.shape[2] - 1
   tmp = data[:,coords[0],coords[1]]
   #m1 = tmp < -999
   #m2 = tmp > -999
   tmp = pcts*tmp
   meteorology[data_var][:,hru] = np.sum(tmp,axis=1)
   #print data_var,np.unique(meteorology[data_var][:,:])

  #Write the meteorology to the netcdf file (single chunk for now...)
  grp = hydroblocks_info['input_fp'].groups['meteorology']
  grp.createVariable(var,'f4',('time','hru'))#,zlib=True)
  grp.variables[data_var][:] = meteorology[data_var][:]


 #Add time information
 dates = []
 date = idate
 while date <= fdate:
  dates.append(date)
  date = date + datetime.timedelta(seconds=dt)
 dates = np.array(dates)
 var = grp.createVariable('time','f8',('time',))
 var.units = 'hours since %4d-01-01' % idate.year
 var.calendar = 'standard'
 dates = nc.date2num(dates,units=var.units,calendar=var.calendar)
 var[:] = dates[:]

 return


def Prepare_Water_Use_Semidistributed(workspace,wbd,OUTPUT,input_dir,info,hydroblocks_info):

 #Define the mapping directory
 mapping_info = {}
 
 #Calculate the fine to coarse scale mapping
 for data_var in wbd['files_water_use']:

  #Define the variable name
  var = data_var#data_var.split('_')[1]
  mapping_info[var] = {}

  #Read in the coarse and fine mapping
  file_coarse = '%s/%s_latlon_coarse.tif' % (workspace,data_var)
  file_fine = '%s/%s_ea_fine.tif' % (workspace,data_var)
  #print(file_coarse)
  mask_coarse = gdal_tools.read_raster(file_coarse)
  mask_fine = gdal_tools.read_raster(file_fine)
  md = gdal_tools.retrieve_metadata(file_fine)
  md['nodata'] = -9999.0
  nlat = mask_coarse.shape[0]
  nlon = mask_coarse.shape[1]

  # NOAH Land Cover code for each water use sector
  water_use_land_cover = {'industrial':[13],'domestic':[6,7,8,9,10,13],'livestock':[6,7,8,9,10], "agriculture":[12,14]}
  
  # 1. Identify location of each type of water use
  # HRU lc map
  hrus_lc = np.copy(OUTPUT['hru_map'])
  for hru in np.arange(hydroblocks_info['nhru']): 
   idx = OUTPUT['hru_map'] == hru
   hrus_lc[idx]= OUTPUT['hru']['land_cover'][hru]
  m = hrus_lc == -9999.0
  lc = gdal_tools.read_raster('%s/lc_ea.tif' % (workspace))
  hrus_lc[m] = lc[m]
 
  for l in np.unique(hrus_lc): 
   idx = hrus_lc == l
   if l in water_use_land_cover[data_var]:
    hrus_lc[idx]=1.0
   else:
    hrus_lc[idx]=0.0
  
  wuse_lc_ea_file = '%s/%s_lc_ea.tif' % (workspace,data_var)
  gdal_tools.write_raster(wuse_lc_ea_file,md,hrus_lc)
  fine_size = hrus_lc.shape
  fine_res = abs(md['resx'])
  
  # Get the coarse water use info and regrid the fine lc to coarser lc
  wuse_lc_coarse_file = '%s/%s_latlon_coarse.tif' % (workspace,data_var)
  md = gdal_tools.retrieve_metadata(wuse_lc_coarse_file)
  minx = md['minx']
  miny = md['miny']
  maxx = md['maxx']
  maxy = md['maxy']
  res  = abs(md['resx'])
  lproj = md['proj4']+' +datum=WGS84'
  file_in = wuse_lc_ea_file
  file_out = '%s/%s_area_latlon_coarse.tif' % (workspace,data_var)
  os.system('gdalwarp -overwrite -t_srs \'%s\' -ot Float32 -dstnodata -9999 -tr %f %f -te %f %f %f %f -r average -q %s %s ' % (lproj,res,res,minx,miny,maxx,maxy,file_in,file_out))

  # Calculate the equivalent area of each grid
  data = gdal_tools.read_raster(file_out)
  md['nodata'] = -9999.0
  data[ data == md['nodata'] ] = 0.0
  coarse_size = data.shape
  #print 'res_fine', fine_res  #Should be in meters
  data_grid_area = fine_res*fine_res*(fine_size[0]/float(coarse_size[0]))*(fine_size[1]/float(coarse_size[1]))
  data = data*data_grid_area
  gdal_tools.write_raster(file_out,md,data)

  #hru_map = np.copy(OUTPUT['hru_map'])
  #hru_map[hru_map<0]=np.nan
  #plt.imshow(hru_map); plt.show()
  #plt.imshow(hrus_lc); plt.show()
  #plt.imshow(data); plt.show()
  
  #Compute the mapping for each hru
  for hru in np.arange(hydroblocks_info['nhru']):
   idx = OUTPUT['hru_map'] == hru
   #print "Catch:",hydroblocks_info['icatch'], "HRU: ", mask_fine[idx].astype(np.int)
   icells = np.unique(mask_fine[idx][mask_fine[idx] != -9999.0].astype(np.int))   # Add != -9999 for unique and bicount - Noemi
   counts = np.bincount(mask_fine[idx][mask_fine[idx] != -9999.0].astype(np.int))
   coords,pcts = [],[]
   for icell in icells:
    ilat = int(np.floor(icell/mask_coarse.shape[1]))
    jlat = icell - ilat*mask_coarse.shape[1]
    #ilat = int(mask_coarse.shape[0] - ilat - 1) #CAREFUL
    pct = float(counts[icell])/float(np.sum(counts))
    coords.append([ilat,jlat])
    pcts.append(pct)
   pcts = np.array(pcts)
   coords = list(np.array(coords).T)
   mapping_info[var][hru] = {'pcts':pcts,'coords':coords}

 #Iterate through variable creating water use product per HSU
 idate = info['time_info']['startdate']
 fdate = info['time_info']['enddate']
 dt = info['time_info']['dt']
 nt = int(3600*24/dt)*((fdate - idate).days+1)

 #Create structured array
 water_use = {}
 for data_var in wbd['files_water_use']:
  water_use[data_var] = np.zeros((nt,hydroblocks_info['nhru']))

 #Load data into structured array
 for data_var in wbd['files_water_use']:
  var = data_var
  date = idate
  file = wbd['files_water_use'][data_var]
  fp = nc.Dataset(file)
  #Determine the time steps to retrieve
  #fidate = ' '.join(fp.variables['t'].units.split(' ')[2::])
  #dates = nc.num2date(fp.variables['t'][:],units='hours since %s' % fidate)
  #mask_dates = (dates >= idate) & (dates <= fdate)
  nc_step = int(fp.variables['t'].units.split(' ')[0].split('h')[0])
  nc_idate = np.array(fp.variables['t'].units.split(' ')[2].split('-'))
  nc_nt = len(fp.variables['t'][:])
  dates = [datetime.datetime(int(nc_idate[0]),int(nc_idate[1]),int(nc_idate[2]))]
  for it in range(1,nc_nt): dates.append(dates[0] + datetime.timedelta(hours=it*nc_step))
  dates=np.array(dates)
  startdate = info['time_info']['startdate']
  enddate  = info['time_info']['enddate']
  mask_dates = (dates >= startdate) & (dates <= enddate)
  data = np.ma.getdata(fp.variables[var][mask_dates,:,:])
  fp.close()
  
  # convert water use volume from m3 to m3/m2
  file_out = '%s/%s_area_latlon_coarse.tif' % (workspace,data_var)
  wuse_area = gdal_tools.read_raster(file_out)
  m = ( wuse_area == 0.0 )
  data[:,m] = 0.0
  wuse_area[m] = 1.0
  data = data/wuse_area
 

  #Assing to hrus
  for hru in mapping_info[var]:
   if OUTPUT['hru']['land_cover'][hru] in water_use_land_cover[data_var]:
    #print data_var,data, data.shape, hru,mapping_info[var][hru]['pcts'],mapping_info[var][hru]['coords'],
    pcts = mapping_info[var][hru]['pcts']
    coords = mapping_info[var][hru]['coords']
    coords[0][coords[0] >= data.shape[1]] = data.shape[1] - 1
    coords[1][coords[1] >= data.shape[2]] = data.shape[2] - 1
    tmp = data[:,coords[0],coords[1]]
    tmp = pcts*tmp
    water_use[data_var][:,hru] = np.sum(tmp,axis=1)  # final variable m3/m2/s --> m/s of water demand
    #print hru, data_var, OUTPUT['hru']['land_cover'][hru], water_use[data_var][:,hru]
   else:
    water_use[data_var][:,hru] = 0.0

  #Write the water use the netcdf file (single chunk for now...)
  grp = hydroblocks_info['input_fp'].groups['water_use']
  grp.createVariable(var,'f4',('time','hru'))#,zlib=True)
  grp.variables[data_var][:] = water_use[data_var][:]

 if hydroblocks_info['water_management']['hwu_flag']:
  if len(wbd['files_water_use']) > 1 :
   #Add time information
   dates = []
   date = idate
   while date <= fdate:
    dates.append(date)
    date = date + datetime.timedelta(seconds=dt)
   dates = np.array(dates)
   var = grp.createVariable('time','f8',('time',))
   var.units = 'hours since %4d-01-01' % idate.year
   var.calendar = 'standard'
   dates = nc.date2num(dates,units=var.units,calendar=var.calendar)
   var[:] = dates[:]

 return


def get_f11_and_qtz_from_soil_type(soil_type):

 mapping = { 
 'sand':1,      
 'loamy_sand':2,
 'sandy_loam':3,
 'silt_loam':4, 
 'silt':5,      
 'loam':6,      
 'sandy_clay_loam':7,                             
 'silty_clay_loam':8,                             
 'clay_loam':9, 
 'sandy_clay':10,                                 
 'silty_clay':11,                                 
 'clay':12,
 'organic_matter':13,
 'water':14,
 'bedrock':15,
 'other':16,
 'playa':17,
 'lava':18,
 'white_sand':19                   
 } 
 
 # return F11 and QTZ
   
 if soil_type == 1:
  return -0.472, 0.92
 if soil_type == 2: 
  return -1.044, 0.82
 if soil_type == 3:
  return -0.569, 0.60
 if soil_type == 4:
  return 0.162, 0.25
 if soil_type == 5:
  return 0.162, 0.10
 if soil_type == 6: 
  return -0.327, 0.40
 if soil_type == 7:
  return -1.491, 0.60
 if soil_type == 8:
  return -1.118, 0.10
 if soil_type == 9:
  return -1.297, 0.35
 if soil_type == 10:
  return -3.209, 0.52
 if soil_type == 11:
  return -1.916, 0.10
 if soil_type == 12:
  return -2.138, 0.25
 if soil_type == 13:
  return -0.327, 0.05
 if soil_type == 14:
  return np.nan, 0.60
 if soil_type == 15:
  return np.nan, 0.07
 if soil_type == 16:
  return -1.044, 0.25
 if soil_type == 17:
  return -10.472, 0.60
 if soil_type == 18:
  return -0.472, 0.52
 if soil_type == 19:
  return -0.472, 0.92

 exit('Missing soil type!')   

 return


def Read_Metadata_File(file):

 import json
 metadata = json.load(open(file))

 return metadata

