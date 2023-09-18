import warnings
warnings.filterwarnings('ignore')
import sys
#import geopandas
import fiona
#sys.path.append('Tools')
import pickle
import datetime
import numpy as np
import scipy.sparse as sparse
import scipy.stats as stats
#import model_tools as mt
import os
import h5py
import netCDF4 as nc
import time
import glob
import numba
from geospatialtools import gdal_tools
from geospatialtools import terrain_tools
import gc
from scipy.interpolate import griddata
import copy
import collections
import shapely.geometry
import rasterio

#dir = os.path.dirname(os.path.abspath(__file__))
#sys.path.append('%s/../HydroBlocks/pyHWU/' % dir )
#import management_funcs as mgmt_funcs

def plot_data(data):

 import matplotlib.pyplot as plt
 data = np.ma.masked_array(data,data==-9999)
 plt.figure(figsize=(10,10))
 plt.imshow(data)
 plt.colorbar()
 plt.savefig('tmp.png')

 return

def Prepare_Model_Input_Data(hydroblocks_info,metadata_file):

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
 input_dir = hydroblocks_info['input_dir']
 os.system('mkdir -p %s' % input_dir)

 #Create soft link to HydroBlocks from within the directory
 HBdir = '%s/model/pyNoahMP' % (("/").join(__file__.split('/')[:-2]))
 HBedir = '%s/pyNoahMP%d' % (input_dir,hydroblocks_info['cid'])
 os.system('ln -s %s %s' % (HBdir,HBedir))

 #Create the dictionary to hold all of the data
 output = {}

 #Create the Latin Hypercube (Clustering)
 nhru = 1#hydroblocks_info['nhru']
 cid = hydroblocks_info['cid']

 #Get metadata
 md = gdal_tools.retrieve_metadata('%s/mask_latlon.tif' % workspace)
 
 #Prepare the input file
 wbd = {}
 wbd['bbox'] = {'minlat':md['miny'],'maxlat':md['maxy'],
                'minlon':md['minx'],'maxlon':md['maxx'],
                'res':abs(md['resx'])}
 wbd['files'] = {
  'WLTSMC':glob.glob('%s/theta1500/*'%workspace), #laura svp
  'TEXTURE_CLASS':'%s/texture_class/texture_class_latlon_2.5cm.tif' % workspace,
  'MAXSMC':glob.glob('%s/thetas/*'%workspace), #laura svp
  'BB':glob.glob('%s/bb/*'%workspace), #laura svp
  'DRYSMC':glob.glob('%s/thetar/*'%workspace), #laura svp
  'QTZ':glob.glob('%s/qtz/*'%workspace), #laura svp
  'SATDW':glob.glob('%s/dsat/*'%workspace), #laura svp
  'REFSMC':glob.glob('%s/theta33/*'%workspace), #laura svp
  'mask':'%s/mask_latlon.tif' % workspace,
  'SATPSI':glob.glob('%s/psisat/*'%workspace), #laura svp
  'lc':'%s/lc_latlon.tif' % workspace,
  'F11':'%s/f11_latlon.tif' % workspace,
  'SATDK':glob.glob('%s/ksat/*'%workspace), #laura svp
  'dem':'%s/dem_latlon.tif' % workspace,
  'acc':'%s/acc_latlon.tif' % workspace,
  'fdir':'%s/fdir_latlon.tif' % workspace,
  'demns':'%s/demns_latlon.tif' % workspace,
  'sand':'%s/sand/sand_latlon_2.5cm.tif' % workspace,
  'clay':'%s/clay/clay_latlon_2.5cm.tif' % workspace,
  'silt':'%s/silt/silt_latlon_2.5cm.tif' % workspace,
  'om':'%s/om/om_latlon_2.5cm.tif' % workspace,
  'bare30':'%s/bare30_latlon.tif' % workspace,
  'water30':'%s/water30_latlon.tif' % workspace,
  'tree30':'%s/tree30_latlon.tif' % workspace,
  'irrig_land':'%s/irrig_land_latlon.tif' % workspace,
  'dbedrock':'%s/dbedrock_latlon.tif' % workspace,
  'lstmean':'%s/lstmean_latlon.tif' % workspace,
  'lststd':'%s/lststd_latlon.tif' % workspace
  }
 #if hydroblocks_info['water_management']['hwu_agric_flag']:
 #  wbd['files']['irrig_land'] = '%s/irrig_land_latlon.tif' % workspace
 #  wbd['files']['start_growing_season'] = '%s/start_growing_season_latlon.tif' % workspace
 #  wbd['files']['end_growing_season']   = '%s/end_growing_season_latlon.tif' % workspace

 wbd['files_meteorology'] = {
  'lwdown':'%s/lwdown.nc' % workspace,
  'swdown':'%s/swdown.nc' % workspace,
  'tair':'%s/tair.nc' % workspace,
  'precip':'%s/precip.nc' % workspace,
  'psurf':'%s/psurf.nc' % workspace,
  'wind':'%s/wind.nc' % workspace,
  'spfh':'%s/spfh.nc' % workspace,
  }

 #if hydroblocks_info['water_management']['hwu_flag'] == True:
 # wbd['files_water_use'] = {}
 # if hydroblocks_info['water_management']['hwu_domest_flag']:
 #  wbd['files_water_use']['domestic']   = '%s/domestic.nc' % workspace
 # if hydroblocks_info['water_management']['hwu_indust_flag']:
 #  wbd['files_water_use']['industrial'] = '%s/industrial.nc' % workspace
 # if hydroblocks_info['water_management']['hwu_lstock_flag']:
 #  wbd['files_water_use']['livestock']  = '%s/livestock.nc' % workspace

 #Create the clusters and their connections
 (output,covariates) = Create_Clusters_And_Connections(workspace,wbd,output,input_dir,nhru,info,hydroblocks_info)

 #Extract the meteorological forcing
 print("Preparing the meteorology",flush=True)
 Prepare_Meteorology_Semidistributed(workspace,wbd,output,input_dir,info,hydroblocks_info,covariates)

 #Extract the water use demands
 #print("Preparing the water use",flush=True)
 #if hydroblocks_info['water_management']['hwu_flag'] == True:
 # Prepare_Water_Use_Semidistributed(workspace,wbd,output,input_dir,info,hydroblocks_info)

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
 grp.dx = 90.0#26.0#25.0#metadata['resx'] #UPDATE WITH DEM!

 #Write out the mapping
 hru_map = np.copy(output['hru_map'])
 hru_map[np.isnan(hru_map) == 1] = -9999.0
 file_ca = '%s/hru_mapping_latlon.tif' % input_dir
 metadata['nodata'] = -9999.0
 gdal_tools.write_raster(file_ca,metadata,hru_map)

 #Write out the hand map
 hand_map = np.copy(output['hand_map'])
 hand_map[np.isnan(hand_map) == 1] = -9999.0
 file_ca = '%s/hand_latlon.tif' % input_dir
 metadata['nodata'] = -9999.0
 gdal_tools.write_raster(file_ca,metadata,hand_map)

 #Write out the basin map
 basin_map = np.copy(output['basin_map'])
 basin_map[np.isnan(basin_map) == 1] = -9999.0
 file_ca = '%s/basins_latlon.tif' % input_dir
 metadata['nodata'] = -9999.0
 gdal_tools.write_raster(file_ca,metadata,basin_map)

 #Write out the basin cluster map
 basin_clusters_map = np.copy(output['basin_clusters_map'])
 basin_clusters_map[np.isnan(basin_clusters_map) == 1] = -9999.0
 file_ca = '%s/basin_clusters_latlon.tif' % input_dir
 metadata['nodata'] = -9999.0
 gdal_tools.write_raster(file_ca,metadata,basin_clusters_map)

 #If fully-distributed, save number of basins in metadata file, laura
 #if (hydroblocks_info['fully_distributed']==True) and (hydroblocks_info['connection_matrix_hbands']==True): 
  #with open(metadata_file,'r') as f:
   #import json
   #json_data=json.load(f)
  #json_data['hmc_parameters']['number_of_characteristic_subbasins_CID_%s'%cid]=len(np.unique(basin_clusters_map))-1
  #with open(metadata_file,'w') as f:
   #json.dump(json_data,f,indent=2)
  #hydroblocks_info=Read_Metadata_File(metadata_file)#Re-read metadata file

 #Write out the hand org map
 hand_org_map = np.copy(output['hand_org_map'])
 hand_org_map[np.isnan(hand_org_map) == 1] = -9999.0
 file_ca = '%s/hand_org_latlon.tif' % input_dir
 metadata['nodata'] = -9999.0
 gdal_tools.write_raster(file_ca,metadata,hand_org_map)

 #Write out the height band id map
 hband_map = np.copy(output['hband_map'])
 hband_map[np.isnan(hband_map) == 1] = -9999.0
 file_ca = '%s/hband_latlon.tif' % input_dir
 metadata['nodata'] = -9999.0
 gdal_tools.write_raster(file_ca,metadata,hband_map)

 #Write out the channels
 channel_map = np.copy(output['channel_map'])
 channel_map[np.isnan(channel_map) == 1] = -9999.0
 file_ca = '%s/channel_mapping_latlon.tif' % input_dir
 metadata['nodata'] = -9999.0
 gdal_tools.write_raster(file_ca,metadata,channel_map)

 #Write the connection matrices
 #width
 #laura's modification start
 if (hydroblocks_info['connection_matrix_hbands']==False):
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
 elif (hydroblocks_info['connection_matrix_hbands']==True): #and (hydroblocks_info['fully_distributed']==False):
  for i in range(1,(int(hydroblocks_info['hmc_parameters']["number_of_characteristic_subbasins"]+1))):
   text='wmatrix_Basin%s' %int(i)
   wmatrix=output['cmatrix_Basin%s' %int(i)]['width']
   nconnections = wmatrix.data.size
   grp = fp.createGroup(text)
   grp.createDimension('connections_columns',wmatrix.indices.size)
   grp.createDimension('connections_rows',wmatrix.indptr.size)
   grp.createVariable('data','f4',('connections_columns',))
   grp.createVariable('indices','f4',('connections_columns',))
   grp.createVariable('indptr','f4',('connections_rows',))
   grp.variables['data'][:] = wmatrix.data
   grp.variables['indices'][:] = wmatrix.indices
   grp.variables['indptr'][:] = wmatrix.indptr
 #elif (hydroblocks_info['connection_matrix_hbands']==True) and (hydroblocks_info['fully_distributed']==True):
  #for i in range(1,(int(hydroblocks_info['hmc_parameters']["number_of_characteristic_subbasins_CID_%s"%cid]+1))):
   #text='wmatrix_Basin%s' %int(i)
   #wmatrix=output['cmatrix_Basin%s' %int(i)]['width']
   #nconnections = wmatrix.data.size
   #grp = fp.createGroup(text)
   #grp.createDimension('connections_columns',wmatrix.indices.size)
   #grp.createDimension('connections_rows',wmatrix.indptr.size)
   #grp.createVariable('data','f4',('connections_columns',))
   #grp.createVariable('indices','f4',('connections_columns',))
   #grp.createVariable('indptr','f4',('connections_rows',))
   #grp.variables['data'][:] = wmatrix.data
   #grp.variables['indices'][:] = wmatrix.indices
   #grp.variables['indptr'][:] = wmatrix.indptr
   #end of laura's modification

 #Write the model parameters
 grp = fp.createGroup('parameters')
 vars = ['slope','area_pct','land_cover','channel',
        'dem','soil_texture_class','carea','area',
        'BB','F11','SATPSI','SATDW','QTZ','clay',
        'WLTSMC','MAXSMC','DRYSMC','REFSMC','SATDK',
        'm','hand','y_aspect','x_aspect','hru','hband',
        'lats','lons']

 #if hydroblocks_info['water_management']['hwu_agric_flag']:
 # for var in ['centroid_lats', 'centroid_lons', 'irrig_land', 'start_growing_season', 'end_growing_season']:
 #   vars.append(var)

 for var in vars:
  if var in ['slope','area_pct','land_cover','channel','dem','soil_texture_class','ti','carea','area','F11','clay','m','hand','y_aspect','x_aspect','hru','hband','lats','lons']: #laura svp
   grp.createVariable(var,'f4',('hru',))#,zlib=True)
   grp.variables[var][:] = data['parameters']['hru'][var] #laura svp
  else: #laura svp
   grp.createVariable(var,'f4',('hru','nsoil'))#,zlib=True) #laura svp
   grp.variables[var][:] = data['soil_properties_model']['hru'][var] #laura svp

 #if hydroblocks_info['water_management']['hwu_flag']:
 # grp.createVariable('hru_min_dist','f4',('hru','hru'))#,zlib=True)
 # grp.variables['hru_min_dist'][:] = data['hru']['hru_min_dist']

 #Write out the stream network info
 grp = fp.createGroup('stream_network')
 grp.createDimension('nc',data['stream_network']['slope'].size)
 for var in data['stream_network']:
  grp.createVariable(var,'f4',('nc'))
  grp.variables[var][:] = data['stream_network'][var][:]

 #Remove info from output
 del output['hru']

 #Add in the catchment info
 output['wbd'] = wbd

 #Close the file
 fp.close()

 return output

def Compute_HRUs_Semidistributed_HMC(covariates,mask,hydroblocks_info,wbd,eares,input_dir):

 #Define the parameters for the hierarchical multivariate clustering
 ncatchments = hydroblocks_info['hmc_parameters']['number_of_characteristic_subbasins']
 dh = hydroblocks_info['hmc_parameters']['average_height_difference_between_bands']
 nclusters = hydroblocks_info['hmc_parameters']['number_of_intraband_clusters']

 #Bring out the mask_all
 mask_all = covariates['mask_all']

 #Bring out the flow direction (Convert flow direction from int to 2d approach)
 fdir = terrain_tools.transform_arcgis_fdir(covariates['fdir'])

 #Pre-process DEM
 dem = covariates['dem']
 demns = np.copy(dem)
 covariates['demns'] = demns
 area_all = covariates['acc']*10**6 #km2->m2
 area_all_cp = np.copy(area_all)
  
 #Calculate slope and aspect
 print("Calculating slope and aspect",flush=True)
 res_array = np.copy(demns)
 res_array[:] = eares
 (slope,aspect) = terrain_tools.ttf.calculate_slope_and_aspect(np.flipud(demns),res_array,res_array)
 slope = np.flipud(slope)
 aspect = np.flipud(aspect)

 #Compute accumulated area
 m2 = np.copy(mask_all)
 m2[m2 > 0] = 1
 mall = np.copy(m2)
 mall[m2 <= 0] = 0
 mall = mall.astype(np.bool)
 print("Calculating accumulated area",flush=True)
 #area = terrain_tools.ttf.calculate_d8_acc_pfdir(demns,m2,eares,fdir)
 area = area_all

 #Calculate channel initiation points (2 parameters)
 C = area/eares*slope**2
 cthrs = hydroblocks_info['channel_initiation']["athrs"]#laura #10**6
 ipoints = ((area > cthrs)).astype(np.int32)
 ipoints[ipoints == 0] = -9999

 #Create area for channel delineation
 ac = np.copy(area_all)
 ac[mask == 0] = -9999 #used to calculate channels within the subdomain only
 fdc = fdir
 ac_all = area_all
 fdc_all = fdir

 #Compute the channels
 print("Defining channels",flush=True)
 (channels,channels_wob,channel_topology,tmp1,crds,channel_outlet_id,channel_target_mp,channel_target_crds,channel_inlet_id,channel_inlet_target_mp,channel_inlet_target_crds) = terrain_tools.ttf.calculate_channels_wocean_wprop_wcrds(ac,ac_all,cthrs,cthrs,fdc,mask,mask_all,np.flipud(covariates['lats']),covariates['lons'])

 #Curate list output
 channel_topology = channel_topology[channel_topology != -9999]
 m = channel_outlet_id != 0
 channel_outlet_id = channel_outlet_id[m]
 channel_target_mp = channel_target_mp[m]
 channel_target_crds = channel_target_crds[m,:]
 crds = crds[crds[:,0,0] != -9999,:,:] 
 m = channel_inlet_id != 0
 channel_inlet_id = channel_inlet_id[m]
 channel_inlet_target_mp = channel_inlet_target_mp[m,:]
 channel_inlet_target_mp[channel_inlet_target_mp == 0] = -9999
 channel_inlet_target_crds = channel_inlet_target_crds[m,:,:]
 #Convert channel ids to start from 0 (instead of 1)
 channel_outlet_id[channel_outlet_id>0] = channel_outlet_id[channel_outlet_id>0] - 1
 channel_inlet_id[channel_inlet_id>0] = channel_inlet_id[channel_inlet_id>0] - 1
 
 #If the dem is undefined then set to undefined
 channels[dem == -9999] = -9999

 #Determine inlets/outlets
 db_routing = {}
 db_routing['mp_connectivity'] = {'channel_outlet_id':channel_outlet_id,
                                  'channel_target_mp':channel_target_mp,
                                  'channel_target_crds':channel_target_crds,
                                  'channel_crds':crds,
                                  'channel_inlet_id':channel_inlet_id,
                                  'channel_inlet_target_mp':channel_inlet_target_mp,
                                  'channel_inlet_target_crds':channel_inlet_target_crds}
 db_routing['i/o'] = terrain_tools.calculate_inlets_oulets(channels_wob,fdir,area_all,mask,np.flipud(covariates['lats']),covariates['lons'],mask_all,area_all)

 #Compute and output the list of the channel positions
 lst_crds = []
 for icrd in range(crds.shape[0]):
   mcrd = crds[icrd,:,0] != -9999
   if (np.sum(mcrd) == 0):break
   crds_i = crds[icrd,mcrd,:]
   if crds_i.shape[0] > 1:
       lst_crds.append(shapely.geometry.LineString(np.fliplr(crds_i)))
   else:
       lst_crds.append(shapely.geometry.Point(np.flipud(crds_i[0,:])))
 db_routing['crds'] = lst_crds

 #Compute the basins
 print("Defining basins",flush=True)
 #basins = terrain_tools.ttf.delineate_basins(channels,m2,fdir)
 basins_wob = terrain_tools.ttf.delineate_basins(channels_wob,mask,fdir)
 basins = basins_wob
 
 #Compute channel properties
 db_channels = terrain_tools.calculate_channel_properties(channels_wob,channel_topology,slope,eares,mask,area_all,area_all_cp,basins_wob,hydroblocks_info['parameter_scaling'])

 #Calculate the height above nearest drainage area
 print("Computing height above nearest drainage area",flush=True)
 hand = terrain_tools.ttf.calculate_depth2channel(channels_wob,basins_wob,fdir,demns)

 #Fill in hand that is undefined (probably flow direction issues)
 hand[(hand == -9999) & (basins_wob!=-9999)] = 0.0

 # cleanup
 slope[mask != 1] = -9999
 aspect[mask != 1] = -9999
 area[mask != 1] = -9999
 channels[mask != 1] = -9999
 basins[mask != 1] = -9999

 # save covariates
 covariates['slope'] = slope
 covariates['aspect'] = aspect
 covariates['x_aspect'] = np.sin(aspect)
 covariates['y_aspect'] = np.cos(aspect)
 covariates['carea'] = area_all_cp#area
 covariates['carea_log10'] = np.log10(area_all_cp)#area
 covariates['hand'] = hand

 #Calculate the subbasin properties
 print("Assembling the subbasin properties",flush=True)
 vars1 = hydroblocks_info['hmc_parameters']['subbasin_clustering_covariates']
 vars = []
 for var in vars1:
  if var not in ['width','bankfull','length','area']:
   vars.append(var)
 hp_in = terrain_tools.calculate_basin_properties_updated(basins_wob,eares,covariates,vars)
 #sort hp_in (should go in geospatialtools)
 argsort = np.argsort(hp_in['bid'])
 for var in hp_in:
  hp_in[var] = hp_in[var][argsort]
 #bring in channel variables
 for var in ['width','bankfull','length','area']:
  hp_in[var] = db_channels[var]

 #Clustering the basins
 print("Clustering the basins",flush=True)

 #Flag fully distributed simulation
 #flag_fd=hydroblocks_info['fully_distributed']
 #if flag_fd==True: #laura
  #basin_clusters=np.copy(basins_wob) #laura
  #nhru=len(np.unique(basin_clusters))-1
  #if np.min(basin_clusters[basin_clusters!=-9999])==0:
   #basin_clusters[basin_clusters!=-9999]=basin_clusters[basin_clusters!=-9999]+1 
   #hydroblocks_info['hmc_parameters']["number_of_characteristic_subbasins"]=len(np.unique(basin_clusters))-1 #laura
  #print(np.unique(basin_clusters),flush=True)
 #else:
  #Set the ncatchments to be at least the number of basins
 ncatchments = min(ncatchments,np.unique(basins_wob)[1:].size)
 subbasin_clustering_cov=hydroblocks_info['hmc_parameters']['subbasin_clustering_covariates']#laura
 #dissaggregate land cover if it is in covariates
 if 'lc' in subbasin_clustering_cov:#laura
  subbasin_clustering_cov.remove('lc') #laura
  subbasin_clustering_cov=subbasin_clustering_cov+['lc_w_now','lc_urb_nourb','lc_grass_forest'] #laura, divide land cover in water_vs_no_water, urban_vs_no_urban, and grass_vs_forest (including grass and shrubs as intermediate values) #laura

  #Assemble input data
 cvs = {}
 for var in subbasin_clustering_cov: #laura
  if var in ['lc_w_now','lc_urb_nourb','lc_grass_forest']: #laura
   if var=='lc_w_now': #laura
    lc_mask=covariates['lc_17'] #laura
   elif var=='lc_urb_nourb': #laura
    lc_mask=covariates['lc_13'] #laura
   elif var=='lc_grass_forest': #laura
    lc_mask=covariates['lc_4'] #deciduous_forest #laura
    if 'lc_2' in covariates: #laura
     lc_mask[covariates['lc_2']==1]=1 #evergreen_forest #laura
    if 'lc_5' in covariates: #laura
     lc_mask[covariates['lc_5']==1]=1 #mixed_forest #laura
    if 'lc_6' in covariates: #laura
     lc_mask[covariates['lc_6']==1]=0.66 #shrub/scrub #laura
    if 'lc_11' in covariates: #laura
     lc_mask[covariates['lc_11']==1]=0.66 #wetlands #laura
    if 'lc_12' in covariates: #laura
     lc_mask[covariates['lc_12']==1]=0.66 #pasture/hay/cultivated_crops #laura
    if 'lc_10' in covariates: #laura
     lc_mask[covariates['lc_10']==1]=0.33 #grassland #laura
    if 'lc_16' in covariates: #laura
     lc_mask[covariates['lc_16']==1]=0.01 #barren_land #laura
  
   cvs[var] = {'min':0, #laura
               'max':1, #laura
               't':-9999, #laura
               'd':lc_mask} #laura
  else: #laura
   tmp = np.copy(hp_in[var])
   cvs[var] = {'min':np.min(tmp),
               'max':np.max(tmp),
               't':-9999,
               'd':tmp}

 (basin_clusters,) = terrain_tools.cluster_basins_updated(basins_wob,cvs,hp_in,ncatchments)
 #Calculate average bankfull depth per basin cluster
 ubcs = np.unique(basin_clusters)
 ubcs = ubcs[ubcs != -9999]
 for ubc in ubcs:
  ubs = np.unique(basins_wob[basin_clusters == ubc])
  ubs = ubs[ubs != -9999]
  #Compute mean width and bankfull depth
  db_channels['width'][ubs-1] = np.mean(db_channels['width'][ubs-1])
  db_channels['bankfull'][ubs-1] = np.mean(db_channels['bankfull'][ubs-1])
 
 #Divide each subbasin into height bands
 print("Discretizing clusters of basins (hbands)",flush=True) #laura
 n_binning = dh #HACK 
 max_nbins = 100
 (tiles,new_hand,tile_position) = terrain_tools.create_basin_tiles_updated(basin_clusters,hand,basins_wob,n_binning,hydroblocks_info['cid'],max_nbins)

 #Assemble river/hillslope database for routing/two-way connectivity
 (db_routing,area_adj,new_hand2) = Build_Hillslope_River_Database(channels_wob,mask,fdir,eares,tiles,hand,basins_wob,basin_clusters,new_hand,db_routing,ubcs,tile_position,db_channels)

 #Disagregate land cover
 intraband_clust_vars = hydroblocks_info['hmc_parameters']['intraband_clustering_covariates']
 if 'lc' in intraband_clust_vars: 
  intraband_clust_vars.remove('lc')
  ##disag = [i for i in covariates.keys() if 'lc_' in i] #laura, commented out so lc not overwhelms clustering
  ##intraband_clust_vars = intraband_clust_vars + disag #laura, commented out so lc not overwhelms clustering
  intraband_clust_vars=intraband_clust_vars+['lc_w_now','lc_urb_nourb','lc_grass_forest'] #laura, divide land cover in water_vs_no_water, urban_vs_no_urban, and grass_vs_forest (including grass and shrubs as intermediate values)

 #Calculate the hrus (kmeans on each tile of each basin)
 cvs = {}
 for var in intraband_clust_vars:
  if var in ['lc_w_now','lc_urb_nourb','lc_grass_forest']:
   if var=='lc_w_now':
    lc_mask=covariates['lc_17']
   elif var=='lc_urb_nourb':
    lc_mask=covariates['lc_13']
   elif var=='lc_grass_forest':
    lc_mask=covariates['lc_4'] #deciduous_forest
    if 'lc_2' in covariates:
     lc_mask[covariates['lc_2']==1]=1 #evergreen_forest
    if 'lc_5' in covariates:
     lc_mask[covariates['lc_5']==1]=1 #mixed_forest
    if 'lc_6' in covariates:
     lc_mask[covariates['lc_6']==1]=0.66 #shrub/scrub
#   lc_mask[covariates['lc_7']==1]=0.66 #dwarf/scrub, Alaska only
    if 'lc_11' in covariates:
     lc_mask[covariates['lc_11']==1]=0.66 #wetlands
    if 'lc_12' in covariates:
     lc_mask[covariates['lc_12']==1]=0.66 #pasture/hay/cultivated_crops
    if 'lc_10' in covariates:
     lc_mask[covariates['lc_10']==1]=0.33 #grassland
#   lc_mask[covariates['lc_19']==1]=0.33 #moss/sedge/lichens, Alaska only
    if 'lc_16' in covariates:
     lc_mask[covariates['lc_16']==1]=0.01 #barren_land  

   cvs[var] = {'min':0,
               'max':1,
               't':-9999,
               'd':lc_mask}
  else:
   cvs[var] = {'min':np.min(covariates[var][covariates[var]!=-9999]),
               'max':np.max(covariates[var][covariates[var]!=-9999]),
               't':-9999,
               'd':covariates[var]}
 
 print("Clustering the height bands into clusters", flush=True)
 #A.Ensure match between basin cluster map and tiles map
 m = (basin_clusters == -9999) | (tiles == -9999)
 basin_clusters[m] = -9999
 tiles[m] = -9999
 
 hrus = terrain_tools.create_hrus_hydroblocks(basin_clusters,tiles,cvs,nclusters,hydroblocks_info['cid']) #laura
 hrus[hrus!=-9999] = hrus[hrus!=-9999] - 1
 nhru = np.unique(hrus[hrus!=-9999]).size
 #print(' CID',hydroblocks_info['cid'],'#HRUs          ',nhru,flush=True)
 #print(' CID',hydroblocks_info['cid'],'#Total pixels  ',np.sum(basin_clusters!=-9999))

 #Save the channel info
 pickle.dump(db_routing,open('%s/routing_info.pck' % input_dir,'wb'))
 pickle.dump(db_routing['i/o'],open('%s/routing_io.pck' % input_dir,'wb'))
 pickle.dump(db_routing['mp_connectivity'],open('%s/routing_mp_connectivity.pck' % input_dir,'wb'))

 #Construct HMC info for creating connections matrix
 HMC_info = {}
 HMC_info['basins'] = basins
 HMC_info['tile_position'] = tile_position
 HMC_info['channel_map'] = channels_wob

 #return (hrus.astype(np.float32),nhru,new_hand,HMC_info,covariates,db_channels,hand,
 return (hrus.astype(np.float32),nhru,new_hand,HMC_info,covariates,db_channels,new_hand2,
         basins,basin_clusters,hand,tiles,area_adj,tile_position)

def Build_Hillslope_River_Database(channels_wob,mask,fdir,eares,tiles,hand,basins_wob,
    basin_clusters,new_hand,db_routing,ubcs,tile_position,db_channels):

 #Calculate histogram of travel distances per height band
 t2c = terrain_tools.ttf.calculate_distance2channel(channels_wob,mask,fdir,eares)
 uhbands = np.unique(tiles)
 uhbands = uhbands[uhbands != -9999]
 bins = np.linspace(0,100,101)
 uhs = []
 for hband in uhbands:
   m = tiles == hband
   hist = np.histogram(t2c[m]/0.1/3600.0,bins=bins,density=True)[0]
   uhs.append(hist)
 uhs = {'data':np.array(uhs),'bins':bins[:]}
 db_routing['uh_per_hband'] = uhs

 #Calculate average hand per height band
 new_hand2 = np.copy(hand)
 for hband in uhbands:
   mhand = tiles == hband
   new_hand2[mhand] = np.mean(new_hand[mhand])

 #Burn the average bankfull depth into newhand2
 for ubc in ubcs:
  ubs = np.unique(basins_wob[basin_clusters == ubc])
  ubs = ubs[ubs != -9999]
  mnw = (basin_clusters == ubc) & (tile_position != 0)
  mnw1 = (basin_clusters == ubc) & (tile_position == 1)
  new_hand2[mnw] = new_hand2[mnw] - np.mean(new_hand2[mnw1]) + np.mean(db_channels['bankfull'][ubs-1])
 new_hand2[np.isnan(new_hand2) == 1] = 0.0
   
 #Compute the areal coverage of each hand value within the basin
 db_routing['reach_hand_area'] = {}
 db_routing['reach_hand_hband'] = {}
 db_routing['reach_hband_area'] = {}
 for i in range(basins_wob.shape[0]):
  for j in range(basins_wob.shape[1]):
   basin = basins_wob[i,j]
   h = new_hand2[i,j]
   hband = tiles[i,j]
   if basin <= 0:continue
   if basin not in db_routing['reach_hand_area']:db_routing['reach_hand_area'][basin] = collections.OrderedDict()
   if h not in db_routing['reach_hand_area'][basin]: db_routing['reach_hand_area'][basin][h] = 0.0
   if basin not in db_routing['reach_hand_hband']:db_routing['reach_hand_hband'][basin] = collections.OrderedDict()
   db_routing['reach_hand_area'][basin][h] += eares**2
   db_routing['reach_hand_hband'][basin][h] = hband
 
 #Compute channel cross section information
 odb = {'Ac':0.0*np.ones((len(db_routing['reach_hand_area'].keys()),100)),
       'Pc':0.0*np.ones((len(db_routing['reach_hand_area'].keys()),100)),
       'Af':0.0*np.ones((len(db_routing['reach_hand_area'].keys()),100)),
       'Pf':0.0*np.ones((len(db_routing['reach_hand_area'].keys()),100)),
       'W':0.0*np.ones((len(db_routing['reach_hand_area'].keys()),100)),
       'M':0.0*np.ones((len(db_routing['reach_hand_area'].keys()),100)),
       'hand':0.0*np.ones((len(db_routing['reach_hand_area'].keys()),100)),
       'hband':-9999*np.ones((len(db_routing['reach_hand_area'].keys()),100)).astype(np.int32)}
 for b in db_routing['reach_hand_area']:
  #Define reach length
  c_length = db_channels['length'][b-1]
  #Sort from lowest to highest hand
  c_hand = np.array(list(db_routing['reach_hand_area'][b].keys()))
  c_area = np.array(list(db_routing['reach_hand_area'][b].values()))
  c_hband = np.array(list(db_routing['reach_hand_hband'][b].values()))
  argsort = np.argsort(c_hand)
  c_hand = c_hand[argsort]
  c_hband = c_hband[argsort]
  odb['hband'][b-1,0:c_hband.size] = c_hband[:]
  #Burn in a channel depth
  if c_hand.size > 1:
   #1.first remove existing difference between channel and adjacent hand value
   c_hand[1:] = c_hand[1:] - (c_hand[1] - c_hand[0])
   #2.then burn in the channel bankfull depth
   c_hand[1:] = c_hand[1:] + db_channels['bankfull'][b-1] #m
  
  c_area = c_area[argsort]
  #Calculate widths of each HRU/height band
  c_width = c_area/c_length
  if c_width.size > 1:
   #Correct channel width using provided estimates
   c_width_diff = db_channels['width'][b-1] - c_width[0]
   #Ensure that the change of width doesn't cause negative values
   if (c_width_diff > 0.9*c_width[1]):
    c_width_diff = 0.9*c_width[1]
   #Update the channel width
   c_width[0] = c_width[0] + c_width_diff 
   #Add the difference to the adjacent HRU
   c_width[1] = c_width[1] - c_width_diff
  #Calculate slope
  c_slope = np.zeros(c_width.size)
  if c_slope.size > 1:
   c_slope[1:-1] = (c_hand[2:] - c_hand[1:-1])/(c_width[1:-1]/2)
   c_slope[-1] = c_slope[-2]
  #Add the channel depth
  odb['M'][b-1,0:c_slope.size] = c_slope[:]
    
  #Adjust the areal coverage of all the HRUs/bands
  c_area = c_length*c_width
  if (np.unique(c_area)[0] <= 0):
   print(c_width)
   print(c_area)
   print(c_length)
   exit()
  #Update values in dictionary (due to correcting for channel info)
  db_routing['reach_hand_area'][b] = collections.OrderedDict()
  db_routing['reach_hand_hband'][b] = collections.OrderedDict()
  db_routing['reach_hband_area'][b] = collections.OrderedDict()
  for ih in range(c_hand.size):
    db_routing['reach_hand_area'][b][c_hand[ih]] = c_area[ih]
    db_routing['reach_hand_hband'][b][c_hand[ih]] = c_hband[ih]
    db_routing['reach_hband_area'][b][c_hband[ih]] = c_area[ih]
  #Update the channel depth
  odb['hand'][b-1,0:c_hand.size] = c_hand[:]
  #Calculate width
  odb['W'][b-1,0:c_width.size] = c_width[:]
  #Calculate wetted perimeter at each stage
  dPc = []
  dPf = []
  for iseg in range(c_width.size-1):
   if iseg == 0:
    dPc.append(c_width[0] + 2*(c_hand[iseg+1]-c_hand[iseg]))
    dPf.append(0.0)
   else:
    dPc.append(0.0)
    dPf.append(c_width[iseg-1] + 2*(c_width[iseg]/2**2 + (c_hand[iseg+1]-c_hand[iseg])**2)**0.5)
  #Compute perimieters for channel and floodplain
  dPc = np.array(dPc)
  dPf = np.array(dPf)
  Pc = np.cumsum(dPc)
  Pf = np.cumsum(dPf)
  odb['Pc'][b-1,0] = 0.0
  odb['Pc'][b-1,1:Pc.size+1] = Pc[:]
  odb['Pf'][b-1,0] = 0.0
  odb['Pf'][b-1,1:Pf.size+1] = Pf[:]
  #Calculate wetted cross sectional area at each stage
  dAc = []
  dAf = []
  dA = []
  for iseg in range(c_width.size-1):
   if iseg == 0:
    dAc.append(c_width[0]*(c_hand[iseg+1]-c_hand[iseg]))
    dAf.append(0.0)
    dA.append(c_width[0]*(c_hand[iseg+1]-c_hand[iseg]))
   else:
    dAc.append(c_width[0]*(c_hand[iseg+1]-c_hand[iseg]))
    pt1 = np.sum(c_width[1:iseg]*(c_hand[iseg+1]-c_hand[iseg]))
    pt2 = 2*c_width[iseg]/2*(c_hand[iseg+1]-c_hand[iseg])/2
    tmp = pt1+pt2
    dAf.append(tmp)
  #Compute cross sectional areas for channel and floodplain
  dAc = np.array(dAc)
  dAf = np.array(dAf)
  Ac = np.cumsum(dAc)
  Af = np.cumsum(dAf)
  odb['Ac'][b-1,0] = 0.0
  odb['Ac'][b-1,1:Ac.size+1] = Ac[:]
  odb['Af'][b-1,0] = 0.0
  odb['Af'][b-1,1:Af.size+1] = Af[:]

 #Calculate inundation height at each stage
 db_routing['reach_cross_section'] = copy.deepcopy(odb)

 #Create array of areas per reach/hband
 reach2hband = np.zeros((np.unique(list(db_routing['reach_hband_area'].keys())).size,uhbands.size))
 for reach in db_routing['reach_hband_area']:
  for hband in db_routing['reach_hband_area'][reach]:
   tmp = db_routing['reach_hband_area'][reach][hband]
   reach2hband[reach-1,hband] = db_routing['reach_hband_area'][reach][hband]

 #Correct the area per grid cell array (and then apply to construct database)
 hband_areas = np.array(np.sum(reach2hband,axis=0))
 area_adj = np.zeros(new_hand.shape)
 area_adj[:] = -9999.0
 for hband in uhbands:
  m = tiles == hband
  area_adj[m] = hband_areas[hband]/np.sum(m)

 return (db_routing,area_adj,new_hand2)

def Assign_Parameters_Semidistributed_svp(covariates,metadata,hydroblocks_info,OUTPUT,cluster_ids,mask,hbands,area_adj,dz_data,dz_model):

 nhru = hydroblocks_info['nhru']
 #Initialize the arrays
 vars = ['area','area_pct','F11','slope','dem','carea','channel',
         'land_cover','soil_texture_class','clay','sand','silt',
         'm','hand','x_aspect','y_aspect','hru','hband','lats','lons'] #laura svp

 vars_s = ['BB','DRYSMC','MAXSMC','REFSMC','SATPSI','SATDK','SATDW','WLTSMC',                 'QTZ'] #laura svp

 #if hydroblocks_info['water_management']['hwu_agric_flag']:
 # for var in ['centroid_lats', 'centroid_lons', 'irrig_land', 'start_growing_season', 'end_growing_season']:
 #   vars.append(var)

 OUTPUT['parameters']={} #laura svp
 OUTPUT['parameters']['hru'] = {} #laura svp
 OUTPUT['soil_properties_model']={} #laura svp
 OUTPUT['soil_properties_model']['hru'] = {} #laura svp
 OUTPUT['soil_properties_data']={} #laura svp
 OUTPUT['soil_properties_data']['hru'] = {} #laura svp

 #if hydroblocks_info['water_management']['hwu_flag']: OUTPUT['hru']['hru_min_dist'] = np.zeros((nhru,nhru))

 for var in vars:
   OUTPUT['parameters']['hru'][var] = np.zeros(nhru)

 for var in vars_s:
   OUTPUT['soil_properties_model']['hru'][var] = np.zeros([nhru,len(dz_model)])
   OUTPUT['soil_properties_data']['hru'][var] = np.zeros([nhru,len(dz_data[var])])
  
 #Metadata
 for hru in np.arange(nhru):
  #Set indices
  idx = np.where(cluster_ids == hru)
  #Define hru
  OUTPUT['parameters']['hru']['hru'][hru] = hru
  #Define height band id
  OUTPUT['parameters']['hru']['hband'][hru] = np.mean(hbands[idx])
  #Calculate area per hru
  OUTPUT['parameters']['hru']['area'][hru] = np.sum(area_adj[idx])
  #Calculate area percentage per hru
  OUTPUT['parameters']['hru']['area_pct'][hru] = 100*OUTPUT['parameters']['hru']['area'][hru]/(np.sum(area_adj[area_adj != -9999]))

  #Constant Soil properties laura svp
  for var in ['F11','clay','sand','silt']:
   OUTPUT['parameters']['hru'][var][hru] = np.mean(covariates[var][idx])

  #Average Slope
  OUTPUT['parameters']['hru']['slope'][hru] = np.nanmean(covariates['slope'][idx])
  #DEM
  OUTPUT['parameters']['hru']['dem'][hru] = np.nanmean(covariates['dem'][idx])
  #HAND
  OUTPUT['parameters']['hru']['hand'][hru] = np.nanmean(covariates['hand'][idx])
  #Average Catchment Area
  OUTPUT['parameters']['hru']['carea'][hru] = np.nanmean(covariates['carea'][idx])
  OUTPUT['parameters']['hru']['x_aspect'][hru] = np.nanmean(covariates['x_aspect'][idx])
  OUTPUT['parameters']['hru']['y_aspect'][hru] = np.nanmean(covariates['y_aspect'][idx])
  #Average geographic coordinates
  OUTPUT['parameters']['hru']['lats'][hru] = np.nanmean(covariates['lats'][idx])
  OUTPUT['parameters']['hru']['lons'][hru] = np.nanmean(covariates['lons'][idx])
  #Land cover type 
  tmp = covariates['lc'][idx]
  tmp = tmp[tmp>=1]
  if len(tmp) >= 1 :
   OUTPUT['parameters']['hru']['land_cover'][hru] = stats.mode(tmp)[0][0]
  else:
   OUTPUT['parameters']['hru']['land_cover'][hru] = 17  # if there is no valid value, set to water #Noemi

  #Soil texture class constant in vertical laura svp
  OUTPUT['parameters']['hru']['soil_texture_class'][hru] = stats.mode(covariates['TEXTURE_CLASS'][idx])[0][0]

  #Define the estimate for the model parameters
  OUTPUT['parameters']['hru']['m'][hru] = np.nanmean(covariates['dbedrock'][idx]) #0.1 #Form of the exponential decline in conductivity (0.01-1.0)

  #Vertically variable Soil properties laura svp
  for var in ['BB','DRYSMC','MAXSMC','SATPSI','SATDK','SATDW','QTZ']:
   #print(var,np.unique(covariates[var][idx]))
   if var in ['SATDK','SATDW']:
    i=0
    for depth in covariates[var]:
     try:
      OUTPUT['soil_properties_data']['hru'][var][hru,i] = stats.mstats.hmean(covariates[var][depth][idx])/3600.0/1000.0 #mm/hr -> m/s
     except:
      OUTPUT['soil_properties_data']['hru'][var][hru,i] = 1.41E-4
     i=i+1
   else:
    i=0
    for depth in covariates[var]:
     OUTPUT['soil_properties_data']['hru'][var][hru,i] = np.mean(covariates[var][depth][idx])
     i=i+1

  i=0
  for depth in covariates[var]:
   OUTPUT['soil_properties_data']['hru']['WLTSMC'][hru,i] = OUTPUT['soil_properties_data']['hru']['MAXSMC'][hru,i]*(OUTPUT['soil_properties_data']['hru']['SATPSI'][hru,i]/150)**(1/OUTPUT['soil_properties_data']['hru']['BB'][hru,i])
   OUTPUT['soil_properties_data']['hru']['REFSMC'][hru,i] = OUTPUT['soil_properties_data']['hru']['MAXSMC'][hru,i]*(OUTPUT['soil_properties_data']['hru']['SATPSI'][hru,i]/3.3)**(1/OUTPUT['soil_properties_data']['hru']['BB'][hru,i])
   i=i+1

 #Sort data depths and soil properties for vertical interpolation laura svp 
 for var in vars_s:
  ind=np.argsort(dz_data[var])
  OUTPUT['soil_properties_data']['hru'][var]=OUTPUT['soil_properties_data']['hru'][var][:,ind[:]]
  dz_data[var]=np.sort(dz_data[var])

 #Vertical interpolation laura svp
 for hru in np.arange(nhru):
  for var in vars_s:
   fp=OUTPUT['soil_properties_data']['hru'][var][hru,:]
   xp=np.array(dz_data[var])
   x=dz_model
   if np.sum(fp==-9999.0)>0 and var=='SATDK':
    fp[fp==-9999.0]=10**-10
   elif np.sum(fp==-9999.0)>0 and var=='BB':
    fp[fp==-9999.0]=11.55
   elif np.sum(fp==-9999.0)>0 and var=='DRYSMC':
    fp[fp==-9999.0]=0.138
   elif np.sum(fp==-9999.0)>0 and var=='MAXSMC':
    fp[fp==-9999.0]=0.468
   elif np.sum(fp==-9999.0)>0 and var=='REFSMC':
    fp[fp==-9999.0]=0.412
   elif np.sum(fp==-9999.0)>0 and var=='SATPSI':
    fp[fp==-9999.0]=0.468
   elif np.sum(fp==-9999.0)>0 and var=='SATDW':
    fp[fp==-9999.0]=1.12E-5
   elif np.sum(fp==-9999.0)>0 and var=='WLTSMC':
    fp[fp==-9999.0]=0.030
   elif np.sum(fp==-9999.0)>0 and var=='QTZ':
    fp[fp==-9999.0]=0.25
   OUTPUT['soil_properties_model']['hru'][var][hru,:]=np.interp(x,xp,fp,left=0)

 return OUTPUT

@numba.jit(nopython=True,cache=True)
def Determine_HMC_Connectivity(h1,h2,b1,b2,tp1,tp2,ivc,irc,ibc):

 if (h2 == -9999):return False
 if (h1 == h2):return True
 if ((tp1 == tp2) & (tp1 == 0) & (b1 != b2) & (ivc)):return True
 if (b1 != b2) & (irc == False):return False
 if (np.abs(tp1 - tp2) != 1) & (ibc == False):return False

 return True

def Calculate_HRU_Connections_Matrix_HMC(covariates,cluster_ids,nhru,dx,HMC_info,hydroblocks_info):

 #Add pointers for simplicity
 tile_position = HMC_info['tile_position']
 basins = HMC_info['basins']
 ivc = hydroblocks_info['hmc_parameters']['intervalley_connectivity']
 irc = hydroblocks_info['hmc_parameters']['interridge_connectivity']
 ibc = hydroblocks_info['hmc_parameters']['intraband_connectivity']
 
 #Perform the work
 (hdst,horg) = Calculate_HRU_Connections_Matrix_HMC_workhorse(cluster_ids,dx,tile_position,
               basins,ivc,irc,ibc)

 #Prepare the sparse matrix
 cmatrix = sparse.coo_matrix((np.ones(hdst.size),(horg,hdst)),shape=(nhru,nhru),dtype=np.float32)
 cmatrix = cmatrix.tocsr()

 #Prepare length, width, and ksat matrices
 wmatrix = cmatrix.copy()
 wmatrix.multiply(dx) #wmatrix[:] = dx*wmatrix[:]

 #Prepare output dictionary
 cdata = {'width':wmatrix.T,}

 return cdata

def Calculate_HRU_Connections_Matrix_HMC_hbands(hbands,dx,HMC_info,hydroblocks_info):
#Removed covariates and cluster ids from parameters, replace nhrus for nhbands, laura
 #Add pointers for simplicity
 tile_position = HMC_info['tile_position']
 basins = HMC_info['basins']
 ivc = hydroblocks_info['hmc_parameters']['intervalley_connectivity']
 irc = hydroblocks_info['hmc_parameters']['interridge_connectivity']
 ibc = hydroblocks_info['hmc_parameters']['intraband_connectivity']

 #Perform the work
 (hdst,horg) = Calculate_HRU_Connections_Matrix_HMC_workhorse(hbands,dx,tile_position,
               basins,ivc,irc,ibc) #laura, nhrus replaced with nhbands

 #Prepare the sparse matrix
 cmatrix = sparse.coo_matrix((np.ones(hdst.size),(horg,hdst)),shape=(int(np.unique(hbands).shape[0]-1),int(np.unique(hbands).shape[0]-1)),dtype=np.float32) #laura, nhrus replaced with hbands
 cmatrix = cmatrix.tocsr()

 #Prepare length, width, and ksat matrices
 wmatrix = cmatrix.copy()
 wmatrix.multiply(dx) #wmatrix[:] = dx*wmatrix[:]

 #Prepare output dictionary
 cdata = {'width':wmatrix.T,}

 return cdata

@numba.jit(nopython=True,cache=True)
def Calculate_HRU_Connections_Matrix_HMC_workhorse(cluster_ids,dx,tile_position,basins,
    ivc,irc,ibc): #laura, removed parameter nhrus 

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
   if h1 == -9999:continue
   #up
   if (i+1) < cluster_ids.shape[0]:
    h2 = cluster_ids[i+1,j]
    b2 = basins[i+1,j]
    tp2 = tile_position[i+1,j]
    if Determine_HMC_Connectivity(h1,h2,b1,b2,tp1,tp2,ivc,irc,ibc):
     horg.append(h1)
     hdst.append(h2)
   #down
   if (i-1) > 0:
    h2 = cluster_ids[i-1,j]
    b2 = basins[i-1,j]
    tp2 = tile_position[i-1,j]
    if Determine_HMC_Connectivity(h1,h2,b1,b2,tp1,tp2,ivc,irc,ibc):
     horg.append(h1)
     hdst.append(h2)
   #left
   if (j-1) > 0:
    h2 = cluster_ids[i,j-1]
    b2 = basins[i,j-1]
    tp2 = tile_position[i,j-1]
    if Determine_HMC_Connectivity(h1,h2,b1,b2,tp1,tp2,ivc,irc,ibc):
     horg.append(h1)
     hdst.append(cluster_ids[i,j-1])
   #right
   if (j+1) < cluster_ids.shape[1]:
    h2 = cluster_ids[i,j+1]
    b2 = basins[i,j+1]
    tp2 = tile_position[i,j+1]
    if Determine_HMC_Connectivity(h1,h2,b1,b2,tp1,tp2,ivc,irc,ibc):
     horg.append(h1)
     hdst.append(cluster_ids[i,j+1])
 horg = np.array(horg)
 hdst = np.array(hdst)

 return (hdst,horg)

def Create_and_Curate_Covariates_svp(wbd,hydroblocks_info):

 covariates = {}
 depths={} #laura svp
 #Read in and curate all the covariates
 #Read in soil vertical properties
 for file in wbd['files']: #laura svp
  if file in ['WLTSMC','MAXSMC','BB','DRYSMC','QTZ','SATDW','REFSMC','SATPSI','SATDK']: #laura svp
   covariates[file]={} #laura svp
   d=[] #laura svp
   for layer in range(0,len(wbd['files'][file])): #laura svp
    d.append(float((wbd['files'][file][layer].split('latlon_')[1]).split('cm')[0])) #laura svp
    covariates[file]['d_%scm'%((wbd['files'][file][layer].split('latlon_')[1]).split('cm')[0])]=gdal_tools.read_data(wbd['files'][file][layer]).data #laura svp
   depths[file]=d #laura svp
  else: #laura svp 
   if os.path.isfile(wbd['files'][file]): 
    covariates[file] = gdal_tools.read_data(wbd['files'][file]).data

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
 mask = np.copy(covariates['mask']).astype(np.int64)
 mask_all = np.copy(mask)
 mask[mask != hydroblocks_info['cid']] = 0
 mask = mask.astype(np.bool)
 
 #Set all nans to the mean
 for var in covariates:
  if var in hydroblocks_info['hmc_parameters']['subbasin_clustering_covariates']:continue
  if var in ['WLTSMC','MAXSMC','BB','DRYSMC','QTZ','SATDW','REFSMC','SATPSI','SATDK']: #laura svp
   for depth in covariates[var]: #laura svp
    mask1 = (np.isinf(covariates[var][depth]) == 0) & (np.isnan(covariates[var][depth]) == 0) #laura svp
    mask0 = (np.isinf(covariates[var][depth]) == 1) | (np.isnan(covariates[var][depth]) == 1) #laura svp
    covariates[var][depth][mask0] = -9999.0 #laura svp
  else: #laura svp
   #covariates[var][mask <= 0] = -9999.0
   mask1 = (np.isinf(covariates[var]) == 0) & (np.isnan(covariates[var]) == 0) 
   mask0 = (np.isinf(covariates[var]) == 1) | (np.isnan(covariates[var]) == 1)
   covariates[var][mask0] = -9999.0# stats.mode(covariates[var][mask1])[0][0]

 #Set everything that is -9999 to the mean
 for var in covariates:
  if var in ['WLTSMC','MAXSMC','BB','DRYSMC','QTZ','SATDW','REFSMC','SATPSI','SATDK']: #laura svp
   for depth in covariates[var]: #laura svp
    m2 = ( mask > 0 ) & (covariates[var][depth] != -9999.0) #laura svp
    missing_ratio = 1.0 - np.sum(m2)/float(np.sum(mask)) #laura svp
    if missing_ratio > 0.99 : #laura svp
     print("Warning: Covariate %s in layer %s in catchment %s has %.2f %% of nan's" % (var,depth,hydroblocks_info['cid'],100*missing_ratio)) # laura svp
    if var not in ['mask',]:
     covariates[var][depth][covariates[var][depth] == -9999.0] = np.mean(covariates[var][depth][covariates[var][depth] != -9999.0])

  else:
   m2 = ( mask > 0 ) & (covariates[var] != -9999.0)
   missing_ratio = 1.0 - np.sum(m2)/float(np.sum(mask))
   if missing_ratio > 0.99 : 
    print("Warning: Covariate %s in catchment %s has %.2f %% of nan's" % (var,hydroblocks_info['cid'],100*missing_ratio)) # Noemi insert
    if var == 'lc': 
     mlc = (covariates[var] == -9999) & mask
     covariates[var][mlc] = 17  # Water
    if var in ['dem','fdir','sand','clay','silt','TEXTURE_CLASS','dbedrock']:
     exit('Error_clustering: %s_full_of_nans %s' % (var,hydroblocks_info['cid']))
   if var not in ['mask',]:
    if var in ['nlcd','TEXTURE_CLASS','lc','irrig_land','bare30','water30','tree30','start_growing_season','end_growing_season']: 
     covariates[var][covariates[var] == -9999.0] = stats.mode(covariates[var][covariates[var] != -9999.0])[0][0]
    else:
     covariates[var][covariates[var] == -9999.0] = np.mean(covariates[var][covariates[var] != -9999.0])

 #Set everything outside of the mask to -9999
 for var in covariates:
  if var in hydroblocks_info['hmc_parameters']['subbasin_clustering_covariates']:continue
  if var in ['dem','fdir','acc']:continue 
  if var in ['WLTSMC','MAXSMC','BB','DRYSMC','QTZ','SATDW','REFSMC','SATPSI','SATDK']: #laura svp
   for depth in covariates[var]: #laura svp
    covariates[var][depth][mask<0]=-9999.0 #laura svp
  else: #laura svp
   covariates[var][mask <= 0] = -9999.0

 #Add the mask_all to the covariates
 covariates['mask_all'] = np.copy(mask_all)
 
 return (covariates,mask,depths) #laura svp returns depths for dataset svp

def Create_Clusters_And_Connections(workspace,wbd,output,input_dir,nhru,info,hydroblocks_info):
 
 dz=hydroblocks_info['dz'] #laura svp
 #Retrieve some metadata
 metadata = gdal_tools.retrieve_metadata(wbd['files']['mask'])
 resx = 90.0#670.0**0.5#26.0

 print("Creating and curating the covariates",flush=True)
 (covariates,mask,z_data)=Create_and_Curate_Covariates_svp(wbd,hydroblocks_info)
 
 #Determine the HRUs (clustering if semidistributed; grid cell if fully distributed)
 print("Computing the HRUs",flush=True)
 (cluster_ids,nhru,new_hand,HMC_info,covariates,dbc,hand,basins,basin_clusters,hand_org,hbands,area_adj,tile_position) = Compute_HRUs_Semidistributed_HMC(covariates,mask,hydroblocks_info,wbd,resx,input_dir)
 #covariates['hand'] = new_hand
 covariates['hand'] = hand
 hydroblocks_info['nhru'] = nhru
  
 #Create the netcdf file
 file_netcdf = '%s/input_file.nc' % hydroblocks_info['input_dir']#hydroblocks_info['input_file']
 hydroblocks_info['input_fp'] = nc.Dataset(file_netcdf, 'w', format='NETCDF4')

 #Create the dimensions (netcdf)
 idate = hydroblocks_info['idate']
 fdate = hydroblocks_info['fdate']
 dt = hydroblocks_info['dt']
 ntime = 24*3600*((fdate - idate).days+1)/dt
 nhru = hydroblocks_info['nhru']
 nsoil= len(hydroblocks_info['dz']) #laura svp
 hydroblocks_info['input_fp'].createDimension('hru',nhru)
 hydroblocks_info['input_fp'].createDimension('time',ntime)
 hydroblocks_info['input_fp'].createDimension('nsoil',nsoil) #laura svp 

 #Create the groups (netcdf)
 hydroblocks_info['input_fp'].createGroup('meteorology')
 hydroblocks_info['input_fp'].createGroup('water_use')

 #Prepare the hru connections matrix (darcy clusters) with laura's modification
 print("Calculating the connections between HRUs",flush=True)
 if (hydroblocks_info['connection_matrix_hbands']==False):
  cmatrix = Calculate_HRU_Connections_Matrix_HMC(covariates,cluster_ids,nhru,resx,HMC_info,hydroblocks_info)
  #Define the metadata
  metadata = gdal_tools.retrieve_metadata(wbd['files']['dem'])
  #Make the output dictionary for the basin
  OUTPUT = {'hru':{},'metadata':metadata,'mask':mask,'cmatrix':cmatrix}
 
 else:
  #Define the metadata
  metadata = gdal_tools.retrieve_metadata(wbd['files']['dem'])
  #Make the output dictionary for the basin
  OUTPUT = {'hru':{},'metadata':metadata,'mask':mask}
  #Create connection matrix per cluster of watersheds, laura
  bcu=np.unique(basin_clusters)
  bcu=bcu[bcu>0]
  for bc in bcu:
   masked_hband=np.empty(hbands.shape)
   if bc==1:
    masked_hband[basin_clusters==int(bc)]=hbands[basin_clusters==bc]
   else:
    masked_hband[basin_clusters==int(bc)]=hbands[basin_clusters==bc]
    masked_hband=masked_hband-(np.min(hbands[basin_clusters==bc]))
   masked_hband[~(basin_clusters==int(bc))]=int(-9999)
   group_name='cmatrix_Basin%s' %bc
   shape=int(((np.unique(masked_hband)).shape[0])-1)
   cmatrix=np.empty([shape,shape])
   cmatrix=Calculate_HRU_Connections_Matrix_HMC_hbands(masked_hband,resx,HMC_info,hydroblocks_info) #laura: removed covariates from the function
   OUTPUT[group_name]=cmatrix #end of laura's modification
 
 #Remember the map of hrus
 OUTPUT['hru_map'] = cluster_ids
 OUTPUT['channel_map'] = HMC_info['channel_map']
 OUTPUT['hand_map'] = hand
 OUTPUT['basin_map'] = basins
 OUTPUT['basin_clusters_map'] = basin_clusters
 OUTPUT['hand_org_map'] = hand_org
 OUTPUT['hband_map'] = hbands

 #Assign the model parameters
 print("Assigning the model parameters",flush=True)
 #Acumulate soil depths and convert to meters laura svp
 dz2=[]
 for elt in dz:
  if len(dz2)>0:
   dz2.append(dz2[-1]+elt)
  else:
   dz2.append(elt)
 z_model=np.array(dz2)*100

 OUTPUT = Assign_Parameters_Semidistributed_svp(covariates,metadata,hydroblocks_info,OUTPUT,cluster_ids,mask,hbands,area_adj,z_data,z_model) #laura svp

 #Add the new number of clusters
 OUTPUT['nhru'] = nhru
 OUTPUT['mask'] = mask
 OUTPUT['stream_network'] = dbc

 return (OUTPUT,covariates)

def Prepare_Meteorology_Semidistributed(workspace,wbd,OUTPUT,input_dir,info,hydroblocks_info,covariates):

 #Define the mapping directory
 mapping_info = {}
 #Calculate the fine to coarse scale mapping
 for data_var in wbd['files_meteorology']:
  
  #Define the variable name
  var = data_var#data_var.split('_')[1]
  mapping_info[var] = {}

  #Read in the coarse and fine mapping
  file_coarse = '%s/%s_latlon_coarse.tif' % (workspace,data_var)
  file_fine = '%s/%s_latlon_fine.tif' % (workspace,data_var)
  mask_coarse = gdal_tools.read_raster(file_coarse)
  mask_fine = gdal_tools.read_raster(file_fine)
  nlat = mask_coarse.shape[0]
  nlon = mask_coarse.shape[1]

  #Compute the mapping for each hru
  for hru in np.arange(hydroblocks_info['nhru']):
   idx = OUTPUT['hru_map'] == hru
   icells = np.unique(mask_fine[idx][mask_fine[idx] != -9999.0].astype(np.int))   # Add != -9999 for unique and bicount - Noemi
   counts = np.bincount(mask_fine[idx][mask_fine[idx] != -9999.0].astype(np.int))
   coords,pcts,dem_coarse = [],[],[] #dem for downscaling
   for icell in icells:
    ilat = int(np.floor(icell/mask_coarse.shape[1]))
    jlat = icell - ilat*mask_coarse.shape[1]
    pct = float(counts[icell])/float(np.sum(counts))
    coords.append([ilat,jlat])
    pcts.append(pct)
    if var == 'tair':dem_coarse.append(np.mean(covariates['dem'][mask_fine == icell]))
   pcts = np.array(pcts)
   coords = list(np.array(coords).T)
   if var == 'tair':
    dem_fine = np.mean(covariates['dem'][idx])
    dem_coarse = np.array(dem_coarse)
    mapping_info[var][hru] = {'pcts':pcts,'coords':coords,'dem_coarse':dem_coarse,'dem_fine':dem_fine}
   else:
    mapping_info[var][hru] = {'pcts':pcts,'coords':coords}

 #Iterate through variable creating forcing product per HSU
 #R
 idate = info['time_info']['startdate']
 fdate = info['time_info']['enddate']
 dt = info['time_info']['dt']
 nt = int(3600*24/dt)*((fdate - idate).days+1)
 #Create structured array
 meteorology = {}
 for data_var in wbd['files_meteorology']:
  meteorology[data_var] = np.zeros((nt,hydroblocks_info['nhru']))
 #Load data into structured array
 db_data = {}
 for data_var in wbd['files_meteorology']:
  var = data_var#data_var.split('_')[1]
  date = idate
  file = wbd['files_meteorology'][data_var]
  fp = nc.Dataset(file)
  
  #Determine the time steps to retrieve
  nc_step = int(60*float(fp.variables['t'].units.split(' ')[0].split('h')[0]))
  nc_idate = np.array(fp.variables['t'].units.split(' ')[2].split('-'))
  nc_nt = len(fp.variables['t'][:])
  dates = [datetime.datetime(int(nc_idate[0]),int(nc_idate[1]),int(nc_idate[2]))]
  #for it in range(1,nc_nt): dates.append(dates[0] + datetime.timedelta(hours=it*nc_step))
  for it in range(1,nc_nt): dates.append(dates[0] + datetime.timedelta(minutes=it*nc_step))
  dates=np.array(dates)
  startdate = info['time_info']['startdate']
  enddate = info['time_info']['enddate']
  mask_dates = (dates >= startdate) & (dates <= enddate)
  db_data[var] = np.ma.getdata(fp.variables[var][mask_dates,:,:])
  fp.close()
 
 #Downscale the variables
 flag_downscale = False
 if flag_downscale == True:db_downscaled_data = Downscale_Meteorology(db_data,mapping_info)

 #Finalize data
 for var in db_data:
  for hru in mapping_info[var]:
   pcts = mapping_info[var][hru]['pcts']
   if flag_downscale == False:
    coords = mapping_info[var][hru]['coords']
    coords[0][coords[0] >= db_data[var].shape[1]] = db_data[var].shape[1] - 1
    coords[1][coords[1] >= db_data[var].shape[2]] = db_data[var].shape[2] - 1
    tmp = db_data[var][:,coords[0],coords[1]]
   else:
    tmp = db_downscaled_data[hru][var]
   tmp = pcts*tmp
   meteorology[var][:,hru] = np.sum(tmp,axis=1)

  #Write the meteorology to the netcdf file (single chunk for now...)
  grp = hydroblocks_info['input_fp'].groups['meteorology']
  grp.createVariable(var,'f4',('time','hru'))#,zlib=True)
  grp.variables[var][:] = meteorology[var][:]

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

def Downscale_Meteorology(db_data,mapping_info):
 
 #Iterate per hru
 db_org = {}
 db_ds = {}
 for hru in mapping_info['tair']:
  db_org[hru] = {}
  db_ds[hru] = {}
  #Collect the data
  for var in db_data:
   pcts = mapping_info[var][hru]['pcts']
   coords = mapping_info[var][hru]['coords']
   coords[0][coords[0] >= db_data[var].shape[1]] = db_data[var].shape[1] - 1
   coords[1][coords[1] >= db_data[var].shape[2]] = db_data[var].shape[2] - 1
   db_org[hru][var] = db_data[var][:,coords[0],coords[1]]
  df = mapping_info['tair'][hru]['dem_fine']
  dc = mapping_info['tair'][hru]['dem_coarse']
  #A.Downscale temperature
  dT = -6.0*10**-3*(df - dc)
  db_ds[hru]['tair'] = dT[np.newaxis,:] + db_org[hru]['tair']
  #db_ds[hru]['tair'] = db_org[hru]['tair'][:]
  #B.Downscale longwave
  #0.Compute radiative temperature 
  sigma = 5.67*10**-8
  emis = 1.0
  trad = (db_org[hru]['lwdown']/sigma/emis)**0.25
  #1.Apply lapse rate to trad
  trad = dT[np.newaxis,:] + trad
  #2.Compute longwave with new radiative tempearture
  db_ds[hru]['lwdown'] = emis*sigma*trad**4
  #db_ds[hru]['lwdown'] = db_org[hru]['lwdown'][:]
  #C.Downscale pressure
  psurf = db_org[hru]['psurf'][:]*np.exp(-10**-3*(df-dc)/7.2)
  db_ds[hru]['psurf'] = psurf[:]
  #D.Downscale specific humidity
  #db_ds[hru]['spfh'] = db_org[hru]['spfh'][:]
  #Convert to vapor pressure
  e = db_org[hru]['psurf'][:]*db_org[hru]['spfh'][:]/0.622 #Pa
  esat = 1000*saturated_vapor_pressure(db_org[hru]['tair'][:] - 273.15) #Pa
  rh = e/esat
  esat = 1000*saturated_vapor_pressure(db_ds[hru]['tair'][:] - 273.15) #Pa
  e = rh*esat
  q = 0.622*e/db_ds[hru]['psurf']
  db_ds[hru]['spfh'] = q[:]
  #E.Downscale shortwave radiation
  db_ds[hru]['swdown'] = db_org[hru]['swdown'][:]
  #F.Downscale wind speed
  db_ds[hru]['wind'] = db_org[hru]['wind'][:]
  #G.Downscale precipitation
  db_ds[hru]['precip'] = db_org[hru]['precip'][:]

 return db_ds

def saturated_vapor_pressure(T):
    es = 0.6112*np.exp(17.67*T/(T + 243.5))
    return es

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
  
  wuse_lc_ea_file = '%s/%s_lc_ea.tif' % (input_dir,data_var)
  gdal_tools.write_raster(wuse_lc_ea_file,md,hrus_lc)
  fine_size = hrus_lc.shape
  #fine_res = abs(md['resx']) #NEED TO UPDATE
  
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
  file_out = '%s/%s_area_latlon_coarse.tif' % (input_dir,data_var)
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
  file_out = '%s/%s_area_latlon_coarse.tif' % (input_dir,data_var)
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
  '''grp = hydroblocks_info['input_fp'].groups['water_use']
  grp.createVariable(var,'f4',('time','hru'))#,zlib=True)
  grp.variables[data_var][:] = water_use[data_var][:]'''

 '''if hydroblocks_info['water_management']['hwu_flag']:
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
   var[:] = dates[:]'''

 return

def driver(comm,metadata_file):

 size = comm.Get_size()
 rank = comm.Get_rank()
 #Read in the metadata
 #metadata_file = '%s/metadata.json' % edir
 metadata = Read_Metadata_File(metadata_file)
 info = metadata
 info['covariates'] = {'lats':'n','lons':'n','lc':'n'}
 info['idate'] = datetime.datetime(metadata['startdate']['year'],
                           metadata['startdate']['month'],
                           metadata['startdate']['day'],0)
 info['fdate'] = datetime.datetime(metadata['enddate']['year'],
                           metadata['enddate']['month'],
                           metadata['enddate']['day'],0) + datetime.timedelta(days=1) - datetime.timedelta(seconds=info['dt'])
 rdir = metadata['rdir']
 edir = '%s/experiments/simulations/%s' % (rdir,metadata['experiment'])
 #Split up the processing across cores
 dfile = '%s/data/shp/domain.shp' % rdir
 fp = fiona.open(dfile,'r')
 cids = np.array(range(1,len(list(fp))+1))
 fp.close()
 for cid in cids[rank::size]:
  print(rank,size,cid)
  metadata['cid'] = cid
  metadata['input_dir'] = "%s/%d" % (edir,cid)
  metadata['workspace'] = "%s/data/cids/%d" % (rdir,cid)
  #Prepare model data
  tic = time.time()
  Prepare_Model_Input_Data(metadata,metadata_file)
  print("Elapsed time: ",time.time() - tic)
 comm.Barrier()

 #Create enhanced input data file
 #Connect_Cell_Networks(rank,size,cids,edir)
 Connect_Cell_Networks_v2(rank,size,cids,edir)
 comm.Barrier()

 #Create downstream channel database for particle tracker routing scheme
 Create_Downstream_Channels_Database(edir,rank,size,cids,comm)
 comm.Barrier()

 #Wait until they are all done
 workspace = '%s/workspace' % (edir)
 os.system('mkdir -p %s' % workspace)
 Finalize_River_Network_Database(rdir,edir,cids,workspace,comm,rank,size)
 comm.Barrier()
 
 #Create the soft link to the NoahMP code from the directory
 #os.system('ln -s 

 #Postprocess the model input 
 if rank == 0:Postprocess_Input(rdir,edir,cids)

 return

def Postprocess_Input(rdir,edir,cids):

 sdir = '%s/postprocess' % (edir)
 ddir = '%s/data/cids' % rdir
 os.system('rm -rf %s' % sdir)
 #Create cid, hru, and channel maps
 vars = ['cids','cids_org','dem','hrus','channels','hand','basins','basin_clusters']
 for var in vars:
  os.system('mkdir -p %s/postprocess/%s' % (edir,var))
 for cid in cids:
  print('Copying files for vrt',cid,flush=True)
  dir = '%s/%s' % (edir,cid)
  #hru
  ifile = '%s/hru_mapping_latlon.tif' % dir
  ofile = '%s/hrus/%d.tif' % (sdir,cid)
  os.system('cp %s %s' % (ifile,ofile))
  #channels
  ifile = '%s/channel_mapping_latlon.tif' % dir
  ofile = '%s/channels/%d.tif' % (sdir,cid)
  os.system('cp %s %s' % (ifile,ofile))
  #cid
  ifile = '%s/%d/mask_latlon.tif' % (ddir,cid)
  ofile = '%s/cids/%d.tif' % (sdir,cid)
  fpi = rasterio.open(ifile)
  profile = fpi.profile
  data = fpi.read(1)
  data[data!=cid] = -9999.0
  fpo = rasterio.open(ofile,'w',**profile)
  fpo.write(data,1)
  fpi.close()
  fpo.close()
  #cid
  ifile = '%s/%d/mask_org_latlon.tif' % (ddir,cid)
  ofile = '%s/cids_org/%d.tif' % (sdir,cid)
  os.system('cp %s %s' % (ifile,ofile))
  #dem
  ifile = '%s/%d/dem_latlon.tif' % (ddir,cid)
  ofile = '%s/dem/%d.tif' % (sdir,cid)
  os.system('cp %s %s' % (ifile,ofile))
  #hand
  ifile = '%s/hand_latlon.tif' % dir
  ofile = '%s/hand/%d.tif' % (sdir,cid)
  os.system('cp %s %s' % (ifile,ofile))
  #basins
  ifile = '%s/basins_latlon.tif' % dir
  ofile = '%s/basins/%d.tif' % (sdir,cid)
  os.system('cp %s %s' % (ifile,ofile))
  #basin clusters
  ifile = '%s/basin_clusters_latlon.tif' % dir
  ofile = '%s/basin_clusters/%d.tif' % (sdir,cid)
  os.system('cp %s %s' % (ifile,ofile))

 #Create vrts
 for var in vars:
  print("creating virtual raster: %s" % var,flush=True)
  os.system('gdalbuildvrt %s/%s.vrt %s/%s/*.tif' % (sdir,var,sdir,var))

 #Create shapefiles
 #os.system('gdal_polygonize.py -f "ESRI Shapefile" -8 %s/basins.vrt %s/basins_shp' % (sdir,sdir))
 #os.system('gdal_polygonize.py -f "ESRI Shapefile" -8 %s/basin_clusters.vrt %s/basin_clusters_shp' % (sdir,sdir))
 #os.system('gdal_polygonize.py -f "ESRI Shapefile" -8 %s/cids.vrt %s/cids_shp' % (sdir,sdir))

 return

def Read_Metadata_File(file):

 import json
 metadata = json.load(open(file))

 return metadata

def Connect_Cell_Networks_v2(rank,size,cids,edir):

 for cid in cids[rank::size]:

  #Change to integer
  cid1 = int(cid)

  #Read in the routing interconnectivitiy dictionary for the given cid
  file = '%s/%s/routing_mp_connectivity.pck' % (edir,cid1)
  db = pickle.load(open(file,'rb'))
    
  #Open input_file.nc for cid in append mode
  file = '%s/%s/input_file.nc' % (edir,cid1)
  fp = h5py.File(file,'a') 

  #Iterate through the outlets to determine the channel id in the target subdomain 
  db2 = {}
 
  #Create the output array
  output_array = -9999*np.ones((db['channel_target_mp'].size,4),dtype=np.int32)
  output_array[:,0] = cid1
  output_array[:,1] = db['channel_outlet_id'][:]
  output_array[:,2] = db['channel_target_mp'][:]
  for ic in range(db['channel_outlet_id'].size):
    cid2 = db['channel_target_mp'][ic]
    if cid2 != -9999:
     #Read in the database for the target cid
     if cid2 not in db2:
      file2 = '%s/%d/routing_mp_connectivity.pck' % (edir,cid2)
      db2[cid2] = pickle.load(open(file2,'rb'))
     #channel_#Determine the channel lat/lon that is closest to create a link 
     lat1 = db['channel_target_crds'][ic][0]
     lon1 = db['channel_target_crds'][ic][1]
     lats2 = db2[cid2]['channel_crds'][:,:,0]
     lons2 = db2[cid2]['channel_crds'][:,:,1]
     dist = ((lats2-lat1)**2 + (lons2-lon1)**2)**0.5
     icd = np.where(dist == np.min(dist))[0][0]
     output_array[ic,3] = icd
        
  #Add array to file
  fp['stream_network']['outlets'] = output_array[:] #out cid, out channel id, in cid, in channel id

  #Iterate through the outlets to determine the channel id in the target subdomain 
  db2 = {}

  #Create the output array
  inlet_array = -9999*np.ones((db['channel_inlet_id'].size,10),dtype=np.int32)
  inlet_array[:,0] = cid1
  inlet_array[:,1] = db['channel_inlet_id'][:]
  inlet_array[:,2:6] = db['channel_inlet_target_mp'][:]
  for ic in range(db['channel_inlet_id'].size):
   for j in range(db['channel_inlet_target_mp'].shape[1]):
    if db['channel_inlet_target_mp'][ic,j] == -9999:break
    cid2 = db['channel_inlet_target_mp'][ic,j]
    #Read in the database for the target cid
    if cid2 not in db2:
     file2 = '%s/%d/routing_mp_connectivity.pck' % (edir,cid2)
     db2[cid2] = pickle.load(open(file2,'rb'))
    #channel_#Determine the channel lat/lon that is closest to create a link 
    lat1 = db['channel_inlet_target_crds'][ic,j,0]
    lon1 = db['channel_inlet_target_crds'][ic,j,1]
    lats2 = db2[cid2]['channel_crds'][:,:,0]
    lons2 = db2[cid2]['channel_crds'][:,:,1]
    dist = ((lats2-lat1)**2 + (lons2-lon1)**2)**0.5
    icd = np.where(dist == np.min(dist))[0][0]
    inlet_array[ic,6+j] = icd 

  #Add array to file
  fp['stream_network']['inlets'] = inlet_array[:]
        
  #Close ammended file
  fp.close()
 
 return

def create_enhanced_topology(topology,outlets,cid):
    
 #create enhanced topology by adding outlet information (channel id and cid)
 topology_enhanced = -1*np.ones((topology.size,2),dtype=np.int32)
 topology_enhanced[:,0] = topology
 topology_enhanced[topology != -1,1] = cid
 topology_enhanced[outlets[:,1],0] = outlets[:,3] 
 topology_enhanced[outlets[:,1],1] = outlets[:,2]
 topology_enhanced[topology_enhanced==-9999] = -1

 return topology_enhanced

def read_channel_database(cid,edir):
    
 db = {}

 #Open input_file.nc for cid in append mode
 file = '%s/%s/input_file.nc' % (edir,cid)
 fp = nc.Dataset(file,'r')
 
 #create enhanced topology by adding outlet information (channel id and cid)
 db['topology_enhanced'] = create_enhanced_topology(fp['stream_network']['topology'][:],fp['stream_network']['outlets'][:],cid)
 db['length'] = fp['stream_network']['length'][:]

 fp.close()
 
 return db

def Create_Downstream_Channels_Database(edir,rank,size,cids,comm):

 #Determine the total distance that can be covered
 dt = 3600 #sec #This should be defined by the dt_routing parameter
 maxu = 2 #m/s #parameter
 maxd = maxu*dt #m
 ncmax = 20 #parameter

 #Iterate per catchment
 for cid in cids[rank::size]:

  #Change to integer
  cid = int(cid)

  #Initialize dictionary where information will be held
  db = {cid:{}}
    
  #create enhanced topology by adding outlet information (channel id and cid)
  db[cid] = read_channel_database(cid,edir)

  #Initialize downstream channel array (channel id, cid)
  downstream_channels = -9999*np.ones((db[cid]['topology_enhanced'].shape[0],2,ncmax),dtype=np.int32)
                                   
  #Iterate through each channel
  for ic in range(db[cid]['topology_enhanced'].shape[0]):
   d = maxd
   cid0 = cid
   ic0 = ic
   count = -1
   while (d > 0) & (count < ncmax):
    count += 1
    ic1 = db[cid0]['topology_enhanced'][ic0,0]
    cid1 = db[cid0]['topology_enhanced'][ic0,1]
    if ic1 == -1:
     downstream_channels[ic,0,count] = -1
     downstream_channels[ic,1,count] = -1
     break
    if cid1 not in db:
     db[cid1]=read_channel_database(cid1,edir)
    #save information
    downstream_channels[ic,0,count] = ic1
    downstream_channels[ic,1,count] = cid1
    #subtract distance
    d = d - db[cid1]['length'][ic1]
    #update ids
    cid0 = cid1
    ic0 = ic1

 comm.Barrier()
 #Add downstream_channels array to input_file.nc
 file = '%s/%s/input_file.nc' % (edir,cid)
 fp = h5py.File(file,'a')
 fp['stream_network']['downstream_channels'] = downstream_channels[:]
 fp.close()

 return

def Finalize_River_Network_Database(rdir,edir,cids,workspace,comm,rank,size):

 core_cids = np.array(cids)[rank::size]
 debug_level = 0

 #Prepare data for cids
 comm.Barrier()
 for cid in core_cids:
  if debug_level >= 0:print(cid,"Assembling the input/output",flush=True)
  db = prepare_data(rank,cid,edir,debug_level,workspace,cids)
  cdir = '%s/%d' % (edir,cid)
  pickle.dump(db,open('%s/octopy.pck' % (cdir,),'wb'),pickle.HIGHEST_PROTOCOL)

 return

def prepare_data(rank,cid,edir,debug_level,workspace,cids):

 #Read in the stream network information
 if debug_level >= 1:print(rank,cid,"Reading in the stream network information",flush=True)
 file = '%s/%d/input_file.nc' % (edir,cid)
 fp = nc.Dataset(file)
 nhband = np.unique(fp['parameters']['hband'][:]).size
 grp = fp['stream_network']
 dbc = {}
 for var in grp.variables:
  dbc[var] = grp[var][:]
 fp.close()

 #Read in reach/height band area
 if debug_level >= 1:print(rank,cid,"Reading in reach/height band relationship",flush=True)
 file = '%s/%d/routing_info.pck' % (edir,cid)
 db = pickle.load(open(file,'rb'))['reach_hband_area']

 #Read in the reach/hand database
 if debug_level >= 1:print(rank,cid,"Read reach/hand database for current cell",flush=True)
 file = '%s/%d/routing_info.pck' % (edir,cid)
 #file = '%s/input_data/domain/%d/routing_info.pck' % (rdir,cid)
 hdb = pickle.load(open(file,'rb'))['reach_cross_section']

 #Read in the unit hydrograph per height band database
 if debug_level >= 1:print(rank,cid,"Read unit hydrograph per height band database",flush=True)
 file = '%s/%d/routing_info.pck' % (edir,cid)
 #file = '%s/input_data/domain/%d/routing_info.pck' % (rdir,cid)
 uhs = pickle.load(open(file,'rb'))['uh_per_hband']

 #Create a reach2hband matrix that describes their relationship
 if debug_level >= 1:print(rank,cid,"Creating reach/hband matrix",flush=True)
 #HERE -> Need a nhband parameter
 #Compute hru average per height band -> Feed into routing (Need the mapping of hru to hband)
 #Apply hband inundation to all hru
 #reach2hband = np.zeros((np.sum(odbc[cid]['cid']==cid),nhband))
 reach2hband = np.zeros((dbc['topology'].size,nhband))
 for reach in db:
  for hband in db[reach]:
   tmp = db[reach][hband]
   if tmp == 0:print(reach,hband,db[reach][hband])
   reach2hband[reach-1,hband] = db[reach][hband]
 reach2hband = sparse.csr_matrix(reach2hband)

 #Initialize arrays
 c_length = dbc['length'][:]
 c_slope = dbc['slope'][:]
 c_width = dbc['width'][:]
 c_bankfull = dbc['bankfull'][:]
 c_n = dbc['manning_channel'][:]
 fp_n = dbc['manning_floodplain'][:]
 Ainit = np.zeros(c_length.size)
 Ainit[:] = 10**-5#0.1
 A0 = np.copy(Ainit)
 A1 = np.copy(Ainit)
 u0 = np.zeros(Ainit.size)
 bcs = np.zeros(c_length.size)
 Qinit = np.zeros(c_length.size)
 Q0 = Qinit[:]
 qin = np.zeros(c_length.size)
 qout = np.zeros(c_length.size)
 dA = np.zeros(c_length.size)

 #Initialize diagnostics
 tsolve = 0.0
 tcount = 0.0

 #Assemble database for simulation stage
 db = {
       'u0':copy.deepcopy(u0),
       'A0':copy.deepcopy(A0),
       'qin':copy.deepcopy(qin),
       'qout':copy.deepcopy(qout),
       'bcs':copy.deepcopy(bcs),
       'hdb':copy.deepcopy(hdb),
       'c_slope':copy.deepcopy(c_slope),
       'c_n':copy.deepcopy(c_n),
       'fp_n':copy.deepcopy(fp_n),
       'c_length':copy.deepcopy(c_length),
       'c_width':copy.deepcopy(c_width),
       'c_bankfull':copy.deepcopy(c_bankfull),
       'tsolve':tsolve,'tcount':tcount,
       'Q0':copy.deepcopy(Q0),
       'u0':copy.deepcopy(u0),
       'dA':copy.deepcopy(dA),
       'uhs':copy.deepcopy(uhs['data']),
       'uh_travel_time':copy.deepcopy(uhs['bins']),
       'reach2hband':copy.deepcopy(reach2hband)
      }

 return db
