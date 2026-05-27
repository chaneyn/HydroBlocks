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
import shutil
from pathlib import Path
import networkx as nx
from sklearn.decomposition import PCA

# Disable HDF5 file locking to avoid MPI file-lock issues on some clusters
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

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
 input_dir = hydroblocks_info['input_dir']
 os.system('mkdir -p %s' % input_dir)

 #Create soft link to HydroBlocks from within the directory
 HBdir = '%s/model/pyNoahMP' % (("/").join(__file__.split('/')[:-2]))
 HBedir = '%s/pyNoahMP%d' % (input_dir,hydroblocks_info['cid'])
 if os.path.exists(HBedir) == False:
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
 (output,covariates,z_data) = Create_Clusters_And_Connections(workspace,wbd,output,input_dir,nhru,info,hydroblocks_info)

 #Determine whether modified HMC / network abstraction behavior is requested
 network_abst_cfg = hydroblocks_info.get('network_abstraction')
 network_abst_flag = bool(network_abst_cfg.get('flag', False)) if isinstance(network_abst_cfg, dict) else False
 second_pass_flag = bool(hydroblocks_info.get('second_pass', False))

 #If network abstraction or modified HMC is active, persist first-pass covariates/output
 if network_abst_flag:
   try:
     cov_pickle = {'covariates': covariates, 'output': output, 'z_data': z_data}
     pickle.dump(cov_pickle, open('%s/covariates.pck' % input_dir, 'wb'))
   except Exception as e:
     print('Warning: could not save covariates.pck: %s' % str(e), flush=True)

 #Extract meteorological forcing immediately for one-pass runs, or during second pass in two-pass runs
 if (not network_abst_flag) or second_pass_flag:
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
 mask_object = gdal_tools.read_data(wbd['files']['mask'])
 terrain_tools.calculate_area(mask_object)
 grp.dx = np.mean(mask_object.area**0.5)

 #Write out the basin cluster map (defer when running two-pass abstraction/mod-HMC)

 #Write out the mapping and hand map
 if (not network_abst_flag) or second_pass_flag:
  hru_map = np.copy(output['hru_map'])
  hru_map[np.isnan(hru_map) == 1] = -9999.0
  file_ca = '%s/hru_mapping_latlon.tif' % input_dir
  metadata['nodata'] = -9999.0
  gdal_tools.write_raster(file_ca,metadata,hru_map)

  hand_map = np.copy(output['hand_map'])
  hand_map[np.isnan(hand_map) == 1] = -9999.0
  file_ca = '%s/hand_latlon.tif' % input_dir
  metadata['nodata'] = -9999.0
  gdal_tools.write_raster(file_ca,metadata,hand_map)
 else:
  print('Deferring write of hru_mapping_latlon.tif and hand_latlon.tif (abstraction/mod-HMC active)',flush=True)

 #Write out the basin map
 if (not network_abst_flag) or second_pass_flag:
  basin_map = np.copy(output['basin_map'])
  basin_map[np.isnan(basin_map) == 1] = -9999.0
  file_ca = '%s/basins_latlon.tif' % input_dir
  metadata['nodata'] = -9999.0
  gdal_tools.write_raster(file_ca,metadata,basin_map)
 else:
  print('Deferring write of basins_latlon.tif (abstraction/mod-HMC active)',flush=True)

 if (not network_abst_flag) or second_pass_flag:
  basin_clusters_map = np.copy(output['basin_clusters_map'])
  basin_clusters_map[np.isnan(basin_clusters_map) == 1] = -9999.0
  file_ca = '%s/basin_clusters_latlon.tif' % input_dir
  metadata['nodata'] = -9999.0
  gdal_tools.write_raster(file_ca,metadata,basin_clusters_map)
  n_cluster_basins = int(len(np.unique(basin_clusters_map)))
 else:
  print('Deferring write of basin_clusters_latlon.tif (abstraction/mod-HMC active)',flush=True)
  n_cluster_basins = int(hydroblocks_info['hmc_parameters']["number_of_characteristic_subbasins"]+1)

 #Write out the hand org and height-band maps
 if (not network_abst_flag) or second_pass_flag:
  hand_org_map = np.copy(output['hand_org_map'])
  hand_org_map[np.isnan(hand_org_map) == 1] = -9999.0
  file_ca = '%s/hand_org_latlon.tif' % input_dir
  metadata['nodata'] = -9999.0
  gdal_tools.write_raster(file_ca,metadata,hand_org_map)

  hband_map = np.copy(output['hband_map'])
  hband_map[np.isnan(hband_map) == 1] = -9999.0
  file_ca = '%s/hband_latlon.tif' % input_dir
  metadata['nodata'] = -9999.0
  gdal_tools.write_raster(file_ca,metadata,hband_map)
 else:
  print('Deferring write of hand_org_latlon.tif and hband_latlon.tif (abstraction/mod-HMC active)',flush=True)

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
  for i in range(1,n_cluster_basins):
   text='wmatrix_Basin%s' %int(i)
   print('Basins in cid %s: %s' %(cid, text),flush=True)
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
 #res_array = terrain_tools.calculate_area(mask_object)
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
 #ipoints = ((area > cthrs)).astype(np.int32)
 #ipoints[ipoints == 0] = -9999

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
 ###
 '''tcid = int(input_dir.split('/')[-1])
 if tcid == 1:
  for i in range(channel_inlet_id.size):
   print(tcid,channel_inlet_id[i],channel_inlet_target_mp[i,:])
 exit()'''
 
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
 #db_routing['i/o'] = terrain_tools.calculate_inlets_oulets(channels_wob,fdir,area_all,mask,np.flipud(covariates['lats']),covariates['lons'],mask_all,area_all)
 #print("got here 3",flush=True)
 #exit()
 #db_routing['i/o'] = terrain_tools.calculate_inlets_oulets(channels,fdir,area_all,mask,np.flipud(covariates['lats']),covariates['lons'],mask_all,area_all)

 #Compute and output the list of the channel positions
 '''lst_crds = []
 for icrd in range(crds.shape[0]):
   mcrd = crds[icrd,:,0] != -9999
   if (np.sum(mcrd) == 0):break
   crds_i = crds[icrd,mcrd,:]
   if crds_i.shape[0] > 1:
       lst_crds.append(shapely.geometry.LineString(np.fliplr(crds_i)))
   else:
       lst_crds.append(shapely.geometry.Point(np.flipud(crds_i[0,:])))
 db_routing['crds'] = lst_crds'''

 #Compute the basins
 print("Defining basins",flush=True)
 #basins = terrain_tools.ttf.delineate_basins(channels,m2,fdir)
 basins_wob = terrain_tools.ttf.delineate_basins(channels_wob,mask,fdir)
 basins = basins_wob
 
 #Compute channel properties
 db_channels = terrain_tools.calculate_channel_properties(channels_wob,channel_topology,slope,eares,mask,area_all,area_all_cp,basins_wob,hydroblocks_info['parameter_scaling'])

 #Compute Shreve order per macroscale polygon
 shreve = np.copy(channel_topology)
 shreve[:] = 0
 #Assign order 1 to all streams that are not in the topology target list
 for i in range(0, len(channel_topology)):
  if i not in np.unique(channel_topology):
   shreve[i] = 1
 #Propagate downstream until reaching domain outlet (topology == -1)
 for first_order in list(np.where(shreve == 1)[0]):
  shreve = go_downstream_shreve(first_order, channel_topology, shreve)
 db_channels['shreve'] = shreve

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

 flag_subgrid = hydroblocks_info.get('channel_initiation', {}).get('flag_subgrid', False)
 if flag_subgrid:
  thr_var = hydroblocks_info.get('channel_initiation', {}).get('var_pca', 'slope')
  cthrs_sg = hydroblocks_info.get('channel_initiation', {}).get('athrs_subgrid', cthrs)
  if cthrs_sg == cthrs:
   cthrs_sg = cthrs_sg - 50000
  dict_sg = {}
  (channels_sg,channels_wob_sg,channel_topology_sg,
   tmp1_sg,crds_sg,
   channel_outlet_id_sg,
   channel_target_mp_sg,
   channel_target_crds_sg,
   channel_inlet_id_sg,
   channel_inlet_target_mp_sg,
   channel_inlet_target_crds_sg) = terrain_tools.ttf.calculate_channels_wocean_wprop_wcrds(ac,ac_all,cthrs_sg,cthrs_sg,fdc,mask,mask_all,np.flipud(covariates['lats']),covariates['lons'])
  channel_topology_sg = channel_topology_sg[channel_topology_sg != -9999]
  dict_sg['channels_wob_sg'] = channels_wob_sg
  dict_sg['channel_topology_sg'] = channel_topology_sg
  basins_wob_sg = terrain_tools.ttf.delineate_basins(channels_wob_sg,mask,fdir)
  db_channels_sg = terrain_tools.calculate_channel_properties(channels_wob_sg,channel_topology_sg,slope,eares,mask,area_all,area_all_cp,basins_wob_sg,hydroblocks_info['parameter_scaling'])
  dict_sg['db_channels_sg'] = db_channels_sg
  principal_components = Subgrid_Indices(channels_wob_sg,channel_topology_sg,basins_wob,db_channels_sg,thr_var,hydroblocks_info['cid'])
  dict_sg['principal_components'] = principal_components
  pickle.dump(dict_sg,open('%s/pca_subgrid_basins.pck' % input_dir,'wb'))

 #Calculate the subbasin properties
 print("Assembling the subbasin properties",flush=True)
 vars1 = hydroblocks_info['hmc_parameters']['subbasin_clustering_covariates']
 vars = []
 for var in vars1:
  if var not in ['width','bankfull','length','area','shreve']:
   vars.append(var)
 hp_in = terrain_tools.calculate_basin_properties_updated(basins_wob,eares,covariates,vars)
 #sort hp_in (should go in geospatialtools)
 argsort = np.argsort(hp_in['bid'])
 for var in hp_in:
  hp_in[var] = hp_in[var][argsort]
 #bring in channel variables
 for var in ['width','bankfull','length','area','shreve']:
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

 if flag_subgrid:
  pca_path = '%s/pca_subgrid_basins.pck' % input_dir
  if os.path.isfile(pca_path):
   principal_components = pickle.load(open(pca_path,'rb'))['principal_components']
   for pc in range(0,principal_components.shape[1]):
    v = 'pc_%s' % (pc+1)
    covariates[v] = principal_components[:,pc]
    subbasin_clustering_cov.append(v)
  else:
   print('Warning: flag_subgrid is enabled but pca_subgrid_basins.pck is missing', flush=True)

  #Assemble input data
 cvs = {}
 for var in subbasin_clustering_cov: #laura
  if var in ['lc_w_now','lc_urb_nourb','lc_grass_forest']: #laura
   lc_mask=np.copy(dem)
   lc_mask[:] = 0.0
   if var=='lc_w_now': #laura
    if 'lc_17' in covariates:
     lc_mask=covariates['lc_17'] #laura
   elif var=='lc_urb_nourb': #laura
    if 'lc_13' in covariates:
     lc_mask=covariates['lc_13'] #laura
   elif var=='lc_grass_forest': #laura
    if 'lc_4' in covariates:
     lc_mask[covariates['lc_4']==1]=1 #deciduous forest
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

 keys = list(cvs.keys())
 if flag_subgrid:
  n_pc = 0
  n_no_pc = 0
  for var in keys:
   if 'pc_' in var:
    n_pc += 1
   else:
    n_no_pc += 1
  for var in keys:
   if 'pc_' in var:
    cvs[var]['w'] = (1/(n_no_pc+1))/n_pc
   else:
    cvs[var]['w'] = (1/(n_no_pc+1))
 else:
  for var in keys:
   cvs[var]['w'] = 1

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
   lc_mask=np.copy(dem)
   lc_mask[:] = 0.0
   if var=='lc_w_now':
    if 'lc_17' in covariates:
     lc_mask=covariates['lc_17']
   elif var=='lc_urb_nourb':
    if 'lc_13' in covariates:
     lc_mask=covariates['lc_13'] #laura
   elif var=='lc_grass_forest':
    if 'lc_4' in covariates:
     lc_mask[covariates['lc_4']==1]=1 #deciduous forest
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
 #pickle.dump(db_routing['i/o'],open('%s/routing_io.pck' % input_dir,'wb'))
 pickle.dump(db_routing['mp_connectivity'],open('%s/routing_mp_connectivity.pck' % input_dir,'wb'))

 #Construct HMC info for creating connections matrix
 HMC_info = {}
 HMC_info['basins'] = basins
 HMC_info['tile_position'] = tile_position
 HMC_info['channel_map'] = channels_wob

 #return (hrus.astype(np.float32),nhru,new_hand,HMC_info,covariates,db_channels,hand,
 return (hrus.astype(np.float32),nhru,new_hand,HMC_info,covariates,db_channels,new_hand2,
         basins,basin_clusters,hand,tiles,area_adj,tile_position)

def Compute_HRUs_Semidistributed_HMC2(hydroblocks_info,eares,input_dir):

 #Define the parameters for the hierarchical multivariate clustering
 ncatchments = hydroblocks_info['hmc_parameters']['number_of_characteristic_subbasins']
 ncatchments_main = hydroblocks_info['network_abstraction']['number_of_characteristic_main_subbasins']
 ncatchments_abst = hydroblocks_info['network_abstraction']['number_of_characteristic_secondary_subbasins']

 main_subbasin_clustering_cov = list(hydroblocks_info['network_abstraction']['main_subbasin_clustering_covariates'])
 abst_subbasin_clustering_cov = list(hydroblocks_info['network_abstraction']['secondary_subbasin_clustering_covariates'])
 subbasin_clustering_cov = list(hydroblocks_info['hmc_parameters']['subbasin_clustering_covariates'])

 dh = hydroblocks_info['hmc_parameters']['average_height_difference_between_bands']
 nclusters = hydroblocks_info['hmc_parameters']['number_of_intraband_clusters']
 vars1 = hydroblocks_info['hmc_parameters']['subbasin_clustering_covariates']
 vars2 = hydroblocks_info['network_abstraction']['main_subbasin_clustering_covariates']
 vars3 = hydroblocks_info['network_abstraction']['secondary_subbasin_clustering_covariates']
 vars4 = list(set(vars1 + vars2 + vars3))
 vars = []
 for var in vars4:
  if var not in ['width','bankfull','length','area','shreve','large_scale_basins']:
   vars.append(var)

 cid = int(hydroblocks_info['cid'])
 input_path = '%s/input_file.nc' % input_dir
 input_file = nc.Dataset(input_path)
 cov = pickle.load(open('%s/covariates.pck' % input_dir,'rb'))
 covariates = cov['covariates']
 z_data = cov.get('z_data', {})
 basins_wob = np.copy(cov['output']['basin_map'])
 if np.min(basins_wob[basins_wob!=-9999]) != 0:
  basins_wob[basins_wob!=-9999] = basins_wob[basins_wob!=-9999] - np.min(basins_wob[basins_wob!=-9999])
 basins = basins_wob
 db_channels = cov['output']['stream_network']
 channels_wob = cov['output']['channel_map']
 db_routing = pickle.load(open('%s/routing_info.pck' % input_dir,'rb'))
 dem = covariates['dem']
 mask = covariates['mask']
 fdir = terrain_tools.transform_arcgis_fdir(covariates['fdir'])

 hp_in = terrain_tools.calculate_basin_properties_updated(basins_wob,eares,covariates,vars)
 argsort = np.argsort(hp_in['bid'])
 for var in hp_in:
  hp_in[var] = hp_in[var][argsort]

 for var in ['width','bankfull','length','area','shreve','large_scale_basins']:
  if var in ['width','bankfull','length','area','shreve']:
   if len(db_channels['length']) < len(np.unique(basins_wob[basins_wob!=-9999])):
    hp_in[var] = [0]
   else:
    hp_in[var] = []
   hp_in[var].extend(input_file['stream_network'][var][:])
   hp_in[var] = np.array(hp_in[var])
  if var in ['large_scale_basins']:
    lsb_lats = []
    lsb_lons = []
    for basin_id in hp_in['bid']:
     basin_mask = basins_wob == basin_id
     lsb_lats.append(np.nanmean(covariates['lats'][basin_mask]))
     lsb_lons.append(np.nanmean(covariates['lons'][basin_mask]))
    hp_in['lsb_lats'] = np.array(lsb_lats)
    hp_in['lsb_lons'] = np.array(lsb_lons)

 if hydroblocks_info['network_abstraction']['flag'] == True: # 2-Step-HMC
  if len(db_channels['length']) < len(np.unique(basins_wob[basins_wob!=-9999])):
   m_main = np.array(np.unique(basins_wob[basins_wob!=-9999]),dtype=bool)
   m_main[:] = 0
   m_main[1:] = input_file['stream_network']['explicit_reach'][:] == 1
   m_abst = np.array(np.unique(basins_wob[basins_wob!=-9999]),dtype=bool)
   m_abst[:] = 1
   m_abst[1:] = input_file['stream_network']['explicit_reach'][:] == 0
  else:
   m_main = np.array(hp_in['shreve'],dtype=bool)
   m_main[:] = 0
   m_main[:] = input_file['stream_network']['explicit_reach'][:] == 1
   m_abst = np.array(hp_in['shreve'],dtype=bool)
   m_abst[:] = 1
   m_abst[:] = input_file['stream_network']['explicit_reach'][:] == 0
 input_file.close()

 if hydroblocks_info['channel_initiation']['flag_subgrid'] == True:
  pca_path = '%s/pca_subgrid_basins.pck' % input_dir
  if os.path.isfile(pca_path):
   y = pickle.load(open(pca_path,'rb'))['principal_components']
   for pc in range(0,y.shape[1]):
    v = 'pc_%s' % (pc+1)
    hp_in[v] = y[:,pc]
    subbasin_clustering_cov.append(v)
    main_subbasin_clustering_cov.append(v)
    abst_subbasin_clustering_cov.append(v)
  else:
   print('Warning: flag_subgrid is enabled but pca_subgrid_basins.pck is missing', flush=True)

 print('Clustering the basins',flush=True)
 if hydroblocks_info['network_abstraction']['flag'] == True: # 2-Step-HMC
  if (np.sum(m_main) != 0) and (np.sum(m_abst) != 0):
   hp_in_main = {}
   hp_in_abst = {}
   for key in list(hp_in.keys()):
    hp_in_main[key] = []
    hp_in_abst[key] = []
    for i in range(0,m_main.shape[0]):
     if m_main[i] == True:
      hp_in_main[key].append(hp_in[key][i])
     elif m_abst[i] == True:
      hp_in_abst[key].append(hp_in[key][i])
    hp_in_main[key] = np.array(hp_in_main[key])
    hp_in_abst[key] = np.array(hp_in_abst[key])
   basins_main = np.copy(basins_wob)
   basins_abst = np.copy(basins_wob)
   basins_main[:] = -9999
   basins_abst[:] = -9999
   for main in hp_in_main['bid']:
    basins_main[basins_wob == main] = 1
   for abst in hp_in_abst['bid']:
    basins_abst[basins_wob == abst] = 1
  elif (np.sum(m_main) != 0) and (np.sum(m_abst) == 0):
   hp_in_main = {}
   for key in list(hp_in.keys()):
    hp_in_main[key] = []
    for i in range(0,m_main.shape[0]):
     if m_main[i] == True:
      hp_in_main[key].append(hp_in[key][i])
    hp_in_main[key] = np.array(hp_in_main[key])
   basins_main = np.copy(basins_wob)
   basins_main[:] = -9999
   for main in hp_in_main['bid']:
    basins_main[basins_wob == main] = 1
  elif (np.sum(m_main) == 0) and (np.sum(m_abst) != 0):
   hp_in_abst = {}
   for key in list(hp_in.keys()):
    hp_in_abst[key] = []
    for i in range(0,m_abst.shape[0]):
     if m_abst[i] == True:
      hp_in_abst[key].append(hp_in[key][i])
    hp_in_abst[key] = np.array(hp_in_abst[key])
   basins_abst = np.copy(basins_wob)
   basins_abst[:] = -9999
   for abst in hp_in_abst['bid']:
    basins_abst[basins_wob == abst] = 1

  ncatchments_main = min(ncatchments_main,np.sum(m_main))
  ncatchments_abst = min(ncatchments_abst,np.sum(m_abst))
  ncatchments = ncatchments_main + ncatchments_abst

  if 'large_scale_basins' in main_subbasin_clustering_cov:
   main_subbasin_clustering_cov.remove('large_scale_basins')
   main_subbasin_clustering_cov = main_subbasin_clustering_cov + ['lsb_lats','lsb_lons']
  if 'large_scale_basins' in abst_subbasin_clustering_cov:
   abst_subbasin_clustering_cov.remove('large_scale_basins')
   abst_subbasin_clustering_cov = abst_subbasin_clustering_cov + ['lsb_lats','lsb_lons']

  if (np.sum(m_main) != 0):
   cvs1 = {}
   for var in main_subbasin_clustering_cov:
    if var in ['lc_w_now','lc_urb_nourb','lc_grass_forest']:
     lc_mask = np.copy(dem)
     lc_mask[:] = 0.0
     if var == 'lc_w_now':
      lc_mask = covariates['lc_17']
     elif var == 'lc_urb_nourb':
      lc_mask = covariates['lc_13']
     elif var == 'lc_grass_forest':
      if 'lc_4' in covariates:
       lc_mask[covariates['lc_4'] == 1] = 1
      if 'lc_2' in covariates:
       lc_mask[covariates['lc_2'] == 1] = 1
      if 'lc_5' in covariates:
       lc_mask[covariates['lc_5'] == 1] = 1
      if 'lc_6' in covariates:
       lc_mask[covariates['lc_6'] == 1] = 0.66
      if 'lc_11' in covariates:
       lc_mask[covariates['lc_11'] == 1] = 0.66
      if 'lc_12' in covariates:
       lc_mask[covariates['lc_12'] == 1] = 0.66
      if 'lc_10' in covariates:
       lc_mask[covariates['lc_10'] == 1] = 0.33
      if 'lc_16' in covariates:
       lc_mask[covariates['lc_16'] == 1] = 0.01
     cvs1[var] = {'min':0,'max':1,'t':-9999,'d':lc_mask}
    else:
     tmp1 = np.copy(hp_in_main[var])
     cvs1[var] = {'min':np.min(tmp1),'max':np.max(tmp1),'t':-9999,'d':tmp1}

  if (np.sum(m_abst) != 0):
   cvs2 = {}
   for var in abst_subbasin_clustering_cov:
    if var in ['lc_w_now','lc_urb_nourb','lc_grass_forest']:
     lc_mask = np.copy(dem)
     lc_mask[:] = 0.0
     if var == 'lc_w_now':
      lc_mask = covariates['lc_17']
     elif var == 'lc_urb_nourb':
      lc_mask = covariates['lc_13']
     elif var == 'lc_grass_forest':
      if 'lc_4' in covariates:
       lc_mask[covariates['lc_4'] == 1] = 1
      if 'lc_2' in covariates:
       lc_mask[covariates['lc_2'] == 1] = 1
      if 'lc_5' in covariates:
       lc_mask[covariates['lc_5'] == 1] = 1
      if 'lc_6' in covariates:
       lc_mask[covariates['lc_6'] == 1] = 0.66
      if 'lc_11' in covariates:
       lc_mask[covariates['lc_11'] == 1] = 0.66
      if 'lc_12' in covariates:
       lc_mask[covariates['lc_12'] == 1] = 0.66
      if 'lc_10' in covariates:
       lc_mask[covariates['lc_10'] == 1] = 0.33
      if 'lc_16' in covariates:
       lc_mask[covariates['lc_16'] == 1] = 0.01
     cvs2[var] = {'min':0,'max':1,'t':-9999,'d':lc_mask}
    else:
     tmp2 = np.copy(hp_in_abst[var])
     cvs2[var] = {'min':np.min(tmp2),'max':np.max(tmp2),'t':-9999,'d':tmp2}

  if (np.sum(m_abst) != 0):
   keys2 = list(cvs2.keys())
  if (np.sum(m_main) != 0):
   keys1 = list(cvs1.keys())

  if hydroblocks_info['channel_initiation']['flag_subgrid'] == True:
   if (np.sum(m_main) != 0):
    for var in keys1:
     cvs1[var]['w'] = 1
   if (np.sum(m_abst) != 0):
    n_pc = 0
    n_no_pc = 0
    for var in keys2:
     if 'pc_' in var:
      n_pc += 1
     else:
      n_no_pc += 1
    for var in keys2:
     if 'pc_' in var:
      cvs2[var]['w'] = (1/(n_no_pc+1))/n_pc
     else:
      cvs2[var]['w'] = (1/(n_no_pc+1))
  else:
   if (np.sum(m_main) != 0):
    for var in keys1:
     cvs1[var]['w'] = 1
   if (np.sum(m_abst) != 0):
    for var in keys2:
     cvs2[var]['w'] = 1

  if (np.sum(m_main) != 0):
   (basin_clusters_main,) = terrain_tools.cluster_basins_hmc_2(basins_wob,cvs1,hp_in_main,ncatchments_main,1)
  if (np.sum(m_abst) != 0):
   (basin_clusters_abst,) = terrain_tools.cluster_basins_hmc_2(basins_wob,cvs2,hp_in_abst,ncatchments_abst,ncatchments_main+1)

  if (np.sum(m_abst) != 0) and (np.sum(m_main) != 0):
   basin_clusters = np.copy(basin_clusters_abst)
   basin_clusters[basin_clusters_main!=-9999] = basin_clusters_main[basin_clusters_main!=-9999]
  elif (np.sum(m_abst) != 0) and (np.sum(m_main) == 0):
   basin_clusters = np.copy(basin_clusters_abst)
  else:
   basin_clusters = np.copy(basin_clusters_main)

 else:
  ncatchments = min(ncatchments,np.unique(basins_wob)[1:].size)
  if 'large_scale_basins' in subbasin_clustering_cov:
   subbasin_clustering_cov.remove('large_scale_basins')
   subbasin_clustering_cov = subbasin_clustering_cov + ['lsb_lats','lsb_lons']
  if 'lc' in subbasin_clustering_cov:
   subbasin_clustering_cov.remove('lc')
   subbasin_clustering_cov = subbasin_clustering_cov + ['lc_w_now','lc_urb_nourb','lc_grass_forest']
  cvs = {}
  for var in subbasin_clustering_cov:
   if var in ['lc_w_now','lc_urb_nourb','lc_grass_forest']:
    lc_mask = np.copy(dem)
    lc_mask[:] = 0.0
    if var == 'lc_w_now':
     lc_mask = covariates['lc_17']
    elif var == 'lc_urb_nourb':
     if 'lc_13' in covariates:
      lc_mask = covariates['lc_13']
    elif var == 'lc_grass_forest':
     if 'lc_4' in covariates:
      lc_mask[covariates['lc_4'] == 1] = 1
     if 'lc_2' in covariates:
      lc_mask[covariates['lc_2'] == 1] = 1
     if 'lc_5' in covariates:
      lc_mask[covariates['lc_5'] == 1] = 1
     if 'lc_6' in covariates:
      lc_mask[covariates['lc_6'] == 1] = 0.66
     if 'lc_11' in covariates:
      lc_mask[covariates['lc_11'] == 1] = 0.66
     if 'lc_12' in covariates:
      lc_mask[covariates['lc_12'] == 1] = 0.66
     if 'lc_10' in covariates:
      lc_mask[covariates['lc_10'] == 1] = 0.33
     if 'lc_16' in covariates:
      lc_mask[covariates['lc_16'] == 1] = 0.01
    cvs[var] = {'min':0,'max':1,'t':-9999,'d':lc_mask}
   else:
    tmp = np.copy(hp_in[var])
    cvs[var] = {'min':np.min(tmp),'max':np.max(tmp),'t':-9999,'d':tmp}
  keys = list(cvs.keys())
  if hydroblocks_info['channel_initiation']['flag_subgrid'] == True:
   n_pc = 0
   n_no_pc = 0
   for var in keys:
    if 'pc_' in var:
     n_pc += 1
    else:
     n_no_pc += 1
   for var in keys:
    if 'pc_' in var:
     cvs[var]['w'] = (1/(n_no_pc+1))/n_pc
    else:
     cvs[var]['w'] = (1/(n_no_pc+1))
  else:
   for var in keys:
    cvs[var]['w'] = 1

  (basin_clusters,) = terrain_tools.cluster_basins_hmc_2(basins_wob,cvs,hp_in,ncatchments,1)

 if np.min(basin_clusters[basin_clusters!=-9999]) == 0:
  basin_clusters[basin_clusters!=-9999] = basin_clusters[basin_clusters!=-9999] + 1

 if len(db_channels['length']) == len(np.unique(basins_wob[basins_wob!=-9999])):
  basins_wob[basins_wob!=-9999] = basins_wob[basins_wob!=-9999] + 1

 if len(np.unique(basin_clusters[basin_clusters!=-9999])) < ncatchments:
  nbcl = 1
  for bcl in np.unique(basin_clusters[basin_clusters!=-9999]):
   basin_clusters[basin_clusters == bcl] = nbcl
   nbcl += 1

 print('Computing height above nearest drainage area',flush=True)
 hand = terrain_tools.ttf.calculate_depth2channel(channels_wob,basins_wob,fdir,dem)
 hand[(hand == -9999) & (basins_wob!=-9999)] = 0.0

 ubcs = np.unique(basin_clusters)
 ubcs = ubcs[ubcs != -9999]
 for ubc in ubcs:
  ubs = np.unique(basins_wob[basin_clusters == ubc])
  ubs = ubs[ubs != -9999]
  db_channels['width'][ubs-1] = np.mean(db_channels['width'][ubs-1])
  db_channels['bankfull'][ubs-1] = np.mean(db_channels['bankfull'][ubs-1])

 print('%s Discretizing clusters of basins (hbands)' % cid,flush=True)
 print(cid,np.unique(basin_clusters),len(np.unique(basin_clusters[basin_clusters!=-9999])),flush=True)
 n_binning = dh
 max_nbins = 100
 (tiles,new_hand,tile_position) = terrain_tools.create_basin_tiles_updated(basin_clusters,hand,basins_wob,n_binning,cid,max_nbins)

 dict_tiling = {}
 dict_tiling['tiles'] = tiles
 dict_tiling['tile_position'] = tile_position
 dict_tiling['basin_clusters'] = basin_clusters
 dict_tiling['hand'] = hand
 dict_tiling['basins_wob'] = basins_wob
 dict_tiling['n_binning'] = n_binning
 dict_tiling['max_nbins'] = max_nbins
 pickle.dump(dict_tiling,open('%s/tiling_params_%s.pck' % (hydroblocks_info['input_dir'],cid),'wb'))

 (db_routing,area_adj,new_hand2) = Build_Hillslope_River_Database(channels_wob,mask,fdir,eares,tiles,hand,basins_wob,basin_clusters,new_hand,db_routing,ubcs,tile_position,db_channels)

 intraband_clust_vars = hydroblocks_info['hmc_parameters']['intraband_clustering_covariates']
 if 'lc' in intraband_clust_vars:
  intraband_clust_vars.remove('lc')
  intraband_clust_vars = intraband_clust_vars + ['lc_w_now','lc_urb_nourb','lc_grass_forest']

 cvs = {}
 for var in intraband_clust_vars:
  if var in ['lc_w_now','lc_urb_nourb','lc_grass_forest']:
   lc_mask = np.copy(dem)
   lc_mask[:] = 0.0
   if var == 'lc_w_now':
    lc_mask = covariates['lc_17']
   elif var == 'lc_urb_nourb':
    if 'lc_13' in covariates:
     lc_mask = covariates['lc_13']
   elif var == 'lc_grass_forest':
    if 'lc_4' in covariates:
     lc_mask[covariates['lc_4'] == 1] = 1
    if 'lc_2' in covariates:
     lc_mask[covariates['lc_2'] == 1] = 1
    if 'lc_5' in covariates:
     lc_mask[covariates['lc_5'] == 1] = 1
    if 'lc_6' in covariates:
     lc_mask[covariates['lc_6'] == 1] = 0.66
    if 'lc_11' in covariates:
     lc_mask[covariates['lc_11'] == 1] = 0.66
    if 'lc_12' in covariates:
     lc_mask[covariates['lc_12'] == 1] = 0.66
    if 'lc_10' in covariates:
     lc_mask[covariates['lc_10'] == 1] = 0.33
    if 'lc_16' in covariates:
     lc_mask[covariates['lc_16'] == 1] = 0.01
   cvs[var] = {'min':0,'max':1,'t':-9999,'d':lc_mask}
  else:
   cvs[var] = {'min':np.min(covariates[var][covariates[var]!=-9999]),
               'max':np.max(covariates[var][covariates[var]!=-9999]),
               't':-9999,
               'd':covariates[var]}

 print('Clustering the height bands into clusters', flush=True)
 m = (basin_clusters == -9999) | (tiles == -9999)
 basin_clusters[m] = -9999
 tiles[m] = -9999

 hrus = terrain_tools.create_hrus_hydroblocks(basin_clusters,tiles,cvs,nclusters,cid)
 hrus[hrus!=-9999] = hrus[hrus!=-9999] - 1
 nhru = np.unique(hrus[hrus!=-9999]).size

 os.system('rm %s/routing_info.pck' % input_dir)
 os.system('rm %s/routing_mp_connectivity.pck' % input_dir)
 pickle.dump(db_routing,open('%s/routing_info.pck' % input_dir,'wb'))
 pickle.dump(db_routing['mp_connectivity'],open('%s/routing_mp_connectivity.pck' % input_dir,'wb'))

 HMC_info = {}
 HMC_info['basins'] = basins
 HMC_info['tile_position'] = tile_position
 HMC_info['channel_map'] = channels_wob

 return (hrus.astype(np.float32),nhru,new_hand,HMC_info,covariates,db_channels,new_hand2,
   basins,basin_clusters,hand,tiles,area_adj,tile_position,z_data)

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


def go_downstream_shreve(channel, topo, shreve):
 if topo[channel] != -1:
  shreve[topo[channel]] = shreve[topo[channel]] + 1
  shreve = go_downstream_shreve(topo[channel], topo, shreve)
 return shreve


def Subgrid_Indices(channels_wob_sg,topology_sg,basins_wob,db_channels_sg,thr_var,cid):
 length_sg = db_channels_sg['length']
 width_sg = db_channels_sg['width']
 slope_sg = db_channels_sg['slope']
 acc_sg = db_channels_sg['acc']
 bankfull_sg = db_channels_sg['bankfull']

 list_nchannels=[]
 list_avrg_id_ilngth=[]
 list_avrg_id_wdth=[]
 list_avrg_id_slp=[]
 list_avrg_shrtst_pth_wdth=[]
 list_avrg_shrtst_pth_slp=[]
 list_avrg_shrtst_pth_lngth=[]
 list_grc_uw=[]
 list_g_uw=[]
 list_g_ilngth=[]
 list_g_wdth=[]
 list_g_slp=[]
 list_g_eff=[]
 list_total_lngth=[]
 list_total_acc=[]
 list_avrg_lngth=[]
 list_avrg_acc=[]
 list_avrg_wdth=[]
 list_avrg_slpe=[]
 list_avrg_bnkfll=[]
 list_drng_dnsty=[]

 for b in np.unique(basins_wob):
  if b == -9999:
   continue
  G = nx.DiGraph()
  nds = np.unique(channels_wob_sg[basins_wob == b])
  nds = nds[nds != 0] - 1
  init_out_nd = -1
  for nd in nds:
   if topology_sg[nd] == -1:
    G.add_node(nd)
    G.add_edge(nd, init_out_nd, length=length_sg[nd], width=width_sg[nd], ilength=1/length_sg[nd], slope=slope_sg[nd], acc=acc_sg[nd], bf=bankfull_sg[nd])
    init_out_nd = init_out_nd - 1
   else:
    G.add_node(nd)
    G.add_node(topology_sg[nd])
    G.add_edge(nd, topology_sg[nd], length=length_sg[nd], width=width_sg[nd], ilength=1/length_sg[nd], slope=slope_sg[nd], acc=acc_sg[nd], bf=bankfull_sg[nd])

  if G.number_of_nodes() == 0:
   list_nchannels.append(-100)
   list_avrg_id_ilngth.append(-100)
   list_avrg_id_wdth.append(-100)
   list_avrg_id_slp.append(-100)
   list_avrg_shrtst_pth_wdth.append(-100)
   list_avrg_shrtst_pth_slp.append(-100)
   list_avrg_shrtst_pth_lngth.append(-100)
   list_grc_uw.append(-100)
   list_g_uw.append(-100)
   list_g_ilngth.append(-100)
   list_g_wdth.append(-100)
   list_g_slp.append(-100)
   list_g_eff.append(-100)
   list_total_lngth.append(-100)
   list_total_acc.append(-100)
   list_avrg_lngth.append(-100)
   list_avrg_acc.append(-100)
   list_avrg_wdth.append(-100)
   list_avrg_slpe.append(-100)
   list_avrg_bnkfll.append(-100)
   list_drng_dnsty.append(-100)
  else:
   list_nchannels.append(len(G.edges))
   id_ilngth = G.degree(weight='ilength')
   id_width = G.degree(weight='width')
   id_slope = G.degree(weight='slope')

   list_id_ilngth = [val for (node, val) in id_ilngth]
   list_id_width = [val for (node, val) in id_width]
   list_id_slope = [val for (node, val) in id_slope]

   list_avrg_id_ilngth.append(sum(list_id_ilngth)/G.number_of_nodes())
   list_avrg_id_wdth.append(sum(list_id_width)/G.number_of_nodes())
   list_avrg_id_slp.append(sum(list_id_slope)/G.number_of_nodes())

   if nx.is_weakly_connected(G):
    list_avrg_shrtst_pth_wdth.append(nx.average_shortest_path_length(G, weight='width'))
    list_avrg_shrtst_pth_slp.append(nx.average_shortest_path_length(G, weight='slope'))
    list_avrg_shrtst_pth_lngth.append(nx.average_shortest_path_length(G, weight='lenght'))
   else:
    shrtst_pth_lnght_w = []
    shrtst_pth_lnght_s = []
    shrtst_pth_lnght_l = []
    sub_graphs = nx.weakly_connected_components(G)
    for sg in sub_graphs:
     SG = G.subgraph(sg).copy()
     shrtst_pth_lnght_w.append(nx.average_shortest_path_length(SG, weight='width'))
     shrtst_pth_lnght_s.append(nx.average_shortest_path_length(SG, weight='slope'))
     shrtst_pth_lnght_l.append(nx.average_shortest_path_length(SG, weight='lenght'))
    list_avrg_shrtst_pth_wdth.append(np.mean(shrtst_pth_lnght_w))
    list_avrg_shrtst_pth_slp.append(np.mean(shrtst_pth_lnght_s))
    list_avrg_shrtst_pth_lngth.append(np.mean(shrtst_pth_lnght_l))

   list_grc_uw.append(nx.global_reaching_centrality(G))
   Gud = G.to_undirected()
   spctrm_uw = nx.adjacency_spectrum(Gud)
   spctrm_ilngth = nx.adjacency_spectrum(Gud, weight='ilength')
   spctrm_wdth = nx.adjacency_spectrum(Gud, weight='width')
   spctrm_slp = nx.adjacency_spectrum(Gud, weight='slope')

   list_g_uw.append(abs(np.max(spctrm_uw)))
   list_g_ilngth.append(abs(np.max(spctrm_ilngth)))
   list_g_wdth.append(abs(np.max(spctrm_wdth)))
   list_g_slp.append(abs(np.max(spctrm_slp)))
   list_g_eff.append(nx.global_efficiency(Gud))
   list_total_lngth.append(G.size(weight='length'))
   list_total_acc.append(G.size(weight='acc'))
   list_avrg_lngth.append(list_total_lngth[-1]/len(G.edges))
   list_avrg_acc.append(list_total_acc[-1]/len(G.edges))
   list_avrg_wdth.append(G.size(weight='width')/len(G.edges))
   list_avrg_slpe.append(G.size(weight='slope')/len(G.edges))
   list_avrg_bnkfll.append(G.size(weight='bf')/len(G.edges))
   list_drng_dnsty.append(list_total_lngth[-1]/list_total_acc[-1])

 list_nchannels = np.array(list_nchannels)
 list_nchannels[list_nchannels == -100] = np.mean(list_nchannels[list_nchannels != -100])
 list_avrg_id_ilngth = np.array(list_avrg_id_ilngth)
 list_avrg_id_ilngth[list_avrg_id_ilngth == -100] = np.mean(list_avrg_id_ilngth[list_avrg_id_ilngth != -100])
 list_avrg_id_wdth = np.array(list_avrg_id_wdth)
 list_avrg_id_wdth[list_avrg_id_wdth == -100] = np.mean(list_avrg_id_wdth[list_avrg_id_wdth != -100])
 list_avrg_id_slp = np.array(list_avrg_id_slp)
 list_avrg_id_slp[list_avrg_id_slp == -100] = np.mean(list_avrg_id_slp[list_avrg_id_slp != -100])
 list_avrg_shrtst_pth_wdth = np.array(list_avrg_shrtst_pth_wdth)
 list_avrg_shrtst_pth_wdth[list_avrg_shrtst_pth_wdth == -100] = np.mean(list_avrg_shrtst_pth_wdth[list_avrg_shrtst_pth_wdth != -100])
 list_avrg_shrtst_pth_slp = np.array(list_avrg_shrtst_pth_slp)
 list_avrg_shrtst_pth_slp[list_avrg_shrtst_pth_slp == -100] = np.mean(list_avrg_shrtst_pth_slp[list_avrg_shrtst_pth_slp != -100])
 list_avrg_shrtst_pth_lngth = np.array(list_avrg_shrtst_pth_lngth)
 list_avrg_shrtst_pth_lngth[list_avrg_shrtst_pth_lngth == -100] = np.mean(list_avrg_shrtst_pth_lngth[list_avrg_shrtst_pth_lngth != -100])
 list_grc_uw = np.array(list_grc_uw)
 list_grc_uw[list_grc_uw == -100] = np.mean(list_grc_uw[list_grc_uw != -100])
 list_g_uw = np.array(list_g_uw)
 list_g_uw[list_g_uw == -100] = np.mean(list_g_uw[list_g_uw != -100])
 list_g_ilngth = np.array(list_g_ilngth)
 list_g_ilngth[list_g_ilngth == -100] = np.mean(list_g_ilngth[list_g_ilngth != -100])
 list_g_wdth = np.array(list_g_wdth)
 list_g_wdth[list_g_wdth == -100] = np.mean(list_g_wdth[list_g_wdth != -100])
 list_g_slp = np.array(list_g_slp)
 list_g_slp[list_g_slp == -100] = np.mean(list_g_slp[list_g_slp != -100])
 list_g_eff = np.array(list_g_eff)
 list_g_eff[list_g_eff == -100] = np.mean(list_g_eff[list_g_eff != -100])
 list_total_lngth = np.array(list_total_lngth)
 list_total_lngth[list_total_lngth == -100] = np.mean(list_total_lngth[list_total_lngth != -100])
 list_total_acc = np.array(list_total_acc)
 list_total_acc[list_total_acc == -100] = np.mean(list_total_acc[list_total_acc != -100])
 list_avrg_lngth = np.array(list_avrg_lngth)
 list_avrg_lngth[list_avrg_lngth == -100] = np.mean(list_avrg_lngth[list_avrg_lngth != -100])
 list_avrg_acc = np.array(list_avrg_acc)
 list_avrg_acc[list_avrg_acc == -100] = np.mean(list_avrg_acc[list_avrg_acc != -100])
 list_avrg_wdth = np.array(list_avrg_wdth)
 list_avrg_wdth[list_avrg_wdth == -100] = np.mean(list_avrg_wdth[list_avrg_wdth != -100])
 list_avrg_slpe = np.array(list_avrg_slpe)
 list_avrg_slpe[list_avrg_slpe == -100] = np.mean(list_avrg_slpe[list_avrg_slpe != -100])
 list_avrg_bnkfll = np.array(list_avrg_bnkfll)
 list_avrg_bnkfll[list_avrg_bnkfll == -100] = np.mean(list_avrg_bnkfll[list_avrg_bnkfll != -100])
 list_drng_dnsty = np.array(list_drng_dnsty)
 list_drng_dnsty[list_drng_dnsty == -100] = np.mean(list_drng_dnsty[list_drng_dnsty != -100])

 metrics = np.zeros((len(list_nchannels), 21))
 metrics[:,0] = list_nchannels
 metrics[:,1] = list_avrg_id_ilngth
 metrics[:,2] = list_avrg_id_wdth
 metrics[:,3] = list_avrg_id_slp
 metrics[:,4] = list_avrg_shrtst_pth_wdth
 metrics[:,5] = list_avrg_shrtst_pth_slp
 metrics[:,6] = list_avrg_shrtst_pth_lngth
 metrics[:,7] = list_grc_uw
 metrics[:,8] = list_g_uw
 metrics[:,9] = list_g_ilngth
 metrics[:,10] = list_g_wdth
 metrics[:,11] = list_g_slp
 metrics[:,12] = list_g_eff
 metrics[:,13] = list_total_lngth
 metrics[:,14] = list_total_acc
 metrics[:,15] = list_avrg_lngth
 metrics[:,16] = list_avrg_acc
 metrics[:,17] = list_avrg_wdth
 metrics[:,18] = list_avrg_slpe
 metrics[:,19] = list_avrg_bnkfll
 metrics[:,20] = list_drng_dnsty

 n_components = min(metrics.shape[0], metrics.shape[1])
 X_std = (metrics - np.mean(metrics, axis=0)) / np.std(metrics, axis=0)
 pca = PCA(n_components=n_components)
 pca.fit(X_std)
 explained_variance = np.cumsum(pca.explained_variance_ / np.sum(pca.explained_variance_))
 n_comp = n_components
 for i in range(explained_variance.shape[0]):
  if explained_variance[i] > thr_var:
   n_comp = i + 1
   break
 pca = PCA(n_components=n_comp)
 pca.fit(X_std)
 return pca.transform(X_std)

def Calculate_HRU_Connections_Matrix_HMC(cluster_ids,nhru,dx,HMC_info,hydroblocks_info):

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

 #If there're not lateral connections (just diagonal) create a single "fake" connection, laura
 if hdst.size == 0:
  hdst = np.array([0])
  horg = np.array([0])

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
 mask_object = gdal_tools.read_data('%s/mask_latlon.tif' % workspace)
 terrain_tools.calculate_area(mask_object)
 resx = np.mean(mask_object.area**0.5) #all pixels in the subdomain have the same resolution in x and y; still not ideal and needs to be revisited, but much better than resx = 90...

 second_pass_flag = hydroblocks_info.get('second_pass', False) # 2-Step-HMC
 network_abst_cfg = hydroblocks_info.get('network_abstraction') # 2-Step-HMC
 network_abst_flag = bool(network_abst_cfg.get('flag', False)) if isinstance(network_abst_cfg, dict) else False

 #Determine the HRUs
 if second_pass_flag and network_abst_flag and os.path.isfile('%s/covariates.pck' % input_dir):
  print("Computing the HRUs (second pass)",flush=True)
  (cluster_ids,nhru,new_hand,HMC_info,covariates,dbc,hand,basins,basin_clusters,hand_org,hbands,area_adj,tile_position,z_data) = Compute_HRUs_Semidistributed_HMC2(hydroblocks_info,resx,input_dir)
  mask = covariates['mask']
 else:
  print("Creating and curating the covariates",flush=True)
  (covariates,mask,z_data)=Create_and_Curate_Covariates_svp(wbd,hydroblocks_info)
  print("Computing the HRUs",flush=True)
  (cluster_ids,nhru,new_hand,HMC_info,covariates,dbc,hand,basins,basin_clusters,hand_org,hbands,area_adj,tile_position) = Compute_HRUs_Semidistributed_HMC(covariates,mask,hydroblocks_info,wbd,resx,input_dir)
 #covariates['hand'] = new_hand
 covariates['hand'] = hand
 hydroblocks_info['nhru'] = nhru
  
 #Create the netcdf file
 file_netcdf = '%s/input_file2.nc' % hydroblocks_info['input_dir'] if second_pass_flag else '%s/input_file.nc' % hydroblocks_info['input_dir']
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
  cmatrix = Calculate_HRU_Connections_Matrix_HMC(cluster_ids,nhru,resx,HMC_info,hydroblocks_info)
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

 return (OUTPUT,covariates,z_data)

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

'''def Prepare_Water_Use_Semidistributed(workspace,wbd,OUTPUT,input_dir,info,hydroblocks_info):

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

 return'''

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
  #for cid in [509,]:
  print(rank,size,cid)
  metadata['cid'] = cid
  metadata['input_dir'] = "%s/%d" % (edir,cid)
  metadata['workspace'] = "%s/data/cids/%d" % (rdir,cid)
  #Prepare model data
  tic = time.time()
  Prepare_Model_Input_Data(metadata)
  print("Elapsed time: ",time.time() - tic)
 comm.Barrier()

 # Additional domain-wide preprocessing for network abstraction.
 network_abst_cfg = metadata.get('network_abstraction')
 network_abst_flag = bool(network_abst_cfg.get('flag', False)) if isinstance(network_abst_cfg, dict) else False
 metadata['network_abstraction'] = network_abst_cfg if isinstance(network_abst_cfg, dict) else {'flag': network_abst_flag}
 
 #Create enhanced input data file
 print('Connect cell networks',flush=True)
 Connect_Cell_Networks_v2(rank,size,cids,edir)
 comm.Barrier() #Wait until they are all done
 
 #Create workspace for intermediate files
 workspace = '%s/workspace' % (edir)
 os.system('mkdir -p %s' % workspace)
 if network_abst_flag:
   # Build connected topology across CIDs
   print('Connect topology',flush=True)
   Topology_Connected(rank,size,cids,edir,comm)
   comm.Barrier()
   # Create trees assigned to domain-spanning outlets
   print('Compute large-scale watersheds',flush=True)
   Create_Trees(rank,size,cids,edir,comm,False)
   comm.Barrier()
   # Correct Shreve order across CIDs
   print('Correct Shreve order',flush=True)
   Correct_Shreve(rank,size,cids,edir,comm)
   comm.Barrier()
   # Compute explicit/abstract reach mask
   Network_Abstraction(rank,size,cids,edir,comm,metadata)
   comm.Barrier()
   # Second pass per-CID (recompute with abstraction-aware clustering)
   for cid in cids[rank::size]:
    print('Second-pass HMC Prepare for CID', cid, flush=True)
    metadata['cid'] = cid
    metadata['input_dir'] = "%s/%d" % (edir,cid)
    metadata['workspace'] = "%s/data/cids/%d" % (rdir,cid)
    # Preserve the first-pass file so the replacement step can recover
    if os.path.isfile('%s/input_file.nc' % metadata['input_dir']):
     os.system('cp %s/input_file.nc %s/input_file3.nc' % (metadata['input_dir'], metadata['input_dir']))
     # Indicate this is the second pass so Prepare_Model_Input_Data will write final TIFFs
     metadata['second_pass'] = True
     Prepare_Model_Input_Data(metadata)
     Replace_Stream_Network(metadata)
     # Clean up the second_pass flag to avoid side effects
     metadata.pop('second_pass', None)
   comm.Barrier()

 #Create downstream channel database for particle tracker routing scheme
 Create_Downstream_Channels_Database(edir,rank,size,cids,comm)
 comm.Barrier()


 Finalize_River_Network_Database(rdir,edir,cids,workspace,comm,rank,size)
 comm.Barrier()
 
 #Postprocess the model input 
 Postprocess_Input(rdir,edir,cids,rank,size,comm)
 comm.Barrier()


def Topology_Connected(rank,size,cids,edir,comm):
 for cid in cids[rank::size]:
  #Connected topology part
  topo_new=[]
  fp = h5py.File('%s/%s/input_file.nc' % (edir,cid),'r')
  topology=fp['stream_network']['topology'][:]
  outlets=fp['stream_network']['outlets'][:]
  fp.close()
  indices = [i for i, item in enumerate(topology) if item == -1]
  for i in range(0,topology.shape[0]):
   if i not in indices:
    topo_new.append('%s-%s'%(int(cid),int(topology[i])))
   else:
    outlets_real=outlets[outlets[:,2]!=-9999]
    if i in list(outlets_real[:,1]):
     for indx in range(0,outlets_real.shape[0]):
      if outlets_real[indx,1]==i:
       topo_new.append('%s-%s'%(int(outlets_real[indx,2]),int(outlets_real[indx,3])))
       break
    else:
     topo_new.append('-1')
  pickle.dump(topo_new,open('%s/workspace/topology_connected_%s.pck'%(edir,cid),'wb'))
 comm.Barrier()
  
 return


def Create_Trees(rank,size,cids,edir,comm,flag_network_abst=False):
  # Load connected topology for all CIDs
  topo_all = []
  for cid in cids:
    topo_file = '%s/workspace/topology_connected_%s.pck' % (edir,cid)
    data = pickle.load(open(topo_file,'rb'))
    topo_all.append(data)

  # For each CID assigned to this rank, build tree assignment per channel
  for cid in cids[rank::size]:
    tree_cid = [-9999] * len(topo_all[int(cid)-1])
    for i in range(0,len(tree_cid)):
      lista = []
      fchid = '%s-%s' % (cid,i)
      fcid = fchid.split('-')[0]
      chid = fchid.split('-')[1]
      lista.append(fchid)
      while topo_all[int(fcid)-1][int(chid)] != '-1':
        fchid = topo_all[int(fcid)-1][int(chid)]
        lista.append(fchid)
        fcid = fchid.split('-')[0]
        chid = fchid.split('-')[1]
      a = lista[-1]
      tree_cid[i] = a
    pickle.dump(tree_cid,open('%s/workspace/trees_%s.pck' %(edir,cid),'wb'))

  # Number trees (global across CIDs)
  trees = []
  list_cids = glob.glob('%s/workspace/trees_*' % edir)
  for p in list_cids:
    data = pickle.load(open(p,'rb'))
    trees.extend(data)
  un_trees = np.unique(trees)
  n_trees = np.linspace(1, len(un_trees), len(un_trees))

  # Assign tree numbers to input_file per CID (parallel across ranks)
  for cid in cids[rank::size]:
    fp = h5py.File('%s/%s/input_file.nc' % (edir,cid), 'a')
    data_trees = pickle.load(open('%s/workspace/trees_%s.pck' % (edir,cid), 'rb'))
    trees_mp = []
    for i in range(0, fp['stream_network']['topology'][:].shape[0]):
      tree = data_trees[i]
      indices = [j for j, item in enumerate(un_trees) if item == tree]
      if len(indices) != 1:
        print('WARNING: Error in unique trees', cid, i, len(indices))
      trees_mp.append(n_trees[indices[0]])
    if 'trees_domain' in fp['stream_network'].keys():
      del fp['stream_network']['trees_domain']
    fp['stream_network']['trees_domain'] = np.array(trees_mp)
    fp.close()

  # Save data used for shreve correction and abstraction per cid
  data = {}
  for cid in cids[rank::size]:
    fp = h5py.File('%s/%s/input_file.nc' % (edir,cid), 'r')
    data['acc'] = fp['stream_network']['acc'][:]
    data['shreve'] = fp['stream_network']['shreve'][:]
    data['length'] = fp['stream_network']['length'][:]
    data['topology'] = fp['stream_network']['topology'][:]
    data['tree'] = fp['stream_network']['trees_domain'][:]
    data['inlets'] = fp['stream_network']['inlets'][:]
    data['outlets'] = fp['stream_network']['outlets'][:]
    if flag_network_abst == True:
      if fp['parameters']['lats'][:].shape[0] == fp['stream_network']['shreve'][:].shape[0]:
        data['lat_basin'] = fp['parameters']['lats'][:]
        data['lon_basin'] = fp['parameters']['lons'][:]
      else:
        data['lat_basin'] = fp['parameters']['lats'][1:]
        data['lon_basin'] = fp['parameters']['lons'][1:]
    fp.close()
    pickle.dump(data,open('%s/workspace/data_channels_%s.pck' % (edir,cid),'wb'))
  comm.Barrier()

  # If modified HMC or abstraction requested, compute mean lat/lon per tree (rank 0)
  if flag_network_abst == True:
    if rank == 0:
      list_cids = glob.glob('%s/workspace/data_channels_*' % edir)
      dict_out = {}
      dict_out['lats'] = []
      dict_out['lons'] = []
      dict_out['tree'] = []
      dict_out['mean_lat'] = []
      dict_out['mean_lon'] = []
      dict_out['index'] = []
      dict_out['cid'] = []
      count = 0
      dict_out['index'].append(count)
      for i in list_cids:
        d = pickle.load(open(i,'rb'))
        cid_i = i.split('/')[-1].split('data_channels_')[-1].split('.pck')[0]
        dict_out['lats'].extend(d['lat_basin'])
        dict_out['lons'].extend(d['lon_basin'])
        dict_out['tree'].extend(d['tree'])
        dict_out['cid'].append(cid_i)
        count += len(d['lat_basin'])
        dict_out['index'].append(count)
      dict_out['lats'] = np.array(dict_out['lats'])
      dict_out['lons'] = np.array(dict_out['lons'])
      dict_out['tree'] = np.array(dict_out['tree'])
      dict_out['mean_lat'] = np.copy(dict_out['lats'])
      dict_out['mean_lon'] = np.copy(dict_out['lons'])
      for t in np.unique(dict_out['tree']):
        mean_lat = np.mean(dict_out['lats'][dict_out['tree'] == t])
        mean_lon = np.mean(dict_out['lons'][dict_out['tree'] == t])
        dict_out['mean_lat'][dict_out['tree'] == t] = mean_lat
        dict_out['mean_lon'][dict_out['tree'] == t] = mean_lon
      dict_out['CID_info'] = {}
      dict_out['CID_info']['lats'] = {}
      dict_out['CID_info']['lons'] = {}
      count2 = 0
      for cid in dict_out['cid']:
        cid_i = int(cid)
        dict_out['CID_info']['lats'][cid_i] = dict_out['mean_lat'][dict_out['index'][count2]:int(dict_out['index'][count2+1])]
        dict_out['CID_info']['lons'][cid_i] = dict_out['mean_lon'][dict_out['index'][count2]:int(dict_out['index'][count2+1])]
        count2 += 1
      pickle.dump(dict_out,open('%s/workspace/mean_lats-lons_trees.pck' % (edir),'wb'))
    comm.Barrier()
  return

def Correct_Shreve(rank,size,cids,edir,comm):
 shreve_all = []
 topo_all = []
 for cid in cids:
  file = '%s/workspace/data_channels_%s.pck' % (edir,cid)
  data = pickle.load(open(file,'rb'))
  shreve_pck = data['shreve']
  file = '%s/workspace/topology_connected_%s.pck' % (edir,cid)
  topo = pickle.load(open(file,'rb'))
  shreve_all.append(shreve_pck)
  topo_all.append(topo)

 for cid in cids[rank::size]:
  shreve_all2 = shreve_all.copy()
  file = '%s/workspace/data_channels_%s.pck' % (edir,cid)
  data = pickle.load(open(file,'rb'))
  outlets = data['outlets']
  out = outlets[outlets[:,2] != -9999]
  for i in cids:
   if i == cid:
    shreve_all2[i-1] = data['shreve']
   else:
    shreve_all2[i-1][:] = 0
  for o in range(0, out.shape[0]):
   lista = []
   add_shreve = shreve_all2[int(cid-1)][out[o,1]]
   fchid2 = '%s-%s' % (int(out[o,2]), int(out[o,3]))
   fcid2 = fchid2.split('-')[0]
   chid2 = fchid2.split('-')[1]
   lista.append(fchid2)
   while topo_all[int(fcid2)-1][int(chid2)] != '-1':
    fchid2 = topo_all[int(fcid2)-1][int(chid2)]
    lista.append(fchid2)
    fcid2 = fchid2.split('-')[0]
    chid2 = fchid2.split('-')[1]
   for s in lista:
    fcid = int(s.split('-')[0])
    chid = int(s.split('-')[1])
    shreve_all2[int(fcid)-1][chid] = shreve_all2[int(fcid)-1][chid] + add_shreve
  pickle.dump(shreve_all2,open('%s/workspace/fix_shreve_%s.pck' % (edir,cid),'wb'))
    
 if rank == 0:
  shreve_final = shreve_all.copy()
  for i in cids:
   shreve_final[i-1][:] = 0
  for cid in cids:
   shreve_cid = pickle.load(open('%s/workspace/fix_shreve_%s.pck' % (edir,cid),'rb'))
   for i in cids:
    shreve_final[i-1][:] = shreve_final[i-1][:] + shreve_cid[i-1][:]
  pickle.dump(shreve_final,open('%s/workspace/fix_shreve_all.pck' % (edir),'wb'))
 comm.Barrier()

 ##Update input_file with corrected Shreve order
 shreve_final = pickle.load(open('%s/workspace/fix_shreve_all.pck' % (edir),'rb'))
 for cid in cids[rank::size]:
  with h5py.File('%s/%s/input_file.nc' % (edir,cid),'r+') as fp3:
   if 'shreve' in fp3['stream_network'].keys():
    del fp3['stream_network']['shreve']
   fp3['stream_network']['shreve'] = shreve_final[cid-1]
   fp3.close()
 
 return


def Network_Abstraction(rank,size,cids,edir,comm,metadata):
 #Evaluates percentiles of selected variable depending on type of abstraction
 var = metadata['network_abstraction']['var']
 if rank == 0:
  list_cids = glob.glob('%s/workspace/data_channels_*' % edir)
  var_values = []
  for i in list_cids:
   d = pickle.load(open(i,'rb'))
   var_values.extend(d[var])

  p = metadata['network_abstraction']['percentile']
  thr = np.percentile(var_values, p)
  pickle.dump(thr, open('%s/workspace/percentile_analysis.pck' % (edir), 'wb'))
 comm.Barrier()

 thr = pickle.load(open('%s/workspace/percentile_analysis.pck' % (edir), 'rb'))
 for cid in cids[rank::size]:
  abst_mask = [] # 0=abstracted, 1=explicit
  fp = h5py.File('%s/%s/input_file.nc' % (edir,cid), 'a')
  if 'explicit_reach' in fp['stream_network'].keys():
   del fp['stream_network']['explicit_reach']
  for i in range(0, fp['stream_network']['shreve'][:].shape[0]):
   if fp['stream_network'][var][i] >= thr:
    abst_mask.append(1)
   else:
    abst_mask.append(0)
  abst_mask = np.array(abst_mask, dtype=int)
  fp['stream_network']['explicit_reach'] = abst_mask
  fp.close()
 return


def Replace_Stream_Network(metadata):
  # Use pathlib and shutil for safer, more explicit file operations
  input_dir = Path(metadata['input_dir'])
  source_path = input_dir / 'input_file.nc'
  backup_path = input_dir / 'input_file3.nc'
  destination_path = input_dir / 'input_file2.nc'

  if not source_path.exists():
    raise FileNotFoundError('Source input_file.nc not found: %s' % str(source_path))

  # In the second-pass flow, destination should already exist because Prepare_Model_Input_Data writes input_file2.nc.
  # Keep a fallback for non-second-pass callers.
  if not destination_path.exists():
    if metadata.get('second_pass', False):
      raise FileNotFoundError('Second-pass input_file2.nc not found: %s' % str(destination_path))
    shutil.copy2(str(source_path), str(destination_path))

  # Read the source file and write its stream-network datasets into the second-pass file.
  explicit = None
  with h5py.File(str(source_path), 'r') as src, h5py.File(str(destination_path), 'a') as dst:
    shreve = src['stream_network']['shreve'][:]
    inlets = src['stream_network']['inlets'][:]
    outlets = src['stream_network']['outlets'][:]
    trees = src['stream_network']['trees_domain'][:]
    if metadata.get('network_abstraction', {}).get('flag', False):
      if 'explicit_reach' in src['stream_network'].keys():
        explicit = src['stream_network']['explicit_reach'][:]

    grp_dst = dst['stream_network']
    # Replace datasets in destination (delete existing then write)
    for name, data in (('shreve', shreve), ('inlets', inlets), ('outlets', outlets), ('trees_domain', trees)):
      if name in grp_dst:
        del grp_dst[name]
      grp_dst[name] = data

    if metadata.get('network_abstraction', {}).get('flag', False) and explicit is not None:
      if 'explicit_reach' in grp_dst:
        del grp_dst['explicit_reach']
      grp_dst['explicit_reach'] = explicit

  # Atomically rotate files: source -> backup, destination -> source
  try:
    os.replace(str(source_path), str(backup_path))
  except Exception:
    shutil.move(str(source_path), str(backup_path))

  try:
    os.replace(str(destination_path), str(source_path))
  except Exception:
    shutil.move(str(destination_path), str(source_path))

  return

def Postprocess_Input(rdir,edir,cids,rank,size,comm):

 sdir = '%s/postprocess' % (edir)
 ddir = '%s/data/cids' % rdir
 os.system('rm -rf %s' % sdir)
 #Create cid, hru, and channel maps
 vars = ['cids','cids_org','dem','hrus','channels','hand','basins','basin_clusters']
 for var in vars:
  os.system('mkdir -p %s/postprocess/%s' % (edir,var))
 for cid in cids[rank::size]:
  print('Copying files for vrt',cid,flush=True)
  dir = '%s/%s' % (edir,cid)
  #hru
  ifile = '%s/hru_mapping_latlon.tif' % dir
  ofile = '%s/hrus/%d.tif' % (sdir,cid)
  os.system('ln -s %s %s' % (ifile,ofile))
  #channels
  ifile = '%s/channel_mapping_latlon.tif' % dir
  ofile = '%s/channels/%d.tif' % (sdir,cid)
  os.system('ln -s %s %s' % (ifile,ofile))
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
  os.system('ln -s %s %s' % (ifile,ofile))
  #dem
  ifile = '%s/%d/dem_latlon.tif' % (ddir,cid)
  ofile = '%s/dem/%d.tif' % (sdir,cid)
  os.system('ln -s %s %s' % (ifile,ofile))
  #hand
  ifile = '%s/hand_latlon.tif' % dir
  ofile = '%s/hand/%d.tif' % (sdir,cid)
  os.system('ln -s %s %s' % (ifile,ofile))
  #basins
  ifile = '%s/basins_latlon.tif' % dir
  ofile = '%s/basins/%d.tif' % (sdir,cid)
  os.system('ln -s %s %s' % (ifile,ofile))
  #basin clusters
  ifile = '%s/basin_clusters_latlon.tif' % dir
  ofile = '%s/basin_clusters/%d.tif' % (sdir,cid)
  os.system('ln -s %s %s' % (ifile,ofile))

 #Create vrts
 comm.Barrier()
 if rank == 0:
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

  #exit()
 
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
 ncmax = 250 #parameter
 dbout = {}

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
  dbout[cid] = np.copy(downstream_channels)

 comm.Barrier()

 for cid in cids[rank::size]:
  #Add downstream_channels array to input_file.nc
  file = '%s/%s/input_file.nc' % (edir,cid)
  fp = h5py.File(file,'a')
  if 'downstream_channels' in fp['stream_network']:del fp['stream_network']['downstream_channels']
  fp['stream_network']['downstream_channels'] = dbout[cid][:]
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
