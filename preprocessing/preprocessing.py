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
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE" #laura
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
import networkx as nx #laura, for topological indices
import sklearn.decomposition #laura, for pca of subgrid indices

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
 if os.path.isdir('%s/experiments/simulations/%s/workspace'%(hydroblocks_info['rdir'],hydroblocks_info['experiment'])) == False:
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
 flag_mod_hmc = hydroblocks_info['flag_mod_hmc']
 (output,covariates,hydroblocks_info,z_data) = Create_Clusters_And_Connections(workspace,wbd,output,input_dir,nhru,info,hydroblocks_info,flag_mod_hmc)

 flag_mod_hmc = hydroblocks_info['flag_mod_hmc']
 if (hydroblocks_info['network_abstraction']['flag']==False) and (flag_mod_hmc == False):#laura
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

 if (hydroblocks_info['network_abstraction']['flag']==False) and (flag_mod_hmc == False):#laura
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

 if (hydroblocks_info['network_abstraction']['flag']==False) and (flag_mod_hmc == False):#laura
  #Write out the basin cluster map
  basin_clusters_map = np.copy(output['basin_clusters_map'])
  basin_clusters_map[np.isnan(basin_clusters_map) == 1] = -9999.0
  file_ca = '%s/basin_clusters_latlon.tif' % input_dir
  metadata['nodata'] = -9999.0
  gdal_tools.write_raster(file_ca,metadata,basin_clusters_map)
  n_cluster_basins = int(len(np.unique(basin_clusters_map)) - 1)

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
 if (hydroblocks_info['network_abstraction']['flag']==False) and (flag_mod_hmc == False):#laura
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
   #for i in range(1,(int(hydroblocks_info['hmc_parameters']["number_of_characteristic_subbasins"]+1))):
   for i in range(1,int(n_cluster_basins+1)):
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
    
 if (hydroblocks_info['network_abstraction']['flag']==True) or (flag_mod_hmc == True):#laura
  dict={}
  dict['z_data'] = z_data
  dict['covariates']=covariates#laura
  dict['output']=output#laura
  pickle.dump(dict,open('%s/covariates.pck'%input_dir,'wb'))#laura

 #Close the file
 fp.close()
  
 return output

def Compute_HRUs_Semidistributed_HMC(covariates,mask,hydroblocks_info,wbd,eares,input_dir,flag_mod_hmc):

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
 cthrs = hydroblocks_info['channel_initiation']["athrs"]#laura
 ipoints = ((area > cthrs)).astype(np.int32)
 ipoints[ipoints == 0] = -9999
    
 #Calculate channel initiation points subgrid, laura
 if hydroblocks_info['channel_initiation']['flag_subgrid']==True:
  cthrs_sg = hydroblocks_info['channel_initiation']["athrs_subgrid"]
  #Ensure that threshold for subgrid network is smaller than main network
  if cthrs_sg == cthrs:
   cthrs_sg = cthrs_sg-50000 #m2

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
 
 #Compute channel_network, basins, and channel_properties for subgrid network, laura   
 if hydroblocks_info['channel_initiation']['flag_subgrid']==True:
  thr_var = hydroblocks_info['channel_initiation']['var_pca']
  dict={}
  #Compute the channels
  (channels_sg,channels_wob_sg,channel_topology_sg,
   tmp1_sg,crds_sg,
   channel_outlet_id_sg,
   channel_target_mp_sg,
   channel_target_crds_sg,
   channel_inlet_id_sg,
   channel_inlet_target_mp_sg,
   channel_inlet_target_crds_sg) = terrain_tools.ttf.calculate_channels_wocean_wprop_wcrds(ac,ac_all,cthrs_sg,cthrs_sg,fdc,mask,mask_all,np.flipud(covariates['lats']),covariates['lons'])
  
  #Curate list output
  channel_topology_sg = channel_topology_sg[channel_topology_sg != -9999]

  dict['channels_wob_sg']=channels_wob_sg
  dict['channel_topology_sg']=channel_topology_sg
  #Compute the basins
  basins_wob_sg = terrain_tools.ttf.delineate_basins(channels_wob_sg,mask,fdir)
  db_channels_sg = terrain_tools.calculate_channel_properties(channels_wob_sg,channel_topology_sg,
                                                              slope,eares,mask,area_all,area_all_cp,
                                                              basins_wob_sg,hydroblocks_info['parameter_scaling'])
  dict['db_channels_sg']=db_channels_sg

  #Function that computes topological, morphologic, graph-theory, and channel-feature indices for subgrid network. It also reduces dimensionality by using PCA that accounts for 95% (hard-coded) of the total variance, laura
  principal_components = Subgrid_Indices(channels_wob_sg,
                                         channel_topology_sg,
                                         basins_wob,
                                         db_channels_sg,
                                         thr_var)
  
  dict['principal_components'] = principal_components
  #Save pickle channel_subgrid
  pickle.dump(dict,open('%s/pca_subgrid_basins.pck' % input_dir,'wb'))

 #Compute channel properties
 db_channels = terrain_tools.calculate_channel_properties(channels_wob,channel_topology,slope,eares,mask,area_all,area_all_cp,basins_wob,hydroblocks_info['parameter_scaling'])

 #Compute Shreve order per macroscale polygon, laura
 shreve=np.copy(channel_topology)
 shreve[:] = 0
 ##Assigns order 1 to all the streams that aren't included in topology
 for i in range(0,len(channel_topology)):
  if i not in np.unique(channel_topology):
   shreve[i]=1
 ##Go dowstream from the order 1 streams until finding a reach draining outside the domain (topology=-1)
 for first_order in list(np.where(shreve==1)[0]):
  (shreve)=go_downstream_shreve(first_order,channel_topology,shreve)
 db_channels['shreve']=shreve

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

 if (hydroblocks_info['network_abstraction']['flag']==False) and (flag_mod_hmc == False): #laura
  (hrus,nhru,new_hand,covariates,db_channels,new_hand2,basin_clusters,hand,tiles,area_adj, tile_position)=regular_HMC(hydroblocks_info,basins_wob,eares,dh,covariates,ncatchments,nclusters,db_channels,channels_wob,mask,fdir,hand,db_routing,dem)

 else: #laura
  basins1 = np.copy(basins)
  basins[basins!=-9999]=basins[basins!=-9999]
  hrus = np.copy(basins1) #laura
  nhru = np.unique(hrus[hrus!=-9999]).size #laura
  new_hand = np.copy(hrus) #laura
  new_hand2 = np.copy(hrus) #laura
  basin_clusters = np.copy(hrus) #laura
  tiles = np.copy(hrus) #laura
  area_adj = np.copy(hrus) #laura
  tile_position = np.copy(hrus) #laura 

 #Save the channel info
 if os.path.isfile('%s/routing_info.pck' % input_dir):
  os.system('rm %s/routing_info.pck' % input_dir)
 if os.path.isfile('%s/routing_io.pck' % input_dir):
  os.system('rm %s/routing_io.pck' % input_dir)
 if os.path.isfile('%s/routing_mp_connectivity' % input_dir):
  os.system('rm %s/routing_mp_connectivity' % input_dir)

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


def Compute_HRUs_Semidistributed_HMC2(hydroblocks_info,eares,input_dir,flag_mod_hmc,):
    
 #Define the parameters for the hierarchical multivariate clustering
 ncatchments = hydroblocks_info['hmc_parameters']['number_of_characteristic_subbasins']
 ncatchments_main = hydroblocks_info['network_abstraction']['number_of_characteristic_main_subbasins']
 ncatchments_abst = hydroblocks_info['network_abstraction']['number_of_characteristic_secondary_subbasins']
    
 main_subbasin_clustering_cov = hydroblocks_info['network_abstraction']['main_subbasin_clustering_covariates']
 abst_subbasin_clustering_cov = hydroblocks_info['network_abstraction']['secondary_subbasin_clustering_covariates']
 subbasin_clustering_cov = hydroblocks_info['hmc_parameters']['subbasin_clustering_covariates']

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
 cov_path = '%s/covariates.pck' % input_dir
 cov = pickle.load(open(cov_path,'rb'))
 z_data = cov['z_data']
 basins_wob = cov['output']['basin_map']
 basins_wob[basins_wob!=-9999] = basins_wob[basins_wob!=-9999] - 1
 basins = basins_wob
 covariates = cov['covariates']
 db_channels = cov['output']['stream_network']
 channels_wob = cov['output']['channel_map']
 dict = pickle.load(open('%s/experiments/simulations/%s/workspace/mean_lats-lons_trees.pck'%(hydroblocks_info['rdir'],hydroblocks_info['experiment']),'rb'))
 db_routing = pickle.load(open('%s/routing_info.pck' % input_dir,'rb'))
 dem = covariates['dem']
 mask = covariates['mask']
 #Bring out the flow direction (Convert flow direction from int to 2d approach)
 fdir = terrain_tools.transform_arcgis_fdir(covariates['fdir'])
    
 hp_in = terrain_tools.calculate_basin_properties_updated(basins_wob,eares,covariates,vars)
 argsort = np.argsort(hp_in['bid'])
 for var in hp_in:
  hp_in[var] = hp_in[var][argsort]
 
 #bring in channel variables
 for var in ['width','bankfull','length','area','shreve','large_scale_basins']:
  if var in ['width','bankfull','length','area','shreve']:
   hp_in[var] = [0]
   hp_in[var].extend(input_file['stream_network'][var][:])
   hp_in[var] = np.array(hp_in[var])
  if var in ['large_scale_basins']:
   hp_in['lsb_lats'] = [input_file['parameters']['lats'][0]]
   hp_in['lsb_lons'] = [input_file['parameters']['lons'][0]]
   hp_in['lsb_lats'].extend(dict['CID_info']['lats'][cid][:])
   hp_in['lsb_lons'].extend(dict['CID_info']['lons'][cid][:])
   hp_in['lsb_lats'] = np.array(hp_in['lsb_lats'])
   hp_in['lsb_lons'] = np.array(hp_in['lsb_lons'])
 
 if hydroblocks_info['network_abstraction']['flag'] == True:
  m_main = np.array(hp_in['shreve'],dtype=bool)
  m_main[:] = 0
  m_main[1:] = input_file['stream_network']['explicit_reach'][:] == 1
  m_abst = np.array(hp_in['shreve'],dtype=bool)
  m_abst[:] = 1
  m_abst[1:] = input_file['stream_network']['explicit_reach'][:] == 0
 input_file.close()
 
 #If subgrid network flag is true, modify hp_in and subbasin_clustering_cov, main_subbasin_cov, and abst_subbasin_cov laura
 if hydroblocks_info['channel_initiation']['flag_subgrid'] == True:
  y = pickle.load(open('%s/pca_subgrid_basins.pck'%(hydroblocks_info['input_dir']),'rb'))['principal_components']
  for pc in range(0,y.shape[1]):
   v = 'pc_%s'%(pc+1)
   hp_in[v] = y[:,pc]
   subbasin_clustering_cov.append(v)
   abst_subbasin_clustering_cov.append(v)
    
 #Clustering the basins
 print("Clustering the basins",flush=True)
 #Mask hp_in for main river basins and secondary basins
 if hydroblocks_info['network_abstraction']['flag'] == True:
  hp_in_main = {}
  hp_in_abst = {}
  for key in list(hp_in.keys()):
   hp_in_main[key] = []
   hp_in_abst[key] = []         
   for i in range(0,m_main.shape[0]):
    if m_main[i]==True:
     hp_in_main[key].append(hp_in[key][i])
    elif m_abst[i]==True:
     hp_in_abst[key].append(hp_in[key][i])
   hp_in_main[key]=np.array(hp_in_main[key])
   hp_in_abst[key]=np.array(hp_in_abst[key])
  basins_main = np.copy(basins_wob)
  basins_abst = np.copy(basins_wob)
  basins_main[:]=-9999
  basins_abst[:]=-9999
  for main in hp_in_main['bid']:
   basins_main[basins_wob==main]=1
  for abst in hp_in_abst['bid']:
   basins_abst[basins_wob==abst]=1

  #Set the ncatchments to be at least the number of basins
  ncatchments_main = min(ncatchments_main,np.sum(m_main))
  ncatchments_abst = min(ncatchments_abst,np.sum(m_abst))
    
  #Assign centroid of large scale basin to covariates for clustering
  if 'large_scale_basins' in main_subbasin_clustering_cov:#laura
   main_subbasin_clustering_cov.remove('large_scale_basins') #laura
   main_subbasin_clustering_cov=main_subbasin_clustering_cov+['lsb_lats','lsb_lons']
  if 'large_scale_basins' in abst_subbasin_clustering_cov:#laura
   abst_subbasin_clustering_cov.remove('large_scale_basins') #laura
   abst_subbasin_clustering_cov=abst_subbasin_clustering_cov+['lsb_lats','lsb_lons']

  #Assemble input data
  cvs1 = {}
  cvs2 = {}
  for var in main_subbasin_clustering_cov: #laura
   if var in ['lc_w_now','lc_urb_nourb','lc_grass_forest']: #laura
    if var=='lc_w_now': #laura
     lc_mask=covariates['lc_17'] #laura
    elif var=='lc_urb_nourb': #laura
     lc_mask=covariates['lc_13'] #laura
    elif var=='lc_grass_forest': #laura
     lc_mask=np.copy(dem)
     lc_mask[:] = 0.0
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
    cvs1[var] = {'min':0, #laura
                 'max':1, #laura
                 't':-9999, #laura
                 'd':lc_mask} #laura
   else: #laura
    tmp1 = np.copy(hp_in_main[var])
    cvs1[var] = {'min':np.min(tmp1),
                 'max':np.max(tmp1),
                 't':-9999,
                 'd':tmp1}
    
  for var in abst_subbasin_clustering_cov: #laura
   if var in ['lc_w_now','lc_urb_nourb','lc_grass_forest']: #laura
    if var=='lc_w_now': #laura
     lc_mask=covariates['lc_17'] #laura
    elif var=='lc_urb_nourb': #laura
     lc_mask=covariates['lc_13'] #laura
    elif var=='lc_grass_forest': #laura
     lc_mask=np.copy(dem)
     lc_mask[:] = 0.0
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
    cvs2[var] = {'min':0, #laura
                 'max':1, #laura
                 't':-9999, #laura
                 'd':lc_mask} #laura
   else: #laura
    tmp2 = np.copy(hp_in_abst[var])
    cvs2[var] = {'min':np.min(tmp2),
                 'max':np.max(tmp2),
                 't':-9999,
                 'd':tmp2}
  
  #If subgrid network flag is true, ensure that topological principal components don't overwhelm clustering of basins, laura
  keys2 = list(cvs2.keys())
  keys1 = list(cvs1.keys())
  if hydroblocks_info['channel_initiation']['flag_subgrid'] == True:
   #cvs1
   for var in keys1:
    cvs1[var]['w'] = 1
   #cvs2
   n_pc = 0
   n_no_pc = 0
   for var in keys2:
    if 'pc_' in var:n_pc += 1
    else:n_no_pc += 1
   for var in keys2:
    if 'pc_' in var:cvs2[var]['w'] = (1/(n_no_pc+1))/n_pc
    else:cvs2[var]['w'] = (1/(n_no_pc+1))
  else:
   for var in keys1:cvs1[var]['w'] = 1
   for var in keys2:cvs2[var]['w'] = 1
    
  (basin_clusters_main,) = terrain_tools.cluster_basins_hmc_2(basins_wob,cvs1,hp_in_main,ncatchments_main,1)
  (basin_clusters_abst,) = terrain_tools.cluster_basins_hmc_2(basins_wob,cvs2,hp_in_abst,ncatchments_abst,ncatchments_main+1)
            
  basin_clusters = np.copy(basin_clusters_abst)
  basin_clusters[basin_clusters_main!=-9999] = basin_clusters_main[basin_clusters_main!=-9999]
  basins_wob[basins_wob!=-9999] = basins_wob[basins_wob!=-9999] + 1
 else:
  #Set the ncatchments to be at least the number of basins
  ncatchments = min(ncatchments,np.unique(basins_wob)[1:].size)
  #Assign centroid of large scale basin to covariates for clustering
  if 'large_scale_basins' in subbasin_clustering_cov:#laura
   subbasin_clustering_cov.remove('large_scale_basins') #laura
   subbasin_clustering_cov=subbasin_clustering_cov+['lsb_lats','lsb_lons']

  #dissaggregate land cover if it is in covariates
  if 'lc' in subbasin_clustering_cov:#laura
   subbasin_clustering_cov.remove('lc') #laura
   subbasin_clustering_cov=subbasin_clustering_cov+['lc_w_now','lc_urb_nourb','lc_grass_forest']
  #Assemble input data
  cvs = {}
  for var in subbasin_clustering_cov: #laura
   if var in ['lc_w_now','lc_urb_nourb','lc_grass_forest']: #laura
    if var=='lc_w_now': #laura
     lc_mask=covariates['lc_17'] #laura
    elif var=='lc_urb_nourb': #laura
     lc_mask=covariates['lc_13'] #laura
    elif var=='lc_grass_forest': #laura
     lc_mask=np.copy(dem)
     lc_mask[:] = 0.0
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
   else:
    tmp = np.copy(hp_in[var])
    cvs[var] = {'min':np.min(tmp),
                'max':np.max(tmp),
                't':-9999,
                'd':tmp}
  #If subgrid network flag is true, ensure that topological principal components don't overwhelm clustering of basins, laura
  keys = list(cvs.keys())
  if hydroblocks_info['channel_initiation']['flag_subgrid'] == True:
   n_pc = 0
   n_no_pc = 0
   for var in keys:
    if 'pc_' in var:n_pc += 1
    else:n_no_pc += 1
   for var in keys:
    if 'pc_' in var:cvs[var]['w'] = (1/(n_no_pc+1))/n_pc
    else:cvs[var]['w'] = (1/(n_no_pc+1))
  else:
   for var in keys:cvs[var]['w'] = 1
    
  (basin_clusters,) = terrain_tools.cluster_basins_hmc_2(basins_wob,cvs,hp_in,ncatchments,1)
  basins_wob[basins_wob!=-9999] = basins_wob[basins_wob!=-9999] + 1
 
 #Calculate the height above nearest drainage area
 print("Computing height above nearest drainage area",flush=True)
 hand = terrain_tools.ttf.calculate_depth2channel(channels_wob,
                                                  basins_wob,
                                                  fdir,
                                                  dem)
 #Fill in hand that is undefined (probably flow direction issues)
 hand[(hand == -9999) & (basins_wob!=-9999)] = 0.0
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
 (tiles,new_hand,tile_position) = terrain_tools.create_basin_tiles_updated(basin_clusters,hand,basins_wob,n_binning,cid,max_nbins)#con HBnew2 entre n_binning y max_nbins va cid
 dict_tiling =  {}
 dict_tiling['tiles']=tiles #laura, oct25
 dict_tiling['tile_position']=tile_position #laura, oct25
 dict_tiling['basin_clusters']=basin_clusters #laura, oct25
 dict_tiling['hand']=hand #laura, oct25
 dict_tiling['basins_wob']=basins_wob #laura, oct25
 dict_tiling['n_binning']=n_binning #laura, oct25
 dict_tiling['max_nbins']=max_nbins #laura, oct25
 pickle.dump(dict_tiling,open('%s/tiling_params_%s.pck'%(hydroblocks_info['input_dir'],cid),'wb'))

 #Assemble river/hillslope database for routing/two-way connectivity
 (db_routing,area_adj,new_hand2) = Build_Hillslope_River_Database(channels_wob,mask,fdir,eares,tiles,hand,basins_wob,basin_clusters,new_hand,db_routing,ubcs,tile_position,db_channels)
 
 #Disagregate land cover
 intraband_clust_vars = hydroblocks_info['hmc_parameters']['intraband_clustering_covariates']
 if 'lc' in intraband_clust_vars:
  intraband_clust_vars.remove('lc')
  intraband_clust_vars=intraband_clust_vars+['lc_w_now','lc_urb_nourb','lc_grass_forest'] #laura
 #Calculate the hrus (kmeans on each tile of each basin)
 cvs = {}
 for var in intraband_clust_vars:
  if var in ['lc_w_now','lc_urb_nourb','lc_grass_forest']:
   if var=='lc_w_now':
    lc_mask=covariates['lc_17']
   elif var=='lc_urb_nourb':
    lc_mask=covariates['lc_13']
   elif var=='lc_grass_forest':
    lc_mask=np.copy(dem)
    lc_mask[:] = 0.0
    if 'lc_4' in covariates:
     lc_mask[covariates['lc_4']==1]=1 #deciduous forest
    if 'lc_2' in covariates:
     lc_mask[covariates['lc_2']==1]=1 #evergreen_forest
    if 'lc_5' in covariates:
     lc_mask[covariates['lc_5']==1]=1 #mixed_forest
    if 'lc_6' in covariates:
     lc_mask[covariates['lc_6']==1]=0.66 #shrub/scrub
    if 'lc_11' in covariates:
     lc_mask[covariates['lc_11']==1]=0.66 #wetlands
    if 'lc_12' in covariates:
     lc_mask[covariates['lc_12']==1]=0.66 #pasture/hay/cultivated_crops
    if 'lc_10' in covariates:
     lc_mask[covariates['lc_10']==1]=0.33 #grassland
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
    
 hrus = terrain_tools.create_hrus_hydroblocks(basin_clusters,tiles,cvs,nclusters,cid) #laura
 hrus[hrus!=-9999] = hrus[hrus!=-9999] - 1
 nhru = np.unique(hrus[hrus!=-9999]).size

 #Save the channel info
 os.system('rm %s/routing_info.pck'%input_dir)
 os.system('rm %s/routing_io.pck'%input_dir)
 os.system('rm %s/routing_mp_connectivity.pck'%input_dir)
 pickle.dump(db_routing,open('%s/routing_info.pck' % input_dir,'wb'))
 pickle.dump(db_routing['i/o'],open('%s/routing_io.pck' % input_dir,'wb'))
 pickle.dump(db_routing['mp_connectivity'],open('%s/routing_mp_connectivity.pck' % input_dir,'wb'))

 #Construct HMC info for creating connections matrix
 HMC_info = {}
 HMC_info['basins'] = basins
 HMC_info['tile_position'] = tile_position
 HMC_info['channel_map'] = channels_wob

 return (hrus.astype(np.float32),nhru,new_hand,HMC_info,covariates,db_channels,new_hand2,
         basins,basin_clusters,hand,tiles,area_adj,tile_position,mask,z_data)

def regular_HMC(hydroblocks_info,basins_wob,eares,dh,covariates,ncatchments,nclusters,db_channels,channels_wob,mask,fdir,hand,db_routing,dem):
 #Clustering the basins
 print("Clustering the basins",flush=True)

 #Set the ncatchments to be at least the number of basins
 ncatchments = min(ncatchments,np.unique(basins_wob)[1:].size)
 subbasin_clustering_cov=hydroblocks_info['hmc_parameters']['subbasin_clustering_covariates']#laura
 
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
 for var in ['width','bankfull','length','area']:
  hp_in[var] = db_channels[var]

 #If subgrid network flag is true, modify hp_in and subbasin_clustering_cov, laura
 if hydroblocks_info['channel_initiation']['flag_subgrid'] == True:
  y = pickle.load(open('%s/pca_subgrid_basins.pck'%(hydroblocks_info['input_dir']),'rb'))['principal_components']
  for pc in range(0,y.shape[1]):
   v = 'pc_%s'%(pc+1)
   hp_in[v] = y[:,pc]
   subbasin_clustering_cov.append(v)
 
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
    lc_mask=np.copy(dem)
    lc_mask[:] = 0.0
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

 #If subgrid network flag is true, ensure that topological principal components don't overwhelm clustering of basins, laura
 keys = list(cvs.keys())
 if hydroblocks_info['channel_initiation']['flag_subgrid'] == True:
  n_pc = 0
  n_no_pc = 0
  for var in keys:
   if 'pc_' in var:n_pc += 1
   else:n_no_pc += 1

  for var in keys:
   if 'pc_' in var:cvs[var]['w'] = (1/(n_no_pc+1))/n_pc
   else:cvs[var]['w'] = (1/(n_no_pc+1))
 else:
  for var in keys:cvs[var]['w'] = 1
 
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
 dict_tiling =  {}
 dict_tiling['tiles']=tiles #laura, oct25
 dict_tiling['tile_position']=tile_position #laura, oct25
 dict_tiling['basin_clusters']=basin_clusters #laura, oct25
 dict_tiling['hand']=hand #laura, oct25
 dict_tiling['basins_wob']=basins_wob #laura, oct25
 dict_tiling['n_binning']=n_binning #laura, oct25
 cid=hydroblocks_info['cid'] #laura, oct25
 dict_tiling['max_nbins']=max_nbins #laura, oct25
 pickle.dump(dict_tiling,open('%s/tiling_params_%s.pck'%(hydroblocks_info['input_dir'],cid),'wb'))

 #Assemble river/hillslope database for routing/two-way connectivity
 (db_routing,area_adj,new_hand2) = Build_Hillslope_River_Database(channels_wob,mask,fdir,eares,tiles,hand,basins_wob,basin_clusters,new_hand,db_routing,ubcs,tile_position,db_channels)

 #Disagregate land cover
 intraband_clust_vars = hydroblocks_info['hmc_parameters']['intraband_clustering_covariates']
 if 'lc' in intraband_clust_vars: 
  intraband_clust_vars.remove('lc')
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
    lc_mask=np.copy(dem)
    lc_mask[:] = 0.0
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
 return (hrus.astype(np.float32),nhru,new_hand,covariates,db_channels,new_hand2,
         basin_clusters,hand,tiles,area_adj,tile_position)

#Function that tracks river network from streams with order 1 to larger streams draining outside of the macroscale polygon (topology=-1), laura
def go_downstream_shreve(channel,topo,shreve):
 if topo[channel]!=-1:
  shreve[topo[channel]]=shreve[topo[channel]]+1
  (shreve)=go_downstream_shreve(topo[channel],topo,shreve)
 return shreve

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

 #If there're not lateral connections (just diagonal) create a single "fake" connection, laura
 if hdst.size == 0:
  hdst = np.array([0])
  horg = np.array([0])
    
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
  if var in ['lats','lons']:continue #laura, ensures that lats and lons don't have to be in covariates for coordinates of channels to be right
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
  if var in ['lats','lons']:continue #laura, ensures that lats and lons don't have to be in covariates for coordinates of channels to be right
  if var in ['dem','fdir','acc']:continue 
  if var in ['WLTSMC','MAXSMC','BB','DRYSMC','QTZ','SATDW','REFSMC','SATPSI','SATDK']: #laura svp
   for depth in covariates[var]: #laura svp
    covariates[var][depth][mask<0]=-9999.0 #laura svp
  else: #laura svp
   covariates[var][mask <= 0] = -9999.0

 #Add the mask_all to the covariates
 covariates['mask_all'] = np.copy(mask_all)
 
 return (covariates,mask,depths) #laura svp returns depths for dataset svp

def Create_Clusters_And_Connections(workspace,wbd,output,input_dir,nhru,info,hydroblocks_info,flag_mod_hmc):
 
 dz=hydroblocks_info['dz'] #laura svp
 #Retrieve some metadata
 metadata = gdal_tools.retrieve_metadata(wbd['files']['mask'])
 resx = 90.0#670.0**0.5#26.0

 if os.path.isdir('%s/experiments/simulations/%s/workspace'%(hydroblocks_info['rdir'],hydroblocks_info['experiment'])) == False:
  print("Creating and curating the covariates",flush=True)
  (covariates,mask,z_data)=Create_and_Curate_Covariates_svp(wbd,hydroblocks_info)
 
 #Determine the HRUs 
 if os.path.isdir('%s/experiments/simulations/%s/workspace'%(hydroblocks_info['rdir'],hydroblocks_info['experiment'])) == False:
  (cluster_ids,nhru,new_hand,HMC_info,covariates,dbc,hand,basins,basin_clusters,hand_org,hbands,area_adj,tile_position) = Compute_HRUs_Semidistributed_HMC(covariates,mask,hydroblocks_info,wbd,resx,input_dir,flag_mod_hmc)

 if os.path.isdir('%s/experiments/simulations/%s/workspace'%(hydroblocks_info['rdir'],hydroblocks_info['experiment'])) == True:
  (cluster_ids,nhru,new_hand,HMC_info,covariates,dbc,hand,basins,basin_clusters,hand_org,hbands,area_adj,tile_position,mask,z_data) = Compute_HRUs_Semidistributed_HMC2(hydroblocks_info,resx,input_dir,flag_mod_hmc)
  hydroblocks_info['network_abstraction']['flag'] = False #laura
  hydroblocks_info['flag_mod_hmc'] = False #laura
    
 flag_mod_hmc = hydroblocks_info['flag_mod_hmc']
 #covariates['hand'] = new_hand
 covariates['hand'] = hand
 hydroblocks_info['nhru'] = nhru
  
 #Create the netcdf file
 if os.path.isdir('%s/experiments/simulations/%s/workspace'%(hydroblocks_info['rdir'],hydroblocks_info['experiment'])) == False:
  file_netcdf = '%s/input_file.nc' % hydroblocks_info['input_dir']#hydroblocks_info['input_file']
  hydroblocks_info['input_fp'] = nc.Dataset(file_netcdf, 'w', format='NETCDF4')
 else:
  file_netcdf = '%s/input_file2.nc' % hydroblocks_info['input_dir']#hydroblocks_info['input_file']
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

 if (hydroblocks_info['network_abstraction']['flag'] == False) and (flag_mod_hmc == False):
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
    cmatrix=Calculate_HRU_Connections_Matrix_HMC_hbands(masked_hband,resx,HMC_info,hydroblocks_info) #laura
    OUTPUT[group_name]=cmatrix #end of laura's modification
 else: #laura
  #Make the output dictionary for the basin
  OUTPUT = {'hru':{},'metadata':metadata,'mask':mask}
    
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

 return (OUTPUT,covariates,hydroblocks_info,z_data)

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
    
 #If connected channel network properties are in covariates for basins, set flag_mod_hmc to True, laura
 if 'shreve' in metadata['hmc_parameters']['subbasin_clustering_covariates']:
  flag_mod_hmc = True
 elif 'large_scale_basins' in metadata['hmc_parameters']['subbasin_clustering_covariates']:
  flag_mod_hmc = True
 #NEED TO ADD CONDITION FOR TOPOLOGICAL INDICES, laura
 else:
  flag_mod_hmc = False

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
  metadata['flag_mod_hmc'] = flag_mod_hmc
  #Prepare model data
  tic = time.time()
  Prepare_Model_Input_Data(metadata,metadata_file)
  flag_network_abst = metadata['network_abstraction']['flag']
  if (flag_network_abst==False) and (flag_mod_hmc==False):
   print("Elapsed time: ",time.time() - tic)
 comm.Barrier()

 #Create enhanced input data file
 #Connect_Cell_Networks(rank,size,cids,edir)
 print('Connect cell networks',flush=True)
 Connect_Cell_Networks_v2(rank,size,cids,edir)
 comm.Barrier()

 #Create downstream channel database for particle tracker routing scheme
 #Create_Downstream_Channels_Database(edir,rank,size,cids,comm)
 #comm.Barrier()

 #Wait until they are all done
 workspace = '%s/workspace' % (edir)
 os.system('mkdir -p %s' % workspace)
    
 #Create topology with connections to other cids, laura
 print('Connect topology',flush=True)
 Topology_Connected(rank,size,cids,edir,comm)
 comm.Barrier()

 #Create self-contained trees of reaches for the domain, laura
 print('Compute large-scale watersheds',flush=True)
 Create_Trees(rank,size,cids,edir,comm,flag_mod_hmc,flag_network_abst)
 comm.Barrier()
    
 #Correct Shreve order (domain-wise), laura
 print('Correct Shreve order',flush=True)
 Correct_Shreve(rank,size,cids,edir)
 comm.Barrier()

 #River network abstraction, laura
 if metadata['network_abstraction']['flag']==True: #laura
  Network_Abstraction(rank,size,cids,edir,comm,metadata)
  comm.Barrier()
        
 #With all the domain-wise variables computed, perform HMC-2step, laura
 if (metadata['network_abstraction']['flag']== True) or (flag_mod_hmc == True): #laura
  flag_replace=True
  for cid in cids[rank::size]:
   metadata['cid'] = cid
   metadata['input_dir'] = "%s/%d" % (edir,cid)
   metadata['workspace'] = "%s/data/cids/%d" % (rdir,cid)
   metadata['flag_mod_hmc'] = flag_mod_hmc
   #Prepare model data
   Prepare_Model_Input_Data(metadata,metadata_file)
   comm.Barrier()
 else:
  flag_replace=False

 #Function that replaces stream_network of input_file2.nc by input_file.nc if HMC2 happened, laura
 if flag_replace == True:
  for cid in cids[rank::size]:
   metadata['cid'] = cid
   metadata['input_dir'] = "%s/%d" % (edir,cid)
   metadata['workspace'] = "%s/data/cids/%d" % (rdir,cid)
   Replace_Stream_Network(metadata)
   print("Elapsed time: ",time.time() - tic)
   comm.Barrier()
 
 Finalize_River_Network_Database(rdir,edir,cids,workspace,comm,rank,size)
 comm.Barrier() 
 
 #Postprocess the model input 
 if rank == 0:Postprocess_Input(rdir,edir,cids)

 return

#Replaces -1 in topology for id of channels in other cids, creates a single hdw array and correct inlets, laura
def Topology_Connected(rank,size,cids,edir,comm):
 for cid in cids[rank::size]:
  #Connected topology part, laura
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

 #Creates a single array of headwaters to use in the tree_creation, laura
 #Outlets part
 list_cids=glob.glob('%s/workspace/topology_connected_*'%edir)
 cids=[]
 for i in list_cids:
  cids.append(int(i.split('/')[-1].split('_')[2].split('.pck')[0]))
 inlets_all=[]
 outlets_all=[]
 for cid in cids:
  fp=h5py.File('%s/%s/input_file.nc' % (edir,cid),'r')
  outlets=np.ravel(fp['stream_network']['outlets'][:])
  fp.close()
  outlets_all.extend(outlets)
 outlets_all=np.array(outlets_all)
 outlets_all=np.reshape(outlets_all,(int(len(outlets_all)/4),4))
 outlets_all=outlets_all[outlets_all[:,2]!=-9999]
 
 #Correction of inlets,laura
 un=[]
 for i in range(0,outlets_all.shape[0]):
  string='%s-%s'%(int(outlets_all[i,2]),int(outlets_all[i,3]))
  if string not in un:
   un.append(string)
 inlets=np.zeros((len(un),10),dtype=int)
 inlets[:]=-9999
 for u in range(0,len(un)):
  inlets[u,0]=int(un[u].split('-')[0])
  inlets[u,1]=int(un[u].split('-')[1])
  for o in range(0,outlets_all.shape[0]):
   strng2='%s-%s'%(int(outlets_all[o,2]),(outlets_all[o,3]))
   if strng2==un[u]:
    for col in range(2,6):
     if inlets[u,col]==-9999:
      inlets[u,col]=int(outlets_all[o,0])
      inlets[u,col+4]=int(outlets_all[o,1])
      break
 dict_hdw={}
 dict_hdw['out']=outlets_all
 dict_hdw['in']=inlets
 pickle.dump(dict_hdw,open('%s/workspace/hdw.pck' %(edir),'wb'))
    
 #Replace wrong inlets with corrected ones based on outlets, laura
 for cid in cids[rank::size]:
  inlets_cid=inlets[inlets[:,0]==int(cid)]
  fp=h5py.File('%s/%s/input_file.nc' % (edir,cid),'a')
  del fp['stream_network']['inlets']
  fp['stream_network']['inlets']=inlets_cid
  fp.close()
  
 return

#Creates trees of reaches draining outside of the domain, laura
def Create_Trees(rank,size,cids,edir,comm,flag_mod_hmc,flag_network_abst):
 dict_trees_domain={}
 for cid in cids[rank::size]:
  hdw=pickle.load(open('%s/workspace/hdw.pck'%edir,'rb'))
  topo=pickle.load(open('%s/workspace/topology_connected_%s.pck'%(edir,cid),'rb'))
  indices = [i for i, item in enumerate(topo) if item == '-1']
  for rid in indices:
   dict_trees_domain['%s-%s'%(cid,rid)]=[]
   dict_trees_domain['%s-%s'%(cid,rid)].append('%s-%s'%(cid,rid))
   (dict_trees_domain['%s-%s'%(cid,rid)])=go_upstream(edir,cid,topo,rid,hdw,dict_trees_domain['%s-%s'%(cid,rid)])

 reaches=[]
 for key in list(dict_trees_domain.keys()):
  reaches.extend(dict_trees_domain[key])
 for hdwo in range(0,hdw['out'].shape[0]):#Check if all outlets are correctly included
  if hdw['out'][hdwo,2]==-9999:
   continue
  string='%s-%s'%(int(hdw['out'][hdwo,2]),int(hdw['out'][hdwo,3]))
  if string in reaches:
   string2='%s-%s'%(int(hdw['out'][hdwo,0]),int(hdw['out'][hdwo,1]))
   if string2 not in reaches:
    cid2=int(hdw['out'][hdwo,0])
    rid2=int(hdw['out'][hdwo,1])
    topo2=pickle.load(open('%s/workspace/topology_connected_%s.pck' % (edir,cid2),'rb'))
    for key2 in list(dict_trees_domain.keys()):
     if string in dict_trees_domain[key2]:
      dict=dict_trees_domain[key2]
      break
    (dict)=go_upstream(edir,cid2,topo2,rid2,hdw,dict)
    dict_trees_domain[key2]=dict
 pickle.dump(dict_trees_domain,open('%s/workspace/trees_%s.pck' %(edir,cid),'wb'))
 comm.Barrier()
 
 #Diagnose if trees need to be corrected further
 reaches=[]
 topo=[]
 tree_files=glob.glob('%s/workspace/trees_*'%edir)
 topo_files=glob.glob('%s/workspace/topology_connected*'%edir)
 for file in tree_files:
  data_tree=pickle.load(open(file,'rb'))
  for key in data_tree:
   reaches.extend(data_tree[key])
 for file in topo_files:
  data_topo=pickle.load(open(file,'rb'))
  topo.extend(data_topo)
 if (len(reaches)!=len(topo)) or (len(np.unique(reaches))!=len(topo)):
  print('WARNING: Trees need to be corrected!!!',len(reaches),len(topo),flush=True)

 #Assign number to trees
 list_cids=glob.glob('%s/workspace/trees_*'%edir)
 cids=[]
 count=1
 dict_numb={}
 for i in list_cids:
  cids.append(int(i.split('/')[-1].split('_')[1].split('.pck')[0]))
  data=pickle.load(open(i,'rb'))
  for key in list(data.keys()):
   dict_numb[count]=data[key]
   count+=1
 comm.Barrier()

 for cid in cids[rank::size]:
  fp=h5py.File('%s/%s/input_file.nc' % (edir,cid),'a')
  trees_mp=[]
  for i in range(0,fp['stream_network']['topology'][:].shape[0]):
   string='%s-%s'%(int(cid),int(i))
   for key in list(dict_numb.keys()):
    if string in dict_numb[key]:
     trees_mp.append(key)
     break
  if 'trees_domain' in fp['stream_network'].keys():
   del fp['stream_network']['trees_domain']
  fp['stream_network']['trees_domain']=np.array(trees_mp)
  fp.close()

 #Save data used for abstraction per cid
 data={}
 for cid in cids[rank::size]:
  fp=h5py.File('%s/%s/input_file.nc' % (edir,cid),'a')
  data['acc']=fp['stream_network']['acc'][:]
  data['shreve']=fp['stream_network']['shreve'][:]
  data['length']=fp['stream_network']['length'][:]
  data['tree']=fp['stream_network']['trees_domain'][:]
  data['lat_basin']=fp['parameters']['lats'][1:]
  data['lon_basin']=fp['parameters']['lons'][1:]
  fp.close()
  pickle.dump(data,open('%s/workspace/data_channels_%s.pck' %(edir,cid),'wb'))
  comm.Barrier()
    
 if (flag_mod_hmc == True) or (flag_network_abst==True):
  #Computes the mean lat and lon per tree to perform clustering instead of tree number
  list_cids=glob.glob('%s/workspace/data_channels_*'%edir)
  dict={}
  dict['lats']=[]
  dict['lons']=[]
  dict['tree']=[]
  dict['mean_lat']=[]
  dict['mean_lon']=[]
  dict['index']=[]
  dict['cid']=[]
  count=0
  dict['index'].append(count)
  for i in list_cids:
   d=pickle.load(open(i,'rb'))
   cid=i.split('/')[-1].split('data_channels_')[-1].split('.pck')[0]
   dict['lats'].extend(d['lat_basin'])
   dict['lons'].extend(d['lon_basin'])
   dict['tree'].extend(d['tree'])
   dict['cid'].append(cid)
   count+=len(d['lat_basin'])
   dict['index'].append(count)
  dict['lats']=np.array(dict['lats'])
  dict['lons']=np.array(dict['lons'])
  dict['tree']=np.array(dict['tree'])
  
  dict['mean_lat']=np.copy(dict['lats'])
  dict['mean_lon']=np.copy(dict['lons'])
  for t in np.unique(dict['tree']):
   mean_lat=np.mean(dict['lats'][dict['tree']==t])
   mean_lon=np.mean(dict['lons'][dict['tree']==t])
   dict['mean_lat'][dict['tree']==t]=mean_lat
   dict['mean_lon'][dict['tree']==t]=mean_lon
   
  dict['CID_info']={}
  dict['CID_info']['lats']={}
  dict['CID_info']['lons']={}
  count2=0
  for cid in dict['cid']:
   cid=int(cid)
   dict['CID_info']['lats'][cid]=dict['mean_lat'][dict['index'][count2]:int(dict['index'][count2+1])]
   dict['CID_info']['lons'][cid]=dict['mean_lon'][dict['index'][count2]:int(dict['index'][count2+1])]
   count2+=1
  pickle.dump(dict,open('%s/workspace/mean_lats-lons_trees.pck' %(edir),'wb'))
 return

#Goes upstream to determine reaches belonging to the same tree
def go_upstream(edir,cid,topo,channel_id,hdw,lista):
 string='%s-%s'%(int(cid),int(channel_id))
 if string not in lista:
  lista.append(string)
 if string in topo: #if analyzed channel has channels draining to it
  indices = [i for i, item in enumerate(topo) if item == string] #finds all the reaches that drain to string
  if indices:#if there're reaches draining to string, continue going upstream one by one
   for ind in indices:
    channel_id=ind
    (lista)=go_upstream(edir,cid,topo,channel_id,hdw,lista)
 for hdwi in range(0,hdw['in'].shape[0]):#check if string is part of headwaters inlets
  string2='%s-%s'%(int(hdw['in'][hdwi,0]),int(hdw['in'][hdwi,1]))
  if string==string2: #If analyzed channel is inlet
   indices_hdw = [i for i, item in enumerate(hdw['in'][hdwi,2:6]) if item != -9999]
   for ind_h in indices_hdw:
    cid2=int(hdw['in'][hdwi,ind_h+2])
    channel_id2=int(hdw['in'][hdwi,ind_h+6])
    topo2=pickle.load(open('%s/workspace/topology_connected_%s.pck' % (edir,cid2),'rb'))
    (lista)=go_upstream(edir,cid2,topo2,channel_id2,hdw,lista)
 return lista

#Function that corrects Shreve order of macroscale polygons accounting for connections, laura
def Correct_Shreve(rank,size,cids,edir):
 for cid in cids[rank::size]:
  #Change to integer
  cid = int(cid)
  fp = h5py.File('%s/%s/input_file.nc' % (edir,cid),'r')
  shreve=fp['stream_network']['shreve'][:]
  topology=fp['stream_network']['topology'][:]
  inlets=fp['stream_network']['inlets'][:]
  fp.close()
  inlets_cids=inlets[:,2:6]
  inlet_cids=np.unique(inlets_cids)[np.unique(inlets_cids)!=-9999]
  for cid2 in inlet_cids:
   cid2=int(cid2)
   fp2 = h5py.File('%s/%s/input_file.nc' % (edir,cid2),'r')
   shreve2=fp2['stream_network']['shreve'][:]
   fp2.close()
   for i in range(0,inlets.shape[0]):
    sequence=[]
    a=np.unique(inlets[i,2:6])[np.unique(inlets[i,2:6])!=-9999]
    if a!=cid2: continue
    sequence.append(inlets[i,1])
    (sequence)=go_downstream_shreve_full(int(inlets[i,1]),topology,sequence)
    id_channel_in=np.unique(inlets[i,6:])[np.unique(inlets[i,6:])!=-9999]
    add_shreve=0
    for idchin in id_channel_in:
     add_shreve+=shreve2[idchin]
    if shreve[sequence[0]]==1:
     for s in sequence:
      shreve[int(s)]=shreve[int(s)]-1+add_shreve
    else:
     for s in sequence:
      shreve[int(s)]=shreve[int(s)]+add_shreve
  ##Update input_file with corrected Shreve order
  with h5py.File('%s/%s/input_file.nc' % (edir,cid),'r+') as fp3:
   if 'shreve' in fp3['stream_network'].keys():
    del fp3['stream_network']['shreve']
   fp3['stream_network']['shreve'] = shreve  
   fp3.close()
 return

#Determines the sequence of channels from upstream all the way down to the end of the path, laura
def go_downstream_shreve_full(channel,topo,sequence):
 if topo[channel]!=-1:
  sequence.append(topo[channel])
  (sequence)=go_downstream_shreve_full(int(topo[channel]),topo,sequence)
 return sequence

#Separates channels into explicit and abstracted based on shreve order or acc. area, laura
def Network_Abstraction(rank,size,cids,edir,comm,metadata):
 #Evaluates percentiles of selected variable depending on type of abstraction
 list_cids=glob.glob('%s/workspace/data_channels_*'%edir)
 dict={}
 var=metadata['network_abstraction']['var']
 dict['var']=[]
 if metadata['network_abstraction']['subreaches']['flag']==True:
  dict['lenght']=[]
 for i in list_cids:
  d=pickle.load(open(i,'rb'))
  dict['var'].extend(d[var])
  if metadata['network_abstraction']['subreaches']['flag']==True:
   dict['lenght'].extend(d['length'])
 n=[]   
 lim_var=[]
 percentiles=[]
 for p in range(60,100,1):
  perc=np.percentile(dict['var'],p)
  lim_var.append(perc)
  percentiles.append(p)
  if metadata['network_abstraction']['subreaches']['flag']==True:
   n_thr=float(metadata['network_abstraction']['subreaches']['max_nmb_subreaches'])
   dx=metadata['network_abstraction']['subreaches']['dx']
   nsr=(np.sum(np.array(dict['lenght'])[np.array(dict['var'])>perc]))/dx
   n.append(nsr)
  elif (metadata['network_abstraction']['nmb_reaches']['flag']==True) or (metadata['network_abstraction']['subreaches']['flag']==False):
   n_thr=float(metadata['network_abstraction']['nmb_reaches']['max_nmb_reaches'])
   nr=np.sum(np.array(dict['var'])>perc)
   n.append(nr)
 dict['percentile']=percentiles
 dict['lim_var']=lim_var
 if metadata['network_abstraction']['subreaches']['flag']==True:
  dict['nmb_subreaches']=n
 else:
  dict['nmb_reaches']=n
 for a in range(0,len(n)):
  if n[a] > n_thr:continue
  else:
   thr=lim_var[a]
   break
 dict['threshold_var']=thr
 pickle.dump(dict,open('%s/workspace/percentile_analysis.pck' %(edir),'wb'))

 abst_mask=[] #0=abstracted, 1=explicit
 for cid in cids[rank::size]:
  fp=h5py.File('%s/%s/input_file.nc' % (edir,cid),'a')
  if 'explicit_reach' in fp['stream_network'].keys():
   del fp['stream_network']['explicit_reach']
  for i in range(0,fp['stream_network']['shreve'][:].shape[0]):
   if fp['stream_network'][var][i]>=thr:
    abst_mask.append(1)
   else:abst_mask.append(0)
  abst_mask=np.array(abst_mask,dtype=int)
  fp['stream_network']['explicit_reach']=abst_mask
  fp.close()    
    
 #This part is determining the drainage network for each reach (needs to be checked), laura
 dict_drainage={}
 for cid in cids[rank::size]:
  hdw=pickle.load(open('%s/workspace/hdw.pck'%edir,'rb'))
  topo=pickle.load(open('%s/workspace/topology_connected_%s.pck'%(edir,cid),'rb'))
  indices = np.linspace(0,len(topo)-1,len(topo))
  for rid in indices:
   dict_drainage['%s-%s'%(int(cid),int(rid))]=[]
   dict_drainage['%s-%s'%(int(cid),int(rid))].append('%s-%s'%(int(cid),int(rid)))
   (dict_drainage['%s-%s'%(int(cid),int(rid))])=go_upstream(edir,cid,topo,rid,hdw,dict_drainage['%s-%s'%(int(cid),int(rid))])

 reaches=[]
 for key in list(dict_drainage.keys()):
  reaches.extend(dict_drainage[key])
 for hdwo in range(0,hdw['out'].shape[0]):#Check if all outlets are correctly included
  if hdw['out'][hdwo,2]==-9999:
   continue
  string='%s-%s'%(int(hdw['out'][hdwo,2]),int(hdw['out'][hdwo,3]))
  if string in reaches:
   string2='%s-%s'%(int(hdw['out'][hdwo,0]),int(hdw['out'][hdwo,1]))
   if string2 not in reaches:
    cid2=int(hdw['out'][hdwo,0])
    rid2=int(hdw['out'][hdwo,1])
    topo2=pickle.load(open('%s/workspace/topology_connected_%s.pck' % (edir,cid2),'rb'))
    for key2 in list(dict_drainage.keys()):
     if string in dict_drainage[key2]:
      dict=dict_drainage[key2]
      break
    (dict)=go_upstream(edir,cid2,topo2,rid2,hdw,dict)
    dict_drainage[key2]=dict
 pickle.dump(dict_drainage,open('%s/workspace/drainage_full_%s.pck' %(edir,cid),'wb'))
 comm.Barrier()
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

#Function that computes topological, morphologic, graph-theory, and channel-feature indices for subgrid network
#Function also reduces dimensionality by using PCA that accounts for user-defined fraction of the total variance
#laura
def Subgrid_Indices(channels_wob_sg,topology_sg,basins_wob,db_channels_sg,thr_var):
 length_sg = db_channels_sg['length']
 width_sg = db_channels_sg['width']
 slope_sg = db_channels_sg['slope']
 acc_sg = db_channels_sg['acc']
 bankfull_sg = db_channels_sg['bankfull']
    
 list_nchannels=[]#0
 list_avrg_id_ilngth=[] #1
 list_avrg_id_wdth=[] #2
 list_avrg_id_slp=[] #3
 list_avrg_shrtst_pth_wdth=[] #4
 list_avrg_shrtst_pth_slp=[] #5
 list_avrg_shrtst_pth_lngth=[] #6
 list_grc_uw=[] #7
 list_g_uw=[] #8
 list_g_ilngth=[] #9
 list_g_wdth=[] #10
 list_g_slp=[] #11
 list_g_eff=[] #12
 list_total_lngth=[] #13
 list_total_acc=[] #14
 list_avrg_lngth=[] #15
 list_avrg_acc=[] #16
 list_avrg_wdth=[] #17
 list_avrg_slpe=[] #18
 list_avrg_bnkfll=[] #19
 list_drng_dnsty=[] #20
    
 for b in np.unique(basins_wob):
  if b==-9999:continue
  else:
   G = nx.DiGraph()
   nds = np.unique(channels_wob_sg[basins_wob==b])[np.unique(channels_wob_sg[basins_wob==b])!=0]-1
   init_out_nd=-1
   for nd in nds:
    if topology_sg[nd]==-1:
     G.add_node(nd)
     G.add_edge(nd,init_out_nd,length=length_sg[nd],width=width_sg[nd],ilength=1/length_sg[nd],slope=slope_sg[nd],acc=acc_sg[nd],bf=bankfull_sg[nd])
     init_out_nd = init_out_nd-1
    else:
     G.add_node(nd)
     G.add_node(topology_sg[nd])
     G.add_edge(nd,topology_sg[nd],length=length_sg[nd],width=width_sg[nd],ilength=1/length_sg[nd],slope=slope_sg[nd],acc=acc_sg[nd],bf=bankfull_sg[nd])
        
   list_nchannels.append(len(G.edges))
   id_ilngth=(G.degree(weight='ilength'))
   id_width=(G.degree(weight='width'))
   id_slope=(G.degree(weight='slope'))
        
   list_id_ilngth=([val for (node, val) in id_ilngth])
   list_id_width=([val for (node, val) in id_width])
   list_id_slope=([val for (node, val) in id_slope])
        
   avrg_id_ilngth=sum(list_id_ilngth)/G.number_of_nodes()
   avrg_id_wdth=sum(list_id_width)/G.number_of_nodes()
   avrg_id_slp=sum(list_id_slope)/G.number_of_nodes()

   list_avrg_id_ilngth.append(avrg_id_ilngth)
   list_avrg_id_wdth.append(avrg_id_wdth)
   list_avrg_id_slp.append(avrg_id_slp)
        
   if nx.is_weakly_connected(G)==True:
    avrg_shrtst_pth_wdth=nx.average_shortest_path_length(G,weight='width')
    avrg_shrtst_pth_slp=nx.average_shortest_path_length(G,weight='slope')
    avrg_shrtst_pth_lngth=nx.average_shortest_path_length(G,weight='lenght')
    list_avrg_shrtst_pth_wdth.append(avrg_shrtst_pth_wdth)
    list_avrg_shrtst_pth_slp.append(avrg_shrtst_pth_slp)
    list_avrg_shrtst_pth_lngth.append(avrg_shrtst_pth_lngth)
   else:
    shrtst_pth_lnght_w = []
    shrtst_pth_lnght_s = []
    shrtst_pth_lnght_l = []
    sub_graphs = nx.weakly_connected_components(G)
    for i, sg in enumerate(sub_graphs):
     SG = G.subgraph(sg).copy()
     shrtst_pth_lnght_w.append(nx.average_shortest_path_length(SG,weight='width'))
     shrtst_pth_lnght_s.append(nx.average_shortest_path_length(SG,weight='slope'))
     shrtst_pth_lnght_l.append(nx.average_shortest_path_length(SG,weight='lenght'))
    list_avrg_shrtst_pth_wdth.append(np.mean(shrtst_pth_lnght_w))
    list_avrg_shrtst_pth_slp.append(np.mean(shrtst_pth_lnght_s))
    list_avrg_shrtst_pth_lngth.append(np.mean(shrtst_pth_lnght_l))
            
   grc_uw=nx.global_reaching_centrality(G)
   list_grc_uw.append(grc_uw)
        
   #convert to undirected graph to compute sprectral properties on a symmetrical matrix              
   Gud=G.to_undirected()
   #compute spectrum of adjacency matrix
   spctrm_uw=nx.adjacency_spectrum(Gud)
   spctrm_ilngth=nx.adjacency_spectrum(Gud, weight='ilength')
   spctrm_wdth=nx.adjacency_spectrum(Gud, weight='width')
   spctrm_slp=nx.adjacency_spectrum(Gud, weight='slope')
                
   #Spectral properties
   ##Spectral gap
   g_uw=abs(np.max(spctrm_uw))
   g_ilngth=abs(np.max(spctrm_ilngth))
   g_wdth=abs(np.max(spctrm_wdth))
   g_slp=abs(np.max(spctrm_slp))
        
   list_g_uw.append(g_uw)
   list_g_ilngth.append(g_ilngth)
   list_g_wdth.append(g_wdth)
   list_g_slp.append(g_slp)
                
   #Non-spectral properties
   ##Global efficiency
   g_eff=nx.global_efficiency(Gud)
   list_g_eff.append(g_eff)
        
   #Channel features
   list_total_lngth.append(G.size(weight="length"))
   list_total_acc.append(G.size(weight="acc"))
   list_avrg_lngth.append(list_total_lngth[-1]/len(G.edges))
   list_avrg_acc.append(list_total_acc[-1]/len(G.edges))
   list_avrg_wdth.append((G.size(weight="width"))/len(G.edges))
   list_avrg_slpe.append((G.size(weight="slope"))/len(G.edges))
   list_avrg_bnkfll.append((G.size(weight="bf"))/len(G.edges))
        
   #River Network morphology
   list_drng_dnsty.append(list_total_lngth[-1]/list_total_acc[-1])
    
   list_nchannels[0] = np.mean(list_nchannels[1:])
   list_avrg_id_ilngth[0] = np.mean(list_avrg_id_ilngth[1:])
   list_avrg_id_wdth[0] = np.mean(list_avrg_id_wdth[1:])
   list_avrg_id_slp[0] = np.mean(list_avrg_id_slp[1:])
   list_avrg_shrtst_pth_wdth[0] = np.mean(list_avrg_shrtst_pth_wdth[1:])
   list_avrg_shrtst_pth_slp[0] = np.mean(list_avrg_shrtst_pth_slp[1:])
   list_avrg_shrtst_pth_lngth[0] = np.mean(list_avrg_shrtst_pth_lngth[1:])
   list_grc_uw[0] = np.mean(list_grc_uw[1:])
   list_g_uw[0] = np.mean(list_g_uw[1:])
   list_g_ilngth[0] = np.mean(list_g_ilngth[1:])
   list_g_wdth[0] = np.mean(list_g_wdth[1:])
   list_g_slp[0] = np.mean(list_g_slp[1:])
   list_g_eff[0] = np.mean(list_g_eff[1:])
   list_total_lngth[0] = np.mean(list_total_lngth[1:])
   list_total_acc[0] = np.mean(list_total_acc[1:])
   list_avrg_lngth[0] = np.mean(list_avrg_lngth[1:])
   list_avrg_acc[0] = np.mean(list_avrg_acc[1:])
   list_avrg_wdth[0] = np.mean(list_avrg_wdth[1:])
   list_avrg_slpe[0] = np.mean(list_avrg_slpe[1:])
   list_avrg_bnkfll[0] = np.mean(list_avrg_bnkfll[1:])
   list_drng_dnsty[0] = np.mean(list_drng_dnsty[1:])

 #Create array of inidices
 metrics=np.zeros((len(list_nchannels),21))
 metrics[:,0]=list_nchannels
 metrics[:,1]=list_avrg_id_ilngth
 metrics[:,2]=list_avrg_id_wdth
 metrics[:,3]=list_avrg_id_slp
 metrics[:,4]=list_avrg_shrtst_pth_wdth
 metrics[:,5]=list_avrg_shrtst_pth_slp
 metrics[:,6]=list_avrg_shrtst_pth_lngth
 metrics[:,7]=list_grc_uw
 metrics[:,8]=list_g_uw
 metrics[:,9]=list_g_ilngth
 metrics[:,10]=list_g_wdth
 metrics[:,11]=list_g_slp
 metrics[:,12]=list_g_eff
 metrics[:,13]=list_total_lngth
 metrics[:,14]=list_total_acc
 metrics[:,15]=list_avrg_lngth
 metrics[:,16]=list_avrg_acc
 metrics[:,17]=list_avrg_wdth
 metrics[:,18]=list_avrg_slpe
 metrics[:,19]=list_avrg_bnkfll
 metrics[:,20]=list_drng_dnsty

 #Standardize the data
 X_std = (metrics - np.mean(metrics,axis=0))/np.std(metrics,axis=0)
    
 #Define the parameters for PCA
 pca = sklearn.decomposition.PCA(n_components=21)
 #Fit the model
 pca.fit(X_std)

 explained_variance = np.cumsum(pca.explained_variance_/np.sum(pca.explained_variance_))
    
 threshold_variance = thr_var
 for i in range(0,explained_variance.shape[0]):
  if explained_variance[i] > threshold_variance:
   n_comp = i+1
   break
            
 #Define the parameters for PCA with n_comp accounting for threshold of variance
 pca = sklearn.decomposition.PCA(n_components=n_comp)
 #Fit the model
 pca.fit(X_std)
 #Transform the data
 Y = pca.transform(X_std)
    
 return Y

def Replace_Stream_Network(metadata):
 # Open the source and destination files
 source_file = h5py.File('%s/input_file.nc'%metadata['input_dir'], 'r')
 destination_file = h5py.File('%s/input_file2.nc'%metadata['input_dir'], 'a')
 shreve = source_file['stream_network']['shreve'][:]
 inlets = source_file['stream_network']['inlets'][:]
 outlets = source_file['stream_network']['outlets'][:]
 trees = source_file['stream_network']['trees_domain'][:]
 if metadata['network_abstraction']['flag']==True:
  explicit = source_file['stream_network']['explicit_reach'][:]
    
 del destination_file['stream_network']['shreve']

 destination_file['stream_network']['shreve'] = shreve
 destination_file['stream_network']['inlets'] = inlets
 destination_file['stream_network']['outlets'] = outlets
 destination_file['stream_network']['trees_domain'] = trees
 if metadata['network_abstraction']['flag']==True:
  destination_file['stream_network']['explicit_reach'] = explicit

 source_file.close()
 destination_file.close()
    
 os.system('mv %s/input_file.nc %s/input_file3.nc'%(metadata['input_dir'],metadata['input_dir']))
 os.system('mv %s/input_file2.nc %s/input_file.nc'%(metadata['input_dir'],metadata['input_dir']))
    
 return
