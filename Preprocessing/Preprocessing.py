import warnings
warnings.filterwarnings('ignore')
import sys
sys.path.append('Tools')
import pickle
import datetime
import numpy as np
import scipy.sparse as sparse
import scipy.stats as stats
import model_tools as mt
import os
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
 input_dir = workspace

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
  'WLTSMC':'%s/theta1500_latlon.tif' % workspace,
  'TEXTURE_CLASS':'%s/texture_class_latlon.tif' % workspace,
  'MAXSMC':'%s/thetas_latlon.tif' % workspace,
  'BB':'%s/bb_latlon.tif' % workspace,
  'DRYSMC':'%s/thetar_latlon.tif' % workspace,
  'QTZ':'%s/qtz_latlon.tif' % workspace,
  'SATDW':'%s/dsat_latlon.tif' % workspace,
  'REFSMC':'%s/theta33_latlon.tif' % workspace,
  'mask':'%s/mask_latlon.tif' % workspace,
  'SATPSI':'%s/psisat_latlon.tif' % workspace,
  'lc':'%s/lc_latlon.tif' % workspace,
  'F11':'%s/f11_latlon.tif' % workspace,
  'SATDK':'%s/ksat_latlon.tif' % workspace,
  'dem':'%s/dem_latlon.tif' % workspace,
  'acc':'%s/acc_latlon.tif' % workspace,
  'demns':'%s/demns_latlon.tif' % workspace,
  'sand':'%s/sand_latlon.tif' % workspace,
  'clay':'%s/clay_latlon.tif' % workspace,
  'silt':'%s/silt_latlon.tif' % workspace,
  'om':'%s/om_latlon.tif' % workspace,
  'bare30':'%s/bare30_latlon.tif' % workspace,
  'water30':'%s/water30_latlon.tif' % workspace,
  'tree30':'%s/tree30_latlon.tif' % workspace,
  'irrig_land':'%s/irrig_land_latlon.tif' % workspace,
  'dbedrock':'%s/dbedrock_latlon.tif' % workspace,
  'lstmean':'%s/lstmean_latlon.tif' % workspace,
  'lststd':'%s/lststd_latlon.tif' % workspace
  }
 if hydroblocks_info['water_management']['hwu_agric_flag']:
   wbd['files']['irrig_land'] = '%s/irrig_land_latlon.tif' % workspace
   wbd['files']['start_growing_season'] = '%s/start_growing_season_latlon.tif' % workspace
   wbd['files']['end_growing_season']   = '%s/end_growing_season_latlon.tif' % workspace

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
 (output,covariates) = Create_Clusters_And_Connections(workspace,wbd,output,input_dir,nhru,info,hydroblocks_info)

 #Extract the meteorological forcing
 print("Preparing the meteorology",flush=True)
 Prepare_Meteorology_Semidistributed(workspace,wbd,output,input_dir,info,hydroblocks_info,covariates)

 #Extract the water use demands
 print("Preparing the water use",flush=True)
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
 grp.dx = 26.0#25.0#metadata['resx'] #UPDATE WITH DEM!

 #Write out the mapping
 hru_map = np.copy(output['hru_map'])
 hru_map[np.isnan(hru_map) == 1] = -9999.0
 file_ca = '%s/hru_mapping_latlon.tif' % workspace
 metadata['nodata'] = -9999.0
 gdal_tools.write_raster(file_ca,metadata,hru_map)

 #Write out the hand map
 hand_map = np.copy(output['hand_map'])
 hand_map[np.isnan(hand_map) == 1] = -9999.0
 file_ca = '%s/hand_latlon.tif' % workspace
 metadata['nodata'] = -9999.0
 gdal_tools.write_raster(file_ca,metadata,hand_map)

 #Write out the basin map
 basin_map = np.copy(output['basin_map'])
 basin_map[np.isnan(basin_map) == 1] = -9999.0
 file_ca = '%s/basins_latlon.tif' % workspace
 metadata['nodata'] = -9999.0
 gdal_tools.write_raster(file_ca,metadata,basin_map)

 #Write out the channels
 channel_map = np.copy(output['channel_map'])
 channel_map[np.isnan(channel_map) == 1] = -9999.0
 file_ca = '%s/channel_mapping_latlon.tif' % workspace
 metadata['nodata'] = -9999.0
 gdal_tools.write_raster(file_ca,metadata,channel_map)

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
        'm','hand','y_aspect','x_aspect']

 if hydroblocks_info['water_management']['hwu_agric_flag']:
  for var in ['centroid_lats', 'centroid_lons', 'irrig_land', 'start_growing_season', 'end_growing_season']:
    vars.append(var)

 for var in vars:
  grp.createVariable(var,'f4',('hru',))#,zlib=True)
  grp.variables[var][:] = data['hru'][var]

 if hydroblocks_info['water_management']['hwu_flag']:
  grp.createVariable('hru_min_dist','f4',('hru','hru'))#,zlib=True)
  grp.variables['hru_min_dist'][:] = data['hru']['hru_min_dist']

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

def Compute_HRUs_Semidistributed_HMC(covariates,mask,hydroblocks_info,wbd,eares):

 #Define the parameters for the hierarchical multivariate clustering
 ncatchments = hydroblocks_info['hmc_parameters']['number_of_characteristic_subbasins']
 dh = hydroblocks_info['hmc_parameters']['average_height_difference_between_bands']
 nclusters = hydroblocks_info['hmc_parameters']['number_of_intraband_clusters']

 #Bring out the mask_all
 mask_all = covariates['mask_all']

 #Pre-process DEM
 dem = covariates['dem']
 #Remove pits in dem
 print("Removing pits in dem",flush=True)
 demns = terrain_tools.ttf.remove_pits_planchon(dem,eares)
 covariates['demns'] = demns
 area_all = covariates['acc']
  
 #Calculate slope and aspect
 print("Calculating slope and aspect",flush=True)
 res_array = np.copy(demns)
 res_array[:] = eares
 (slope,aspect) = terrain_tools.ttf.calculate_slope_and_aspect(np.flipud(demns),res_array,res_array)
 slope = np.flipud(slope)
 aspect = np.flipud(aspect)

 #Compute accumulated area
 #m2 = np.copy(demns)
 m2 = np.copy(mask_all)
 m2[m2 > 0] = 1
 mall = np.copy(m2)
 mall[m2 <= 0] = 0
 mall = mall.astype(np.bool)
 #m2[:] = 1
 print("Calculating accumulated area",flush=True)
 (area,fdir) = terrain_tools.ttf.calculate_d8_acc(demns,mask_all,eares)
 #(area,fdir) = terrain_tools.ttf.calculate_d8_acc(demns,mall,eares)
 area_all = area[:] #This could be defined for entire domain instead

 #Calculate channel initiation points (2 parameters)
 C = area/eares*slope**2
 #ipoints = ((C > 200) & (area > 10**5)).astype(np.int32)
 #ipoints = ((C > 100) & (area > 10**5)).astype(np.int32)
 #ipoints = ((C > 100) & (area > 10**4)).astype(np.int32)
 ipoints = ((C > 25) & (area > 10**4)).astype(np.int32)
 #ipoints = (C > 100).astype(np.int32)
 ipoints[ipoints == 0] = -9999

 #Create area for channel delineation
 (ac,fdc) = terrain_tools.ttf.calculate_d8_acc_wipoints(demns,mask,ipoints,eares)
 (ac_all,fdc_all) = terrain_tools.ttf.calculate_d8_acc_wipoints(demns,mall,ipoints,eares)
 ac[(ac != 0) & (mask == 1)] = area[(ac != 0) & (mask == 1)]
 ac_all[(ac_all != 0) & (mall == 1)] = area[(ac_all != 0) & (mall == 1)]
 ac_all = np.ma.masked_array(ac_all,ac_all<=0)

 #Compute the channels
 print(np.unique(mask_all[mask]))
 print("Defining channels",flush=True)
 #(channels,channels_wob,channel_topology,tmp1) = terrain_tools.ttf.calculate_channels_wocean_wprop(ac,10**4,10**4,fdir,mask)
 (channels,channels_wob,channel_topology,tmp1,crds) = terrain_tools.ttf.calculate_channels_wocean_wprop_wcrds(ac,10**4,10**4,fdc,mask,np.flipud(covariates['lats']),covariates['lons'])
 #(tmp1,tmp2,tmp3,shreve_order) = terrain_tools.ttf.calculate_channels_wocean_wprop(ac_all,10**4,10**4,fdc,mask_all)
 (tmp1,tmp2,tmp3,shreve_order) = terrain_tools.ttf.calculate_channels_wocean_wprop(ac_all,10**4,10**4,fdc,mask_all)
 #Curate channel_topology
 channel_topology = channel_topology[channel_topology != -9999]
 
 #Compute channel properties
 #db_channels = terrain_tools.calculate_channel_properties(channels_wob,channel_topology,slope,eares,mask,area_all)
 #import matplotlib.pyplot as plt
 #channels_wob = np.ma.masked_array(channels_wob,channels_wob<=0)
 #plt.imshow(channels_wob)
 #shreve_order = np.ma.masked_array(shreve_order,shreve_order<=0)
 #plt.imshow(np.log(ac_all))
 #plt.imshow(channels_wob)
 #plt.imshow(shreve_order)
 #plt.show()
 #exit()

 #If the dem is undefined then set to undefined
 channels[dem == -9999] = -9999

 #Determine inlets/outlets
 os.system('mkdir -p routing')
 #print(np.unique(mask))
 #print(np.unique(area))
 #print(np.unique(mask_all))
 #print(np.unique(channels_wob))
 db_routing = {}
 #db_routing['i/o'] = terrain_tools.calculate_inlets_oulets(channels_wob,fdir,area,mask,np.flipud(covariates['lats']),covariates['lons'],mask_all,area_all)
 db_routing['i/o'] = terrain_tools.calculate_inlets_oulets(channels_wob,fdir,area,mask,np.flipud(covariates['lats']),covariates['lons'],mask_all,area_all)

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
 basins = terrain_tools.ttf.delineate_basins(channels,m2,fdir)
 basins_wob = terrain_tools.ttf.delineate_basins(channels_wob,mask,fdir)
 
 #Compute channel properties
 db_channels = terrain_tools.calculate_channel_properties(channels_wob,channel_topology,slope,eares,mask,area_all,basins_wob,shreve_order)

 #Burn bankfull depth into DEM (the values are adjusted to ensure it works with the DEM res)
 demns_adj = np.copy(demns)
 for i in range(channels_wob.shape[0]):
    for j in range(channels_wob.shape[1]):
        if channels_wob[i,j] > 0:
            ic = channels_wob[i,j]
            cwidth = db_channels['width'][ic-1]
            cbankfull = db_channels['bankfull'][ic-1]
            cA = cwidth*cbankfull
            #Determine bankfull to burn in to ensure that A holds
            demns_adj[i,j]=demns[i,j]-cA/eares

 #Calculate the height above nearest drainage area
 print("Computing height above nearest drainage area",flush=True)
 hand = terrain_tools.ttf.calculate_depth2channel(channels_wob,basins_wob,fdir,demns_adj)

 '''#Compute the areal coverage of each hand value within the basin
 db_routing['reach_hand_area'] = {}
 for i in range(basins_wob.shape[0]):
  for j in range(basins_wob.shape[1]):
   basin = basins_wob[i,j]
   h = hand[i,j]
   if basin <= 0:continue
   if basin not in db_routing['reach_hand_area']:db_routing['reach_hand_area'][basin] = {}
   if h not in db_routing['reach_hand_area'][basin]: db_routing['reach_hand_area'][basin][h] = 0.0
   db_routing['reach_hand_area'][basin][h] += eares**2
 
 #Compute channel cross section information
 odb = {'A':0.0*np.ones((len(db_routing['reach_hand_area'].keys()),1000)),
       'P':0.0*np.ones((len(db_routing['reach_hand_area'].keys()),1000)),
       'W':0.0*np.ones((len(db_routing['reach_hand_area'].keys()),1000)),
       'hand':0.0*np.ones((len(db_routing['reach_hand_area'].keys()),1000))}
 for b in db_routing['reach_hand_area']:
  #Define reach length
  c_length = db_channels['length'][b-1]
  #Sort from lowest to highest hand
  c_hand = np.array(list(db_routing['reach_hand_area'][b].keys()))
  c_area = np.array(list(db_routing['reach_hand_area'][b].values()))
  argsort = np.argsort(c_hand)
  c_hand = c_hand[argsort]
  c_area = c_area[argsort]
  #Need to reduce the size of the c_hand array by binning by percentiles (only if necessary)
  if c_hand.size > 1000:
   argsort = np.argsort(c_hand)
   pcts = np.copy(c_hand)
   pcts[argsort] = np.linspace(0,1,c_hand.size)
   (hist, bin_edges) = np.histogram(pcts,bins=100)
   #Update c_hand and c_area
   tmp_hand = []
   tmp_area = []
   for ib in range(bin_edges.size-1):
    if ib == 0:m = (pcts >= bin_edges[ib]) & (pcts <= bin_edges[ib+1])
    else:m = (pcts > bin_edges[ib]) & (pcts <= bin_edges[ib+1])
    tmp_hand.append(np.mean(c_hand[m])) #Compute the mean hand for within the bins
    tmp_area.append(np.sum(c_area[m])) #Compute the area sum for within the bins
   c_hand = np.array(tmp_hand)
   c_area = np.array(tmp_area)
  odb['hand'][b-1,0:c_hand.size] = c_hand[:]
  #Calculate width
  c_width = c_area/c_length
  odb['W'][b-1,0:c_width.size] = c_width[:]
  #Calculate wetted perimeter at each stage
  dP = c_width[0:-1] + 2*np.diff(c_hand)
  P = np.cumsum(dP)
  odb['P'][b-1,0] = 0.0
  odb['P'][b-1,1:P.size+1] = P[:]
  #Calculate wetted cross sectional area at each stage
  dA = c_width[0:-1]*np.diff(c_hand)
  A = np.cumsum(dA)
  odb['A'][b-1,0] = 0.0
  odb['A'][b-1,1:A.size+1] = A[:]
 db_routing['reach_cross_section'] = copy.deepcopy(odb)'''

 #Calculate topographic index
 print("Computing topographic index",flush=True)
 ti = np.copy(area)
 m = (area != -9999) & (slope != -9999) & (slope != 0.0)
 ti[m] = np.log(area[m]/eares/slope[m])
 ti[slope == 0] = 15.0

 # cleanup
 slope[mask != 1] = -9999
 aspect[mask != 1] = -9999
 area[mask != 1] = -9999
 channels[mask != 1] = -9999
 basins[mask != 1] = -9999
 #hand[mask != 1] = -9999
 ti[mask != 1] = -9999
 shreve_order[mask != 1] = -9999

 covariates['slope'] = slope
 covariates['aspect'] = aspect
 covariates['x_aspect'] = np.sin(aspect)
 covariates['y_aspect'] = np.cos(aspect)
 covariates['carea'] = area
 covariates['shreve_order'] = shreve_order
 covariates['hand'] = hand
 #covariates['fdir'] = fdir
 covariates['ti'] = ti

 #Calculate the subbasin properties
 print("Assembling the subbasin properties",flush=True)
 #print(eares,np.unique(basins))
 hp_in = terrain_tools.calculate_basin_properties_updated(basins,eares,covariates,hydroblocks_info['hmc_parameters']['subbasin_clustering_covariates'])

 #Clustering the basins
 print("Clustering the basins",flush=True)
 #Assemble input data
 cvs = {}
 for var in hydroblocks_info['hmc_parameters']['subbasin_clustering_covariates']:
  tmp = np.copy(hp_in[var])
  cvs[var] = {'min':np.min(tmp),
              'max':np.max(tmp),
              't':-9999,
              'd':tmp}
 
 (basin_clusters,) = terrain_tools.cluster_basins_updated(basins,cvs,hp_in,ncatchments)
  
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
  if missing_ratio > 0.9 : 
    basin_clusters[m] = -9999
    basins[m] = -9999
 if len(basin_clusters[basin_clusters != -9999]) < 1 : 
  exit('Error_basin_clustering: hand_full_of_nans %s' % hydroblocks_info['cid']) 
 
 #Divide each subbasin into height bands
 (tiles,new_hand,tile_position) = terrain_tools.create_basin_tiles(basin_clusters,hand,basins,dh)

 #Disagregate land cover
 intraband_clust_vars = hydroblocks_info['hmc_parameters']['intraband_clustering_covariates']
 if 'lc' in intraband_clust_vars: 
  intraband_clust_vars.remove('lc')
  disag = [i for i in covariates.keys() if 'lc_' in i]
  intraband_clust_vars = intraband_clust_vars + disag

 #Calculate the hrus (kmeans on each tile of each basin)
 cvs = {}
 #print(intraband_clust_vars)
 #print(covariates.keys())
 #exit()
 for var in intraband_clust_vars:
  cvs[var] = {'min':np.min(covariates[var][covariates[var]!=-9999]),
              'max':np.max(covariates[var][covariates[var]!=-9999]),
              't':-9999,
              'd':covariates[var]}
 #print("Clustering the height bands into %d clusters" % nclusters)
 hrus = terrain_tools.create_hrus_hydroblocks(basin_clusters,tiles,cvs,nclusters)
 hrus[hrus!=-9999] = hrus[hrus!=-9999] - 1
 nhru = np.unique(hrus[hrus!=-9999]).size

 #Calculate histogram of travel distances per hru
 t2c = terrain_tools.ttf.calculate_distance2channel(channels_wob,mask,fdir,eares)
 uhrus = np.unique(hrus)
 uhrus = uhrus[uhrus != -9999]
 bins = np.linspace(0,100,101)
 uhs = []
 for hru in uhrus:
   m = hrus == hru
   hist = np.histogram(t2c[m]/0.1/3600.0,bins=bins,density=True)[0]
   #hist = hist/np.sum(hist)
   uhs.append(hist)
 uhs = {'data':np.array(uhs),'bins':bins[:]}
 db_routing['uh_per_hru'] = uhs

 #Calculate averaged hand per hru per basin
 hand_tmp = {}
 new_hand2 = np.copy(hand)
 for i in range(basins_wob.shape[0]):
  for j in range(basins_wob.shape[1]):
   basin = basins_wob[i,j]
   if basin <= 0:continue
   if basin not in hand_tmp:hand_tmp[basin] = {}
   hru = hrus[i,j]
   if hru < 0:continue
   if hru not in hand_tmp[basin]:
    hand_tmp[basin][hru] = {'sum':0.0,'count':0}
   hand_tmp[basin][hru]['sum'] += hand[i,j]
   hand_tmp[basin][hru]['count'] += 1
 for i in range(basins_wob.shape[0]):
  for j in range(basins_wob.shape[1]):
   basin = basins_wob[i,j]
   if basin <= 0:continue
   hru = hrus[i,j]
   if hru < 0:continue
   new_hand2[i,j] = hand_tmp[basin][hru]['sum']/hand_tmp[basin][hru]['count'] 
   
 #THIS WILL PROBABLY CAUSE PROBLEMS WITH INTRABAND CLUSTERING
 #Compute the areal coverage of each hand value within the basin
 db_routing['reach_hand_area'] = {}
 db_routing['reach_hand_hru'] = {}
 db_routing['reach_hru_area'] = {}
 for i in range(basins_wob.shape[0]):
  for j in range(basins_wob.shape[1]):
   basin = basins_wob[i,j]
   #h = hand[i,j]
   h = new_hand2[i,j]
   hru = hrus[i,j]
   if basin <= 0:continue
   if basin not in db_routing['reach_hand_area']:db_routing['reach_hand_area'][basin] = collections.OrderedDict()
   if h not in db_routing['reach_hand_area'][basin]: db_routing['reach_hand_area'][basin][h] = 0.0
   if basin not in db_routing['reach_hand_hru']:db_routing['reach_hand_hru'][basin] = collections.OrderedDict()
   db_routing['reach_hand_area'][basin][h] += eares**2
   db_routing['reach_hand_hru'][basin][h] = hru
 
 #Compute channel cross section information
 odb = {'A':0.0*np.ones((len(db_routing['reach_hand_area'].keys()),100)),
       'P':0.0*np.ones((len(db_routing['reach_hand_area'].keys()),100)),
       'W':0.0*np.ones((len(db_routing['reach_hand_area'].keys()),100)),
       'hand':0.0*np.ones((len(db_routing['reach_hand_area'].keys()),100)),
       'hru':-9999*np.ones((len(db_routing['reach_hand_area'].keys()),100)).astype(np.int32)}
 for b in db_routing['reach_hand_area']:
  #Define reach length
  c_length = db_channels['length'][b-1]
  #Sort from lowest to highest hand
  c_hand = np.array(list(db_routing['reach_hand_area'][b].keys()))
  c_area = np.array(list(db_routing['reach_hand_area'][b].values()))
  c_hru = np.array(list(db_routing['reach_hand_hru'][b].values()))
  argsort = np.argsort(c_hand)
  c_hand = c_hand[argsort]
  c_hru = c_hru[argsort]
  odb['hru'][b-1,0:c_hru.size] = c_hru[:]
  #Burn in a channel depth
  if c_hand.size > 1:
   #1.first remove existing difference between channel and adjacent hand value
   c_hand[1:] = c_hand[1:] - (c_hand[1] - c_hand[0])
   #2.then burn in the channel bankfull depth
   c_hand[1:] = c_hand[1:] + db_channels['bankfull'][b-1] #m
  c_area = c_area[argsort]
  #Update values in dictionary (due to correcting for channel info)
  db_routing['reach_hand_area'][b] = collections.OrderedDict()
  db_routing['reach_hand_hru'][b] = collections.OrderedDict()
  db_routing['reach_hru_area'][b] = collections.OrderedDict()
  for ih in range(c_hand.size):
    db_routing['reach_hand_area'][b][c_hand[ih]] = c_area[ih]
    db_routing['reach_hand_hru'][b][c_hand[ih]] = c_hru[ih]
    db_routing['reach_hru_area'][b][c_hru[ih]] = c_area[ih] #CAUTION: THIS IS PROBLEMATIC FOR INTRABAND CLUSTERING
  #Calculate widths of each HRU/height band
  c_width = c_area/c_length
  if c_width.size > 1:
   #Correct channel width using provided estimates
   c_width_diff = db_channels['width'][b-1] - c_width[0]
   c_width[0] = c_width[0] + c_width_diff 
   #Add the difference to the adjacent HRU
   c_width[1] = c_width[1] - c_width_diff
  #Adjust the areal coverage of all the HRUs/bands
  c_area = c_length*c_width
  #Update the channel depth
  odb['hand'][b-1,0:c_hand.size] = c_hand[:]
  #Calculate width
  #c_width = c_area/c_length
  odb['W'][b-1,0:c_width.size] = c_width[:]
  #Calculate wetted perimeter at each stage
  dP = c_width[0:-1] + 2*np.diff(c_hand)
  P = np.cumsum(dP)
  odb['P'][b-1,0] = 0.0
  odb['P'][b-1,1:P.size+1] = P[:]
  #Calculate wetted cross sectional area at each stage
  dA = np.cumsum(c_width[0:-1])*np.diff(c_hand)
  A = np.cumsum(dA)
  odb['A'][b-1,0] = 0.0
  odb['A'][b-1,1:A.size+1] = A[:]
  #Calculate inundation height at each stage
 db_routing['reach_cross_section'] = copy.deepcopy(odb)

 #Compute areal coverage of each HRU per basin (or hillslope) (THIS WILL PROBABLE BE WRONG )
 #db_routing['reach_hru_area'] = {}
 '''for i in range(basins_wob.shape[0]):
  for j in range(basins_wob.shape[1]):
   basin = basins_wob[i,j]
   if basin <= 0:continue
   if basin not in db_routing['reach_hru_area']:db_routing['reach_hru_area'][basin] = {}
   if hrus[i,j] not in db_routing['reach_hru_area'][basin]: db_routing['reach_hru_area'][basin][hrus[i,j]] = 0.0
   db_routing['reach_hru_area'][basin][hrus[i,j]] += eares**2 #EARES'''
 pickle.dump(db_routing,open('routing_info.pck','wb'))
 pickle.dump(db_routing['i/o'],open('routing_io.pck','wb'))

 #Construct HMC info for creating connections matrix
 HMC_info = {}
 HMC_info['basins'] = basins
 HMC_info['tile_position'] = tile_position
 HMC_info['channel_map'] = channels_wob

 #return (hrus.astype(np.float32),nhru,new_hand,HMC_info,covariates,db_channels,hand,
 return (hrus.astype(np.float32),nhru,new_hand,HMC_info,covariates,db_channels,new_hand2,
         basins)

def Assign_Parameters_Semidistributed(covariates,metadata,hydroblocks_info,OUTPUT,cluster_ids,mask):

 nhru = hydroblocks_info['nhru']
 #Initialize the arrays
 vars = ['area','area_pct','BB','DRYSMC','F11','MAXSMC','REFSMC','SATPSI',
         'SATDK','SATDW','WLTSMC','QTZ','slope','ti','dem','carea','channel',
         'land_cover','soil_texture_class','clay','sand','silt',
         'm','hand','x_aspect','y_aspect']

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
  #Calculate area per hru
  #OUTPUT['hru']['area'][hru] = metadata['resx']**2*idx[0].size
  OUTPUT['hru']['area'][hru] = 26.0**2*idx[0].size #NEED TO FIX WITH DEM
  #Calculate area percentage per hru
  #OUTPUT['hru']['area_pct'][hru] = 100*OUTPUT['hru']['area'][hru]/(metadata['resx']**2*mask[mask].size)
  OUTPUT['hru']['area_pct'][hru] = 100*OUTPUT['hru']['area'][hru]/(26.0**2*mask[mask].size) #NEED TO FIX WITH DEM
  #Soil properties
  #for var in ['BB','DRYSMC','F11','MAXSMC','REFSMC','SATPSI','SATDK','SATDW','WLTSMC','QTZ','clay','sand','silt']:
  for var in ['BB','DRYSMC','F11','MAXSMC','SATPSI','SATDK','SATDW','QTZ','clay','sand','silt']:
   #print(var,np.unique(covariates[var][idx]))
   if var in ['SATDK','SATDW']:
    try:
     OUTPUT['hru'][var][hru] = stats.mstats.hmean(covariates[var][idx])/3600.0/1000.0 #mm/hr -> m/s
    except:
     OUTPUT['hru'][var][hru] = 1.41E-4
   else:
    OUTPUT['hru'][var][hru] = np.mean(covariates[var][idx])
  OUTPUT['hru']['WLTSMC'][hru] = OUTPUT['hru']['MAXSMC'][hru]*(OUTPUT['hru']['SATPSI'][hru]/150)**(1/OUTPUT['hru']['BB'][hru])
  OUTPUT['hru']['REFSMC'][hru] = OUTPUT['hru']['MAXSMC'][hru]*(OUTPUT['hru']['SATPSI'][hru]/3.3)**(1/OUTPUT['hru']['BB'][hru])
  #Average Slope
  OUTPUT['hru']['slope'][hru] = np.nanmean(covariates['slope'][idx])
  #Topographic index
  OUTPUT['hru']['ti'][hru] = np.nanmean(covariates['ti'][idx])
  #DEM
  OUTPUT['hru']['dem'][hru] = np.nanmean(covariates['dem'][idx])
  #HAND
  OUTPUT['hru']['hand'][hru] = np.nanmean(covariates['hand'][idx])
  #Average Catchment Area
  OUTPUT['hru']['carea'][hru] = np.nanmean(covariates['carea'][idx])
  #Channel?
  #OUTPUT['hru']['channel'][hru] = stats.mode(covariates['channels'][idx])[0]
  #Average aspect
  #OUTPUT['hru']['aspect'][hru] = np.arctan(y_aspect/x_aspect)
  #OUTPUT['hru']['aspect'][hru] = np.arctan(x_aspect/y_aspect)
  OUTPUT['hru']['x_aspect'][hru] = np.nanmean(covariates['x_aspect'][idx])
  OUTPUT['hru']['y_aspect'][hru] = np.nanmean(covariates['y_aspect'][idx])

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
 # for hru in np.arange(nhru):
 #  #HRU distance between the centroids of the hru and all the other hrus 
 #  OUTPUT['hru']['hru_min_dist'][hru,:] = mgmt_funcs.calculate_min_distance(hru, nhru, cluster_ids, covariates['lats'], covariates['lons'], OUTPUT['hru']['centroid_lats'], OUTPUT['hru']['centroid_lons'])
 #  OUTPUT['hru']['hru_min_dist'][hru,hru] = 0.0
   
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
 (hdst,horg) = Calculate_HRU_Connections_Matrix_HMC_workhorse(cluster_ids,nhru,dx,tile_position,
               basins,ivc,irc,ibc)

 #Prepare the sparse matrix
 cmatrix = sparse.coo_matrix((np.ones(hdst.size),(horg,hdst)),shape=(nhru,nhru),dtype=np.float32)
 cmatrix = cmatrix.tocsr()

 #Prepare length, width, and ksat matrices
 wmatrix = cmatrix.copy()
 wmatrix[:] = dx*wmatrix[:]

 #Prepare output dictionary
 cdata = {'width':wmatrix.T,}

 return cdata

@numba.jit(nopython=True,cache=True)
def Calculate_HRU_Connections_Matrix_HMC_workhorse(cluster_ids,nhru,dx,tile_position,basins,
    ivc,irc,ibc):

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

def Create_and_Curate_Covariates(wbd,hydroblocks_info):

 covariates = {}
 #Read in and curate all the covariates
 for file in wbd['files']:
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
  #covariates[var][mask <= 0] = -9999.0
  mask1 = (np.isinf(covariates[var]) == 0) & (np.isnan(covariates[var]) == 0) 
  mask0 = (np.isinf(covariates[var]) == 1) | (np.isnan(covariates[var]) == 1)
  covariates[var][mask0] = -9999.0# stats.mode(covariates[var][mask1])[0][0]

 #Set everything that is -9999 to the mean
 for var in covariates:

  m2 = ( mask > 0 ) & (covariates[var] != -9999.0)
  missing_ratio = 1.0 - np.sum(m2)/float(np.sum(mask))
  if missing_ratio > 0.99 : 
   print("Warning: Covariate %s in catchment %s has %.2f %% of nan's" % (var,hydroblocks_info['cid'],100*missing_ratio)) # Noemi insert
   if var == 'lc': 
    mlc = (covariates[var] == -9999) & mask
    covariates[var][mlc] = 17  # Water
   if var in ['dem','fdir','sand','clay','silt','TEXTURE_CLASS','dbedrock']:
    exit('Error_clustering: %s_full_of_nans %s' % (var,hydroblocks_info['cid']))
   #else: sys.stderr.write("Error_clustering: variable %s has %.2f %% of nan's" % (var,100*missing_ratio))

  #print var
  '''if var not in ['mask']:
   if var in ['fdir','nlcd','TEXTURE_CLASS','lc','irrig_land','bare30','water30','tree30','start_growing_season','end_growing_season']: 
    covariates[var] = spatial_imputation(covariates[var],-9999.0,'nearest')  
   else:
    #covariates[var] = spatial_imputation(covariates[var],-9999.0,'nearest') #faster
    covariates[var] = spatial_imputation(covariates[var],-9999.0,'linear')'''
  if var not in ['mask',]:
   if var in ['nlcd','TEXTURE_CLASS','lc','irrig_land','bare30','water30','tree30','start_growing_season','end_growing_season']: 
    covariates[var][covariates[var] == -9999.0] = stats.mode(covariates[var][covariates[var] != -9999.0])[0][0]
   else:
    covariates[var][covariates[var] == -9999.0] = np.mean(covariates[var][covariates[var] != -9999.0])

 #Set everything outside of the mask to -9999
 for var in covariates:
  if var in hydroblocks_info['hmc_parameters']['subbasin_clustering_covariates']:continue
  if var == 'dem':continue 
  covariates[var][mask <= 0] = -9999.0

 #Add the mask_all to the covariates
 covariates['mask_all'] = np.copy(mask_all)
 
 return (covariates,mask)

'''def spatial_imputation(array,fill_val,method):
 pos_fill = (array==fill_val)
 if np.sum(pos_fill) == 0: return array
 
 pos_mask = np.zeros(array.shape)
 pos_mask[pos_fill] = 1.0
 from skimage.morphology import binary_dilation, square
 mask_good = binary_dilation(pos_mask.astype(np.uint8),square(5))
 
 array[pos_fill]=np.nan
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
 filled  = griddata((x1, y1), newarr.ravel(),(xf, yf), method=method)
 #import matplotlib.pyplot as plt
 #plt.imshow(array)
 #plt.show()
 #print array
 #plot_data(array) 
 array[array.mask] = filled  # replace fill vals
 #plt.imshow(array)
 #plt.show()
 if method in ['linear','cubic']:
  # fill the missing boundaries in nearest
  array = np.ma.masked_invalid(array)
  xf = xx[array.mask] # fill vals
  yf = yy[array.mask] #fill vals
  filled  = griddata((x1, y1), newarr.ravel(),(xf, yf), method='nearest')
  array[array.mask] = filled
 return np.ma.filled(array)'''

def Create_Clusters_And_Connections(workspace,wbd,output,input_dir,nhru,info,hydroblocks_info):
 
 #Retrieve some metadata
 metadata = gdal_tools.retrieve_metadata(wbd['files']['mask'])
 #resx = metadata['resx'] #NEED TO REDEFINE FROM DEM
 #eares = 25.0#30 #meters NEED TO REDEFINE FROM DEM
 resx = 26.0#670.0**0.5#26.0

 print("Creating and curating the covariates",flush=True)
 (covariates,mask) = Create_and_Curate_Covariates(wbd,hydroblocks_info)
 
 #Determine the HRUs (clustering if semidistributed; grid cell if fully distributed)
 print("Computing the HRUs",flush=True)
 (cluster_ids,nhru,new_hand,HMC_info,covariates,dbc,hand,basins) = Compute_HRUs_Semidistributed_HMC(covariates,mask,hydroblocks_info,wbd,resx)
 covariates['hand'] = new_hand
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
 print("Calculating the connections between HRUs",flush=True)
 cmatrix = Calculate_HRU_Connections_Matrix_HMC(covariates,cluster_ids,nhru,resx,HMC_info,hydroblocks_info)

 #Define the metadata
 metadata = gdal_tools.retrieve_metadata(wbd['files']['dem'])

 #Make the output dictionary for the basin
 OUTPUT = {'hru':{},'metadata':metadata,'mask':mask,'cmatrix':cmatrix}

 #Remember the map of hrus
 OUTPUT['hru_map'] = cluster_ids
 OUTPUT['channel_map'] = HMC_info['channel_map']
 OUTPUT['hand_map'] = hand
 OUTPUT['basin_map'] = basins

 #Assign the model parameters
 print("Assigning the model parameters",flush=True)
 OUTPUT = Assign_Parameters_Semidistributed(covariates,metadata,hydroblocks_info,OUTPUT,cluster_ids,mask)

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
  nc_step = int(fp.variables['t'].units.split(' ')[0].split('h')[0])
  nc_idate = np.array(fp.variables['t'].units.split(' ')[2].split('-'))
  nc_nt = len(fp.variables['t'][:])
  dates = [datetime.datetime(int(nc_idate[0]),int(nc_idate[1]),int(nc_idate[2]))]
  for it in range(1,nc_nt): dates.append(dates[0] + datetime.timedelta(hours=it*nc_step))
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
  
  wuse_lc_ea_file = '%s/%s_lc_ea.tif' % (workspace,data_var)
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

def Read_Metadata_File(file):

 import json
 metadata = json.load(open(file))

 return metadata

