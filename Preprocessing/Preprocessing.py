import warnings
warnings.filterwarnings('ignore')
import sys
sys.path.append('Tools')
import cPickle as pickle
import datetime
import gdal_tools
import grads_tools
import numpy as np
import scipy.sparse as sparse
import scipy.stats as stats
import model_tools as mt
import os
import netCDF4 as nc
import time
import glob

def Prepare_Model_Input_Data(hydrobloks_info):

 #Prepare the info dictionary
 info = {}

 #Define the start/end dates
 info['time_info'] = {}
 info['time_info']['startdate'] = hydrobloks_info['idate']
 info['time_info']['enddate'] = hydrobloks_info['fdate']

 #Define the workspace
 workspace = hydrobloks_info['workspace']

 #Define the model input data directory
 input_dir = workspace#'%s/input' % workspace

 #Read in the metadata
 file = '%s/workspace_info.pck' % workspace
 wbd = pickle.load(open(file))

 #Create the dictionary to hold all of the data
 output = {}

 #Create the Latin Hypercube (Clustering)
 nclusters = hydrobloks_info['nclusters']
 ncores = hydrobloks_info['ncores']
 icatch = hydrobloks_info['icatch']

 #Prepare the input file
 wbd['files'] = {
  'WLTSMC':'%s/WLTSMC.tif' % workspace,
  'TEXTURE_CLASS':'%s/TEXTURE_CLASS.tif' % workspace,
  'cslope':'%s/cslope.tif' % workspace,
  'MAXSMC':'%s/MAXSMC.tif' % workspace,
  'BB':'%s/BB.tif' % workspace,
  'DRYSMC':'%s/DRYSMC.tif' % workspace,
  'fdir':'%s/fdir.tif' % workspace,
  'QTZ':'%s/QTZ.tif' % workspace,
  'SATDW':'%s/SATDW.tif' % workspace,
  'REFSMC':'%s/REFSMC.tif' % workspace,
  'mask':'%s/mask.tif' % workspace,
  'channels':'%s/channels.tif' % workspace,
  'SATDW':'%s/SATDW.tif' % workspace,
  'REFSMC':'%s/REFSMC.tif' % workspace,
  'SATPSI':'%s/SATPSI.tif' % workspace,
  'nlcd':'%s/nlcd.tif' % workspace,
  'carea':'%s/carea.tif' % workspace,
  'ti':'%s/ti.tif' % workspace,
  'ndvi':'%s/ndvi.tif' % workspace,
  'F11':'%s/F11.tif' % workspace,
  'SATDK':'%s/SATDK.tif' % workspace,
  'dem':'%s/dem.tif' % workspace,
  #'demns':'%s/workspace/demns.tif' % workspace,
  'strahler':'%s/strahler.tif' % workspace
  }
 wbd['files_meteorology'] = {
  'dlwrf':'%s/nldas/dlwrf/dlwrf.nc' % workspace,
  'dswrf':'%s/nldas/dswrf/dswrf.nc' % workspace,
  'tair':'%s/nldas/tair/tair.nc' % workspace,
  'prec':'%s/nldas/prec/prec.nc' % workspace,
  'pres':'%s/nldas/pres/pres.nc' % workspace,
  'wind':'%s/nldas/wind/wind.nc' % workspace,
  'rh':'%s/nldas/rh/rh.nc' % workspace,
  'apcpsfc':'%s/stageiv/apcpsfc/apcpsfc.nc' % workspace,
  }

 #Create the clusters and their connections
 output = Create_Clusters_And_Connections(workspace,wbd,output,input_dir,nclusters,ncores,info,hydrobloks_info)

 #Extract the meteorological forcing
 print "Preparing the meteorology"
 if hydrobloks_info['model_type'] == 'semi':
  Prepare_Meteorology_Semidistributed(workspace,wbd,output,input_dir,info,hydrobloks_info)
 elif hydrobloks_info['model_type'] == 'full':
  Prepare_Meteorology_Fulldistributed(workspace,wbd,output,input_dir,info,hydrobloks_info)

 #Write out the files to the netcdf file
 fp = hydrobloks_info['input_fp']
 data = output

 #Write out the metadata
 grp = fp.createGroup('metadata')
 grp.latitude = (wbd['bbox']['minlat'] + wbd['bbox']['maxlat'])/2
 grp.longitude = (360.0+(wbd['bbox']['minlon'] + wbd['bbox']['maxlon'])/2)

 #Write the HRU mapping
 #CONUS conus_albers metadata
 metadata = gdal_tools.retrieve_metadata(wbd['files']['mask']) 
 metadata['nodata'] = -9999.0
 #Save the conus_albers metadata
 grp = fp.createGroup('conus_albers_mapping')
 grp.createDimension('nx',metadata['nx'])
 grp.createDimension('ny',metadata['ny'])
 hmca = grp.createVariable('hmca','f4',('ny','nx')) 
 hmca.gt = metadata['gt']
 hmca.projection = metadata['projection']
 hmca.description = 'HSU mapping (conus albers)'
 hmca.nodata = metadata['nodata']
 #Save the conus albers mapping
 hsu_map = np.copy(output['hsu_map'])
 hsu_map[np.isnan(hsu_map) == 1] = metadata['nodata']
 hmca[:] = hsu_map

 if hydrobloks_info['create_mask_flag'] == True:

  #Write out the mapping
  file_ca = '%s/hsu_mapping_conus_albers.tif' % workspace
  gdal_tools.write_raster(file_ca,metadata,hsu_map)

  #Map the mapping to regular lat/lon
  file_ll = '%s/hsu_mapping_latlon.tif' % workspace
  os.system('rm -f %s' % file_ll)
  res = wbd['bbox']['res']
  minlat = wbd['bbox']['minlat']
  minlon = wbd['bbox']['minlon']
  maxlat = wbd['bbox']['maxlat']
  maxlon = wbd['bbox']['maxlon']
  log = '%s/log.txt' % workspace
  os.system('gdalwarp -tr %.16f %.16f -dstnodata %.16f -t_srs EPSG:4326 -s_srs EPSG:102039 -te %.16f %.16f %.16f %.16f %s %s >> %s 2>&1' % (res,res,metadata['nodata'],minlon,minlat,maxlon,maxlat,file_ca,file_ll,log))

 
  #Write a map for the catchment id
  file_icatch = '%s/icatch_latlon.tif' % workspace
  metadata = gdal_tools.retrieve_metadata(file_ll)
  metadata['nodata'] = -9999.0
  tmp = gdal_tools.read_raster(file_ll)
  tmp[tmp >= 0] = hydrobloks_info['icatch']
  gdal_tools.write_raster(file_icatch,metadata,tmp)

  #Add the lat/lon mapping
  #Retrieve the lat/lon metadata
  metadata = gdal_tools.retrieve_metadata(file_ll)
  metadata['nodata'] = -9999.0
  #Save the lat/lon metadata
  grp = fp.createGroup('latlon_mapping')
  grp.createDimension('nlon',metadata['nx'])
  grp.createDimension('nlat',metadata['ny'])
  hmll = grp.createVariable('hmll','f4',('nlat','nlon'))
  hmll.gt = metadata['gt']
  hmll.projection = metadata['projection']
  hmll.description = 'HSU mapping (regular lat/lon)'
  hmll.nodata = metadata['nodata']
  #Save the lat/lon mapping
  hsu_map = np.copy(gdal_tools.read_raster(file_ll))
  hsu_map[np.isnan(hsu_map) == 1] = metadata['nodata']
  hmll[:] = hsu_map

 #Write the flow matrix
 flow_matrix = output['flow_matrix']
 nconnections = flow_matrix.data.size
 grp = fp.createGroup('flow_matrix')
 grp.createDimension('connections_columns',flow_matrix.indices.size)
 grp.createDimension('connections_rows',flow_matrix.indptr.size)
 grp.createVariable('data','f4',('connections_columns',))
 grp.createVariable('indices','f4',('connections_columns',))
 grp.createVariable('indptr','f4',('connections_rows',))
 grp.variables['data'][:] = flow_matrix.data
 grp.variables['indices'][:] = flow_matrix.indices
 grp.variables['indptr'][:] = flow_matrix.indptr

 #Write the outlet information
 outlet = output['outlet']
 grp = fp.createGroup('outlet')
 full = grp.createGroup('full')
 full.createDimension('cell',outlet['full']['hru_org'].size)
 full.createVariable('i','i4',('cell',))
 full.createVariable('j','i4',('cell',))
 full.createVariable('hru_org','i4',('cell',))
 full.createVariable('hru_dst','i4',('cell',))
 full.createVariable('d8','i4',('cell',))
 full.variables['i'][:] = outlet['full']['i']
 full.variables['j'][:] = outlet['full']['j']
 full.variables['hru_org'][:] = outlet['full']['hru_org']
 full.variables['hru_dst'][:] = outlet['full']['hru_dst']
 full.variables['d8'][:] = outlet['full']['d8']
 summary = grp.createGroup('summary')
 summary.createDimension('hru',outlet['summary']['hru_org'].size)
 summary.createVariable('hru_org','i4',('hru',))
 summary.createVariable('hru_dst','i4',('hru',))
 summary.createVariable('counts','i4',('hru',))
 summary.variables['hru_org'][:] = outlet['summary']['hru_org']
 summary.variables['hru_dst'][:] = outlet['summary']['hru_dst']
 summary.variables['counts'][:] = outlet['summary']['counts']
 #outlet = {'full':{'i':outlet_icoord,'j':outlet_jcoord,'hru_org':outlet_hru_org,'hru_dst':outlet_hru_dst,'d8':outlet_d8},
 #          'summary':{'hru_org':outlet_hru_org_summary,'hru_dst':outlet_hru_dst_summary,'counts':counts}}

 #Write the model parameters
 grp = fp.createGroup('parameters')
 vars = ['slope','area_pct','land_cover','channel',
        'dem','soil_texture_class','ti','carea','area',
        'WLTSMC','MAXSMC','DRYSMC','REFSMC','SATDK',
        'mannings','m','psoil','pksat','sdmax']
 for var in vars:
  grp.createVariable(var,'f4',('hsu',))
  grp.variables[var][:] = data['hsu'][var]

 #Write other metadata
 #grp = fp.createGroup('metadata')
 #grp.outlet_hsu = data['outlet']['hsu']

 #Remove info from output
 del output['hsu']

 #Add in the catchment info
 output['wbd'] = wbd

 #Close the file
 fp.close()

 return output

def Compute_HRUs_Fulldistributed(covariates,mask,nclusters):

 cluster_ids = np.empty(covariates['ti'].shape)
 cluster_ids[:] = -9999
 cluster_ids[mask == True] = np.arange(nclusters)

 return (cluster_ids,)

def Compute_HRUs_Semidistributed(covariates,mask,nclusters,hydrobloks_info):

 #Define the number of requested hrus (channel and non-channel)
 nclusters_c = hydrobloks_info['nclusters_c']
 nclusters_nc = hydrobloks_info['nclusters_nc']
 
 #Find the mask for and without channels
 mask_c = (covariates['carea'] >= 30000.0) & (mask == True)
 mask_nc = (covariates['carea'] < 30000.0) & (mask == True)
 #mask_c = (covariates['channels'] >= 1) | (covariates['channels'] == -1)
 #mask_nc = (covariates['channels'] == 0)

 #Define the covariates
 info = {#'area':{'data':covariates['carea'][mask_nc == True],},
        #'slope':{'data':covariates['cslope'][mask == True],},
        #'sms':{'data':covariates['MAXSMC'][mask == True],},
        #'smw':{'data':covariates['WLTSMC'][mask == True],},
        'clay':{'data':covariates['clay'][mask_nc == True],},
        'sand':{'data':covariates['sand'][mask_nc == True],},
        #'ndvi':{'data':covariates['ndvi'][mask_nc ==True],},
        #'nlcd':{'data':covariates['nlcd'][mask_woc ==True],},
        'ti':{'data':covariates['ti'][mask_nc == True],},
        #'dem':{'data':covariates['dem'][mask == True],},
        #'demns':{'data':covariates['dem'][mask == True],},
        'strahler':{'data':covariates['strahler'][mask_nc == True],},
        'lats':{'data':covariates['lats'][mask_nc == True],},
        'lons':{'data':covariates['lons'][mask_nc == True],},
        }

 #Define the covariates for the channels
 info_channels = {
		 'area':{'data':np.log(covariates['carea'][mask_c == True]),},
                 #'strahler':{'data':covariates['strahler'][mask_c == True],},
                 #'lats':{'data':covariates['lats'][mask_c == True],},
                 #'lons':{'data':covariates['lons'][mask_c == True],},
                 }

 #Scale all the variables (Calculate the percentiles
 for var in info:
  #if var in ['strahler','area','ndvi','sms']:
  tmp = info[var]['data']
  tmp = (tmp - np.nanmin(tmp))/(np.nanmax(tmp) - np.nanmin(tmp))
  info[var]['data'] = tmp
  #else:
  # argsort = np.argsort(info[var]['data'])
  # pcts = np.copy(info[var]['data'])
  # pcts[argsort] = np.linspace(0,1,len(info[var]['data']))
  # info[var]['data'] = pcts
 for var in info_channels:
  tmp = info_channels[var]['data']
  tmp = (tmp - np.nanmin(tmp))/(np.nanmax(tmp) - np.nanmin(tmp))
  info_channels[var]['data'] = tmp

 import sklearn.cluster
 #Cluster the non channels regions
 bins,data = [],[]
 X = []
 for id in info:
  #Set all nans to the mean
  info[id]['data'][np.isnan(info[id]['data']) == 1] = np.nanmean(info[id]['data'])
  X.append(info[id]['data'])

 time0 = time.time()
 X = np.array(X).T
 #Subsample the array
 np.random.seed(1)
 minsamples = 10**5
 if X.shape[0] > minsamples:
  Xf = X[np.random.choice(np.arange(X.shape[0]),minsamples),:]
 else:
  Xf = X
 #Initialize all points at the 0.5 point
 init = 0.5*np.ones((nclusters_nc,Xf.shape[1]))
 batch_size = 25*nclusters_nc
 init_size = 3*batch_size
 clf = sklearn.cluster.MiniBatchKMeans(nclusters_nc,random_state=1,init=init,batch_size=batch_size,init_size=init_size)
 clf.fit(Xf)#
 clf_output = clf.predict(X)
 #Reassign the ids
 clf_output_copy = np.copy(clf_output)
 for cid in xrange(len(np.unique(clf_output))):
  clf_output[clf_output_copy == np.unique(clf_output)[cid]] = cid
 cluster_ids = np.empty(covariates['ti'].shape)
 cluster_ids[:] = -9999
 cluster_ids[mask_nc == True] = clf_output
 nclusters_old = nclusters_nc
 nclusters_nc = np.unique(clf_output).size
 #Redefine the number of clusters (We are sampling regions that just don't have data...)
 print 'clustering (non-channels) %d->%d' % (nclusters_old,nclusters_nc),time.time() - time0

 #Cluster the channels
 bins,data = [],[]
 X = []
 for id in info_channels:
  #Set all nans to the mean
  info_channels[id]['data'][np.isnan(info_channels[id]['data']) == 1] = np.nanmean(info_channels[id]['data'])
  X.append(info_channels[id]['data'])

 time0 = time.time()
 X = np.array(X).T
 #Subsample the array
 np.random.seed(1)
 minsamples = 10**5
 if X.shape[0] > minsamples:
  Xf = X[np.random.choice(np.arange(X.shape[0]),minsamples),:]
 else:
  Xf = X
 #Initialize all points at the 0.5 point
 init = 0.5*np.ones((nclusters_c,Xf.shape[1]))
 batch_size = 25*nclusters_c
 init_size = 3*batch_size
 clf = sklearn.cluster.MiniBatchKMeans(nclusters_c,random_state=1,init=init,batch_size=batch_size,init_size=init_size)
 clf.fit(Xf)#
 clf_output = clf.predict(X)
 #Reassign the ids
 clf_output_copy = np.copy(clf_output)
 for cid in xrange(len(np.unique(clf_output))):
  clf_output[clf_output_copy == np.unique(clf_output)[cid]] = cid
 #cluster_ids = np.empty(covariates['ti'].shape)
 #cluster_ids[:] = -9999
 cluster_ids[mask_c == True] = clf_output + np.nanmax(cluster_ids) + 1
 nclusters_old = nclusters_c
 nclusters_c = np.unique(clf_output).size
 #Redefine the number of clusters (We are sampling regions that just don't have data...)
 print 'clustering (channels) %d->%d' % (nclusters_old,nclusters_c),time.time() - time0
 #Add in the channel clusters
 #channels = np.unique(covariates['channels'][covariates['channels'] > 0])
 #for channel in channels:
 # cid = int(np.nanmax(cluster_ids) + 1)
 # idx = np.where(covariates['channels'] == channel)
 # cluster_ids[idx] = cid
 # nclusters = nclusters + 1
 nclusters = nclusters_nc + nclusters_c

 #Reorder according to areas
 areas = []
 for cid in xrange(nclusters):
  msk = cluster_ids == cid
  #areas.append(np.nanmean(covariates['carea'][msk]))
  areas.append(np.nanmax(covariates['carea'][msk]))
 argsort = np.argsort(np.array(areas))
 cluster_ids_new = np.copy(cluster_ids)
 for cid in xrange(nclusters):
  msk = cluster_ids == argsort[cid]
  cluster_ids_new[msk] = cid
 cluster_ids = np.copy(cluster_ids_new)
 areas = []
 for cid in xrange(nclusters):
  msk = cluster_ids == cid
  areas.append(np.nanmean(covariates['carea'][msk]))

 return (cluster_ids,nclusters)

def Assign_Parameters_Fulldistributed(covariates,metadata,hydrobloks_info,OUTPUT,cluster_ids,mask):

 nclusters = hydrobloks_info['nclusters']
 #Initialize the arrays
 vars = ['area','area_pct','BB','DRYSMC','F11','MAXSMC','REFSMC','SATPSI',
         'SATDK','SATDW','WLTSMC','QTZ','slope','ti','dem','carea','channel',
         'land_cover','soil_texture_class',
         'mannings','m','psoil','pksat','sdmax']
 #vars = ['area','area_pct','BB','DRYSMC','F11','MAXSMC','REFSMC','SATPSI',
 #        'SATDK','SATDW','WLTSMC','QTZ','slope','ti','dem','carea','channel',         'vchan','vof','land_cover','soil_texture_class']
 OUTPUT['hsu'] = {}
 for var in vars:
  OUTPUT['hsu'][var] = np.zeros(nclusters) 

 #Metadata
 NLCD2NOAH = {11:17,12:15,21:10,22:10,23:10,24:13,31:16,41:4,42:1,43:5,51:6,52:6,71:10,72:10,73:19,74:19,81:10,82:12,90:11,95:11}

 #Calculate area per hsu
 OUTPUT['hsu']['area'][:] = metadata['resx']**2
 #Calculate area percentage per hsu
 OUTPUT['hsu']['area_pct'][:] = 100*OUTPUT['hsu']['area']/(metadata['resx']**2*nclusters)
 #Soil properties
 for var in ['BB','DRYSMC','F11','MAXSMC','REFSMC','SATPSI','SATDK','SATDW','WLTSMC','QTZ']:
  if var in ['SATDK','SATDW']:
   OUTPUT['hsu'][var][:] = covariates[var][mask]
  else:
   OUTPUT['hsu'][var][:] = covariates[var][mask]
 #Average Slope
 OUTPUT['hsu']['slope'][:] = covariates['cslope'][mask]
 #Topographic index
 OUTPUT['hsu']['ti'][:] = covariates['ti'][mask]
 #DEM
 OUTPUT['hsu']['dem'][:] = covariates['dem'][mask]
 #Average Catchment Area
 OUTPUT['hsu']['carea'][:] = covariates['carea'][mask]
 #Channel?
 OUTPUT['hsu']['channel'][:] = covariates['channels'][mask]
 #Land cover type  
 land_cover = np.copy(covariates['nlcd'])
 for lc in np.unique(land_cover)[1:]:
  land_cover[land_cover == lc] = NLCD2NOAH[lc]
 OUTPUT['hsu']['land_cover'][:] = land_cover[mask]
 #Soil texture class
 OUTPUT['hsu']['soil_texture_class'][:] = covariates['TEXTURE_CLASS'][mask]
 #Define the estimate for the model parameters
 OUTPUT['hsu']['m'][:] = 0.1 #Form of the exponential decline in conductivity (0.01-1.0)
 OUTPUT['hsu']['pksat'][:] = 1.0 #saturated hydraulic conductivity scalar multiplier (0.1-1.0)
 OUTPUT['hsu']['psoil'][:] = 1.0 #soil hydraulic properties (residual,wilting,field capacity, and porosity) (0.1-10.0)
 OUTPUT['hsu']['sdmax'][:] = 5.0 #maximum effective deficit of subsurface saturated zone (0.1-10.0)
 OUTPUT['hsu']['mannings'][covariates['carea'][mask] >= 100000.0] = 0.03 #manning's n for channel flow (0.01-0.1)
 OUTPUT['hsu']['mannings'][covariates['carea'][mask] < 100000.0] = 0.15 #manning's n for overland flow (0.01-0.8)

 return OUTPUT

def Assign_Parameters_Semidistributed(covariates,metadata,hydrobloks_info,OUTPUT,cluster_ids,mask):

 nclusters = hydrobloks_info['nclusters']
 #Initialize the arrays
 vars = ['area','area_pct','BB','DRYSMC','F11','MAXSMC','REFSMC','SATPSI',
         'SATDK','SATDW','WLTSMC','QTZ','slope','ti','dem','carea','channel',
         'land_cover','soil_texture_class',
         'mannings','m','psoil','pksat','sdmax']
 OUTPUT['hsu'] = {}
 for var in vars:
  OUTPUT['hsu'][var] = np.zeros(nclusters)

 #Metadata
 NLCD2NOAH = {11:17,12:15,21:10,22:10,23:10,24:13,31:16,41:4,42:1,43:5,51:6,52:6,71:10,72:10,73:19,74:19,81:10,82:12,90:11,95:11}
 for hsu in np.arange(nclusters):

  #Set indices
  idx = np.where(cluster_ids == hsu)
  #Calculate area per hsu
  OUTPUT['hsu']['area'][hsu] = metadata['resx']**2*idx[0].size
  #Calculate area percentage per hsu
  OUTPUT['hsu']['area_pct'][hsu] = 100*OUTPUT['hsu']['area'][hsu]/(metadata['resx']**2*mask[mask].size)
  #Soil properties
  for var in ['BB','DRYSMC','F11','MAXSMC','REFSMC','SATPSI','SATDK','SATDW','WLTSMC','QTZ']:
   if var in ['SATDK','SATDW']:
    OUTPUT['hsu'][var][hsu] = stats.mstats.hmean(covariates[var][idx])
   else:
    OUTPUT['hsu'][var][hsu] = stats.mstats.gmean(covariates[var][idx])
  #Average Slope
  OUTPUT['hsu']['slope'][hsu] = np.mean(covariates['cslope'][idx])
  #Topographic index
  OUTPUT['hsu']['ti'][hsu] = np.mean(covariates['ti'][idx])
  #DEM
  OUTPUT['hsu']['dem'][hsu] = np.mean(covariates['dem'][idx])
  #Average Catchment Area
  OUTPUT['hsu']['carea'][hsu] = np.mean(covariates['carea'][idx])
  #Channel?
  OUTPUT['hsu']['channel'][hsu] = stats.mode(covariates['channels'][idx])[0]
  #Land cover type  
  OUTPUT['hsu']['land_cover'][hsu] = NLCD2NOAH[stats.mode(covariates['nlcd'][idx])[0][0]]
  #Soil texture class
  OUTPUT['hsu']['soil_texture_class'][hsu] = stats.mode(covariates['TEXTURE_CLASS'][idx])[0][0]
  #Define the estimate for the model parameters
  OUTPUT['hsu']['m'][hsu] = 0.1 #Form of the exponential decline in conductivity (0.01-1.0)
  OUTPUT['hsu']['pksat'][hsu] = 1.0 #saturated hydraulic conductivity scalar multiplier (0.1-1.0)
  OUTPUT['hsu']['psoil'][hsu] = 1.0 #soil hydraulic properties (residual,wilting,field capacity, and porosity) (0.1-10.0)
  OUTPUT['hsu']['sdmax'][hsu] = 5.0 #maximum effective deficit of subsurface saturated zone (0.1-10.0)
  if np.max(covariates['carea'][idx]) >= 100000.0: OUTPUT['hsu']['mannings'][hsu] = 0.03 #manning's n for channel flow (0.01-0.1)
  else: OUTPUT['hsu']['mannings'][hsu] = 0.15 #manning's n for overland flow (0.01-0.8)


 return OUTPUT

def Calculate_Flow_Matrix(covariates,cluster_ids,nclusters):

 #Prepare the flow matrix
 mask1 = covariates['fdir'] < 0
 covariates['fdir'][mask1] = -9999.0
 cluster_ids_copy = np.copy(cluster_ids)
 cluster_ids_copy[mask1] = np.nan
 max_nhru = np.sum(cluster_ids >= 0)
 #tp_matrix = mt.preprocessor.calculate_connections_d8(cluster_ids_copy,covariates['fdir'],nclusters,max_nhru)
 (hrus_dst,hrus_org,outlet_icoord,outlet_jcoord,outlet_hru,outlet_d8) = mt.preprocessor.calculate_connections_d8(cluster_ids_copy,covariates['fdir'],covariates['carea'],nclusters,max_nhru)
 #Only use the non -9999 values
 hrus_dst = hrus_dst[hrus_dst != -9999]-1
 hrus_org = hrus_org[hrus_org != -9999]-1
 outlet_icoord = outlet_icoord[outlet_icoord != -9999]-1
 outlet_jcoord = outlet_jcoord[outlet_jcoord != -9999]-1
 outlet_hru_org = outlet_hru[outlet_hru != -9999] - 1
 outlet_d8 = outlet_d8[outlet_d8 != -9999]

 #Create hrus for the outlets
 outlet_hru_dst_summary = np.arange(nclusters,nclusters+np.unique(outlet_hru_org).size)
 outlet_hru_org_summary = np.unique(outlet_hru_org)
 outlet_hru_dst = np.zeros(outlet_hru_org.size)
 counts = []
 for (hru_org,hru_dst) in zip(outlet_hru_org_summary,outlet_hru_dst_summary):
  idx = outlet_hru_org == hru_org
  counts.append(np.sum(idx))
  outlet_hru_dst[idx] = hru_dst
 counts = np.array(counts)
 
 #Create a dictionary of outlet information
 outlet = {'full':{'i':outlet_icoord,'j':outlet_jcoord,'hru_org':outlet_hru_org,'hru_dst':outlet_hru_dst,'d8':outlet_d8},
           'summary':{'hru_org':outlet_hru_org_summary,'hru_dst':outlet_hru_dst_summary,'counts':counts}}

 #Update the input to create the sparse matrix
 hrus_dst = np.append(hrus_dst,outlet_hru_dst)
 hrus_org = np.append(hrus_org,outlet_hru_org)
 hrus_dst[hrus_dst == -1] = outlet_hru_org[0] #CAREFUL. Designed to push all the extra small outlets to the same exit

 #Prepare the sparse matrix
 flow_matrix = sparse.coo_matrix((np.ones(hrus_dst.size),(hrus_org,hrus_dst)),dtype=np.float32)
 flow_matrix = flow_matrix.tocsr()
 #Normalize the rows (sum to 1)
 fm_sum = flow_matrix.sum(axis=1)
 fm_data = flow_matrix.data
 fm_indptr = flow_matrix.indptr
 for row in xrange(fm_sum.size):
  fm_data[fm_indptr[row]:fm_indptr[row+1]] = fm_data[fm_indptr[row]:fm_indptr[row+1]]/fm_sum[row]
 flow_matrix.data = fm_data

 return (flow_matrix.T,outlet)

def Create_Soils_File(hydrobloks_info,OUTPUT,input_dir):

 #Read in table of NOAH soil parameter values
 dir = os.path.dirname(os.path.abspath(__file__))
 fp = open('%s/../HydroBloks/pyNoahMP/data/SOILPARM.TBL' % dir)
 iline = 0
 soils_data = {'MAXSMC':[],'DRYSMC':[],'REFSMC':[],'ID':[],'BB':[],'F11':[],'SATPSI':[],'SATDK':[]
,'SATDW':[],'WLTSMC':[],'QTZ':[]}
 for line in fp:
  if (iline > 2) & (iline < 15):
   tmp = line.split(',')
   soils_data['ID'].append(float(tmp[0]))
   soils_data['BB'].append(float(tmp[1]))
   soils_data['DRYSMC'].append(float(tmp[2]))
   soils_data['F11'].append(float(tmp[3]))
   soils_data['MAXSMC'].append(float(tmp[4]))
   soils_data['REFSMC'].append(float(tmp[5]))
   soils_data['SATPSI'].append(float(tmp[6]))
   soils_data['SATDK'].append(float(tmp[7]))
   soils_data['SATDW'].append(float(tmp[8]))
   soils_data['WLTSMC'].append(float(tmp[9]))
   soils_data['QTZ'].append(float(tmp[10]))
  iline = iline + 1
 fp.close()

 #Soil properties
 soil_vars = ['BB','DRYSMC','F11','MAXSMC','REFSMC','SATPSI','SATDK','SATDW','WLTSMC','QTZ']
 nhsus = hydrobloks_info['nclusters']
 soilsdir = '%s/soils' % input_dir
 os.system('mkdir -p %s' % soilsdir)
 soils_lookup = hydrobloks_info['soil_file']
 fp = open(soils_lookup,'w')
 fp.write('Soil Parameters\n')
 fp.write('CUST\n')
 fp.write("%d,1   'BB      DRYSMC      F11     MAXSMC   REFSMC   SATPSI  SATDK       SATDW     WLTSMC  QTZ    '\n" % nhsus)
 for hsu in np.arange(nhsus):
  fp.write('%d ' % (hsu+1))
  for var in soil_vars:
   if var in ['DRYSMC','MAXSMC','REFSMC','WLTSMC','SATDK']:
    fp.write(',%.10f ' % OUTPUT['hsu'][var][hsu])
   else:
    idx = soils_data['ID'].index(OUTPUT['hsu']['soil_texture_class'][hsu])
    fp.write(',%.10f ' % soils_data[var][idx])
  fp.write('\n')
 fp.close()

 return soils_lookup

def Create_and_Curate_Covariates(wbd):

 covariates = {}
 #Read in and curate all the covariates
 root = wbd['files']['MAXSMC'][0:-11]
 wbd['files']['sand'] = '%s/dssurgo/sandtotal_r.tif' % root
 wbd['files']['clay'] = '%s/dssurgo/claytotal_r.tif' % root
 for file in wbd['files']:
  covariates[file] = gdal_tools.read_raster(wbd['files'][file])
  if file == 'cslope':
   mask = covariates[file] == 0.0
   covariates[file][mask] = 0.000001

 #Create lat/lon grids
 lats = np.linspace(wbd['bbox']['minlat']+wbd['bbox']['res']/2,wbd['bbox']['maxlat']-wbd['bbox']['res']/2,covariates['ti'].shape[0])
 lons = np.linspace(wbd['bbox']['minlon']+wbd['bbox']['res']/2,wbd['bbox']['maxlon']-wbd['bbox']['res']/2,covariates['ti'].shape[1])
 lats, lons = np.meshgrid(lats, lons)
 covariates['lats'] = lats.T
 covariates['lons'] = lons.T

 #Define the mask
 mask = np.copy(covariates['mask'])
 mask[mask > 0] = 1
 mask[mask < 0] = 0
 mask = mask.astype(np.bool)

 #Set all nans to the mean
 for var in covariates:
  mask1 = (np.isinf(covariates[var]) == 0) & (np.isnan(covariates[var]) == 0)
  mask0 = (np.isinf(covariates[var]) == 1) | (np.isnan(covariates[var]) == 1)
  if var in ['fdir','nlcd']:
   covariates[var][mask0] = stats.mode(covariates[var][mask1])[0][0]
  else:
   covariates[var][mask0] = np.mean(covariates[var][mask1])

 #Set everything that is -9999 to the mean
 for var in covariates:
  if var in ['fdir','nlcd','TEXTURE_CLASS']:
   covariates[var][covariates[var] == -9999.0] = stats.mode(covariates[var][covariates[var] != -9999.0])[0][0]
  else:
   covariates[var][covariates[var] == -9999.0] = np.mean(covariates[var][covariates[var] != -9999.0])

 #Set everything outside of the mask to -9999
 for var in covariates:
  covariates[var][mask <= 0] = -9999.0
 
 return (covariates,mask)

def Create_Clusters_And_Connections(workspace,wbd,output,input_dir,nclusters,ncores,info,hydrobloks_info):
 
 print "Creating and curating the covariates"
 (covariates,mask) = Create_and_Curate_Covariates(wbd)

 #Determine the HRUs (clustering if semidistributed; grid cell if fully distributed)
 print "Computing the HRUs"
 if hydrobloks_info['model_type'] == 'semi':
  (cluster_ids,nclusters) = Compute_HRUs_Semidistributed(covariates,mask,nclusters,hydrobloks_info)
  hydrobloks_info['nclusters'] = nclusters
 elif hydrobloks_info['model_type'] == 'full':
  nclusters = np.sum(mask == True)
  hydrobloks_info['nclusters'] = nclusters
  (cluster_ids,) = Compute_HRUs_Fulldistributed(covariates,mask,nclusters)

 #Create the netcdf file
 file_netcdf = hydrobloks_info['input_file']
 hydrobloks_info['input_fp'] = nc.Dataset(file_netcdf, 'w', format='NETCDF4')

 #Create the dimensions (netcdf)
 idate = hydrobloks_info['idate']
 fdate = hydrobloks_info['fdate']
 dt = hydrobloks_info['dt']
 ntime = 24*3600*((fdate - idate).days+1)/dt
 nhsu = hydrobloks_info['nclusters']
 hydrobloks_info['input_fp'].createDimension('hsu',nhsu)
 hydrobloks_info['input_fp'].createDimension('time',ntime)

 #Create the groups (netcdf)
 hydrobloks_info['input_fp'].createGroup('meteorology')

 #Prepare the flow matrix
 print "Calculating the flow matrix"
 (flow_matrix,outlet) = Calculate_Flow_Matrix(covariates,cluster_ids,nclusters)

 #Define the metadata
 metadata = gdal_tools.retrieve_metadata(wbd['files']['ti'])

 #Make the output dictionary for the basin
 OUTPUT = {'hsu':{},'metadata':metadata,'mask':mask,'flow_matrix':flow_matrix}
 OUTPUT['outlet'] = outlet

 #Determine outlet cell
 #covariates['carea'][mask == False] = np.nan
 #outlet_idx = np.where(covariates['carea'] == np.max(covariates['carea'][np.isnan(covariates['carea']) == 0]))
 #outlet_idx = [int(outlet_idx[0]),int(outlet_idx[1])]
 #OUTPUT['outlet'] = {'idx':outlet_idx,'hsu':cluster_ids[outlet_idx[0],outlet_idx[1]]}

 #Remember the map of hrus
 OUTPUT['hsu_map'] = cluster_ids

 #Assign the model parameters
 print "Assigning the model parameters"
 if hydrobloks_info['model_type'] == 'semi':
  OUTPUT = Assign_Parameters_Semidistributed(covariates,metadata,hydrobloks_info,OUTPUT,cluster_ids,mask)
 elif hydrobloks_info['model_type'] == 'full':
  OUTPUT = Assign_Parameters_Fulldistributed(covariates,metadata,hydrobloks_info,OUTPUT,cluster_ids,mask)

 #Create the soil parameters file
 print "Creating the soil file"
 soils_lookup = Create_Soils_File(hydrobloks_info,OUTPUT,input_dir)

 #Add the new number of clusters
 OUTPUT['nclusters'] = nclusters
 OUTPUT['mask'] = mask

 return OUTPUT

def Prepare_Meteorology_Fulldistributed(workspace,wbd,OUTPUT,input_dir,info,hydrobloks_info):

 #Assign other variables
 mask = OUTPUT['mask']

 #Define the mapping directory
 mapping_dir = '%s/mapping' % workspace
 mapping_info = {}
 #Calculate the fine to coarse scale mapping
 for data_var in wbd['files_meteorology']:

  #Define the variable name
  var = data_var#data_var.split('_')[1]
  mapping_info[var] = {}

  #Read in the coarse and fine mapping
  file_coarse = '%s/%s_coarse.tif' % (mapping_dir,data_var)
  file_fine = '%s/%s_fine.tif' % (mapping_dir,data_var)
  mask_coarse = gdal_tools.read_raster(file_coarse)
  mask_fine = gdal_tools.read_raster(file_fine)
  nlat = mask_coarse.shape[0]
  nlon = mask_coarse.shape[1]

  print data_var,"Creating the i and j mapping"
  #Create maps of i and j for mapping
  coords = {'i':np.zeros(mask_fine.shape,dtype=np.int),
            'j':np.zeros(mask_fine.shape,dtype=np.int)}
  for value in np.unique(mask_coarse):
   if value < 0:continue
   idx_fine = mask_fine == value
   idx_coarse = np.where(mask_coarse == value)
   coords['i'][idx_fine] = nlat - idx_coarse[0] - 1 #Careful
   coords['j'][idx_fine] = idx_coarse[1]
  coords['i'] = coords['i'][mask]
  coords['j'] = coords['j'][mask]

  print data_var,"Extracting all the data"
  #Extract all of the data for that variable
  idate = info['time_info']['startdate']
  fdate = info['time_info']['enddate']
  nt = 24*((fdate - idate).days+1)
  var = data_var#data_var.split('_')[1]
  date = idate
  file = wbd['files_meteorology'][data_var]
  fp = nc.Dataset(file)
  #Determine the time steps to retrieve
  dates = nc.num2date(fp.variables['t'][:],units='hours since %02d-%02d-%02d 00:00:00' % (idate.year,idate.month,idate.day))
  mask_dates = (dates >= idate) & (dates <= fdate)
  data = np.ma.getdata(fp.variables[var][mask_dates])
  fp.close()

  print data_var,"Assigning the data"
  #Create the variable
  grp = hydrobloks_info['input_fp'].groups['meteorology']
  grp.createVariable(var,'f4',('time','hsu'))
  #Place all the data per time step
  tmp = np.zeros(np.sum(mask))
  for itime in np.arange(data.shape[0]):
   tmp[:] = data[itime,coords['i'],coords['j']]
   #Write the timestep
   grp.variables[var][itime,:] = tmp[:]

 return

def Prepare_Meteorology_Semidistributed(workspace,wbd,OUTPUT,input_dir,info,hydrobloks_info):

 #Define the mapping directory
 mapping_dir = '%s/mapping' % workspace
 mapping_info = {}
 #Calculate the fine to coarse scale mapping
 for data_var in wbd['files_meteorology']:
  
  #Define the variable name
  var = data_var#data_var.split('_')[1]
  mapping_info[var] = {}

  #Read in the coarse and fine mapping
  file_coarse = '%s/%s_coarse.tif' % (mapping_dir,data_var)
  file_fine = '%s/%s_fine.tif' % (mapping_dir,data_var)
  mask_coarse = gdal_tools.read_raster(file_coarse)
  mask_fine = gdal_tools.read_raster(file_fine)
  nlat = mask_coarse.shape[0]
  nlon = mask_coarse.shape[1]

  #Compute the mapping for each hsu
  for hsu in np.arange(hydrobloks_info['nclusters']):
   idx = OUTPUT['hsu_map'] == hsu
   icells = np.unique(mask_fine[idx].astype(np.int))
   counts = np.bincount(mask_fine[idx].astype(np.int))
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
   mapping_info[var][hsu] = {'pcts':pcts,'coords':coords}

 #Iterate through variable creating forcing product per HSU
 idate = info['time_info']['startdate']
 fdate = info['time_info']['enddate']
 nt = 24*((fdate - idate).days+1)
 #Create structured array
 meteorology = {}
 for data_var in wbd['files_meteorology']:
  meteorology[data_var] = np.zeros((nt,hydrobloks_info['nclusters']))
 #Load data into structured array
 for data_var in wbd['files_meteorology']:
  var = data_var#data_var.split('_')[1]
  date = idate
  file = wbd['files_meteorology'][data_var]
  fp = nc.Dataset(file)
  #fp = h5py.File(file)
  #Determine the time steps to retrieve
  dates = nc.num2date(fp.variables['t'][:],units='hours since %02d-%02d-%02d 00:00:00' % (idate.year,idate.month,idate.day))
  #dates = nc.num2date(fp['t'][:],units='hours since %02d-%02d-%02d 00:00:00' % (idate.year,idate.month,idate.day))
  mask_dates = (dates >= idate) & (dates <= fdate)
  data = np.ma.getdata(fp.variables[var][mask_dates,:,:])
  fp.close()
  #Assing to hsus
  for hsu in mapping_info[var]:
   pcts = mapping_info[var][hsu]['pcts']
   coords = mapping_info[var][hsu]['coords']
   coords[0][coords[0] >= data.shape[1]] = data.shape[1] - 1
   coords[1][coords[1] >= data.shape[2]] = data.shape[2] - 1
   tmp = pcts*data[:,coords[0],coords[1]]
   #Combine stage iv and nldas here
   if data_var not in ['apcpsfc',]:tmp[tmp < -999] = np.mean(tmp[tmp > -999])
   meteorology[data_var][:,hsu] = np.sum(tmp,axis=1)

  #Write the meteorology to the netcdf file (single chunk for now...)
  grp = hydrobloks_info['input_fp'].groups['meteorology']
  grp.createVariable(var,'f4',('time','hsu'))
  grp.variables[data_var][:] = meteorology[data_var][:]

 return 

def Read_Metadata_File(file):

 import json
 metadata = json.load(open(file))

 return metadata
