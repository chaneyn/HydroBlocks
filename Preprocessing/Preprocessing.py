import warnings
warnings.filterwarnings('ignore')
import sys
sys.path.append('Tools')
import cPickle as pickle
import datetime
#import gdal_tools
import numpy as np
import scipy.sparse as sparse
import scipy.stats as stats
import model_tools as mt
import os
import netCDF4 as nc
import time
import glob
from geospatialtools import gdal_tools
from geospatialtools import terrain_tools
import matplotlib.pyplot as plt

def Prepare_Model_Input_Data(hydrobloks_info):

 #Prepare the info dictionary
 info = {}

 #Define the start/end dates
 info['time_info'] = {}
 info['time_info']['startdate'] = hydrobloks_info['idate']
 info['time_info']['enddate'] = hydrobloks_info['fdate']
 info['time_info']['dt'] = hydrobloks_info['dt']

 #Define the workspace
 workspace = hydrobloks_info['workspace']

 #Define the model input data directory
 input_dir = workspace#'%s/input' % workspace

 #Read in the metadata
 #file = '%s/workspace_info.pck' % workspace
 #wbd = pickle.load(open(file))

 #Create the dictionary to hold all of the data
 output = {}

 #Create the Latin Hypercube (Clustering)
 nclusters = hydrobloks_info['nclusters']
 ncores = hydrobloks_info['ncores']
 icatch = hydrobloks_info['icatch']

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
  'cslope':'%s/cslope_ea.tif' % workspace,
  'MAXSMC':'%s/thetas_ea.tif' % workspace,
  'BB':'%s/bb_ea.tif' % workspace,
  'DRYSMC':'%s/thetar_ea.tif' % workspace,
  'fdir':'%s/fdir_ea.tif' % workspace,
  'QTZ':'%s/qtz_ea.tif' % workspace,
  'SATDW':'%s/dsat_ea.tif' % workspace,
  'REFSMC':'%s/theta33_ea.tif' % workspace,
  'mask':'%s/mask_ea.tif' % workspace,
  'SATPSI':'%s/psisat_ea.tif' % workspace,
  'lc':'%s/lc_ea.tif' % workspace,
  'carea':'%s/carea_ea.tif' % workspace,
  'ti':'%s/ti_ea.tif' % workspace,
  'ndvi':'%s/ndvi_ea.tif' % workspace,
  'F11':'%s/f11_ea.tif' % workspace,
  'SATDK':'%s/ksat_ea.tif' % workspace,
  'dem':'%s/dem_ea.tif' % workspace,
  'demns':'%s/demns_ea.tif' % workspace,
  'sand':'%s/sand_ea.tif' % workspace,
  'clay':'%s/clay_ea.tif' % workspace,
  'silt':'%s/silt_ea.tif' % workspace,
  'om':'%s/om_ea.tif' % workspace,
  'bare30':'%s/bare30_ea.tif' % workspace,
  'water30':'%s/water30_ea.tif' % workspace,
  'tree30':'%s/tree30_ea.tif' % workspace,
  }
 wbd['files_meteorology'] = {
  'lwdown':'%s/lwdown.nc' % workspace,
  'swdown':'%s/swdown.nc' % workspace,
  'tair':'%s/tair.nc' % workspace,
  'precip':'%s/precip.nc' % workspace,
  'psurf':'%s/psurf.nc' % workspace,
  'wind':'%s/wind.nc' % workspace,
  'spfh':'%s/spfh.nc' % workspace,
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
 lon = (wbd['bbox']['minlon'] + wbd['bbox']['maxlon'])/2 
 if lon < 0:lon += 360
 grp.longitude = lon
 #grp.longitude = (360.0+(wbd['bbox']['minlon'] + wbd['bbox']['maxlon'])/2)

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

 #Write out the mapping
 file_ca = '%s/hsu_mapping_ea.tif' % workspace
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
 os.system('gdalwarp -tr %.16f %.16f -dstnodata %.16f -t_srs \'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs \' -te %.16f %.16f %.16f %.16f %s %s >> %s 2>&1' % (res,res,metadata['nodata'],minlon,minlat,maxlon,maxlat,file_ca,file_ll,log))

 #Write a map for the catchment id
 file_icatch = '%s/icatch_latlon.tif' % workspace
 metadata = gdal_tools.retrieve_metadata(file_ll)
 metadata['nodata'] = -9999.0
 tmp = gdal_tools.read_raster(file_ll)
 tmp[tmp >= 0] = hydrobloks_info['icatch']
 gdal_tools.write_raster(file_icatch,metadata,tmp)

 #Retrieve the lat/lon metadata
 #metadata = gdal_tools.retrieve_metadata(file_ll)
 #metadata['nodata'] = -9999.0
 #Save the lat/lon metadata
 #grp = fp.createGroup('latlon_mapping')
 #grp.createDimension('nlon',metadata['nx'])
 #grp.createDimension('nlat',metadata['ny'])
 #hmll = grp.createVariable('hmll','f4',('nlat','nlon'))
 #hmll.gt = metadata['gt']
 #hmll.projection = metadata['projection']
 #hmll.description = 'HSU mapping (regular lat/lon)'
 #hmll.nodata = metadata['nodata']
 #Save the lat/lon mapping
 #hsu_map = np.copy(gdal_tools.read_raster(file_ll))
 #hsu_map[np.isnan(hsu_map) == 1] = metadata['nodata']
 #hmll[:] = hsu_map

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
        'BB','F11','SATPSI','SATDW','QTZ',
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

def Compute_HRUs_Semidistributed_Basin(covariates,mask,nclusters,hydrobloks_info):

 time0 = time.time()
 dem = covariates['demns']
 res = hydrobloks_info['dx']
 basin_info = hydrobloks_info['basin_info']
 nbasins = basin_info['nbasins'] #characteristic hillslope
 ntiles = basin_info['nclusters']
 #Calculate mfd accumulation area
 mfd_area = terrain_tools.ttf.calculate_mfd_acc(dem,res,1)
 mfd_area[mask == 0] = 0.0
 #Calculate mfd topographic index
 mfd_ti = np.log(mfd_area/res/covariates['cslope'])
 #Calculate d8 accumulation area
 (area,fdir) = terrain_tools.ttf.calculate_d8_acc(dem,res)
 area[mask == 0] = 0.0
 #Define the channels
 channels = terrain_tools.ttf.calculate_channels(area,10**6,10**6,fdir)
 #Define the basins
 basins = terrain_tools.ttf.delineate_basins(channels,mask,fdir)
 #Calculate the basin properties
 lats = np.linspace(0,1,basins.shape[0])
 lons = np.linspace(0,1,basins.shape[1])
 (latitude,longitude) = np.meshgrid(lats,lons)
 bp = terrain_tools.calculate_basin_properties(basins,res,latitude.T,longitude.T,fdir)
 #Reduce the number of basins to the desired number
 basins = terrain_tools.reduce_basin_number(basins,bp,nbasins)
 #Define the depth to channel
 depth2channel = terrain_tools.ttf.calculate_depth2channel(channels,basins,fdir,dem)
 #Compute the accumulated values
 (acc_ksat,fdir) = terrain_tools.ttf.calculate_d8_acc_alt(dem,res,covariates['SATDK'])
 (acc_slope,fdir) = terrain_tools.ttf.calculate_d8_acc_alt(dem,res,covariates['cslope'])
 #(acc_ndvi,fdir) = terrain_tools.ttf.calculate_d8_acc_alt(dem,res,covariates['ndvi'])
 #Determine how many neighboring cells converge to this one
 neighbor_count = terrain_tools.ttf.calculate_d8_acc_neighbors(dem,res,np.ones(dem.shape))
 neighbor_sti = terrain_tools.ttf.calculate_d8_acc_neighbors(dem,res,covariates['sti'])
 neighbor_sti[neighbor_count > 0] = neighbor_sti[neighbor_count > 0]/neighbor_count[neighbor_count > 0]
 covariates = {#'depth2channel':{'data':depth2channel},
              #'clay':{'data':covariates['clay']},
              #'sand':{'data':covariates['sand']},
              #'ndvi':{'data':acc_ndvi/area},
              'ndvi':{'data':covariates['ndvi']},
              #'nlcd':{'data':covariates['nlcd']},
	      'ksat_alt':{'data':np.log10(900*acc_ksat/area)},
	      'ksat':{'data':np.log10(covariates['SATDK'])},
	      #'slope':{'data':np.log10(covariates['cslope'])},
	      #'slope_alt':{'data':900*acc_slope/area},
              'neighbor_count':{'data':neighbor_count},
              #'neighbor_sti':{'data':neighbor_sti},
              #'refsmc':{'data':covariates['REFSMC']},
	      #'slope':{'data':covariates['cslope']},
	      #'smcmax':{'data':covariates['MAXSMC']},
              #'d2c':{'data':depth2channel},
              'acc':{'data':np.log10(area)},
              #'sti':{'data':covariates['sti']},
              'lats':{'data':covariates['lats']},
              'lons':{'data':covariates['lons']},
              }
 hrus = terrain_tools.create_tiles_kmeans(basins,covariates,ntiles)
 #hrus = terrain_tools.create_nd_histogram(basins,covariates)-1
 hrus = hrus.astype(np.float32) 
 hrus[hrus < 0] = -9999
 #hrus[channels > 0] = np.max(hrus) + 1
 #nhrus = int(np.max(hrus) + 1)
 nhrus = int(np.max(hrus) + 1)
 print '%d hrus' % (nhrus),time.time() - time0

 return (hrus,nhrus)

def Compute_HRUs_Semidistributed_Hillslope(covariates,mask,nclusters,hydrobloks_info):

 time0 = time.time()
 dem = covariates['demns']
 res = hydrobloks_info['dx']
 hillslope_info = hydrobloks_info['hillslope_info']
 nclusters = hillslope_info['nhillslopes'] #characteristic hillslope
 #Calculate mfd accumulation area
 mfd_area = terrain_tools.ttf.calculate_mfd_acc(dem,res,1)
 mfd_area[mask == 0] = 0.0
 #Calculate mfd topographic index
 mfd_ti = np.log(mfd_area/res/covariates['cslope'])
 #Calculate d8 accumulation area
 (area,fdir) = terrain_tools.ttf.calculate_d8_acc(dem,res)
 area[mask == 0] = 0.0
 #Define the channels
 channels = terrain_tools.ttf.calculate_channels(area,10**6,3*10**7,fdir)
 #Define the basins
 basins = terrain_tools.ttf.delineate_basins(channels,mask,fdir)
 #Define the depth to channel
 depth2channel = terrain_tools.ttf.calculate_depth2channel(channels,basins,fdir,dem)
 #Define the hillslopes
 hillslopes = terrain_tools.ttf.calculate_hillslopesd8(channels,basins,fdir)
 #Calculate the hillslope properties
 latitude = covariates['lats'].T
 longitude = covariates['lons'].T
 hp = terrain_tools.calculate_hillslope_properties(hillslopes,dem,
                       basins,res,latitude.T,longitude.T,depth2channel)
 #Cluster the hillslopes
 (hillslopes_clusters,nhillslopes) = terrain_tools.cluster_hillslopes(hp,hillslopes,nclusters)
 #Create the hrus
 '''covariates = {#'depth2channel':{'data':depth2channel,'nbins':nclusters_tiles},
              #'clay':{'data':covariates['clay'],'nbins':hillslope_info['clay']},
              'ndvi':{'data':covariates['ndvi'],'nbins':hillslope_info['ndvi']},
              #'ti':{'data':mfd_ti,'nbins':hillslope_info['ti']},
              'acc':{'data':area,'nbins':10},
              #'depth2channel':{'data':depth2channel,'nbins':10
              }'''
 covariates = {#'depth2channel':{'data':depth2channel},
              #'clay':{'data':covariates['clay']},
              #'sand':{'data':covariates['sand']},
              #'ndvi':{'data':covariates['ndvi']},
              #'ksat':{'data':np.log10(covariates['SATDK'])},
              #'refsmc':{'data':covariates['REFSMC']},
              #'slope':{'data':covariates['cslope']},
              #'smcmax':{'data':covariates['MAXSMC']},
              #'d2c':{'data':depth2channel},
              #'acc':{'data':np.log10(area)},
              'ti':{'data':mfd_ti},
              }
 hrus = terrain_tools.create_tiles_kmeans(hillslopes_clusters,covariates,100)
 #hrus = terrain_tools.create_nd_histogram(hillslopes_clusters,covariates)-1
 hrus = hrus.astype(np.float32)
 hrus[channels > 0] = np.max(hrus) + 1
 nhrus = int(np.max(hrus) + 1)
 print '%d hrus' % (nhrus),time.time() - time0

 return (hrus,nhrus)

def Compute_HRUs_Semidistributed_Kmeans(covariates,mask,nclusters,hydrobloks_info):

 #Define the number of requested hrus (channel and non-channel)
 nclusters_c = hydrobloks_info['nclusters_c']
 nclusters_nc = hydrobloks_info['nclusters_nc']
 
 #Find the mask for and without channels
 minchannelarea = 3*10**4
 if nclusters_c > 0:
  mask_c = (covariates['carea'] >= minchannelarea) & (mask == True)
  mask_nc = (covariates['carea'] < minchannelarea) & (mask == True)
 else:
  mask_nc = mask

 #Define the covariates
 info_general = {'area':{'data':np.log(covariates['carea'][mask_nc == True]),},
        'slope':{'data':covariates['cslope'][mask_nc == True],},
        'ksat':{'data':covariates['SATDK'][mask_nc == True],},
        'fdir':{'data':covariates['fdir'][mask_nc == True],},
        'bc_b':{'data':covariates['BB'][mask_nc == True],},
        'bc_psisat':{'data':covariates['SATPSI'][mask_nc == True],},
        'sms':{'data':covariates['MAXSMC'][mask == True],},
        'smw':{'data':covariates['WLTSMC'][mask == True],},
        'clay':{'data':covariates['clay'][mask_nc == True],},
        'sand':{'data':covariates['sand'][mask_nc == True],},
        'silt':{'data':covariates['silt'][mask_nc == True],},
        'om':{'data':covariates['om'][mask_nc == True],},
        'lc':{'data':covariates['lc'][mask_nc ==True],},
        'ndvi':{'data':covariates['ndvi'][mask_nc ==True],},
        'ti':{'data':covariates['ti'][mask_nc == True],},
        'dem':{'data':covariates['dem'][mask == True],},
        'demns':{'data':covariates['dem'][mask == True],},
        'bare30':{'data':covariates['bare30'][mask_nc == True],},
        'water30':{'data':covariates['water30'][mask_nc == True],},
        'tree30':{'data':covariates['tree30'][mask_nc == True],},
        'lats':{'data':covariates['lats'][mask_nc == True],},
        'lons':{'data':covariates['lons'][mask_nc == True],},
        }

 #Subset the chosen ones for the non-channels
 info = {}
 for var in hydrobloks_info['covariates']:
  info[var] = info_general[var]

 if nclusters_c > 0:
  info_channels = {}
  #Subset the chosen ones for the channels
  for var in hydrobloks_info['covariates']:
   info_channels[var] = info_general[var]

 #Scale all the variables (Calculate the percentiles
 for var in info:
  #print var, info[var]['data'], np.nanmax(info[var]['data']), np.nanmin(info[var]['data']), np.nanmean(info[var]['data'])
  #if var in ['clay','ndvi','ti','lats','lons']:
  if hydrobloks_info['covariates'][var] == 'p':
   #print hydrobloks_info['covariates'][var] 
   argsort = np.argsort(info[var]['data'])
   pcts = np.copy(info[var]['data'])
   pcts[argsort] = np.linspace(0,1,len(info[var]['data']))
   uniques,counts = np.unique(info[var]['data'],return_counts=True)
   for ival in xrange(uniques.size):
    value = uniques[ival]
    count = counts[ival]
    if count <= 10:continue
    pcts[info[var]['data'] == value] = np.nanmean(pcts[info[var]['data'] == value])
   info[var]['data'][:] = pcts[:]
  else:
   tmp = info[var]['data']
   utmp = np.unique(tmp)
   utmp = utmp[utmp != -9999.0]
   if len(utmp) > 1:
    tmp = (tmp - np.nanmin(tmp))/(np.nanmax(tmp) - np.nanmin(tmp))
    info[var]['data'] = tmp
   else:
    tmp[tmp != -9999.0] = 0.5
    info[var]['data'] = tmp
   #print var, info[var]['data'], np.nanmax(info[var]['data']), np.nanmin(info[var]['data']), np.nanmean(info[var]['data'])
  #argsort = np.argsort(info[var]['data'])
  #pcts = np.copy(info[var]['data'])
  #pcts[argsort] = np.linspace(0,1,len(info[var]['data']))
  #info[var]['data'] = pcts
 if nclusters_c > 0:
  for var in info_channels:
   tmp = info_channels[var]['data']
   utmp = np.unique(tmp)
   utmp = utmp[utmp != -9999.0]
   if len(utmp) > 1:
    tmp = (tmp - np.nanmin(tmp))/(np.nanmax(tmp) - np.nanmin(tmp))
    info_channels[var]['data'] = tmp
   else:
    tmp[tmp != -9999.0] = 0.5
    info_channels[var]['data'] = tmp
   #tmp = (tmp - np.nanmin(tmp))/(np.nanmax(tmp) - np.nanmin(tmp))
   #info_channels[var]['data'] = tmp

 import sklearn.cluster
 #Cluster the non channels regions
 bins,data = [],[]
 X = []
 for id in info:
  #Set all nans to the mean
  
  # Noemi addition
  if np.isnan(np.nanmean(info[id]['data'])) == 1:
    print "ERROR: all values for clustering are Nan, catch:",hydrobloks_info['icatch'],id, info[id]['data']
    info[id]['data'][:] = 0.5

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
 if nclusters_nc > len(info[id]['data']): nclusters_nc=len(info[id]['data']) # Insert Noemi
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

 if nclusters_c > 0:

  #Cluster the channels
  bins,data = [],[]
  X = []
  for id in info_channels:
   #Set all nans to the mean
   if np.isnan(np.nanmean(info_channels[id]['data'])) == 1: 
     print id, hydrobloks_info['icatch'], info_channels[id]['data']
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
 #Define the actual number of clusters
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
 #NLCD2NOAH = {11:17,12:15,21:10,22:10,23:10,24:13,31:16,41:4,42:1,43:5,51:6,52:6,71:10,72:10,73:19,74:19,81:10,82:12,90:11,95:11}

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
 #for lc in np.unique(land_cover)[1:]:
 # land_cover[land_cover == lc] = NLCD2NOAH[lc]
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
 #NLCD2NOAH = {11:17,12:15,21:10,22:10,23:10,24:13,31:16,41:4,42:1,43:5,51:6,52:6,71:10,72:10,73:19,74:19,81:10,82:12,90:11,95:11}
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
  #OUTPUT['hsu']['SATDW'][hsu] = 0.0
  #Average Slope
  OUTPUT['hsu']['slope'][hsu] = np.nanmean(covariates['cslope'][idx])
  #print 'mean',np.mean(np.sin(covariates['cslope'][idx]))
  #print 'arcsin mean',np.arcsin(np.mean(np.sin(covariates['cslope'][idx])))
  #OUTPUT['hsu']['slope'][hsu] = np.arcsin(np.mean(np.sin(covariates['cslope'][idx])))
  #print np.min(covariates['cslope'][idx]),np.mean(covariates['cslope'][idx]),np.max(covariates['cslope'][idx])
  #print OUTPUT['hsu']['slope'][hsu]
  #Topographic index
  OUTPUT['hsu']['ti'][hsu] = np.nanmean(covariates['ti'][idx])
  #DEM
  OUTPUT['hsu']['dem'][hsu] = np.nanmean(covariates['dem'][idx])
  #Average Catchment Area
  OUTPUT['hsu']['carea'][hsu] = np.nanmean(covariates['carea'][idx])
  #Channel?
  #OUTPUT['hsu']['channel'][hsu] = stats.mode(covariates['channels'][idx])[0]
  #Land cover type  
  #OUTPUT['hsu']['land_cover'][hsu] = NLCD2NOAH[95]#stats.mode(covariates['nlcd'][idx])[0][0]]
  OUTPUT['hsu']['land_cover'][hsu] = stats.mode(covariates['lc'][idx])[0][0]
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
 mask2 = cluster_ids_copy < 0 #new
 covariates['fdir'][mask2] = -9999.0 #new
 cluster_ids_copy[mask2] = np.nan #new
 covariates['carea'][mask1] = -9999.0 #new
 covariates['carea'][mask2] = -9999.0 #new
 #plt.imshow(covariates['fdir'])
 #plt.imshow(cluster_ids_copy)
 #plt.show()
 #exit()
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
 #hrus_dst[hrus_dst == -1] = outlet_hru_org[0] #CAREFUL. Designed to push all the extra small outlets to the same exit
 hrus_dst[hrus_dst == -1] = hrus_org[hrus_dst == -1]

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
 #os.system('mkdir -p %s' % soilsdir)
 #soils_lookup = hydrobloks_info['soil_file']
 #fp = open(soils_lookup,'w')
 #fp.write('Soil Parameters\n')
 #fp.write('CUST\n')
 #fp.write("%d,1   'BB      DRYSMC      F11     MAXSMC   REFSMC   SATPSI  SATDK       SATDW     WLTSMC  QTZ    '\n" % nhsus)
 for hsu in np.arange(nhsus):
  #fp.write('%d ' % (hsu+1))
  for var in soil_vars:
   #if var in ['DRYSMC','MAXSMC','REFSMC','WLTSMC','SATDK','SATPSI']:#,'BB','SATPSI']:
   #if var in ['DRYSMC','MAXSMC','REFSMC','WLTSMC','SATDK','SATPSI']:#,'BB','SATPSI']:
   if var in ['ASDASD']:
    donothing = 1
    # fp.write(',%.10f ' % OUTPUT['hsu'][var][hsu])
   #elif var == 'DRYSMC':
   # OUTPUT['hsu'][var][hsu] = 0.0
   #elif var == 'SATDW':
   # OUTPUT['hsu'][var][hsu] = 0.0
   else:
    idx = soils_data['ID'].index(OUTPUT['hsu']['soil_texture_class'][hsu])
    # fp.write(',%.10f ' % soils_data[var][idx])
    OUTPUT['hsu'][var][hsu] = soils_data[var][idx]
  #fp.write('\n')
 #fp.close()

 return #soils_lookup

def Create_and_Curate_Covariates(wbd):

 covariates = {}
 #Read in and curate all the covariates
 for file in wbd['files']:
  covariates[file] = gdal_tools.read_raster(wbd['files'][file])
  if file == 'cslope':
   mask = covariates[file] == 0.0
   covariates[file][mask] = 0.000001

 #Map the NLCD to IGBP
 #NLCD2NOAH = {11:17,12:15,21:10,22:10,23:10,24:13,31:16,41:4,42:1,43:5,51:6,52:6,71:10,72:10,73:19,74:19,81:10,82:12,90:11,95:11}
 #tmp = np.copy(covariates['nlcd'])
 #for lc in np.unique(covariates['nlcd']):
 # tmp[covariates['nlcd'] == lc] = NLCD2NOAH[lc]
 #covariates['nlcd'][:] = tmp[:]


 #print covariates['ndvi'],covariates['lc']
 #Curate the NDVI to match land cover 
 for lc in np.unique(covariates['lc']):
  masklc = covariates['lc'] == lc
  covariates['ndvi'][masklc] = np.nanmean(covariates['ndvi'][masklc])
 #print covariates['ndvi']


 #Create lat/lon grids
 lats = np.linspace(wbd['bbox']['minlat']+wbd['bbox']['res']/2,wbd['bbox']['maxlat']-wbd['bbox']['res']/2,covariates['ti'].shape[0])
 lons = np.linspace(wbd['bbox']['minlon']+wbd['bbox']['res']/2,wbd['bbox']['maxlon']-wbd['bbox']['res']/2,covariates['ti'].shape[1])

 #Constrain the lat/lon grid by the effective resolution (4km...)
 nres = int(np.floor(2000.0/30.0))
 for ilat in xrange(0,lats.size,nres):
  lats[ilat:ilat+nres] = np.mean(lats[ilat:ilat+nres])
 for ilon in xrange(0,lons.size,nres):
  lons[ilon:ilon+nres] = np.mean(lons[ilon:ilon+nres])
 
 #Need to fix so that it doesn't suck up all the clustering:
 lats, lons = np.meshgrid(lats, lons)
 covariates['lats'] = lats.T
 covariates['lons'] = lons.T
 
 #Add sti
 sti = covariates['ti'] - np.log(0.1*covariates['SATDK']) + np.log(np.mean(0.1*covariates['SATDK']))
 covariates['sti'] = sti

 #Define the mask
 mask = np.copy(covariates['mask'])
 mask[mask >= 0] = 1
 mask[mask < 0] = 0
 mask = mask.astype(np.bool)

 #Set all nans to the mean
 for var in covariates:
  mask1 = (np.isinf(covariates[var]) == 0) & (np.isnan(covariates[var]) == 0) 
  mask0 = (np.isinf(covariates[var]) == 1) | (np.isnan(covariates[var]) == 1)

  if var in ['fdir','nlcd','TEXTURE_CLASS','lc']:
   if len(covariates[var][mask1]) < 1: print "Full of Nan's Error: 1", var, covariates[var],  wbd['files'][var] # Noemi insert
   covariates[var][mask0] = -9999.0 #stats.mode(covariates[var][mask1])[0][0]
  else:
   covariates[var][mask0] = -9999.0 #np.mean(covariates[var][mask1])

 #Set everything that is -9999 to the mean
 for var in covariates:
  if var in ['fdir','nlcd','TEXTURE_CLASS','lc']:  
   if len(covariates[var][covariates[var] != -9999.0])  < 1: print "Full of Nan's Error: ", var, covariates[var],  wbd['files'][var] # Noemi insert
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
  if hydrobloks_info['clustering_type'] == 'kmeans':
   (cluster_ids,nclusters) = Compute_HRUs_Semidistributed_Kmeans(covariates,mask,nclusters,hydrobloks_info)
  elif hydrobloks_info['clustering_type'] == 'hillslope':
   (cluster_ids,nclusters) = Compute_HRUs_Semidistributed_Hillslope(covariates,mask,nclusters,hydrobloks_info)
  elif hydrobloks_info['clustering_type'] == 'basin':
   (cluster_ids,nclusters) = Compute_HRUs_Semidistributed_Basin(covariates,mask,nclusters,hydrobloks_info)
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
 #print "Creating the soil file"
 Create_Soils_File(hydrobloks_info,OUTPUT,input_dir)
 #soils_lookup = Create_Soils_File(hydrobloks_info,OUTPUT,input_dir)

 #Add the new number of clusters
 OUTPUT['nclusters'] = nclusters
 OUTPUT['mask'] = mask

 return OUTPUT

def Prepare_Meteorology_Fulldistributed(workspace,wbd,OUTPUT,input_dir,info,hydrobloks_info):

 #Assign other variables
 mask = OUTPUT['mask']

 mapping_info = {}
 #Calculate the fine to coarse scale mapping
 for data_var in wbd['files_meteorology']:

  #Define the variable name
  var = data_var#data_var.split('_')[1]
  mapping_info[var] = {}

  #Read in the coarse and fine mapping
  file_coarse = '%s/%s_ea_coarse.tif' % (workspace,var)
  file_fine = '%s/%s_ea_fine.tif' % (workspace,var)
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
 mapping_info = {}
 #Calculate the fine to coarse scale mapping
 for data_var in wbd['files_meteorology']:
  
  #Define the variable name
  var = data_var#data_var.split('_')[1]
  mapping_info[var] = {}

  #Read in the coarse and fine mapping
  file_coarse = '%s/%s_latlon_coarse.tif' % (workspace,data_var)
  file_fine = '%s/%s_ea_fine.tif' % (workspace,data_var)
  mask_coarse = gdal_tools.read_raster(file_coarse)
  mask_fine = gdal_tools.read_raster(file_fine)
  nlat = mask_coarse.shape[0]
  nlon = mask_coarse.shape[1]

  #Compute the mapping for each hsu
  for hsu in np.arange(hydrobloks_info['nclusters']):
   idx = OUTPUT['hsu_map'] == hsu
   #print "Catch:",hydrobloks_info['icatch'], "HRU: ", mask_fine[idx].astype(np.int)
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
   mapping_info[var][hsu] = {'pcts':pcts,'coords':coords}

 #Iterate through variable creating forcing product per HSU
 idate = info['time_info']['startdate']
 fdate = info['time_info']['enddate']
 dt = info['time_info']['dt']
 nt = int(3600*24/dt)*((fdate - idate).days+1)
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
  fidate = ' '.join(fp.variables['t'].units.split(' ')[2::])
  #dates = nc.num2date(fp.variables['t'][:],units='hours since %02d-%02d-%02d 00:00:00' % (idate.year,idate.month,idate.day))
  dates = nc.num2date(fp.variables['t'][:],units='hours since %s' % fidate)
  #dates = nc.num2date(fp['t'][:],units='hours since %02d-%02d-%02d 00:00:00' % (idate.year,idate.month,idate.day))
  mask_dates = (dates >= idate) & (dates <= fdate)
  data = np.ma.getdata(fp.variables[var][mask_dates,:,:])
  fp.close()
  #Assing to hsus
  for hsu in mapping_info[var]:
   #print data_var,data, data.shape, hsu,mapping_info[var][hsu]['pcts'],mapping_info[var][hsu]['coords'],
   pcts = mapping_info[var][hsu]['pcts']
   coords = mapping_info[var][hsu]['coords']
   coords[0][coords[0] >= data.shape[1]] = data.shape[1] - 1
   coords[1][coords[1] >= data.shape[2]] = data.shape[2] - 1
   tmp = data[:,coords[0],coords[1]]
   #m1 = tmp < -999
   #m2 = tmp > -999
   tmp = pcts*tmp
   meteorology[data_var][:,hsu] = np.sum(tmp,axis=1)
   #print data_var,np.unique(meteorology[data_var][:,:])

  #Write the meteorology to the netcdf file (single chunk for now...)
  grp = hydrobloks_info['input_fp'].groups['meteorology']
  grp.createVariable(var,'f4',('time','hsu'))
  grp.variables[data_var][:] = meteorology[data_var][:]

 return 

def Read_Metadata_File(file):

 import json
 metadata = json.load(open(file))

 return metadata
