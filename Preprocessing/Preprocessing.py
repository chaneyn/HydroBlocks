import warnings
warnings.filterwarnings('ignore')
import sys
sys.path.append('Tools')
#import cPickle as pickle
import pickle
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
#import matplotlib.pyplot as plt
#import random
import gc

dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append('%s/../HydroBlocks/pyHWU/' % dir )
import management_funcs as mgmt_funcs

def plot_data(data):

 import matplotlib.pyplot as plt
 data = np.ma.masked_array(data,data==-9999)
 plt.imshow(data)
 plt.colorbar()
 plt.show()

 return


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
 nhru = hydrobloks_info['nhru']
 #ncores = hydrobloks_info['ncores']
 icatch = hydrobloks_info['icatch']

 #Get metadata
 print('%s/mask_latlon.tif' % workspace)
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
  'irrig_land':'%s/irrig_land_ea.tif' % workspace,
  'dbedrock':'%s/dbedrock_ea.tif' % workspace
  }
 if hydrobloks_info['water_management']['hwu_agric_flag']:
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

 if hydrobloks_info['water_management']['hwu_flag'] == True:
  wbd['files_water_use'] = {}
  if hydrobloks_info['water_management']['hwu_domest_flag']:
   wbd['files_water_use']['domestic']   = '%s/domestic.nc' % workspace
  if hydrobloks_info['water_management']['hwu_indust_flag']:
   wbd['files_water_use']['industrial'] = '%s/industrial.nc' % workspace
  if hydrobloks_info['water_management']['hwu_lstock_flag']:
   wbd['files_water_use']['livestock']  = '%s/livestock.nc' % workspace
  

 #Create the clusters and their connections
 output = Create_Clusters_And_Connections(workspace,wbd,output,input_dir,nhru,info,hydrobloks_info)

 #Extract the meteorological forcing
 print("Preparing the meteorology")
 if hydrobloks_info['model_type'] == 'semi':
  Prepare_Meteorology_Semidistributed(workspace,wbd,output,input_dir,info,hydrobloks_info)
 elif hydrobloks_info['model_type'] == 'full':
  Prepare_Meteorology_Fulldistributed(workspace,wbd,output,input_dir,info,hydrobloks_info)

 #Extract the water use demands
 print("Preparing the water use")
 if hydrobloks_info['model_type'] == 'semi':
  if hydrobloks_info['water_management']['hwu_flag'] == True:
   Prepare_Water_Use_Semidistributed(workspace,wbd,output,input_dir,info,hydrobloks_info) 
 elif hydrobloks_info['model_type'] == 'full':
   exit('Error: Human Water Management not implemented in the fully distributed mode yet.')

 #Write out the files to the netcdf file
 fp = hydrobloks_info['input_fp']
 data = output

 #Write out the metadata
 grp = fp.createGroup('metadata')
 grp.latitude = (wbd['bbox']['minlat'] + wbd['bbox']['maxlat'])/2
 lon = (wbd['bbox']['minlon'] + wbd['bbox']['maxlon'])/2 
 if lon < 0:lon += 360
 grp.longitude = lon
 metadata = gdal_tools.retrieve_metadata(wbd['files']['mask']) 
 grp.dx = metadata['resx']
 #grp.longitude = (360.0+(wbd['bbox']['minlon'] + wbd['bbox']['maxlon'])/2)

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
 os.system('gdalwarp -ot Int16 -tr %.16f %.16f -dstnodata %.16f -t_srs \'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs \' -te %.16f %.16f %.16f %.16f %s %s >> %s 2>&1' % (res,res,metadata['nodata'],minlon,minlat,maxlon,maxlat,file_ca,file_ll,log))

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
        'BB','F11','SATPSI','SATDW','QTZ','clay',
        'WLTSMC','MAXSMC','DRYSMC','REFSMC','SATDK',
        'mannings','m','psoil','pksat','sdmax']

 if hydrobloks_info['water_management']['hwu_agric_flag']:
  for var in ['centroid_lats', 'centroid_lons', 'irrig_land', 'start_growing_season', 'end_growing_season']:
    vars.append(var)

 for var in vars:
  grp.createVariable(var,'f4',('hsu',))#,zlib=True)
  grp.variables[var][:] = data['hsu'][var]

 if hydrobloks_info['water_management']['hwu_flag']:
  grp.createVariable('hru_min_dist','f4',('hsu','hsu'))#,zlib=True)
  grp.variables['hru_min_dist'][:] = data['hsu']['hru_min_dist']


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

def Compute_HRUs_Fulldistributed(covariates,mask,nhru):

 cluster_ids = np.empty(covariates['ti'].shape)
 cluster_ids[:] = -9999
 cluster_ids[mask == True] = np.arange(nhru)

 return (cluster_ids,)

def Compute_HRUs_Semidistributed_Kmeans(covariates,mask,nhru,hydrobloks_info,wbd):

 #Define the number of requested hrus (channel and non-channel)
 nhru = hydrobloks_info['nhru']
 
 #Find the mask for and without channels
 mask_nc = mask
 
 #Mask regions without soil properties
 tmp_clay = gdal_tools.read_raster(wbd['files']['clay'])
 tmp_lc = gdal_tools.read_raster(wbd['files']['lc'])
 nodata = (tmp_clay == -9999) & (tmp_lc == 17)
 mask_nc[nodata]=False
 del tmp_clay,tmp_lc
 gc.collect()

 
 #Define the covariates
 info_general = {}
 for var in covariates:
   if var == 'carea':    info_general['area'] =  {'data':np.log(covariates['carea'][mask_nc == True]),}
   elif var == 'cslope': info_general['slope'] = {'data':covariates['cslope'][mask_nc == True],}
   elif var == 'SATDK':  info_general['ksat'] =  {'data':covariates['SATDK'][mask_nc == True],}
   elif var == 'BB':     info_general['bc_b'] =  {'data':covariates['BB'][mask_nc == True],}
   elif var == 'SATPSI': info_general['bc_psisat'] = {'data':covariates['SATPSI'][mask_nc == True],}
   elif var == 'MAXSMC': info_general['sms'] =   {'data':covariates['MAXSMC'][mask == True],}
   elif var == 'WLTSMC': info_general['smw'] =   {'data':covariates['WLTSMC'][mask == True],}
   else: info_general[var] =   {'data':covariates[var][mask == True],}

 # check if lc is a covariates, and disagregate it in classes
 if 'lc' in hydrobloks_info['covariates']:
  for lc in np.unique(covariates['lc'][mask == True]):
   if lc >= 0 :
    vnam = 'lc_%i' % lc
    masklc = covariates['lc'] == lc
    covariates[vnam] = np.zeros(covariates['lc'].shape)
    covariates[vnam][masklc] = 1.0
    hydrobloks_info['covariates'][vnam] = 'n'
    info_general[vnam] = {'data':covariates[vnam][mask_nc == True],}
  del hydrobloks_info['covariates']['lc']    

 #Subset the chosen ones for the non-channels
 info = {}
 for var in hydrobloks_info['covariates']:
  info[var] = info_general[var]

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
   for ival in range(uniques.size):
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



 import sklearn.cluster
 #Cluster the non channels regions
 bins,data = [],[]
 X = []
 for id in info:
  #Set all nans to the mean
  
  # Noemi addition
  if np.isnan(np.nanmean(info[id]['data'])) == 1:
    print("ERROR: all values for clustering are Nan, catch:",hydrobloks_info['icatch'],id, info[id]['data'])
    info[id]['data'][:] = 0.5

  info[id]['data'][np.isnan(info[id]['data']) == 1] = np.nanmean(info[id]['data'])
  X.append(info[id]['data'])


 time0 = time.time()
 X = np.array(X).T
 #Subsample the array
 np.random.seed(1)
 #minsamples = 3*10**5
 minsamples = X.shape[0] #Noemi
 if X.shape[0] > minsamples:
  Xf = X[np.random.choice(np.arange(X.shape[0]),minsamples),:]
 else:
  Xf = X

 #Initialize all points at the 0.5 point
 if nhru > len(info[id]['data'])/10:
   nhru=int(len(info[id]['data'])/10)  # Noemi
 init = 0.5*np.ones((nhru,Xf.shape[1]))
 batch_size = 25*nhru
 init_size = 3*batch_size
 clf = sklearn.cluster.MiniBatchKMeans(nhru,random_state=1,init=init,batch_size=batch_size,init_size=init_size)
 clf.fit(Xf)#
 clf_output = clf.predict(X)

 #Reassign the ids
 clf_output_copy = np.copy(clf_output)
 for cid in range(len(np.unique(clf_output))):
  clf_output[clf_output_copy == np.unique(clf_output)[cid]] = cid
 cluster_ids = np.empty(covariates['ti'].shape)
 cluster_ids[:] = -9999
 cluster_ids[mask_nc == True] = clf_output
 nhru_old = nhru
 nhru = np.unique(cluster_ids).size-1

 #Redefine the number of clusters (We are sampling regions that just don't have data...)
 print('clustering %d->%d' % (nhru_old,nhru),time.time() - time0)

 #Reorder according to areas
 areas = []
 for cid in range(nhru):
  msk = cluster_ids == cid
  #areas.append(np.nanmean(covariates['carea'][msk]))
  areas.append(np.nanmax(covariates['carea'][msk]))
 argsort = np.argsort(np.array(areas))
 cluster_ids_new = np.copy(cluster_ids)
 for cid in range(nhru):
  msk = cluster_ids == argsort[cid]
  cluster_ids_new[msk] = cid
 cluster_ids = np.copy(cluster_ids_new)
 areas = []
 for cid in range(nhru):
  msk = cluster_ids == cid
  areas.append(np.nanmean(covariates['carea'][msk]))

 return (cluster_ids,nhru)

def Assign_Parameters_Fulldistributed(covariates,metadata,hydrobloks_info,OUTPUT,cluster_ids,mask):

 nhru = hydrobloks_info['nhru']
 #Initialize the arrays
 vars = ['area','area_pct','BB','DRYSMC','F11','MAXSMC','REFSMC','SATPSI',
         'SATDK','SATDW','WLTSMC','QTZ','slope','ti','dem','carea','channel',
         'land_cover','soil_texture_class',
         'mannings','m','psoil','pksat','sdmax']
 #vars = ['area','area_pct','BB','DRYSMC','F11','MAXSMC','REFSMC','SATPSI',
 #        'SATDK','SATDW','WLTSMC','QTZ','slope','ti','dem','carea','channel',         'vchan','vof','land_cover','soil_texture_class']
 OUTPUT['hsu'] = {}
 for var in vars:
  OUTPUT['hsu'][var] = np.zeros(nhru) 

 #Calculate area per hsu
 OUTPUT['hsu']['area'][:] = metadata['resx']**2
 #Calculate area percentage per hsu
 OUTPUT['hsu']['area_pct'][:] = 100*OUTPUT['hsu']['area']/(metadata['resx']**2*nhru)
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
 OUTPUT['hsu']['m'][:] = covariates['dbedrock'][mask] #Form of the exponential decline in conductivity (0.01-1.0)
 OUTPUT['hsu']['pksat'][:] = 1.0 #saturated hydraulic conductivity scalar multiplier (0.1-1.0)
 OUTPUT['hsu']['psoil'][:] = 1.0 #soil hydraulic properties (residual,wilting,field capacity, and porosity) (0.1-10.0)
 OUTPUT['hsu']['sdmax'][:] = 5.0 #maximum effective deficit of subsurface saturated zone (0.1-10.0)
 OUTPUT['hsu']['mannings'][covariates['carea'][mask] >= 100000.0] = 0.03 #manning's n for channel flow (0.01-0.1)
 OUTPUT['hsu']['mannings'][covariates['carea'][mask] < 100000.0] = 0.15 #manning's n for overland flow (0.01-0.8)

 return OUTPUT

def Assign_Parameters_Semidistributed(covariates,metadata,hydrobloks_info,OUTPUT,cluster_ids,mask):

 nhru = hydrobloks_info['nhru']
 #Initialize the arrays
 vars = ['area','area_pct','BB','DRYSMC','F11','MAXSMC','REFSMC','SATPSI',
         'SATDK','SATDW','WLTSMC','QTZ','slope','ti','dem','carea','channel',
         'land_cover','soil_texture_class','clay',
         'mannings','m','psoil','pksat','sdmax']

 if hydrobloks_info['water_management']['hwu_agric_flag']:
  for var in ['centroid_lats', 'centroid_lons', 'irrig_land', 'start_growing_season', 'end_growing_season']:
    vars.append(var)
 #import numpy as np
 OUTPUT['hsu'] = {}
 if hydrobloks_info['water_management']['hwu_flag']: OUTPUT['hsu']['hru_min_dist'] = np.zeros((nhru,nhru))
 for var in vars:
   OUTPUT['hsu'][var] = np.zeros(nhru)


 #Metadata
 for hsu in np.arange(nhru):
  #Set indices
  idx = np.where(cluster_ids == hsu)
  #Calculate area per hsu
  OUTPUT['hsu']['area'][hsu] = metadata['resx']**2*idx[0].size
  #Calculate area percentage per hsu
  OUTPUT['hsu']['area_pct'][hsu] = 100*OUTPUT['hsu']['area'][hsu]/(metadata['resx']**2*mask[mask].size)
  #Soil properties
  for var in ['BB','DRYSMC','F11','MAXSMC','REFSMC','SATPSI','SATDK','SATDW','WLTSMC','QTZ','clay']:
   #if var in ['SATDK']:
   # OUTPUT['hsu'][var][hsu] = stats.mstats.hmean(covariates[var][idx])/3600.0/1000.0 #mm/hr -> m/s
   #elif var in ['SATDW']:
   # OUTPUT['hsu'][var][hsu] = stats.mstats.hmean(covariates[var][idx])/3600.0/(1000.0*1000.0)  #mm2/hr -> m2/s #Noemi
   #elif var in ['SATPSI']:
   # OUTPUT['hsu'][var][hsu] = stats.mstats.gmean(covariates[var][idx])/1000.0 #mm -> m  # Noemi
   if var in ['SATDK','SATDW']:
    OUTPUT['hsu'][var][hsu] = stats.mstats.hmean(covariates[var][idx])/3600.0/1000.0 #mm/hr -> m/s
   else:
    OUTPUT['hsu'][var][hsu] = stats.mstats.gmean(covariates[var][idx])
  #Average Slope
  OUTPUT['hsu']['slope'][hsu] = np.nanmean(covariates['cslope'][idx])
  #Topographic index
  OUTPUT['hsu']['ti'][hsu] = np.nanmean(covariates['ti'][idx])
  #DEM
  OUTPUT['hsu']['dem'][hsu] = np.nanmean(covariates['dem'][idx])
  #Average Catchment Area
  OUTPUT['hsu']['carea'][hsu] = np.nanmean(covariates['carea'][idx])
  #Channel?
  #OUTPUT['hsu']['channel'][hsu] = stats.mode(covariates['channels'][idx])[0]
  #Land cover type 
  tmp = covariates['lc'][idx]
  tmp = tmp[tmp>=1]
  if len(tmp) >= 1 :
   OUTPUT['hsu']['land_cover'][hsu] = stats.mode(tmp)[0][0]
  else:
   OUTPUT['hsu']['land_cover'][hsu] = 17  # if there is no valid value, set to water #Noemi
      
  #Soil texture class
  OUTPUT['hsu']['soil_texture_class'][hsu] = stats.mode(covariates['TEXTURE_CLASS'][idx])[0][0]
  #Define the estimate for the model parameters
  OUTPUT['hsu']['m'][hsu] = np.nanmean(covariates['dbedrock'][idx]) #0.1 #Form of the exponential decline in conductivity (0.01-1.0)
  OUTPUT['hsu']['pksat'][hsu] = 1.0 #saturated hydraulic conductivity scalar multiplier (0.1-1.0)
  OUTPUT['hsu']['psoil'][hsu] = 1.0 #soil hydraulic properties (residual,wilting,field capacity, and porosity) (0.1-10.0)
  OUTPUT['hsu']['sdmax'][hsu] = 5.0 #maximum effective deficit of subsurface saturated zone (0.1-10.0)
  if np.max(covariates['carea'][idx]) >= 100000.0: OUTPUT['hsu']['mannings'][hsu] = 0.03 #manning's n for channel flow (0.01-0.1)
  else: OUTPUT['hsu']['mannings'][hsu] = 0.15 #manning's n for overland flow (0.01-0.8)

  # Water Management Variables
  if hydrobloks_info['water_management']['hwu_agric_flag'] == True:
   # Irrigation: 1 Irrigated, 2 paddy crop
   irrig_vec = covariates['irrig_land'][idx].flatten()
   nii = np.nansum([x==1 for x in irrig_vec])
   npi = np.nansum([x==2 for x in irrig_vec])
   ttt = np.nansum([x>=0 for x in irrig_vec])
   iratio = 0
   if ttt > 0: iratio = float((nii+npi))/ttt
   if   iratio >= 0.50 and nii > npi: OUTPUT['hsu']['irrig_land'][hsu] = 1
   elif iratio >= 0.50 and nii < npi: OUTPUT['hsu']['irrig_land'][hsu] = 2
   else: OUTPUT['hsu']['irrig_land'][hsu] = 0
 
    # Crop Calendar
   OUTPUT['hsu']['start_growing_season'][hsu] = int(stats.mode(covariates['start_growing_season'][idx])[0][0])
   OUTPUT['hsu']['end_growing_season'][hsu] = int(stats.mode(covariates['end_growing_season'][idx])[0][0])
  
  if hydrobloks_info['water_management']['hwu_flag'] == True:
   #HRU Centroids for water management
   OUTPUT['hsu']['centroid_lons'][hsu] = np.nanmean(covariates['lons'][idx])
   OUTPUT['hsu']['centroid_lats'][hsu] = np.nanmean(covariates['lats'][idx])
   #print OUTPUT['hsu']['centroid_lons']  

 if hydrobloks_info['water_management']['hwu_flag'] == True: 
  for hsu in np.arange(nhru):
   #HRU distance between the centroids of the hru and all the other hrus 
   OUTPUT['hsu']['hru_min_dist'][hsu,:] = mgmt_funcs.calculate_min_distance(hsu, nhru, cluster_ids, covariates['lats'], covariates['lons'], OUTPUT['hsu']['centroid_lats'], OUTPUT['hsu']['centroid_lons'])
   OUTPUT['hsu']['hru_min_dist'][hsu,hsu] = 0.0
   
 return OUTPUT


def Calculate_Flow_Matrix(covariates,cluster_ids,nhru):

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
 max_nhru = np.sum(cluster_ids >= 0)
 (hrus_dst,hrus_org,outlet_icoord,outlet_jcoord,outlet_hru,outlet_d8) = mt.preprocessor.calculate_connections_d8(cluster_ids_copy,covariates['fdir'],covariates['carea'],nhru,max_nhru)
 #Only use the non -9999 values
 hrus_dst = hrus_dst[hrus_dst != -9999]-1
 hrus_org = hrus_org[hrus_org != -9999]-1
 outlet_icoord = outlet_icoord[outlet_icoord != -9999]-1
 outlet_jcoord = outlet_jcoord[outlet_jcoord != -9999]-1
 outlet_hru_org = outlet_hru[outlet_hru != -9999] - 1
 outlet_d8 = outlet_d8[outlet_d8 != -9999]

 #Create hrus for the outlets
 outlet_hru_dst_summary = np.arange(nhru,nhru+np.unique(outlet_hru_org).size)
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
 hrus_dst[hrus_dst == -1] = hrus_org[hrus_dst == -1]

 #Prepare the sparse matrix
 flow_matrix = sparse.coo_matrix((np.ones(hrus_dst.size),(hrus_org,hrus_dst)),dtype=np.float32)
 flow_matrix = flow_matrix.tocsr()
 #Normalize the rows (sum to 1)
 fm_sum = flow_matrix.sum(axis=1)
 fm_data = flow_matrix.data
 fm_indptr = flow_matrix.indptr
 for row in range(fm_sum.size):
  fm_data[fm_indptr[row]:fm_indptr[row+1]] = fm_data[fm_indptr[row]:fm_indptr[row+1]]/fm_sum[row]
 flow_matrix.data = fm_data

 return (flow_matrix.T,outlet)

def Calculate_HRU_Connections_Matrix(covariates,cluster_ids,nhru,dx):

 #Define spatial resolution
 res = dx
 
 horg = []
 hdst = []
 #Count the connections
 for i in range(cluster_ids.shape[0]):
  for j in range(cluster_ids.shape[1]):
   h1 = cluster_ids[i,j]
   if h1 == -9999:continue
   #up
   if (i+1) < cluster_ids.shape[0]:
    h2 = cluster_ids[i+1,j]
    if h2 != -9999:
     horg.append(h1)
     hdst.append(h2)
   #down
   if (i-1) > 0:
    h2 = cluster_ids[i-1,j]
    if h2 != -9999:
     horg.append(h1)
     hdst.append(h2)
   #left
   if (j-1) > 0:
    h2 = cluster_ids[i,j-1]
    if h2 != -9999:
     horg.append(h1)
     hdst.append(cluster_ids[i,j-1])
   #right
   if (j+1) < cluster_ids.shape[1]:
    h2 = cluster_ids[i,j+1]
    if h2 != -9999:
     horg.append(h1)
     hdst.append(cluster_ids[i,j+1])
 horg = np.array(horg)
 hdst = np.array(hdst)

 #Prepare the sparse matrix
 cmatrix = sparse.coo_matrix((np.ones(hdst.size),(horg,hdst)),dtype=np.float32)
 cmatrix = cmatrix.tocsr()

 #Prepare length, width, and ksat matrices
 wmatrix = cmatrix.copy()
 wmatrix[:] = res*wmatrix[:]

 #Prepare output dictionary
 cdata = {'width':wmatrix.T,}

 return cdata

def Create_and_Curate_Covariates(wbd,hydrobloks_info):

 covariates = {}
 #Read in and curate all the covariates
 for file in wbd['files']:
  if os.path.isfile(wbd['files'][file]): 
   covariates[file] = gdal_tools.read_raster(wbd['files'][file])
   if file == 'cslope':
    mask = covariates[file] == 0.0
    covariates[file][mask] = 0.000001
  
 #Curate the NDVI to match LC # Noemi
 #tmp = covariates['ndvi']
 #tmp[tmp== -9999.0] = np.nan
 #for lc in np.unique(covariates['lc']):
 # masklc = covariates['lc'] == lc
 # covariates['ndvi'][masklc] = np.nanmean(tmp[masklc])
 
 # Create sub-covariates from the land cover

 #Create lat/lon grids
 lats = np.linspace(wbd['bbox']['minlat']+wbd['bbox']['res']/2,wbd['bbox']['maxlat']-wbd['bbox']['res']/2,covariates['ti'].shape[0])
 lons = np.linspace(wbd['bbox']['minlon']+wbd['bbox']['res']/2,wbd['bbox']['maxlon']-wbd['bbox']['res']/2,covariates['ti'].shape[1])

 #Constrain the lat/lon grid by the effective resolution (4km...)
 nres = int(np.floor(90.0/30.0))
 for ilat in range(0,lats.size,nres):
  lats[ilat:ilat+nres] = np.mean(lats[ilat:ilat+nres])
 for ilon in range(0,lons.size,nres):
  lons[ilon:ilon+nres] = np.mean(lons[ilon:ilon+nres])
 
 #Need to fix so that it doesn't suck up all the clustering:
 lats, lons = np.meshgrid(lats, lons)
 covariates['lats'] = lats.T
 covariates['lons'] = lons.T
 
 #Add sti
 satdk = covariates['SATDK']
 satdk[satdk==-9999.0] = np.nan
 sti = covariates['ti'] - np.log(0.1*satdk) + np.log(np.nanmean(0.1*satdk))
 sti[sti==np.nan]=-9999.0
 covariates['sti'] = sti

 #Define the mask
 mask = np.copy(covariates['mask'])
 mask[mask >= 0] = 1
 mask[mask < 0] = 0
 mask = mask.astype(np.bool)

 #Set all nans to the mean
 for var in covariates:
  covariates[var][mask <= 0] = -9999.0
  mask1 = (np.isinf(covariates[var]) == 0) & (np.isnan(covariates[var]) == 0) 
  mask0 = (np.isinf(covariates[var]) == 1) | (np.isnan(covariates[var]) == 1)
  covariates[var][mask0] = -9999.0# stats.mode(covariates[var][mask1])[0][0]

 #Set everything that is -9999 to the mean
 for var in covariates:
  if len(covariates[var][covariates[var] != -9999.0]) < 1: 
   print("Error: Covariate %s is full of Nan's" % (var), covariates[var], wbd['files'][var]) # Noemi insert
   if var == 'lc':
    covariates[var][mask == 1] = 17  # Water
   else:
    covariates[var][mask == 1] = 0.5
    print("Error: Covariate was set to 0.5")
  if var in ['fdir','nlcd','TEXTURE_CLASS','lc','irrig_land','bare30','water30','tree30']: 
   covariates[var][covariates[var] == -9999.0] = stats.mode(covariates[var][covariates[var] != -9999.0])[0][0]  
  else:
   covariates[var][covariates[var] == -9999.0] = np.mean(covariates[var][covariates[var] != -9999.0])

 #Set everything outside of the mask to -9999
 for var in covariates:
  covariates[var][mask <= 0] = -9999.0
 
 return (covariates,mask)

def Create_Clusters_And_Connections(workspace,wbd,output,input_dir,nhru,info,hydrobloks_info):
 
 print("Creating and curating the covariates")
 (covariates,mask) = Create_and_Curate_Covariates(wbd,hydrobloks_info)

 #Determine the HRUs (clustering if semidistributed; grid cell if fully distributed)
 print("Computing the HRUs")
 if hydrobloks_info['model_type'] == 'semi':
  (cluster_ids,nhru) = Compute_HRUs_Semidistributed_Kmeans(covariates,mask,nhru,hydrobloks_info,wbd)
 elif hydrobloks_info['model_type'] == 'full':
  nhru = np.sum(mask == True)
  (cluster_ids,) = Compute_HRUs_Fulldistributed(covariates,mask,nhru)
 hydrobloks_info['nhru'] = nhru
  
 #Create the netcdf file
 file_netcdf = hydrobloks_info['input_file']
 hydrobloks_info['input_fp'] = nc.Dataset(file_netcdf, 'w', format='NETCDF4')

 #Retrieve some metadata
 metadata = gdal_tools.retrieve_metadata(wbd['files']['mask'])
 resx = metadata['resx']

 #Create the dimensions (netcdf)
 idate = hydrobloks_info['idate']
 fdate = hydrobloks_info['fdate']
 dt = hydrobloks_info['dt']
 ntime = 24*3600*((fdate - idate).days+1)/dt
 nhsu = hydrobloks_info['nhru']
 hydrobloks_info['input_fp'].createDimension('hsu',nhsu)
 hydrobloks_info['input_fp'].createDimension('time',ntime)
 
 #Create the groups (netcdf)
 hydrobloks_info['input_fp'].createGroup('meteorology')
 hydrobloks_info['input_fp'].createGroup('water_use')

 #Prepare the flow matrix (dynamic topmodel)
 print("Calculating the flow matrix")
 (flow_matrix,outlet) = Calculate_Flow_Matrix(covariates,cluster_ids,nhru)

 #Prepare the hru connections matrix (darcy clusters)
 cmatrix = Calculate_HRU_Connections_Matrix(covariates,cluster_ids,nhru,resx)

 #Define the metadata
 metadata = gdal_tools.retrieve_metadata(wbd['files']['ti'])

 #Make the output dictionary for the basin
 OUTPUT = {'hsu':{},'metadata':metadata,'mask':mask,'flow_matrix':flow_matrix,'cmatrix':cmatrix}
 OUTPUT['outlet'] = outlet

 #Remember the map of hrus
 OUTPUT['hsu_map'] = cluster_ids

 #Assign the model parameters
 print("Assigning the model parameters")
 if hydrobloks_info['model_type'] == 'semi':
  OUTPUT = Assign_Parameters_Semidistributed(covariates,metadata,hydrobloks_info,OUTPUT,cluster_ids,mask)
 elif hydrobloks_info['model_type'] == 'full':
  OUTPUT = Assign_Parameters_Fulldistributed(covariates,metadata,hydrobloks_info,OUTPUT,cluster_ids,mask)

 #Add the new number of clusters
 OUTPUT['nhru'] = nhru
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

  print(data_var,"Creating the i and j mapping")
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

  print(data_var,"Extracting all the data")
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

  print(data_var,"Assigning the data")
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
  for hsu in np.arange(hydrobloks_info['nhru']):
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
  meteorology[data_var] = np.zeros((nt,hydrobloks_info['nhru']))
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
  dates = [datetime.datetime(int(nc_idate[0]),int(nc_idate[1]),1)]
  for it in range(0,nc_nt-1): dates.append( dates[it] + datetime.timedelta(hours=nc_step) )
  dates=np.array(dates)
  startdate = datetime.datetime(idate.year,idate.month,1,0)
  enddate = datetime.datetime(fdate.year,fdate.month,31,23)
  mask_dates = (dates >= startdate) & (dates <= enddate)
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
  grp.createVariable(var,'f4',('time','hsu'))#,zlib=True)
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


def Prepare_Water_Use_Semidistributed(workspace,wbd,OUTPUT,input_dir,info,hydrobloks_info):

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
  hrus_lc = np.copy(OUTPUT['hsu_map'])
  for hsu in np.arange(hydrobloks_info['nhru']): 
   idx = OUTPUT['hsu_map'] == hsu
   hrus_lc[idx]= OUTPUT['hsu']['land_cover'][hsu]
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

  #hru_map = np.copy(OUTPUT['hsu_map'])
  #hru_map[hru_map<0]=np.nan
  #plt.imshow(hru_map); plt.show()
  #plt.imshow(hrus_lc); plt.show()
  #plt.imshow(data); plt.show()
  
  #Compute the mapping for each hsu
  for hsu in np.arange(hydrobloks_info['nhru']):
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

 #Iterate through variable creating water use product per HSU
 idate = info['time_info']['startdate']
 fdate = info['time_info']['enddate']
 dt = info['time_info']['dt']
 nt = int(3600*24/dt)*((fdate - idate).days+1)

 #Create structured array
 water_use = {}
 for data_var in wbd['files_water_use']:
  water_use[data_var] = np.zeros((nt,hydrobloks_info['nhru']))

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
  dates = [datetime.datetime(int(nc_idate[0]),int(nc_idate[1]),1)]
  for it in range(0,nc_nt-1): dates.append( dates[it] + datetime.timedelta(hours=nc_step) )
  dates=np.array(dates)
  startdate = datetime.datetime(idate.year,idate.month,1,0)
  enddate = datetime.datetime(fdate.year,fdate.month,31,23)
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
 

  #Assing to hsus
  for hsu in mapping_info[var]:
   if OUTPUT['hsu']['land_cover'][hsu] in water_use_land_cover[data_var]:
    #print data_var,data, data.shape, hsu,mapping_info[var][hsu]['pcts'],mapping_info[var][hsu]['coords'],
    pcts = mapping_info[var][hsu]['pcts']
    coords = mapping_info[var][hsu]['coords']
    coords[0][coords[0] >= data.shape[1]] = data.shape[1] - 1
    coords[1][coords[1] >= data.shape[2]] = data.shape[2] - 1
    tmp = data[:,coords[0],coords[1]]
    tmp = pcts*tmp
    water_use[data_var][:,hsu] = np.sum(tmp,axis=1)  # final variable m3/m2/s --> m/s of water demand
    #print hsu, data_var, OUTPUT['hsu']['land_cover'][hsu], water_use[data_var][:,hsu]
   else:
    water_use[data_var][:,hsu] = 0.0

  #Write the water use the netcdf file (single chunk for now...)
  grp = hydrobloks_info['input_fp'].groups['water_use']
  grp.createVariable(var,'f4',('time','hsu'))#,zlib=True)
  grp.variables[data_var][:] = water_use[data_var][:]

 if hydrobloks_info['water_management']['hwu_flag']:
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

