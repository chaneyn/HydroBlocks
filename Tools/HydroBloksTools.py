import sys
sys.path.append('Model')
import HydroBloks as HB
import matplotlib
matplotlib.use('Agg')
import cPickle as pickle
import datetime
import gdal_tools
import grads_tools
import numpy as np
import scipy.stats as stats
import model_tools as mt
import matplotlib.pyplot as plt
import os
import netCDF4 as nc
import time

def Deterministic(info):

 #Read in the catchment database
 wbd = pickle.load(open(info['wbd']))

 #Define the dates
 idate = datetime.datetime(2000,1,1,0)
 fdate = datetime.datetime(2000,1,31,23)

 #Iterate through all the catchments until done
 for icatch in [3637,]:#len(wbd.keys()):

  dir = info['dir']
  #Define the parameters
  parameters = {}
  parameters['log10m'] = -2.582977995297425888e+00
  parameters['lnTe'] = -1.963648774068431635e-01
  parameters['log10soil'] = 1.389834359162560144e-02
  parameters['sdmax'] = 1.938762117265730334e+00

  #Define the info
  hydrobloks_info = {
        'input':'%s/input/data.pck' % dir,
        'dt':3600.,
        'nsoil':20,
        'wbd':wbd[icatch],
        'ncores':1,
        'idate':idate,
        'fdate':fdate,
        'parameters':parameters,
	'dir':dir,
        'nbins':{'area':10,'slope':1,'sms':1,'ndvi':1,'ti':1,'dem':1,'channels':2}
        }

 #Cluster the data
 Prepare_Model_Input_Data(hydrobloks_info)

 #Run the model
 output = HB.run_model(hydrobloks_info)

 return

def Convergence_Analysis(info):

 #Define the rank and size
 rank = info['rank']
 size = info['size']
 ncores = info['ncores']

 #Read in the catchment database
 wbd = pickle.load(open(info['wbd']))

 #Define the dates
 idate = datetime.datetime(2000,1,1,0)
 fdate = datetime.datetime(2000,1,31,23)

 #Initialize the element count
 ielement = 0
 nens = 100
 elements = {}

 #Create a dictionary of information
 for icatch in [8756,500,3637]:#len(wbd.keys()):

  dir = info['dir']
  #Define the parameters
  parameters = {}
  parameters['log10m'] = -2.582977995297425888e+00
  parameters['lnTe'] = -1.963648774068431635e-01
  parameters['log10soil'] = 1.389834359162560144e-02
  parameters['sdmax'] = 1.938762117265730334e+00

  #Cycle through the ensemble of clusters
  for iens in xrange(nens):
  
   #Define the number of bins
   nclusters = np.random.randint(1,2500)

   #Add the info to the dictionary
   elements[ielement] = {
		'parameters':parameters,
		'nclusters':nclusters,
		'icatch':icatch,
		'iens':iens,
		 } 

   #Update the element
   ielement += 1

  #Initialize metrics dictionary
  metrics = {'icatch':[],'dt':[],'nclusters':[],'vars':{}}

  #Add the output variables
  vars = ['lh','sh','smc1','prcp','qexcess','qsurface','swe']
  for var in vars:
   metrics['vars'][var] = {'mean':[],'std':[]}

  #Iterate through the dictionary elements
  for ielement in np.arange(len(elements.keys()))[rank::size]:

   #Define the info
   element = elements[ielement]

   #Print where we are at
   print 'Catchment %d, Ensemble %d' % (element['icatch'],element['iens']),element['nclusters']

   #Define the info
   hydrobloks_info = {
        'input':'%s/input/data.pck' % dir,
        'dt':3600.,
        'nsoil':20,
        'wbd':wbd[icatch],
        'ncores':ncores,
        'idate':idate,
        'fdate':fdate,
        'parameters':parameters,
        'dir':dir,
	'nclusters':element['nclusters']
        }

   #Cluster the data
   #time0 = time.time()
   input = Prepare_Model_Input_Data(hydrobloks_info)
   #dt = time.time() - time0
   #nclusters = 1
   #for var in nbins:
   # nclusters = nclusters*nbins[var]
   #pickle.dump(input,open('data.pck','wb')) 
   #exit()

   #Run the model
   time0 = time.time()
   output = HB.run_model(hydrobloks_info,input)
   dt = time.time() - time0

   #Compute heterogeneity metrics
   pcts = output['misc']['pct']
   nclusters = len(pcts)
   metrics['icatch'].append(icatch)
   metrics['nclusters'].append(nclusters)
   metrics['dt'].append(dt)
   for var in vars:
    output['variables'][var] = np.array(output['variables'][var])
    #Compute the mean
    mean = np.sum(pcts*output['variables'][var],axis=1)
    #Compute the standard deviation
    std = np.sum(pcts*(output['variables'][var] - mean[:,np.newaxis])**2,axis=1)**0.5
    #Save the time series for the mean and standard deviation
    metrics['vars'][var]['mean'].append(mean)
    metrics['vars'][var]['std'].append(std)
    #mean
    #percentiles = []
    #for percentile in [1,10,25,50,75,90,99]:
    # percentiles.append(np.percentile(mean,percentile))
    #metrics['vars'][var]['mean'].append(percentiles)
    #std
    #percentiles = []
    #for percentile in [1,10,25,50,75,90,99]:
    # percentiles.append(np.percentile(std,percentile))
    #metrics['vars'][var]['std'].append(percentiles)

   #Save time info and metrics to file
   file = 'Output/%d.pck' % rank
   pickle.dump(metrics,open(file,'wb'))

 return

def Latin_Hypercube_Sample(info):

 print info

 return

def Prepare_Model_Input_Data(hydrobloks_info):

 #Prepare the info dictionary
 info = {}

 #Define the binaries
 info['binaries'] = {}
 info['binaries']['grads'] = '/u/sciteam/nchaney/libraries/opengrads/Contents/grads'
 info['binaries']['gdalwarp'] = '/u/sciteam/nchaney/local/bin/gdalwarp'

 #Define the start/end dates
 info['time_info'] = {}
 info['time_info']['startdate'] = hydrobloks_info['idate']
 info['time_info']['enddate'] = hydrobloks_info['fdate']

 #Define the workspace
 workspace = '%s/workspace' % hydrobloks_info['dir']

 #Define the model input data directory
 input_dir = '%s/input' % hydrobloks_info['dir']

 #Read in the metadata
 file = '%s/workspace_info.pck' % workspace
 wbd = pickle.load(open(file))

 #Create the dictionary to hold all of the data
 output = {}

 #Create the Latin Hypercube (Clustering)
 nclusters = hydrobloks_info['nclusters']
 ncores = hydrobloks_info['ncores']
 output = Create_Clusters_And_Connections(workspace,wbd,output,input_dir,nclusters,ncores)

 #Extract the meteorological forcing
 output = Prepare_HSU_Meteorology(workspace,wbd,output,input_dir,info)

 #Add in the catchment info
 output['wbd'] = wbd

 return output

def Create_Clusters_And_Connections(workspace,wbd,output,input_dir,nclusters,ncores):

 covariates = {}
 #Read in all the covariates
 for file in wbd['files']:
  original = '/scratch/sciteam/nchaney/data/CONUS_SIMULATIONS_HUC10/catchments/catch_3637' #HERE
  final = '/u/sciteam/nchaney/projects/HydroBloks/ReynoldsCreek' #HERE
  wbd['files'][file] = wbd['files'][file].replace(original,final) #HERE
  covariates[file] = gdal_tools.read_raster(wbd['files'][file])
  if file == 'carea': covariates[file] = np.log(covariates[file])
  if file == 'cslope':
   mask = covariates[file] == 0.0
   covariates[file][mask] = 0.000001

 #Create lat/lon grids
 lats = np.linspace(wbd['bbox']['minlat']+wbd['bbox']['res']/2,wbd['bbox']['maxlat']-wbd['bbox']['res']/2,covariates['ti'].shape[0])
 lons = np.linspace(wbd['bbox']['minlon']+wbd['bbox']['res']/2,wbd['bbox']['maxlon']-wbd['bbox']['res']/2,covariates['ti'].shape[1])
 lats, lons = np.meshgrid(lats, lons)
 covariates['lats'] = lats.T
 covariates['lons'] = lons.T

 #Clean up the covariates
 covariates['ti'][covariates['ti'] > 14] = 14
 #covariates['channels'][covariates['channels'] > 1] = 1
 #covariates['channels'][covariates['channels'] < 0] = 0

 #Define the mask
 mask = np.copy(covariates['mask'])
 mask[mask > 0] = 1
 mask[mask < 0] = 0
 #channels_original[mask == 0] = -9999
 #mask[covariates['channels'] >= 1] = 0
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
  covariates[var][mask <= 0] = -9999

 #Create channels mask
 mask_woc = np.copy(mask)
 mask_woc[covariates['channels'] > 0] = 0
 mask_wc = np.copy(mask)
 mask_wc[covariates['channels'] < 0] = 0

 #Define the covariates
 info = {'area':{'data':covariates['carea'][mask_woc == True],},
        #'slope':{'data':covariates['cslope'][mask_woc == True],},
        'sms':{'data':covariates['MAXSMC'][mask_woc == True],},
        'ndvi':{'data':covariates['ndvi'][mask_woc ==True],},
        'ti':{'data':covariates['ti'][mask_woc == True],},
        'dem':{'data':covariates['dem'][mask_woc == True],},
        'lats':{'data':covariates['lats'][mask_woc == True],},
        'lons':{'data':covariates['lons'][mask_woc == True],},
        }
 
 #Create the LHS bins
 import sklearn.cluster
 bins,data = [],[]
 X = []
 for id in info:
  #Set all nans to the mean
  info[id]['data'][np.isnan(info[id]['data']) == 1] = np.nanmean(info[id]['data'])
  X.append(info[id]['data'])

 X = np.array(X).T
 #Subsample the array
 Xf = X[np.random.choice(np.arange(X.shape[0]),10000),:]
 clf = sklearn.cluster.KMeans(nclusters,n_jobs=ncores)
 clf.fit(Xf)
 clf_output = clf.predict(X)
 cluster_ids = np.empty(covariates['ti'].shape)
 cluster_ids[:] = -9999
 cluster_ids[mask_woc == True] = clf_output

 #Create a dictionary of class info
 clusters = {}
 #Hfl = H.flat
 for cid in xrange(nclusters):
  #Determine the percentage coverage
  pct = float(np.sum(cluster_ids == cid))/float(np.sum(mask))
  clusters[cid] = {'pct':pct}
  idx = np.where(cluster_ids == cid)
  clusters[cid]['idx'] = idx
  
 #Add in the channel clusters
 channels = np.unique(covariates['channels'][covariates['channels'] > 0])
 for channel in channels:
  cid = int(np.nanmax(cluster_ids) + 1)
  pct = float(np.sum(covariates['channels'] == channel))/float(np.sum(mask))
  clusters[cid] = {'pct':pct}
  idx = np.where(covariates['channels'] == channel)
  clusters[cid]['idx'] = idx
  cluster_ids[idx] = cid

 #Determine the links between clusters
 mask1 = covariates['fdir'] < 0
 covariates['fdir'][mask1] = -9999.0
 cluster_ids[mask1] = np.nan
 nclusters = len(clusters.keys())
 tp_matrix = mt.preprocessor.calculate_connections_d8(cluster_ids,covariates['fdir'],nclusters)

 #Define the metadata
 metadata = gdal_tools.retrieve_metadata(wbd['files']['ti'])

 #Make the output dictionary for the basin
 OUTPUT = {'hsu':{},'tp':tp_matrix,'metadata':metadata,'mask':mask}

 #Determine outlet cell
 covariates['carea'][mask == False] = np.nan
 outlet_idx = covariates['carea'] == np.max(covariates['carea'][np.isnan(covariates['carea']) == 0])
 OUTPUT['outlet'] = {'idx':outlet_idx,'hsu':cluster_ids[outlet_idx]}
 OUTPUT['hsu_map'] = cluster_ids

 for hsu in clusters:

  #Set indices
  OUTPUT['hsu'][hsu] = {'idx':clusters[hsu]['idx']}
  #Calculate area per hsu
  OUTPUT['hsu'][hsu]['area'] = metadata['resx']*OUTPUT['hsu'][hsu]['idx'][0].size
  #Calculate area percentage per hsu
  OUTPUT['hsu'][hsu]['area_pct'] = 100*OUTPUT['hsu'][hsu]['area']/(metadata['resx']*mask[mask].size)
  #SOIL
  OUTPUT['hsu'][hsu]['soil_parameters'] = {}
  for var in ['BB','DRYSMC','F11','MAXSMC','REFSMC','SATPSI','SATDK','SATDW','WLTSMC','QTZ']:
   OUTPUT['hsu'][hsu]['soil_parameters'][var] = np.mean(covariates[var][OUTPUT['hsu'][hsu]['idx']])
  #Average Slope
  OUTPUT['hsu'][hsu]['slope'] = np.mean(covariates['cslope'][OUTPUT['hsu'][hsu]['idx']])
  #Topographic index
  OUTPUT['hsu'][hsu]['ti'] = np.mean(covariates['ti'][OUTPUT['hsu'][hsu]['idx']])
  #DEM
  OUTPUT['hsu'][hsu]['dem'] = np.mean(covariates['dem'][OUTPUT['hsu'][hsu]['idx']])
  #Average Catchment Area
  OUTPUT['hsu'][hsu]['carea'] = np.mean(covariates['carea'][OUTPUT['hsu'][hsu]['idx']])
  #Channel?
  OUTPUT['hsu'][hsu]['channel'] = np.mean(covariates['channels'][OUTPUT['hsu'][hsu]['idx']])
  #Vchan
  OUTPUT['hsu'][hsu]['vchan'] = 1000 #m/hr
  #Vof
  OUTPUT['hsu'][hsu]['vof'] = 100 #m/hr

 #Vegetation
 NLCD2NOAH = {11:17,12:15,21:10,22:10,23:10,24:13,31:16,41:4,42:1,43:5,51:6,52:6,71:10,72:10,73:19,74:19,81:10,82:12,90:11,95:11}
 #Determine the most frequent vegetation class per hsu and assign type
 for hsu in OUTPUT['hsu']:
  idx = OUTPUT['hsu'][hsu]['idx']
  #if stats.mode(covariates['nlcd'][idx])[0][0] == 0.0:return
  OUTPUT['hsu'][hsu]['land_cover'] = NLCD2NOAH[stats.mode(covariates['nlcd'][idx])[0][0]]

 #Soil
 #Determine the most frequent soil texture class per hsu and assign type
 for hsu in OUTPUT['hsu']:
  idx = OUTPUT['hsu'][hsu]['idx']
  OUTPUT['hsu'][hsu]['soil_texture_class'] = stats.mode(covariates['TEXTURE_CLASS'][idx])[0][0]

 return OUTPUT

def Prepare_HSU_Meteorology(workspace,wbd,OUTPUT,input_dir,info):

 #Define the mapping directory
 mapping_dir = '%s/mapping' % workspace
 #Calculate the fine to coarse scale mapping
 for data_var in wbd['files_meteorology']:
  
  #Define the variable name
  var = data_var.split('_')[1]

  #Read in the coarse and fine mapping
  file_coarse = '%s/%s_coarse.tif' % (mapping_dir,data_var)
  file_fine = '%s/%s_fine.tif' % (mapping_dir,data_var)
  mask_coarse = gdal_tools.read_raster(file_coarse)
  mask_fine = gdal_tools.read_raster(file_fine)

  #Compute the mapping for each hsu
  for hsu in OUTPUT['hsu']:
   idx = OUTPUT['hsu'][hsu]['idx']
   icells = np.unique(mask_fine[idx].astype(np.int))
   counts = np.bincount(mask_fine[idx].astype(np.int))
   coords,pcts = [],[]
   for icell in icells:
    ilat = int(np.floor(icell/mask_coarse.shape[1]))
    jlat = icell - ilat*mask_coarse.shape[1]
    pct = float(counts[icell])/float(np.sum(counts))
    coords.append([ilat,jlat])
    pcts.append(pct)
   pcts = np.array(pcts)
   coords = list(np.array(coords).T)
   OUTPUT['hsu'][hsu][var] = {'pcts':pcts,'coords':coords}

 #Iterate through variable creating forcing product per HSU
 idate = info['time_info']['startdate']
 fdate = info['time_info']['enddate']
 nt = 24*((fdate - idate).days+1)
 #Create structured array
 meteorology = {}
 formats,names = [],[]
 for name in OUTPUT['hsu'].keys():
  formats.append('f8')
  names.append(str(name))
 for data_var in wbd['files_meteorology']:
  meteorology[data_var] = np.zeros((nt,len(names)))
 #Load data into structured array
 for data_var in wbd['files_meteorology']:
  var = data_var.split('_')[1]
  date = idate
  ctl = wbd['files_meteorology'][data_var]
  dir = ctl[0:-(len(var)+5)]
  file = '%s/%s/%s.nc' % (dir,var,var)
  fp = nc.Dataset(file)
  #Determine the time steps to retrieve
  dates = nc.num2date(fp.variables['t'][:],units='hours since 2000-01-01 00:00:00')
  mask_dates = (dates >= idate) & (dates <= fdate)
  data = np.ma.getdata(fp.variables[var][mask_dates])
  fp.close()
  #Assing to hsus
  for hsu in OUTPUT['hsu']:
   pcts = OUTPUT['hsu'][hsu][var]['pcts']
   coords = OUTPUT['hsu'][hsu][var]['coords']
   tmp = pcts*data[:,coords[0],coords[1]]
   #Combine stage iv and nldas here
   if data_var not in ['stageiv_prec',]:tmp[tmp < -999] = np.mean(tmp[tmp > -999])
   meteorology[data_var][:,hsu] = np.sum(tmp,axis=1)

 #Append the meteorology to the output dictionary
 OUTPUT['meteorology'] = meteorology

 return OUTPUT
