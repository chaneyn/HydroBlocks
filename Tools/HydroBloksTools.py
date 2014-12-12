import warnings
warnings.filterwarnings('ignore')
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
import glob

def Deterministic(info):

 #Read in the catchment database
 wbd = pickle.load(open(info['wbd']))

 #Define the dates
 idate = datetime.datetime(2000,1,1,0)
 fdate = datetime.datetime(2000,1,31,23)

 #Iterate through all the catchments until done
 for icatch in [500,]:#[3637,]:#len(wbd.keys()):

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

 #Read all the catchments
 catchments = glob.glob('%s/*' % info['dir'])
 icatchs = []
 for dir in catchments:
  icatchs.append(int(dir.split('/')[-1].split('_')[-1]))
 #Shuffle the catchments
 icatchs = np.array(icatchs)
 np.random.seed(1)
 np.random.shuffle(icatchs)
 icatchs = [8072,3637,8756,500]#icatchs[0:1]#1000]#1000]#1000] #SUBSET
 #icatchs = icatchs[0:5000]#1000]#1000] #SUBSET

 #Define the dates
 idate = datetime.datetime(2004,1,1,0)
 #fdate = datetime.datetime(2004,1,31,23)
 fdate = datetime.datetime(2006,12,31,23)

 #Initialize the element count
 ielement = 0
 nens = 200#12#50#400
 elements = {}

 #Create a dictionary of information
 for icatch in icatchs:#xrange(800):#[501,502,503,504]:#[8756,500,3637]:#len(wbd.keys()):

  dir = info['dir']
  #Define the parameters
  parameters = {}
  parameters['log10m'] = -2.0#-2.582977995297425888e+00
  parameters['lnTe'] = -3.12#-1.963648774068431635e-01
  parameters['log10soil'] = np.log10(1.0)#1.389834359162560144e-02
  parameters['sdmax'] = 1.0#1.938762117265730334e+00

  #Cycle through the ensemble of clusters
  for iens in xrange(nens):
  
   #Define the number of bins
   #nclusters = int(np.linspace(2,1000,nens)[iens])#np.random.randint(1,1000)
   #nclusters = int(np.logspace(2,1000,nens)[iens])#np.random.randint(1,1000)
   nclusters = np.logspace(np.log(2),np.log(5000),nens,base=np.exp(1))
   np.random.seed(1)
   np.random.shuffle(nclusters)
   nclusters = np.int(np.ceil(nclusters[iens]))
   #nclusters = np.ceil(int(np.logspace(np.log(2),np.log(2000),nens,base=np.exp(1))[iens]))
   #nclusters = 200#250#100#250#25#2500

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
 #vars = ['lh','sh','smc1','prcp','qexcess','qsurface','swe']
 vars = ['smc1',]
 for var in vars:
  metrics['vars'][var] = {'mean':[],'std':[]}

 #Setup output redirection
 fout,ferr = os.open('/dev/null',os.O_RDWR|os.O_CREAT),os.open('/dev/null',os.O_RDWR|os.O_CREAT)
 so,se = os.dup(1),os.dup(2)

 #Randomize the data to minimize cores repeating the hard tasks
 #np.random.seed(1)
 ielements = np.arange(len(elements.keys()))
 #np.random.shuffle(ielements)

 #Iterate through the dictionary elements
 for ielement in ielements[rank::size]:

  #Define the info
  element = elements[ielement]

  #Print where we are at
  print 'Catchment %d, Ensemble %d' % (element['icatch'],element['iens']),element['nclusters']

  #Define the info
  hydrobloks_info = {
        'input':'%s/input/data.pck' % dir,
        'dt':3600.,
        'nsoil':20,
        'wbd':wbd[element['icatch']],
        'ncores':ncores,
        'idate':idate,
        'fdate':fdate,
        'parameters':parameters,
        'dir':'%s/catch_%d' % (dir,element['icatch']),
	'nclusters':element['nclusters']
        }

  try:
   #if 1 == 1:

   #Cluster the data
   time0 = time.time()
   input = Prepare_Model_Input_Data(hydrobloks_info)
   elements['nclusters'] = input['nclusters']
   dt = time.time() - time0
   print "time to prepare the data",element['nclusters'],dt
   #pickle.dump(input,open('data.pck','wb')) 
   #exit()

   #Flush out the output
   sys.stdout.flush()
   #sys.stderr.flush()

   #Redirect output
   #os.dup2(fout, 1),os.dup2(ferr,2)
   os.dup2(fout,1)

   #Run the model
   time0 = time.time()
   output = HB.run_model(hydrobloks_info,input,output_type='Summary')
   dt = time.time() - time0

   #Flush out the output
   sys.stdout.flush()
   #sys.stderr.flush()

   #Redirect the output back to the terminal 
   #os.dup2(so, 1),os.dup2(se,2)
   os.dup2(so,1)
   print 'time to run HydroBloks',element['nclusters'],dt

   #Compute heterogeneity metrics
   pcts = output['misc']['pct']
   nclusters = len(pcts)
   metrics['icatch'].append(element['icatch'])
   metrics['nclusters'].append(element['nclusters'])
   metrics['dt'].append(dt)
   for var in vars:
    metrics['vars'][var]['mean'].append(np.array(output['variables'][var]['mean']))
    metrics['vars'][var]['std'].append(np.array(output['variables'][var]['std']))

  except:
   #else:

   print "catchment %d Failed" % element['icatch']

 #Save time info and metrics to file
 #file = '/u/sciteam/nchaney/scratch/data/CONUS_SIMULATIONS_HUC10/output/%d.pck' % rank
 file = '/u/sciteam/nchaney/projects/HydroBloks/Output/output/%d.pck' % rank
 #file = '/u/sciteam/nchaney/scratch/data/CONUS_SIMULATIONS_HUC10/output/%d.pck' % rank
 pickle.dump(metrics,open(file,'wb'))

 return

def Latin_Hypercube_Sample(info):

 sys.path.append('Tools')
 from SALib.sample import latin_hypercube
 from SALib.util import scale_samples
 import random as rd
 import netCDF4 as nc

 #Define the rank and size
 rank = info['rank']
 size = info['size']
 ncores = info['ncores']
 nens = 2#10

 #Read in the catchment database
 wbd = pickle.load(open(info['wbd']))

 #Read all the catchments
 catchments = glob.glob('%s/*' % info['dir'])
 icatchs = []
 for dir in catchments:
  icatchs.append(int(dir.split('/')[-1].split('_')[-1]))
 #Shuffle the catchments
 icatchs = np.array(icatchs)
 np.random.seed(1)
 np.random.shuffle(icatchs)
 icatchs = [8072,3637,8756,500]#icatchs[0:1] #SUBSET
 #icatchs = [500,]#3637,]

 #Define the dates
 idate = datetime.datetime(2005,1,1,0)
 #fdate = datetime.datetime(2004,12,31,23)
 fdate = datetime.datetime(2007,12,31,23)

 #Obtain Latin Hypercube Sample 
 #1.Set random seed
 np.random.seed(1)
 rd.seed(1)

 #Define parameters and ranges
 parameters = []
 parameters.append(['log10m',np.log10(0.001),np.log10(0.1)])
 parameters.append(['lnTe',np.log(np.exp(-8.0)/3600.0),np.log(np.exp(8.0)/3600.0)])
 parameters.append(['log10soil',np.log10(1.0),np.log10(2.00)])
 parameters.append(['sdmax',0.1,2.0])

 #Generate samples (choose method here)
 param_values = latin_hypercube.sample(nens,len(parameters))
 bounds = []
 for iparam in xrange(len(parameters)):
  bounds.append([parameters[iparam][1],parameters[iparam][2]])
 scale_samples(param_values,np.array(bounds))

 #Initialize the element count
 ielement = 0
 elements = {}

 #Create a dictionary of information
 for icatch in icatchs:

  #Define the directory
  dir = info['dir']

  #Define the number of clusters (Change to catchment)
  nclusters = 250

  #Cycle through the ensemble of clusters
  for iens in xrange(nens):

   #Define the parameters
   parameters = {}
   parameters['log10m'] = param_values[iens,0]
   parameters['lnTe'] = param_values[iens,1]
   parameters['log10soil'] = param_values[iens,2]
   parameters['sdmax'] = param_values[iens,3]

   #Add the info to the dictionary
   elements[ielement] = {
		'parameters':parameters,
		'nclusters':nclusters,
		'icatch':icatch,
		'iens':iens,
                'cdir':'%s/catch_%d' % (dir,icatch)
		 } 

   #Update the element
   ielement += 1

 #Setup output redirection
 fout,ferr = os.open('/dev/null',os.O_RDWR|os.O_CREAT),os.open('/dev/null',os.O_RDWR|os.O_CREAT)
 so,se = os.dup(1),os.dup(2)

 #Iterate through the dictionary elements
 icatch = -9999
 for ielement in np.arange(len(elements.keys()))[rank::size]:

  #Define the info
  element = elements[ielement]

  #Print where we are at
  print 'Catchment %d, Ensemble %d' % (element['icatch'],element['iens']),element['nclusters']#,element['parameters']

  #Define the info
  hydrobloks_info = {
        'input':'%s/input/data.pck' % dir,
        'dt':3600.,
        'nsoil':20,
        'wbd':wbd[element['icatch']],
        'ncores':ncores,
        'idate':idate,
        'fdate':fdate,
        'parameters':element['parameters'],
        'dir':'%s/catch_%d' % (dir,element['icatch']),
        'nclusters':element['nclusters']
        }

  #If a new catchment then prepare the data
  if element['icatch'] != icatch:
   time0 = time.time()
   input = Prepare_Model_Input_Data(hydrobloks_info)
   elements['nclusters'] = input['nclusters']
   dt = time.time() - time0
   print "time to prepare the data",element['nclusters'],dt
   #pickle.dump(input,open('data.pck','wb')) 
   #exit()

  #Flush out the output
  sys.stdout.flush()
  #sys.stderr.flush()

  #Redirect output
  #os.dup2(fout, 1),os.dup2(ferr,2)
  #os.dup2(fout,1)

  #Run the model
  time0 = time.time()
  output = HB.run_model(hydrobloks_info,input,output_type='Full')
  dt = time.time() - time0
  print 'time to run HydroBloks',element['nclusters'],dt

  #Flush out the output
  sys.stdout.flush()
  #sys.stderr.flush()

  #Redirect the output back to the terminal 
  #os.dup2(so, 1),os.dup2(se,2)
  #os.dup2(so,1)

  #Save info to netcdf
  print 'Catchment number: %d, Ensemble number: %d (Writing to netcdf)' % (element['icatch'],element['iens'])
  file_netcdf = 'LHSoutput/%d.nc' % element['icatch']
  update_netcdf(element['cdir'],element['iens'],nens,parameters,file_netcdf,input,output)

  #Remember the catchment number
  icatch = element['icatch']

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

 #time0 = time.time()
 #Create the Latin Hypercube (Clustering)
 nclusters = hydrobloks_info['nclusters']
 ncores = hydrobloks_info['ncores']
 output = Create_Clusters_And_Connections(workspace,wbd,output,input_dir,nclusters,ncores)
 #print time.time() - time0

 #Extract the meteorological forcing
 time0 = time.time()
 output = Prepare_HSU_Meteorology(workspace,wbd,output,input_dir,info)
 #print time.time() - time0

 #Add in the catchment info
 output['wbd'] = wbd

 return output

def Create_Clusters_And_Connections(workspace,wbd,output,input_dir,nclusters,ncores):

 covariates = {}
 #Read in all the covariates
 for file in wbd['files']:
  #original = '/scratch/sciteam/nchaney/data/CONUS_SIMULATIONS_HUC10/catchments/catch_3637' #HERE
  #final = '/u/sciteam/nchaney/projects/HydroBloks/ReynoldsCreek' #HERE
  #wbd['files'][file] = wbd['files'][file].replace(original,final) #HERE
  covariates[file] = gdal_tools.read_raster(wbd['files'][file])
  #if file == 'carea': 
  # #covariates[file][covariates[file] == 0.0] = 1.0
  # #covariates[file] = np.log(covariates[file])
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
 #covariates['ti'][covariates['ti'] > 14] = 14
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
  covariates[var][mask <= 0] = -9999.0

 #Create channels mask
 mask_woc = np.copy(mask)
 #mask_woc[covariates['channels'] > 0] = 0
 mask_wc = np.copy(mask)
 mask_wc[covariates['channels'] <= 0] = 0

 #Define the covariates
 info = {'area':{'data':covariates['carea'][mask_woc == True],},
        'slope':{'data':covariates['cslope'][mask_woc == True],},
        'sms':{'data':covariates['MAXSMC'][mask_woc == True],},
        'ndvi':{'data':covariates['ndvi'][mask_woc ==True],},
        #'nlcd':{'data':covariates['nlcd'][mask_woc ==True],},
        'ti':{'data':covariates['ti'][mask_woc == True],},
        'dem':{'data':covariates['dem'][mask_woc == True],},
        'lats':{'data':covariates['lats'][mask_woc == True],},
        'lons':{'data':covariates['lons'][mask_woc == True],},
        }
 
 #Scale all the variables (Calculate the percentiles
 for var in info:
  #info[var]['data'] = (info[var]['data'] - np.min(info[var]['data']))/(np.max(info[var]['data']) - np.min(info[var]['data']))
  #if var == 'area':info[var]['data'] = 10*info[var]['data']
  argsort = np.argsort(info[var]['data'])
  pcts = np.copy(info[var]['data'])
  pcts[argsort] = np.linspace(0,1,len(info[var]['data']))
  #if np.unique(info[var]['data']).size < 1000:
  # for value in np.unique(info[var]['data']):
  #  pcts[pcts == value] = np.mean(pcts[pcts == value])
  info[var]['data'] = pcts

 #Create the LHS bins
 import sklearn.cluster
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
 minsamples = 2.5*10**4
 if X.shape[0] > minsamples:
  Xf = X[np.random.choice(np.arange(X.shape[0]),minsamples),:]
  #Make sure we have the extremes
  #for ivar in xrange(X.shape[1]):
  # imin = np.argmin(X[:,ivar])
  # imax = np.argmax(X[:,ivar])
  # Xf = np.append(Xf,X[imin,:][np.newaxis,:],axis=0)
  # Xf = np.append(Xf,X[imax,:][np.newaxis,:],axis=0)
 else:
  Xf = X
 #Construct grid to fit the data
 #npoints = np.ceil(minsamples**(1./len(info.keys())))
 #if npoints < 2: npoints = 2
 #for id in xrange(len(info.keys())):
 # if id == 0: string = "np.meshgrid(np.linspace(0,1,%d)" % npoints 
 # else: string = string + ",np.linspace(0,1,%d)" % npoints
 #string = string + ")"
 #tmp = eval(string)
 #Xf = []
 #Reshape data for the fitting
 #for id in xrange(len(info.keys())):
 # Xf.append(np.reshape(tmp[id],tmp[id].size))
 #Xf = np.array(Xf).T
 #print Xf
 #idx = eval('np.where(%s)' % string)
 #Xf = np.mgrid(
 #Construct a features array (0,1)
 #Xf = np.random.rand(minsamples,X.shape[1])
 #Initialize all points at the 0.5 point
 init = 0.5*np.ones((nclusters,Xf.shape[1]))
 #clf = sklearn.cluster.KMeans(nclusters,n_jobs=ncores,n_init=1,init=init,tol=1e-4,max_iter=300)
 batch_size = 25*nclusters
 init_size = 3*batch_size
 #clf = sklearn.cluster.MiniBatchKMeans(nclusters,random_state=1,init=init,batch_size=batch_size,init_size=init_size)
 clf = sklearn.cluster.MiniBatchKMeans(nclusters,init=init,batch_size=batch_size,init_size=init_size)
 #clf = sklearn.cluster.AgglomerativeClustering(nclusters)
 #clf = sklearn.cluster.DBSCAN(eps=0.3, min_samples=10)#clusters)
 clf.fit(Xf)#
 clf_output = clf.predict(X)
 #Reassign the ids
 clf_output_copy = np.copy(clf_output)
 for cid in xrange(len(np.unique(clf_output))):
  clf_output[clf_output_copy == np.unique(clf_output)[cid]] = cid
 cluster_ids = np.empty(covariates['ti'].shape)
 cluster_ids[:] = -9999
 cluster_ids[mask_woc == True] = clf_output
 nclusters_old = nclusters
 nclusters = np.unique(clf_output).size
 #Redefine the number of clusters (We are sampling regions that just don't have data...)
 print 'clustering %d->%d' % (nclusters_old,nclusters),time.time() - time0

 #Add in the channel clusters
 #channels = np.unique(covariates['channels'][covariates['channels'] > 0])
 #for channel in channels:
 # cid = int(np.nanmax(cluster_ids) + 1)
 # idx = np.where(covariates['channels'] == channel)
 # cluster_ids[idx] = cid
 # nclusters = nclusters + 1

 #Reorder according to areas
 areas = []
 for cid in xrange(nclusters):
  msk = cluster_ids == cid
  areas.append(np.nanmean(covariates['carea'][msk]))
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

 #Create a dictionary of class info
 clusters = {}
 areas = []
 for cid in xrange(nclusters):
  #Determine the percentage coverage
  pct = float(np.sum(cluster_ids == cid))/float(np.sum(mask))
  clusters[cid] = {'pct':pct}
  idx = np.where(cluster_ids == cid)
  clusters[cid]['idx'] = idx

 #Determine the links between clusters
 mask1 = covariates['fdir'] < 0
 covariates['fdir'][mask1] = -9999.0
 cluster_ids_copy = np.copy(cluster_ids)
 cluster_ids_copy[mask1] = np.nan
 nclusters = len(clusters.keys())
 tp_matrix = mt.preprocessor.calculate_connections_d8(cluster_ids_copy,covariates['fdir'],nclusters)

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

 #Create the image of the covariates
 '''covariates['hsu_map'] = OUTPUT['hsu_map']
 vars= ['carea','cslope','MAXSMC','nlcd','ti','dem','hsu_map']
 ns = int(np.ceil(float(len(vars))**0.5))
 plt.figure(figsize=(25,25))
 for var in vars:
  #Apply the mask
  plotdata = np.copy(covariates[var])
  plotdata[mask <= 0] = np.nan
  #Plot the data
  plt.subplot(ns,ns,vars.index(var) + 1)
  plt.title(var,fontsize=45)
  plt.imshow(plotdata,interpolation='nearest')
  plt.axis('off')

 #Save the figure
 plt.tight_layout()
 figure = 'covariates.png'
 plt.savefig(figure)'''

 #Add the new number of clusters
 OUTPUT['nclusters'] = nclusters

 return OUTPUT

def Prepare_HSU_Meteorology(workspace,wbd,OUTPUT,input_dir,info):

 #Define the mapping directory
 mapping_dir = '%s/mapping' % workspace
 #Calculate the fine to coarse scale mapping
 for data_var in wbd['files_meteorology']:
  
  #Define the variable name
  var = data_var#data_var.split('_')[1]

  #Read in the coarse and fine mapping
  file_coarse = '%s/%s_coarse.tif' % (mapping_dir,data_var)
  file_fine = '%s/%s_fine.tif' % (mapping_dir,data_var)
  mask_coarse = gdal_tools.read_raster(file_coarse)
  mask_fine = gdal_tools.read_raster(file_fine)
  nlat = mask_coarse.shape[0]
  nlon = mask_coarse.shape[1]

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
  var = data_var#data_var.split('_')[1]
  date = idate
  file = wbd['files_meteorology'][data_var]
  #dir = ctl[0:-(len(var)+5)]
  #file = '%s/%s/%s.nc' % (dir,var,var)
  #if var == 'apcpsfc': file = '%s/%s.nc' % (dir,var)
  fp = nc.Dataset(file)
  #Determine the time steps to retrieve
  dates = nc.num2date(fp.variables['t'][:],units='hours since %02d-%02d-%02d 00:00:00' % (idate.year,idate.month,idate.day))
  mask_dates = (dates >= idate) & (dates <= fdate)
  data = np.ma.getdata(fp.variables[var][mask_dates])
  fp.close()
  #Assing to hsus
  for hsu in OUTPUT['hsu']:
   pcts = OUTPUT['hsu'][hsu][var]['pcts']
   coords = OUTPUT['hsu'][hsu][var]['coords']
   coords[0][coords[0] >= data.shape[1]] = data.shape[1] - 1
   coords[1][coords[1] >= data.shape[2]] = data.shape[2] - 1
   tmp = pcts*data[:,coords[0],coords[1]]
   #Combine stage iv and nldas here
   if data_var not in ['apcpsfc',]:tmp[tmp < -999] = np.mean(tmp[tmp > -999])
   meteorology[data_var][:,hsu] = np.sum(tmp,axis=1)

 #Append the meteorology to the output dictionary
 OUTPUT['meteorology'] = meteorology

 return OUTPUT

def create_netcdf_file(file_netcdf,output,nens,input,cdir):

 #Create the file
 fp = nc.Dataset(file_netcdf, 'w', format='NETCDF4')

 #Create the dimensions
 ntime = output['variables']['lh'].shape[0]/24
 nhsu = output['variables']['lh'].shape[1]
 nsoil = 4
 fp.createDimension('hsu',nhsu)
 fp.createDimension('time',ntime)
 fp.createDimension('ensemble',nens)
 fp.createDimension('soil',nsoil)

 #Create the summary group and its variables
 grp = fp.createGroup('summary')
 for var in output['variables']:
  ncvar = grp.createVariable(var,'f4',('time','ensemble'))
  ncvar.description = output['misc']['metadata'][var]['description']
  ncvar.units = output['misc']['metadata'][var]['units']

 #Create the output
 grp = fp.createGroup('catchment')
 for var in output['variables']:
  grp.createVariable(var,'f4',('time','hsu','ensemble'))
  ncvar.description = output['misc']['metadata'][var]['description']
  ncvar.units = output['misc']['metadata'][var]['units']

 #Create the metadata
 grp = fp.createGroup('metadata')
 #dates
 times = grp.createVariable('dates','f8',('time',))
 dates = []
 for date in output['misc']['dates']:
  if date.hour == 0:dates.append(date)
 times.units = 'days since 1900-01-01'
 times.calendar = 'standard'
 times[:] = nc.date2num(np.array(dates),units=times.units,calendar=times.calendar)
 #Soil depth
 soil = grp.createVariable('soil','f4',('soil',))
 soil[:] = np.array(output['misc']['sldpth'])[0,0:nsoil]
 soil.units = 'meters'
 soil.description = 'soil depth'
 #HSU percentage coverage
 pcts = grp.createVariable('pct','f4',('hsu',))
 pcts[:] = 100*np.array(output['misc']['pct'])
 pcts.description = 'hsu percentage coverage'
 pcts.units = '%'
 #HSU Spatial resolution
 dx = grp.createVariable('dx','f4',('hsu',))
 dx[:] = np.array(output['misc']['dx'])
 dx.units = 'meters'
 #HSU area
 area = grp.createVariable('area','f4',('hsu',))
 area[:] = np.array(output['misc']['area'])
 area.units = 'meters squared'
 #HSU ids
 hsu = grp.createVariable('hsu','i4',('hsu',))
 hsus =[]
 for value in xrange(len(output['misc']['pct'])):hsus.append(value)
 hsu[:] = np.array(hsus)
 hsu.description = 'hsu ids'
 #Add wbd metadata
 grp.HUC = input['wbd']['HUC10']
 #Define outlet hsu
 grp.outlet_hsu = int(output['misc']['outlet_hsu'])
 #Define catchment name
 grp.catchment_name = input['wbd']['Name']
 #Define catchment area
 grp.AreaSqKm = input['wbd']['AreaSqKm']
 #Add HSU transition probability matrix
 tp = grp.createVariable('tpm','f4',('hsu','hsu'))
 tp.description = 'Transition probability matrix between hsus'
 tp[:] =  input['tp']

 #Retrieve the conus_albers metadata
 metadata = gdal_tools.retrieve_metadata(input['wbd']['files']['ti']) 
 metadata['nodata'] = -9999.0
 #Save the conus_albers metadata
 fp.createDimension('nx',metadata['nx'])
 fp.createDimension('ny',metadata['ny'])
 hmca = grp.createVariable('hmca','f4',('ny','nx')) 
 hmca.gt = metadata['gt']
 hmca.projection = metadata['projection']
 hmca.description = 'HSU mapping (conus albers)'
 hmca.nodata = metadata['nodata']
 #Save the conus albers mapping
 hsu_map = np.copy(input['hsu_map'])
 hsu_map[np.isnan(hsu_map) == 1] = metadata['nodata']
 hmca[:] = hsu_map
 #Write out the mapping
 file_ca = '%s/workspace/hsu_mapping_conus_albers.tif' % cdir
 gdal_tools.write_raster(file_ca,metadata,hsu_map)

 #Map the mapping to regular lat/lon
 file_ll = '%s/workspace/hsu_mapping_latlon.tif' % cdir
 os.system('rm -f %s' % file_ll)
 res = input['wbd']['bbox']['res']
 minlat = input['wbd']['bbox']['minlat']
 minlon = input['wbd']['bbox']['minlon']
 maxlat = input['wbd']['bbox']['maxlat']
 maxlon = input['wbd']['bbox']['maxlon']
 log = '%s/workspace/log.txt' % cdir
 os.system('gdalwarp -tr %.16f %.16f -dstnodata %.16f -t_srs EPSG:4326 -s_srs EPSG:102039 -te %.16f %.16f %.16f %.16f %s %s >> %s 2>&1' % (res,res,metadata['nodata'],minlon,minlat,maxlon,maxlat,file_ca,file_ll,log))
 #Add the lat/lon mapping
 #Retrieve the lat/lon metadata
 metadata = gdal_tools.retrieve_metadata(file_ll)
 metadata['nodata'] = -9999.0
 #Save the lat/lon metadata
 fp.createDimension('nlon',metadata['nx'])
 fp.createDimension('nlat',metadata['ny'])
 hmll = grp.createVariable('hmll','f4',('nlat','nlon'))
 hmll.gt = metadata['gt']
 hmll.projection = metadata['projection']
 hmll.description = 'HSU mapping (regular lat/lon)'
 hmll.nodata = metadata['nodata']
 #Save the lat/lon mapping
 hsu_map = np.copy(gdal_tools.read_raster(file_ll))
 hsu_map[np.isnan(hsu_map) == 1] = metadata['nodata']
 hmll[:] = hsu_map
 
 #Close the file 
 fp.close()

 return

def update_netcdf(cdir,iens,nens,parameters,file_netcdf,input,output):
 
 #Convert variables to arrays
 for var in output['variables']:
  output['variables'][var] = np.array(output['variables'][var])

 #Create the netcdf file if necessary
 if iens == 0: create_netcdf_file(file_netcdf,output,nens,input,cdir)

 #Open the file
 fp = nc.Dataset(file_netcdf, 'a')
 
 #Output the data
 for var in output['variables']:
  
  #Compute the daily average
  data = []
  for itime in xrange(output['variables'][var].shape[0]/24):
   data.append(np.mean(output['variables'][var][24*itime:24*(itime+1)],axis=0))
  data = np.array(data)

  #Write the catchment info
  fp.groups['catchment'].variables[var][:,:,iens] = data

  #Compute the catchment summary 
  if var not in ['qout_surface','qout_subsurface']:
   data = np.sum(output['misc']['pct']*data,axis=1)
  else:
   imax = output['misc']['outlet_hsu']
   data = data[:,imax]

  #Write the catchment summary
  fp.groups['summary'].variables[var][:,iens] = data

 #Close the file
 fp.close()

 return
