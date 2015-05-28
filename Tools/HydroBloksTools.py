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
import scipy.sparse as sparse
import scipy.stats as stats
import model_tools as mt
import matplotlib.pyplot as plt
import os
import netCDF4 as nc
import time
import glob

def Deterministic(info):

 #Define the metadata
 nclusters = 3
 rank = 0

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
        'icatch':icatch,
        'rank':rank,
        #'input':'%s/input/data.pck' % dir,
        'input':'%s/catch_%d/input_data.nc' % (dir,icatch),
        'dt':3600.,
        'nsoil':20,
        'wbd':wbd[icatch],
        'ncores':1,
        'idate':idate,
        'fdate':fdate,
        'parameters':parameters,
        'dir':'%s/catch_%d' % (dir,icatch),
        'nclusters':nclusters,
        'model_type':'semi',
        'output_type':'Full',
        'soil_file':'%s/catch_%d/workspace/soils/SOILPARM_%d_%d.TBL' % (dir,icatch,icatch,rank),
        'output':'%s/catch_%d/output_data.nc' % (dir,icatch),
        }

 #Cluster the data
 input = Prepare_Model_Input_Data(hydrobloks_info)
 #pickle.dump(input,open('tmp.pck','wb'),pickle.HIGHEST_PROTOCOL)
 #input = pickle.load(open('workspace/tmp.pck'))

 #Run the model
 #output = HB.run_model(hydrobloks_info,input,output_type='Full')
 output = HB.run_model(hydrobloks_info)

 #Write to netcdf
 file_netcdf = hydrobloks_info['output']
 #Update the netcdf file
 vars = output['variables'].keys()
 update_netcdf(hydrobloks_info['dir'],0,1,parameters,file_netcdf,output,0)

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
 #MISSING CATCHMENTS REPLACE!!!!!!!!!!!!
 #icatchs = pickle.load(open('/u/sciteam/nchaney/projects/HydroBloks/Output/miscellanous/screened.pck'))
 np.random.seed(1)
 np.random.shuffle(icatchs)
 icatchs = [8072,3637,8756,500]#icatchs[0:1]#1000]#1000]#1000] #SUBSET
 #icatchs = icatchs[:]#[0:5000]#1000]#1000] #SUBSET

 #Define the dates
 idate = datetime.datetime(2004,1,1,0)
 #fdate = datetime.datetime(2004,1,31,23)
 fdate = datetime.datetime(2007,12,31,23)

 #Initialize the element count
 ielement = 0
 nens = 1000#12#50#400
 elements = {}

 #Create a dictionary of information
 for icatch in icatchs:#xrange(800):#[501,502,503,504]:#[8756,500,3637]:#len(wbd.keys()):

  dir = info['dir']
  #Define the parameters
  parameters = {'log10m': -1.279675506557842, 'log10soil': 0.17854995314368682, 'sdmax': 0.73625873496330596, 'lnTe': -15.163693541009348}
  #parameters = {}
  #parameters['log10m'] = -2.0#-2.582977995297425888e+00
  #parameters['lnTe'] = -3.12#-1.963648774068431635e-01
  #parameters['log10soil'] = np.log10(1.0)#1.389834359162560144e-02
  #parameters['sdmax'] = 1.0#1.938762117265730334e+00

  #Cycle through the ensemble of clusters
  for iens in xrange(nens):
  
   #Define the number of bins
   #nclusters = int(np.linspace(2,1000,nens)[iens])#np.random.randint(1,1000)
   #nclusters = int(np.logspace(2,1000,nens)[iens])#np.random.randint(1,1000)
   nclusters = np.logspace(np.log(2),np.log(2000),nens,base=np.exp(1))
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
 vars = ['lh','sh','smc1','prcp','qexcess','qsurface','swe']
 #vars = ['smc1',]
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
 #file = '/u/sciteam/nchaney/scratch/data/CONUS_SIMULATIONS_HUC10/missing/%d.pck' % rank
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
 nens = 250

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
 clusters_info = {8072:469,3637:794,8756:379,500:393}
 #clusters_info = {8072:100,3637:100,8756:100,500:100}##393}
 #icatchs = [3637,]#3637,]

 #Define the dates
 idate = datetime.datetime(2004,1,1,0)
 #fdate = datetime.datetime(2005,12,31,23)
 #fdate = datetime.datetime(2006,12,31,23)
 fdate = datetime.datetime(2009,12,31,23)

 #Obtain Latin Hypercube Sample 
 #1.Set random seed
 np.random.seed(1)
 rd.seed(1)

 #Define parameters and ranges
 parameters = []
 parameters.append(['log10m',np.log10(0.001),np.log10(10.0)])
 #parameters.append(['lnTe',np.log(np.exp(-8.0)/3600.0),np.log(np.exp(8.0)/3600.0)])
 #parameters.append(['lnTe',np.log(np.exp(-8.0)/3600.0),np.log(np.exp(8.0)/3600.0)])
 #parameters.append(['log10soil',np.log10(1.0),np.log10(2.00)])
 #parameters.append(['psoil',1.0,1.0])
 parameters.append(['psoil',0.01,1.0])
 parameters.append(['pksat',0.25,4.0])
 parameters.append(['sdmax',0.1,1.0])

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
  nclusters = clusters_info[icatch]

  #Cycle through the ensemble of clusters
  for iens in xrange(nens):

   #Define the parameters
   parameters = {}
   parameters['log10m'] = param_values[iens,0]
   #parameters['lnTe'] = param_values[iens,1]
   #parameters['log10soil'] = param_values[iens,2]
   parameters['psoil'] = param_values[iens,1]
   parameters['pksat'] = param_values[iens,2]
   parameters['sdmax'] = param_values[iens,3]

   #RIGGED
   #parameters = {}
   #parameters['log10m'] = -2.0#-2.582977995297425888e+00
   #parameters['lnTe'] = -3.12#-1.963648774068431635e-01
   #parameters['log10soil'] = np.log10(1.0)#1.389834359162560144e-02
   #parameters['sdmax'] = 1.0#1.938762117265730334e+00

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
 nelements_rank = len(icatchs)*nens/size#int(len(icatchs)*float(nens/float(size)))
 ncatch_rank = int(float(nens)/(float(size)/float(len(icatchs))))
 for ielement in np.arange(len(elements.keys()))[rank*nelements_rank:(rank+1)*nelements_rank]:

  #Define the info
  element = elements[ielement]

  #Print where we are at
  print 'Catchment %d, Ensemble %d' % (element['icatch'],element['iens']),element['nclusters'],element['parameters']

  #Define the info
  hydrobloks_info = {
        'icatch':element['icatch'],
	'rank':rank,
        'input':'%s/input/data.pck' % dir,
        'dt':3600.,
        #'dt':900.,
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
   #pickle.dump(input,open('%d.pck' % icatch,'wb')) 
   #input = pickle.load(open('%d.pck' % icatch)) 
   #exit()
   elements['nclusters'] = input['nclusters']
   dt = time.time() - time0
   print "time to prepare the data",element['nclusters'],dt
   #Memorize the original soil file
   soilfile = input['files']['soils']
   #exit()

  #Update the soils file
  input = Update_Soils(input,soilfile,element['parameters'],element['icatch'],rank)

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
  #Determine the rank of the ensemble
  nfile = element['iens']/ncatch_rank
  file_netcdf = 'LHSoutput/%d_%d.nc' % (element['icatch'],nfile)
  #Set ensemble number
  iens = element['iens']
  if nfile > 0:iens = iens % nfile
  #Update the netcdf file
  vars = output['variables'].keys()
  for var in vars:
   if var not in ['smc1','lh']:del output['variables'][var]
  if element['iens'] < ncatch_rank: update_netcdf(element['cdir'],iens,nens,parameters,file_netcdf,input,output,nfile)
  else: update_netcdf(element['cdir'],iens,ncatch_rank,parameters,file_netcdf,input,output,nfile)

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
 input_dir = workspace#'%s/input' % hydrobloks_info['dir']

 #Read in the metadata
 file = '%s/workspace_info.pck' % workspace
 wbd = pickle.load(open(file))

 #Create the netcdf file
 file_netcdf = hydrobloks_info['input']
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


 #Create the dictionary to hold all of the data
 output = {}

 #Create the Latin Hypercube (Clustering)
 nclusters = hydrobloks_info['nclusters']
 ncores = hydrobloks_info['ncores']
 icatch = hydrobloks_info['icatch']
 rank = hydrobloks_info['rank']

 #Create the clusters and their connections
 output = Create_Clusters_And_Connections(workspace,wbd,output,input_dir,nclusters,ncores,icatch,rank,info,hydrobloks_info)

 #Extract the meteorological forcing
 output = Prepare_HSU_Meteorology(workspace,wbd,output,input_dir,info,hydrobloks_info)

 #Write out the files to the netcdf file
 fp = hydrobloks_info['input_fp']
 data = output
 #Create the dimensions
 #ntime = data['meteorology']['wind'].shape[0]
 #nhsu = hydrobloks_info['nclusters']
 #fp.createDimension('hsu',nhsu)
 #fp.createDimension('time',ntime)

 #Write the meteorology
 #grp = fp.createGroup('meteorology')
 #for var in data['meteorology']:
 # grp.createVariable(var,'f4',('time','hsu'))
 # grp.variables[var][:] = data['meteorology'][var][:]

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

 #Write the model parameters
 grp = fp.createGroup('parameters')
 vars = ['slope','vof','area_pct','land_cover','channel',
        'vchan','dem','soil_texture_class','ti','carea','area',
        'WLTSMC','MAXSMC','DRYSMC','REFSMC','SATDK']
 for var in vars:
  grp.createVariable(var,'f4',('hsu',))
  grp.variables[var][:] = data['hsu'][var]
 #for hsu in data['hsu']:
 # for var in vars:
 #  if var not in ['WLTSMC','MAXSMC','DRYSMC','REFSMC','SATDK']:
 #   grp.variables[var][hsu] = data['hsu'][hsu][var]
 #  else:
 #   grp.variables[var][hsu] = data['hsu'][hsu]['soil_parameters'][var]

 #Write other metadata
 grp = fp.createGroup('metadata')
 grp.outlet_hsu = data['outlet']['hsu']

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

def Compute_HRUs_Semidistributed(covariates,mask,nclusters):

 #Define the covariates
 info = {'area':{'data':covariates['carea'][mask == True],},
        'slope':{'data':covariates['cslope'][mask == True],},
        'sms':{'data':covariates['MAXSMC'][mask == True],},
        'smw':{'data':covariates['WLTSMC'][mask == True],},
        #'clay':{'data':covariates['clay'][mask_woc == True],},
        #'sand':{'data':covariates['sand'][mask_woc == True],},
        'ndvi':{'data':covariates['ndvi'][mask ==True],},
        #'nlcd':{'data':covariates['nlcd'][mask_woc ==True],},
        #'ti':{'data':covariates['ti'][mask_woc == True],},
        'dem':{'data':covariates['dem'][mask == True],},
        'lats':{'data':covariates['lats'][mask == True],},
        'lons':{'data':covariates['lons'][mask == True],},
        }

 #Scale all the variables (Calculate the percentiles
 for var in info:
  argsort = np.argsort(info[var]['data'])
  pcts = np.copy(info[var]['data'])
  pcts[argsort] = np.linspace(0,1,len(info[var]['data']))
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
 minsamples = 10**5
 if X.shape[0] > minsamples:
  Xf = X[np.random.choice(np.arange(X.shape[0]),minsamples),:]
 else:
  Xf = X
 #Initialize all points at the 0.5 point
 init = 0.5*np.ones((nclusters,Xf.shape[1]))
 batch_size = 25*nclusters
 init_size = 3*batch_size
 clf = sklearn.cluster.MiniBatchKMeans(nclusters,random_state=1,init=init,batch_size=batch_size,init_size=init_size)
 clf.fit(Xf)#
 clf_output = clf.predict(X)
 #Reassign the ids
 clf_output_copy = np.copy(clf_output)
 for cid in xrange(len(np.unique(clf_output))):
  clf_output[clf_output_copy == np.unique(clf_output)[cid]] = cid
 cluster_ids = np.empty(covariates['ti'].shape)
 cluster_ids[:] = -9999
 cluster_ids[mask == True] = clf_output
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

 return (cluster_ids,)

def Assign_Parameters_Fulldistributed(covariates,metadata,hydrobloks_info,OUTPUT,cluster_ids,mask):

 nclusters = hydrobloks_info['nclusters']
 #Initialize the arrays
 vars = ['area','area_pct','BB','DRYSMC','F11','MAXSMC','REFSMC','SATPSI',
         'SATDK','SATDW','WLTSMC','QTZ','slope','ti','dem','carea','channel',         'vchan','vof','land_cover','soil_texture_class']
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
  #Vchan
  OUTPUT['hsu']['vchan'][:] = 1000 #m/hr
  #Vof
  OUTPUT['hsu']['vof'][:] = 100 #m/hr
  #Land cover type  
  land_cover = np.copy(covariates['nlcd'])
  for lc in np.unique(land_cover)[1:]:
   land_cover[land_cover == lc] = NLCD2NOAH[lc]
  #OUTPUT['hsu']['land_cover'][:] = NLCD2NOAH[stats.mode(covariates['nlcd'][idx])[0][0]]
  #Soil texture class
  OUTPUT['hsu']['soil_texture_class'][:] = covariates['TEXTURE_CLASS'][mask]

 return OUTPUT

def Assign_Parameters_Semidistributed(covariates,metadata,hydrobloks_info,OUTPUT,cluster_ids,mask):

 nclusters = hydrobloks_info['nclusters']
 #Initialize the arrays
 vars = ['area','area_pct','BB','DRYSMC','F11','MAXSMC','REFSMC','SATPSI',
         'SATDK','SATDW','WLTSMC','QTZ','slope','ti','dem','carea','channel',
         'vchan','vof','land_cover','soil_texture_class']
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
  OUTPUT['hsu']['channel'][hsu] = np.mean(covariates['channels'][idx])
  #Vchan
  OUTPUT['hsu']['vchan'][hsu] = 1000 #m/hr
  #Vof
  OUTPUT['hsu']['vof'][hsu] = 100 #m/hr
  #Land cover type  
  OUTPUT['hsu']['land_cover'][hsu] = NLCD2NOAH[stats.mode(covariates['nlcd'][idx])[0][0]]
  #Soil texture class
  OUTPUT['hsu']['soil_texture_class'][hsu] = stats.mode(covariates['TEXTURE_CLASS'][idx])[0][0]

 return OUTPUT

def Calculate_Flow_Matrix(covariates,cluster_ids,nclusters):

 #Prepare the flow matrix
 mask1 = covariates['fdir'] < 0
 covariates['fdir'][mask1] = -9999.0
 cluster_ids_copy = np.copy(cluster_ids)
 cluster_ids_copy[mask1] = np.nan
 max_nhru = np.sum(mask1)
 #tp_matrix = mt.preprocessor.calculate_connections_d8(cluster_ids_copy,covariates['fdir'],nclusters,max_nhru)
 (hrus_dst,hrus_org) = mt.preprocessor.calculate_connections_d8(cluster_ids_copy,covariates['fdir'],nclusters,max_nhru)
 #Only use the non -9999 values
 hrus_dst = hrus_dst[hrus_dst != -9999]-1
 hrus_org = hrus_org[hrus_org != -9999]-1
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

 return flow_matrix.T

def Create_Soils_File(hydrobloks_info,OUTPUT,input_dir,icatch,rank):

 #Read in table of NOAH soil parameter values
 fp = open('Model/pyNoahMP/data/SOILPARM.TBL')
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
 soils_lookup = '%s/SOILPARM_%d_%d.TBL' % (soilsdir,icatch,rank)
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

def Create_Clusters_And_Connections(workspace,wbd,output,input_dir,nclusters,ncores,icatch,rank,info,hydrobloks_info):
 
 print "Creating and curating the covariates"
 (covariates,mask) = Create_and_Curate_Covariates(wbd)

 #Determine the HRUs (clustering if semidistributed; grid cell if fully distributed)
 print "Computing the HRUs"
 if hydrobloks_info['model_type'] == 'semi':
  (cluster_ids,) = Compute_HRUs_Semidistributed(covariates,mask,nclusters)
 elif hydrobloks_info['model_type'] == 'full':
  nclusters = np.sum(mask == True)
  hydrobloks_info['nclusters'] = nclusters
  (cluster_ids,) = Compute_HRUs_Fulldistributed(covariates,mask,nclusters)

 #Prepare the flow matrix
 print "Calculating the flow matrix"
 flow_matrix = Calculate_Flow_Matrix(covariates,cluster_ids,nclusters)

 #Define the metadata
 metadata = gdal_tools.retrieve_metadata(wbd['files']['ti'])

 #Make the output dictionary for the basin
 OUTPUT = {'hsu':{},'metadata':metadata,'mask':mask,'flow_matrix':flow_matrix}

 #Determine outlet cell
 covariates['carea'][mask == False] = np.nan
 outlet_idx = np.where(covariates['carea'] == np.max(covariates['carea'][np.isnan(covariates['carea']) == 0]))
 outlet_idx = [int(outlet_idx[0]),int(outlet_idx[1])]
 OUTPUT['outlet'] = {'idx':outlet_idx,'hsu':cluster_ids[outlet_idx[0],outlet_idx[1]]}

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
 soils_lookup = Create_Soils_File(hydrobloks_info,OUTPUT,input_dir,icatch,rank)

 #Save the soils file name to the model input
 OUTPUT['files'] = {'soils':soils_lookup,}

 #Add the new number of clusters
 OUTPUT['nclusters'] = nclusters

 return OUTPUT

def Prepare_HSU_Meteorology(workspace,wbd,OUTPUT,input_dir,info,hydrobloks_info):

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
  #Determine the time steps to retrieve
  dates = nc.num2date(fp.variables['t'][:],units='hours since %02d-%02d-%02d 00:00:00' % (idate.year,idate.month,idate.day))
  mask_dates = (dates >= idate) & (dates <= fdate)
  data = np.ma.getdata(fp.variables[var][mask_dates])
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

 #Append the meteorology to the output dictionary
 #OUTPUT['meteorology'] = meteorology

 return OUTPUT

def create_netcdf_file(file_netcdf,output,nens,cdir,rank):

 #Create the file
 fp = nc.Dataset(file_netcdf, 'w', format='NETCDF4')

 #Create the dimensions
 ntime = output['variables']['smc1'].shape[0]/24
 nhsu = output['variables']['smc1'].shape[1]
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

 if rank != 0: 
  fp.close()
  return

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
 #grp.HUC = input['wbd']['HUC10']
 #Define outlet hsu
 grp.outlet_hsu = int(output['misc']['outlet_hsu'])
 #Define catchment name
 #grp.catchment_name = input['wbd']['Name']
 #Define catchment area
 #grp.AreaSqKm = input['wbd']['AreaSqKm']
 #Add HSU transition probability matrix
 #tp = grp.createVariable('tpm','f4',('hsu','hsu'))
 #tp.description = 'Transition probability matrix between hsus'
 #tp[:] =  input['tp']

 #Retrieve the conus_albers metadata
 #metadata = gdal_tools.retrieve_metadata(input['wbd']['files']['ti']) 
 #metadata['nodata'] = -10000.0
 #Save the conus_albers metadata
 #fp.createDimension('nx',metadata['nx'])
 #fp.createDimension('ny',metadata['ny'])
 #hmca = grp.createVariable('hmca','f4',('ny','nx')) 
 #hmca.gt = metadata['gt']
 #hmca.projection = metadata['projection']
 #hmca.description = 'HSU mapping (conus albers)'
 #hmca.nodata = metadata['nodata']
 #Save the conus albers mapping
 #hsu_map = np.copy(input['hsu_map'])
 #hsu_map[np.isnan(hsu_map) == 1] = metadata['nodata']
 #hmca[:] = hsu_map
 '''#Write out the mapping
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
 hmll[:] = hsu_map'''
 
 #Close the file 
 fp.close()

 return

def update_netcdf(cdir,iens,nens,parameters,file_netcdf,output,rank):
 
 #Convert variables to arrays
 for var in output['variables']:
  output['variables'][var] = np.array(output['variables'][var])

 #Create the netcdf file if necessary
 if iens == 0: create_netcdf_file(file_netcdf,output,nens,cdir,rank)

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

def Update_Soils(input,soilfile,parameters,icatch,rank):

 #Read in table of NOAH soil parameter values
 dtype = {'names':('ID','BB','DRYSMC','F11','MAXSMC','REFSMC','SATPSI','SATDK','SATDW','WLTSMC','QTZ'),
          'formats':('d4','f4','f4','f4','f4','f4','f4','f4','f4','f4','f4')}
          
 soils_data = np.loadtxt(soilfile,skiprows=3,delimiter=',',dtype=dtype)
 soilsdir = "/".join(soilfile.split('/')[0:-1])
 soils_lookup = '%s/SOILPARM_ensemble_%d_%d.TBL' % (soilsdir,icatch,rank)
 nhsus = soils_data['ID'].size
 soil_vars = ['BB','DRYSMC','F11','MAXSMC','REFSMC','SATPSI','SATDK','SATDW','WLTSMC','QTZ']
 fp = open(soils_lookup,'w')
 fp.write('Soil Parameters\n')
 fp.write('CUST\n')
 fp.write("%d,1   'BB      DRYSMC      F11     MAXSMC   REFSMC   SATPSI  SATDK       SATDW     WLTSMC  QTZ    '\n" % nhsus)
 for hsu in input['hsu']:
  fp.write('%d, ' % (hsu+1))
  idx = soils_data['ID'] == hsu + 1
  for var in soil_vars:
   #if var in ['DRYSMC','MAXSMC','REFSMC','WLTSMC']:
   if var in ['DRYSMC','WLTSMC','REFSMC']:
    tmp = parameters['psoil']*soils_data[var][idx]
    #tmp = soils_data[var][idx]
    input['hsu'][hsu]['soil_parameters'][var] = tmp
    fp.write('%.10f, ' % tmp)
   elif var in ['SATDK',]:
    tmp = parameters['pksat']*soils_data[var][idx]
    #tmp = soils_data[var][idx]
    input['hsu'][hsu]['soil_parameters'][var] = tmp
    fp.write('%.10f, ' % tmp)
   else:
    fp.write('%.10f, ' % soils_data[var][idx])
  fp.write('\n')
 fp.close()

 #Save the soils file name to the model input
 input['files']['soils'] = soils_lookup

 return input
