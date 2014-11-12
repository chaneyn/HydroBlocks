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

 print info
 #Define the rank and size
 rank = info['rank']
 size = info['size']

 #Read in the catchment database
 wbd = pickle.load(open(info['wbd']))

 #Define the dates
 idate = datetime.datetime(2000,1,1,0)
 fdate = datetime.datetime(2000,1,31,23)

 #Initialize the element count
 ielement = 0
 nens = 10
 elements = {}

 #Create a dictionary of information
 for icatch in [3637,]:#len(wbd.keys()):

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
   nbins={
	'area':np.random.randint(1,100),
	'slope':1,
	'sms':1,
	'ndvi':1,
	'ti':1,
	'dem':1,
	'channels':2
	}
  
   #Add the info to the dictionary
   elements[ielement] = {
		'parameters':parameters,
		'nbins':nbins,
		'icatch':icatch,
		'iens':iens,
		 } 

   #Update the element
   ielement += 1

  #Iterate through the dictionary elements
  for ielement in np.arange(len(elements.keys()))[rank::size]:

   #Define the info
   element = elements[ielement]

   #Print where we are at
   print 'Catchment %d, Ensemble %d' % (element['icatch'],element['iens']),element['nbins']

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
	'nbins':nbins
        }

   #Cluster the data
   Prepare_Model_Input_Data(hydrobloks_info)

   #Run the model
   output = HB.run_model(hydrobloks_info)

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
 print "Clustering the covariates and calculating the connections"
 output = Create_Clusters_And_Connections(workspace,wbd,output,input_dir,hydrobloks_info['nbins'])

 #Extract the meteorological forcing
 print "Preparing the meteorological forcing"
 output = Prepare_HSU_Meteorology(workspace,wbd,output,input_dir,info)

 #Add in the catchment info
 output['wbd'] = wbd

 #Save the data 
 file = '%s/data.pck' % input_dir
 pickle.dump(output,open(file,'wb'),pickle.HIGHEST_PROTOCOL)

 return

def Create_Clusters_And_Connections(workspace,wbd,output,input_dir,nbins):

 covariates = {}
 #Read in all the covariates
 for file in wbd['files']:
  covariates[file] = gdal_tools.read_raster(wbd['files'][file])
  if file == 'carea': covariates[file] = np.log(covariates[file])
  if file == 'cslope':
   mask = covariates[file] == 0.0
   covariates[file][mask] = 0.000001

 #Clean up the covariates
 covariates['ti'][covariates['ti'] > 14] = 14
 covariates['channels'][covariates['channels'] > 1] = 1
 covariates['channels'][covariates['channels'] < 0] = 0

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
  covariates[var][mask <= 0] = -9999

 #Might want to add lat/lon
 #Define the arrays to be used to create the LHS
 arrays = {'slope':covariates['cslope'].astype(np.float32),
          'area':covariates['carea'].astype(np.float32),
          'ndvi':covariates['ndvi'].astype(np.float32),
          'sms':covariates['MAXSMC'].astype(np.float32),
          'channels':covariates['channels'].astype(np.float32),
          'dem':covariates['dem'].astype(np.float32),
          'ti':covariates['ti'].astype(np.float32)}

 #Define the binning
 info = {'area':{'nbins':nbins['area'],'data':covariates['carea'][mask == True]},
        'slope':{'nbins':nbins['slope'],'data':covariates['cslope'][mask == True]},
        'sms':{'nbins':nbins['sms'],'data':covariates['MAXSMC'][mask == True]},
        'ndvi':{'nbins':nbins['ndvi'],'data':covariates['ndvi'][mask==True]},
        'ti':{'nbins':nbins['ti'],'data':covariates['ti'][mask==True]},
        'dem':{'nbins':nbins['dem'],'data':covariates['dem'][mask==True]},
        'channels':{'nbins':nbins['channels'],'data':covariates['channels'][mask==True]}
        }
 
 #Create the LHS bins
 bins,data = [],[]
 for id in info:
  #Set all nans to the mean
  info[id]['data'][np.isnan(info[id]['data']) == 1] = np.nanmean(info[id]['data'])
  #Find the percentiles
  pcts = np.linspace(0,100,info[id]['nbins']+1)
  scores = []
  for pct in pcts:
   scores.append(np.percentile(info[id]['data'],pct))
  bins.append(scores)
  data.append(info[id]['data'])
 #Remove non-unique edges
 for id in xrange(len(info.keys())):
  if info.keys()[id] != 'channels':
   bins[id] = np.unique(np.array(bins[id]))
  else:
   bins[id][1] = 0.5
   bins[id] = np.array(bins[id])
 bins = np.array(bins)
 data = np.array(data).T

 #Create the histogram
 H,edges = np.histogramdd(data,bins=bins)
 H = H/np.sum(H)

 #Create a dictionary of class info
 clusters = {}
 Hfl = H.flat
 icluster = -1
 for i in xrange(H.size):
  coords = Hfl.coords
  if H[coords] > 0:
   icluster = icluster + 1
   clusters[icluster] = {'pct':H[coords]}
   clusters[icluster]['bounds'] = {}  
   for id in info:
    key = info.keys().index(id)   
    clusters[icluster]['bounds'][id] = [edges[key][coords[key]],edges[key][coords[key]+1]]
  Hfl.next()

 #Map the cluster number to the grid
 cluster_ids = np.empty(covariates['ti'].shape)
 cluster_ids[:] = np.nan
 cmask = np.zeros(cluster_ids.shape).astype(np.int8)
 #print "Computing the hsu indices"
 for cid in clusters.keys():
  cmask[:] = 0
  for id in info.keys():
   if info.keys().index(id) == 0: string = "(arrays['%s'] >= clusters[%d]['bounds']['%s'][0]) & (arrays['%s'] <= clusters[%d]['bounds']['%s'][1])" % (id,cid,id,id,cid,id)
   else: string = string +  " & (arrays['%s'] >= clusters[%d]['bounds']['%s'][0]) & (arrays['%s'] <= clusters[%d]['bounds']['%s'][1])" % (id,cid,id,id,cid,id)
  idx = eval('np.where(%s)' % string)
  clusters[cid]['idx'] = idx
  cluster_ids[idx] = cid

 #Determine the links between clusters
 mask1 = covariates['fdir'] < 0
 covariates['fdir'][mask1] = -9999.0
 cluster_ids[mask1] = np.nan
 nclusters = len(clusters.keys())
 tp_matrix = mt.preprocessor.calculate_connections_d8(cluster_ids,covariates['fdir'],nclusters)
 #print tp_matrix

 #Create a plot of the map of clusters and transition probabilities
 plt.figure(figsize=(50,20))
 plt.subplot(121)
 plt.imshow(tp_matrix,interpolation='nearest')
 plt.axis('off')
 cb = plt.colorbar()
 cb.ax.tick_params(labelsize=45)
 plt.subplot(122)
 plt.imshow(cluster_ids,interpolation='nearest')
 plt.axis('off')
 cb = plt.colorbar()
 cb.ax.tick_params(labelsize=45)
 plt.tight_layout()
 file = '%s/cluster_connections.png' % workspace
 plt.savefig(file)

 #exit()
 #Create the input directory
 os.system('mkdir -p %s' % input_dir)

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

 #Soil properties
 soil_vars = ['BB','DRYSMC','F11','MAXSMC','REFSMC','SATPSI','SATDK','SATDW','WLTSMC','QTZ']
 nhsus = len(OUTPUT['hsu'])
 soils_lookup = '%s/SOILPARM.TBL' % input_dir
 fp = open(soils_lookup,'w')
 fp.write('Soil Parameters\n')
 fp.write('CUST\n')
 fp.write("%d,1   'BB      DRYSMC      F11     MAXSMC   REFSMC   SATPSI  SATDK       SATDW     WLTSMC  QTZ    '\n" % nhsus)
 for hsu in OUTPUT['hsu']:
  fp.write('%d, ' % (hsu+1))
  for var in soil_vars:
   fp.write('%.10f, ' % OUTPUT['hsu'][hsu]['soil_parameters'][var])
  fp.write('\n')
 fp.close()

 return OUTPUT

def Prepare_HSU_Meteorology(workspace,wbd,OUTPUT,input_dir,info):

 #Open grads
 ga = grads_tools.open_grads(info['binaries']['grads'])

 #Define gdalwarp
 gdalwarp = info['binaries']['gdalwarp']
 
 #Create the mapping
 mapping_dir = '%s/mapping' % workspace
 os.system('mkdir -p %s' % mapping_dir)
 #Calculate the fine to coarse scale mapping
 for data_var in wbd['files_meteorology']:

  #print data_var
  ctl = wbd['files_meteorology'][data_var]
  var = data_var.split('_')[1]
  ga("xdfopen %s" % ctl)

  #Write out a sample file for each variable
  ga("set gxout geotiff")
  tmp = '%s/tmp.tif' % (mapping_dir,)

  #Define th coarse and fine scale mapping
  file_coarse = '%s/%s_coarse.tif' % (mapping_dir,data_var)
  file_fine = '%s/%s_fine.tif' % (mapping_dir,data_var)
  maskij = ga.expr(var)
  maskij = ga.exp(var)
  for i in xrange(maskij.shape[0]):
   maskij[i,:] = np.arange(i*maskij.shape[1],(i+1)*maskij.shape[1])
  ga.imp('mask',maskij)
  ga("set geotiff %s" % file_coarse)
  ga("d mask")

  #Close access to grads file
  ga("close 1")

  #Regrid and downscale
  log = '%s/log.txt' % workspace
  os.system('%s -t_srs EPSG:102039 -overwrite -dstnodata -9999 -r near -tr %.16f %.16f -te %.16f %.16f %.16f %.16f %s %s >> %s 2>&1' % 
           (gdalwarp,wbd['bbox_albers']['res'],wbd['bbox_albers']['res'],wbd['bbox_albers']['minx'],wbd['bbox_albers']['miny'],
            wbd['bbox_albers']['maxx'],wbd['bbox_albers']['maxy'],file_coarse,file_fine,log))
  #Read in the fine scale version
  mask_fine = gdal_tools.read_raster(file_fine)

  #Compute the mapping for each hsu
  for hsu in OUTPUT['hsu']:
   idx = OUTPUT['hsu'][hsu]['idx']
   icells = np.unique(mask_fine[idx].astype(np.int))
   counts = np.bincount(mask_fine[idx].astype(np.int))
   coords,pcts = [],[]
   for icell in icells:
    ilat = int(np.floor(icell/maskij.shape[1]))
    jlat = icell - ilat*maskij.shape[1]
    pct = float(counts[icell])/float(np.sum(counts))
    coords.append([ilat,jlat])
    pcts.append(pct)
   pcts = np.array(pcts)
   coords = list(np.array(coords).T)
   OUTPUT['hsu'][hsu][var] = {'pcts':pcts,'coords':coords}

 #Reinitialize grads
 ga("reinit")

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

 #Close grads
 del ga

 return OUTPUT
