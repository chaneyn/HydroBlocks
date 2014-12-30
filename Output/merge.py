import cPickle as pickle
import numpy as np
import metrics
import glob
import os
import scipy.special
from scipy.optimize import curve_fit
from mpi4py import MPI
comm = MPI.COMM_WORLD

def Extract_Info(icatch,mapping):

 #Extract the info for the catchment
 data = {'dt':[],'nclusters':[]}
 for i in mapping[icatch]:
  file = '/u/sciteam/nchaney/scratch/data/CONUS_SIMULATIONS_HUC10/output/%d.pck' % i
  if os.path.exists(file) == False: continue
  tmp = pickle.load(open(file))
  if icatch in tmp['icatch']: j = tmp['icatch'].index(icatch)
  else: continue
  if 'vars' not in data:
   data['vars'] = {}
   for var in tmp['vars']:
    data['vars'][var] = {'mean':[],'std':[]}
  #Fill in the data
  flag = 1
  for var in tmp['vars']:
   if tmp['vars'][var]['std'][j].size <= 1: 
    flag = 0
    break
   data['vars'][var]['mean'].append(tmp['vars'][var]['mean'][j])
   data['vars'][var]['std'].append(tmp['vars'][var]['std'][j])
  if flag == 0: continue
  data['dt'].append(tmp['dt'][j])
  data['nclusters'].append(tmp['nclusters'][j])

 #Convert to arrays
 for var in data['vars']:
  data['vars'][var]['mean'] = np.array(data['vars'][var]['mean'])
  data['vars'][var]['std'] = np.array(data['vars'][var]['std'])

 return data

def Compute_Metrics(data):

 data['nclusters'] = np.array(data['nclusters'])
 data['dt'] = np.array(data['dt'])
 nclusters = np.array(data['nclusters'])
 argsort = np.argsort(nclusters)
 nclusters = nclusters[argsort]
 snclusters = np.arange(5000000)
 dt = data['dt'][argsort]
 mean = np.mean(data['vars']['smc1']['std'][:,8500:],axis=1)
 mean = mean[argsort]
 #Calculate the rate of change
 dnclusters = nclusters[1::] - nclusters[0:-1]
 dmean = np.abs(mean[1::] - mean[0:-1])
 dmean = dmean[dnclusters != 0]
 nclusters = nclusters[1::][dnclusters != 0]
 mean = mean[1::][dnclusters != 0]
 dnclusters = dnclusters[dnclusters != 0]
 dmdc = dmean/dnclusters
 #Fit a power law to the rate of change
 p = np.poly1d(np.polyfit(np.log(nclusters),np.log(dmdc),1))
 snclusters = np.arange(2,5000000)
 sdmdc = np.exp(p(np.log(snclusters)))
 nclustersout = []
 for threshold in [10**-4,10**-5,10**-6]:
  mask = sdmdc <= threshold
  nclustersout.append(snclusters[mask][0])

 #If the number of runs is below x then skip this cluster
 if nclusters.size < 10:nclustersout = [np.nan,np.nan,np.nan]
 #Save the output
 output = {
           'nclusters':nclustersout,
          }

 #Add the mean and standard deviation of all the variables
 output['data'] = {}
 for var in data['vars']:
  output['data'][var] = {}
  for metric in ['mean','std']:
   mean = np.mean(data['vars'][var][metric][:,8500:],axis=1)
   mean = mean[argsort][-1]
   output['data'][var][metric] = mean
   
 return output

rank = comm.rank#0
size = comm.size#01

#Read in the list of icatch and core mapping
mapping = pickle.load(open('list/info_conus.pck'))

for icatch in mapping.keys()[rank::size]:

 print icatch
 #if icatch not in [10243,]:continue
 try:
 #if 1 == 1:
  #Extract the info
  data = Extract_Info(icatch,mapping)
  #pickle.dump(data,open('catch.pck','wb'))

  #Compute final metrics
  cametrics = Compute_Metrics(data)
  print cametrics

  #Output the metrics
  pickle.dump(cametrics,open('metrics/%d.pck' % icatch,'wb',pickle.HIGHEST_PROTOCOL))
 except:continue

