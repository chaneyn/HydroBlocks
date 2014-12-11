import cPickle as pickle
import numpy as np
import metrics
import glob
import os
from scipy.optimize import curve_fit
from mpi4py import MPI
comm = MPI.COMM_WORLD

def func(x,a,b):
 return a*x**b

def Extract_Info(icatch,mapping):

 #Extract the info for the catchment
 data = {'dt':[],'nclusters':[]}
 for i in mapping[icatch]:
  #print icatch,i
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
 #args = np.argsort(data['nclusters'])
 #print data['nclusters'][args == 1]
 #exit()
 #argmax = np.argmax(data['nclusters'])
 #kgesall = []
 #vars = ['lh','sh','smc1','prcp','qexcess','qsurface','swe']
 #for var in vars:
 # ivar = vars.index(var) + 1
 # #Compute the metrics
 # kges = []
 #or i in xrange(len(data['nclusters'])):
 #  #kges.append(metrics.KGE(data['vars'][var]['std'][i,:],data['vars'][var]['std'][argmax,:]))
 # kgesall.append(kges)
 #kgesall = np.array(kgesall).T
 nclusters = np.array(data['nclusters'])
 argsort = np.argsort(nclusters)
 nclusters = nclusters[argsort]
 rts = []
 for var in ['smc1']:#data['vars']:
  mean = np.mean(data['vars'][var]['std'][:,8500:],axis=1)
  mean = mean[argsort]
  #mean = (mean - np.min(mean))/(np.max(mean) - np.min(mean))
  #Construct relative tolerance 
  rt = []
  val = 5
  for i in xrange(val,len(mean)):
   #rt.append(np.abs((mean[i] - mean[i-1])/min(mean[i],mean[i-1])))
   #rt.append(np.abs((mean[i] - mean[i-1])))
   rt.append(np.mean(np.abs(mean[i] - mean[i-val:i-1])))
  #rt.append(np.max(np.abs(numerator/denominator)))
  #Scale from 0-1
  rts.append(rt)
 rts = np.array(rts)
 rts = np.mean(rts,axis=0)
 #Fit the curve
 popt, pcov = curve_fit(func,nclusters[val::],rts)
 #Estimate rts for a lot more...
 snclusters = np.arange(2,10**6)
 srts = func(snclusters,popt[0],popt[1])
 thresholds = [0.0001,0.001,0.01]
 nclustersout = []
 for threshold in thresholds:
  mask = srts <= threshold
  tmp = snclusters[mask]
  if len(tmp) > 0:nclustersout.append(tmp[0])
  else: nclustersout.append(np.nan)
 #If the number of runs is below x then skip this cluster
 if nclusters.size < 50:nclustersout = [np.nan,np.nan,np.nan]

 return nclustersout

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

