import sys
import processing
import numpy as np
import datetime
import os
import cPickle as pickle

#Get general info
dir = sys.argv[1]
ncores = int(sys.argv[2])
rank = int(sys.argv[3])
size = int(sys.argv[4])

#Read in the LHS sampled parameter space
data = np.loadtxt('LHS/LHS_sampling.txt',
       dtype={'names':('log10m','lnTe','log10soil','sdmax'),
              'formats':('f4','f4','f4','f4')})

#Determine the number of ensembles
nens = data.size

#Read in the general info dataset
general_info = pickle.load(open('../PREPROCESSING/general_info.pck'))
idate = general_info['time_info']['startdate']#datetime.datetime(2009,1,1,0)
fdate = datetime.datetime(2003,12,31,23) #HERE
#fdate = general_info['time_info']['enddate']#datetime.datetime(2009,12,31,23)

#Read in the wbd dataset
wbd = pickle.load(open('../PREPROCESSING/wbd.pck'))
ncatch = len(wbd)

#Define the core info
#ncores = 16

#Read the scan sites info 
metadata = pickle.load(open('../ANALYSIS/SCAN/scan_metadata.pck'))

#Iterate through all this cores catchments
#for icatch in [14009,]:
#for icatch in [11918,]:
#for icatch in [11705,]:#11531,]:
#for icatch in [11911,]:
hucs = ['031102040101',]#'031102040103','031102040102','031102040104']
#Randomize the catchments 
np.random.seed(1)
catchments = np.arange(ncatch)
np.random.shuffle(catchments)

#for icatch in xrange(rank,ncatch,size):
for icatch in catchments[rank::size]:
 #if wbd[icatch]['HUC12'] != "050500030107":continue
 #if wbd[icatch]['HUC12'] not in hucs:continue
 #if icatch != 62656:continue#50915:continue#63313: continue
 #if icatch != 86189:continue
 if icatch != 3637:continue#43772:continue
 #if icatch != 8756:continue#43772:continue,500
 #if icatch != 500:continue#8756:continue#43772:continue,500
 #icatch = int(metadata[metadata.keys()[rank]]['icatch']) #ONLY FOR SCAN SITES
 #print icatch

 #Define the catchment directory
 cdir = '%s/catchments/catch_%d' % (dir,icatch)

 #Determine if it exists
 if os.path.exists(cdir) == False: continue

 #Make the output directory
 os.system('mkdir -p %s/output' % cdir)

 #Open the redirected output
 stdout = '%s/output/stdout.txt' % cdir
 stderr = '%s/output/stderr.txt' % cdir
 os.system('rm -f %s' % stdout)
 os.system('rm -f %s' % stderr)
 fout = os.open(stdout,os.O_RDWR|os.O_CREAT)
 ferr = os.open(stderr,os.O_RDWR|os.O_CREAT)

 #Run the model for all parameter sets
 for iens in xrange(nens):
 
  #Extract the parameters
  parameters = {}
  for var in ['log10m','lnTe','log10soil','sdmax']:
   parameters[var] = data[var][iens]
 
  print 'Catchment number: %d, Ensemble number: %d (Running model)' % (icatch,iens)

  #save the current file descriptors to a tuple
  save = os.dup(1),os.dup(2)

  #Redirect the output to log
  #os.dup2(fout, 1),os.dup2(ferr,2)
  #print 'Catchment number: %d, Ensemble number: %d (Running model)' % (icatch,iens)

  #Run the model
  processing.run_model(parameters,cdir,idate,fdate,ncores,wbd[icatch])
  
  #Redirect the output to terminal
  os.dup2(save[0], 1),os.dup2(save[1],2)

  #Add output to the ensemble file
  print 'Catchment number: %d, Ensemble number: %d (Writing to netcdf)' % (icatch,iens)
  processing.update_netcdf(cdir,iens,nens,parameters)

 #Close the redirected output
 os.close(fout),os.close(ferr)

 exit()

