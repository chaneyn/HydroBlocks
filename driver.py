import sys
import processing
import numpy as np
import datetime
import os
import cPickle as pickle

#Get general info
ncores = 1
dir = '/scratch/sciteam/nchaney/data/CONUS_SIMULATIONS_HUC10'

#Read in the LHS sampled parameter space
data = np.loadtxt('LHS/LHS_sampling.txt',
       dtype={'names':('log10m','lnTe','log10soil','sdmax'),
              'formats':('f4','f4','f4','f4')})

#Read in the watershed dataset
wbd = pickle.load(open('wbd.pck'))

#Determine the number of ensembles
nens = data.size

#Read in the general info dataset
idate = datetime.datetime(2000,1,1,0)
fdate = datetime.datetime(2003,12,31,23)

#Run HydroBloks on the model
icatch = 3637
catchid = 100

#Define the catchment directory
cdir = '%s/catchments/catch_%d' % (dir,icatch)

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

 #Run the model
 processing.run_model(parameters,cdir,idate,fdate,ncores,wbd[icatch])
  
 #Redirect the output to terminal
 os.dup2(save[0], 1),os.dup2(save[1],2)

 #Add output to the ensemble file
 print 'Catchment number: %d, Ensemble number: %d (Writing to netcdf)' % (icatch,iens)
 processing.update_netcdf(cdir,iens,nens,parameters)

#Close the redirected output
os.close(fout),os.close(ferr)
