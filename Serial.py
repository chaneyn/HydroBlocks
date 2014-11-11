import sys
sys.path.append('Model')
import HydroBloks as HB
#import Tools as processing
import numpy as np
import datetime
import os
import cPickle as pickle

#Get general info
ncores = 1
dir = 'ReynoldsCreek'
        
#Read in the LHS sampled parameter space
data = np.loadtxt('%s/parameters/LHS_sampling.txt' % dir,
       dtype={'names':('log10m','lnTe','log10soil','sdmax'),
              'formats':('f4','f4','f4','f4')})
parameters = {}
iens = 0
for var in ['log10m','lnTe','log10soil','sdmax']:
 parameters[var] = data[var][iens]

#Read in the watershed dataset
wbd = pickle.load(open('%s/wbd.pck' % dir))

#Read in the general info dataset
idate = datetime.datetime(2000,1,1,0)
fdate = datetime.datetime(2000,1,31,23)

#Run HydroBloks on the model
icatch = 3637
hydrobloks_info = {
        'input':'%s/input/data.pck' % dir,
        'dt':3600.,
        'nsoil':20,
        'wbd':wbd[icatch],
        'ncores':1,
        'idate':idate,
        'fdate':fdate,
	'parameters':parameters
        }

#Define the catchment directory
cdir = dir #'%s/catchments/catch_%d' % (dir,icatch)

#Make the output directory
os.system('mkdir -p %s/output' % cdir)

#Run the model
output = HB.run_model(hydrobloks_info)
