import cPickle as pickle

def Deterministic(info):

 #Read in the catchment database
 wbd = pickle.load(open(info['wbd']))

 #Define the dates
 idate = datetime.datetime(2000,1,1,0)
 fdate = datetime.datetime(2000,1,31,23)

 #Iterate through all the catchments until done
 for icatch in [3637,]:#len(wbd.keys()):

  #Define the info
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

  #Define the parameters
  parameters['log10m']
  parameters['lnTe']
  parameters['log10soil']
  parameters['sdmax']

  #Run the model

 return

def Latin_Hypercube_Sample(info):

 print info

 return

def Convergence_Analysis(info):

 print info
 #Create a sample of all the combinations

 #Make way through all of them to find convergence

 return
