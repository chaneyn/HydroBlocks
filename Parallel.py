import os

#Run the model for each catchment
dir = 'test'
njobs = 4
ncores = 8
type = 'Deterministic'
os.system('aprun -n %d -d %d python driver_node.py %s %s' % (njobs,ncores,dir,type))
