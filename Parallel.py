import os

#Run the model for each catchment
dir = 'ReynoldsCreek'
njobs = 4
ncores = 8
type = 'Deterministic'
os.system('aprun -n %d -d %d python Tools/driver_node.py %s %s' % (njobs,ncores,dir,type))
