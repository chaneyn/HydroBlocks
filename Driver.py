import sys
sys.path.append('Tools')
import HydroBloksTools as HBM

#INFO
#parallel_flag = sys.argv[1]
parallel_flag = sys.argv[1]
#dir = '/scratch/02179/chaneyna/catchments'#'../LittleWashitaRegion'#'Test'
#dir = '/scratch/02179/chaneyna/conus30m/CONUS_SIMULATIONS/catchments'
dir = 'Test'
wbd_file = '%s/catchments.pck' % dir#'Test/catchments.pck'
#wbd_file = '../test.pck'
#run_flag = 'Convergence Analysis'
#run_flag = 'Latin Hypercube Sample'
run_flag = 'Deterministic'
ncores = 32

if parallel_flag == 'parallel':
 from mpi4py import MPI
 comm = MPI.COMM_WORLD
 
 info = {
	'rank':comm.rank,
	'size':comm.size,
	'wbd':wbd_file,
	'dir':dir,
	'ncores':ncores,
	}

elif parallel_flag == 'serial':
  info = {
        'rank':0,
        'size':1,
        'wbd':wbd_file,
        'dir':dir,
	'ncores':ncores,
        }

#Deterministic
if run_flag == 'Deterministic': HBM.Deterministic(info)

#Latin Hypercube Sample (CURRENTLY NOT WORKING)
if run_flag == 'Latin Hypercube Sample': HBM.Latin_Hypercube_Sample(info)
 
#Convergence Analysis (CURRENTLY NOT WORKING)
if run_flag == 'Convergence Analysis': HBM.Convergence_Analysis(info)
