import sys
sys.path.append('Tools')
import HydroBloksTools as HBM

#INFO
#parallel_flag = sys.argv[1]
parallel_flag = sys.argv[1]
dir = '/scratch/sciteam/nchaney/data/CONUS_SIMULATIONS_HUC10/catchments'
run_flag = 'Convergence Analysis'
ncores = 16#32

if parallel_flag == 'parallel':
 #from mpi4py import MPI
 #comm = MPI.COMM_WORLD
 
 info = {
	'rank':int(sys.argv[2]),#comm.rank,
	'size':int(sys.argv[3]),#comm.size,
	'wbd':'wbd.pck',
	'dir':dir,
	'ncores':ncores,
	}

elif parallel_flag == 'serial':
  info = {
        'rank':0,
        'size':1,
        'wbd':'wbd.pck',
        'dir':dir,
	'ncores':ncores,
        }

#Deterministic
if run_flag == 'Deterministic': HBM.Deterministic(info)

#Latin Hypercube Sample
if run_flag == 'Latin Hypercube Sample': HBM.Latin_Hypercube_Sample(info)

#Deterministic
if run_flag == 'Convergence Analysis': HBM.Convergence_Analysis(info)
