import sys
sys.path.append('Tools')
import HydroBloksTools as HBM

#INFO
parallel_flag = sys.argv[1]
dir = 'ReynoldsCreek'
run_flag = 'Convergence Analysis'

if parallel_flag == 'parallel':
 from mpi4py import MPI
 comm = MPI.COMM_WORLD
 info = {
	'rank':comm.rank,
	'size':comm.size,
	'wbd':'%s/wbd.pck' % dir,
	'dir':dir,
	}

elif parallel_flag == 'serial':
  info = {
        'rank':0,
        'size':1,
        'wbd':'%s/wbd.pck' % dir,
        'dir':dir,
        }

#Deterministic
if run_flag == 'Deterministic': HBM.Deterministic(info)

#Latin Hypercube Sample
if run_flag == 'Latin Hypercube Sample': HBM.Latin_Hypercube_Sample(info)

#Deterministic
if run_flag == 'Convergence Analysis': HBM.Convergence_Analysis(info)
