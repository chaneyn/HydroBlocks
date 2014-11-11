from mpi4py import MPI
comm = MPI.COMM_WORLD
import HydroBloksMPI as HBM
import sys

#Define the directory and run type
dir = sys.argv[1]
type = sys.argv[2]

#Gather info
info = {'rank':comm.rank,
	'size':comm.size,}

#Deterministic
if type == 'Deterministic': HBM.Deterministic(info)

#Latin Hypercube Sample
if type == 'Latin Hypercube Sample': HBM.Latin_Hypercube_Sample(info)

#Deterministic
if type == 'Convergence Analysis': HBM.Convergence_Analysis(info)
