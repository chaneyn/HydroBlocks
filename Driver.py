import sys
sys.path.append('Tools')
import HydroBloksTools as HBM

#Read in the metadata file
metadata_file = sys.argv[1]
metadata = HBM.Read_Metadata_File(metadata_file)
parallel_flag = metadata['parallel_flag']
ncores = metadata['parallel_ncores']
run_flag = metadata['run_type']
dir = metadata['directory']
wbd_file = '%s/catchments.pck' % dir

if parallel_flag == True:
 from mpi4py import MPI
 comm = MPI.COMM_WORLD
 
 info = {
	'rank':comm.rank,
	'size':comm.size,
	'wbd':wbd_file,
	'dir':dir,
	'ncores':ncores,
	}

elif parallel_flag == False:
  info = {
        'rank':0,
        'size':1,
        'wbd':wbd_file,
        'dir':dir,
	'ncores':ncores,
        }

#Deterministic
if run_flag == 'deterministic': HBM.Deterministic(info)

#Latin Hypercube Sample (CURRENTLY NOT WORKING)
if run_flag == 'Latin Hypercube Sample': HBM.Latin_Hypercube_Sample(info)
 
#Convergence Analysis (CURRENTLY NOT WORKING)
if run_flag == 'Convergence Analysis': HBM.Convergence_Analysis(info)
