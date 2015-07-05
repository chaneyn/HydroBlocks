import os
import sys
sys.path.append('Tools')
from HydroBloksTools import Read_Metadata_File

#Read in the metadata file
metadata_file = sys.argv[1]
metadata = Read_Metadata_File(metadata_file)
Parallel = metadata['parallel_flag']#False#True 
njobs = metadata['parallel_nnodes']#800#1600#800#80000#5000#8000#800#400
ncores = metadata['parallel_ncores']#2
mpi_type = metadata['parallel_mpitype']#'ibrun'
python = '/home1/02179/chaneyna/libraries/python/lib/python2.7/site-packages/mpi4py/bin/python-mpi'

if Parallel == True:
 
 #Run the model using MPI (w/OMP)
 #os.system('aprun -n %d -d %d python Driver.py parallel' % (njobs,ncores))
 if mpi_type == 'ibrun':
  os.system('ibrun -np %d %s Driver.py %s' % (njobs,python,metadata_file))
 if mpi_type == 'mpirun':
  os.system('mpirun -n %d python Driver.py %s' % (njobs,metadata_file))

elif Parallel == False:

 #Run the model without MPI (w/OMP)
 os.system('python Driver.py %s' % metadata_file)
