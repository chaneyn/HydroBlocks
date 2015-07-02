import os
Parallel = True 
njobs = 1#256#800#80000#5000#8000#800#400
ncores = 2
mpi_type = 'ibrun'
python = '/home1/02179/chaneyna/libraries/python/lib/python2.7/site-packages/mpi4py/bin/python-mpi'

if Parallel == True:
 
 #Run the model using MPI (w/OMP)
 #os.system('aprun -n %d -d %d python Driver.py parallel' % (njobs,ncores))
 if mpi_type == 'ibrun':
  os.system('ibrun -np %d %s Driver.py parallel' % (njobs,python))
 if mpi_type == 'mpirun':
  os.system('mpirun -n %d python Driver.py parallel' % (njobs,))

elif Parallel == False:

 #Run the model without MPI (w/OMP)
 os.system('python Driver.py serial')
