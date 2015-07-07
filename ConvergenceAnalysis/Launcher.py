import os
#Initialize all the processes
os.system('rm -rf /scratch/02179/chaneyna/convergence_analysis/workspace')
mpi_type = 'ibrun'
njobs = 10
if mpi_type == 'ibrun':
 python = '/home1/02179/chaneyna/libraries/python/lib/python2.7/site-packages/mpi4py/bin/python-mpi'
 os.system('ibrun -np %d %s Driver.py' % (njobs,python))
if mpi_type == 'mpirun':
 os.system('mpirun -n %d python Driver.py' % (njobs,))
