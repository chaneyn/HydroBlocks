import os
Parallel = True#False#True 
njobs = 1#800#80000#5000#8000#800#400
ncores = 32
mpi_type = 'ibrun'

if Parallel == True:
 
 #Run the model using MPI (w/OMP)
 #os.system('aprun -n %d -d %d python Driver.py parallel' % (njobs,ncores))
 if mpi_type == 'ibrun':
  os.system('ibrun python Driver.py parallel')

elif Parallel == False:

 #Run the model without MPI (w/OMP)
 os.system('python Driver.py serial')
