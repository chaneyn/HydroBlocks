import os
Parallel = False#True 
njobs = 1#800#80000#5000#8000#800#400
ncores = 4

if Parallel == True:
 
 #Run the model using MPI (w/OMP)
 #os.system('aprun -n %d -d %d python Driver.py parallel' % (njobs,ncores))
 os.system('aprun -n %d -d %d run_mpi "python Driver.py parallel"' % (njobs,ncores))

elif Parallel == False:

 #Run the model without MPI (w/OMP)
 os.system('python Driver.py serial')
