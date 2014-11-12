import os
Parallel = True
njobs = 1
ncores = 8

if Parallel == True:
 
 #Run the model using MPI (w/OMP)
 os.system('aprun -n %d -d %d python Driver.py parallel' % (njobs,ncores))

elif Parallel == False:

 #Run the model without MPI (w/OMP)
 os.system('python Driver.py serial')
