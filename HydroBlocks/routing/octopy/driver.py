import os
ncores = 64
#Run the preprocessor
os.system('mpirun -n %d python preprocessor.py' % (ncores,))
#Run the model
os.system('mpirun -n %d python driver_core.py' % (ncores,))
#Initialize each sub-region with the necessary information
#Calculate the lateral inflow/outflow per reach within the targe region
#Communicate with the other sub-regions
#Update the PDE solution
