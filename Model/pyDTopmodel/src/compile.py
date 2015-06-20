import os
os.system('rm -f *.mod')
os.system('rm -f *.pyf')
os.system('rm -f *.o')

#Compile the code to create the module
cmd = 'gfortran -c dynamic_topmodel_tools.f90 -lgomp -liomp5 -lmkl_rt -lpthread -lm -w -O3 -funroll-loops -m64 -Wl,--no-as-needed -openmp -fopenmp -g -fbacktrace -o dynamic_topmodel_tools.o'
os.system(cmd)

#Create a version with the include commented out
cmd = 'tail -n +2 dynamic_topmodel_tools.f90 | head -n -1 > tmp.f90'
os.system(cmd)

#Create the signature file
cmd = 'f2py tmp.f90 -h dynamic_topmodel_tools.pyf -m dynamic_topmodel_tools --overwrite-signature only: update :'
os.system(cmd)

#Create driver library
cmd = 'f2py tmp.f90 -lgomp -liomp5 -lmkl_rt -lpthread -lm -c dynamic_topmodel_tools.pyf --fcompiler=gnu95 --f90flags="-w -O3 -funroll-loops -m64 -Wl,--no-as-needed -openmp -fopenmp -g -fbacktrace"'
os.system(cmd)

#Move to the previous directory
os.system('mv dynamic_topmodel_tools.so ../.')

#Remove the temporary file
os.system('rm -f tmp.f90')
os.system('rm -f *.mod')
os.system('rm -f *.pyf')
os.system('rm -f *.o')
