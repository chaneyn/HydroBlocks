import os
#Compile NOAH
os.chdir('src/Noah')
os.system('make')
#Create signature file
cmd = 'f2py pyNoahMP.f90 -h NoahMP.pyf -m NoahMP --overwrite-signature'
os.system(cmd)
#Create driver library
cmd = 'f2py -lgomp -c NoahMP.pyf *.o ../Utility_routines/*.o --fcompiler=gnu95 --f90flags="-w -fopenmp"'
os.system(cmd)
#Clean up both directories
#os.system('make clean')
#os.system('rm *.f90')
os.system('rm NoahMP.pyf')
os.system('mv NoahMP.so ../../.')
