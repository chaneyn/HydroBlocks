import os
#Compile NOAH
os.chdir('src/Utility_routines')
os.system('make clean')
os.system('make')
os.chdir('../Noah')
os.system('make clean')
os.system('make')
#Create signature file
cmd = 'f2py pyNoahMP.f90 -h NoahMP.pyf -m NoahMP --overwrite-signature'
os.system(cmd)
#Create driver library
#cmd = 'f2py -lgomp -c NoahMP.pyf *.o ../Utility_routines/*.o --fcompiler=gnu95 --f90flags="-w -fopenmp"'
#cmd = 'f2py -lgomp -c NoahMP.pyf *.o ../Utility_routines/*.o --fcompiler=gnu95 --f90flags="-w -fopenmp -Wall -g"'
cmd = 'f2py --debug -lgomp -c NoahMP.pyf *.o ../Utility_routines/*.o --fcompiler=gnu95 --f90flags="-w -fopenmp -g -Werror -fmodule-private -fimplicit-none -fbounds-check"'
os.system(cmd)
#Clean up both directories
#os.system('make clean')
#os.system('rm *.f90')
#os.system('rm NoahMP.pyf')
os.system('mv NoahMP.so ../../.')
