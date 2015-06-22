import os
os.system('rm -f *.mod')
os.system('rm -f *.pyf')
os.system('rm -f *.o')

#User defined info
mkl_include = '/opt/intel/mkl/include'

#Compile the code to create the module
cmd = 'gfortran -c dynamic_topmodel_tools.f90 -lgomp -liomp5 -lmkl_rt -lpthread -lm -w -O3 -funroll-loops -m64 -Wl,--no-as-needed -openmp -fopenmp -g -fbacktrace -fpic -o dynamic_topmodel_tools.o -I%s' % (mkl_include,)
os.system(cmd)

#Create a version with the include commented out
#cmd = 'tail -n +12 dynamic_topmodel_tools.f90 | head -n -1 > tmp.f90'
#os.system(cmd)
cmd = 'f2py wrapper.f90 -h dynamic_topmodel_tools.pyf -m dynamic_topmodel_tools --overwrite-signature'
os.system(cmd)
#cmd = 'f2py -c dynamic_topmodel_tools.pyf wrapper.f90 *.o -lgomp -liomp5 -lmkl_rt -lpthread -lm --fcompiler=gnu95 --f90flags="-w -O3 -funroll-loops -m64 -Wl,--no-as-needed -openmp -fopenmp -g -fbacktrace"'
cmd = 'f2py -c dynamic_topmodel_tools.pyf wrapper.f90 *.o -lgomp -lmkl_rt --fcompiler=gnu95'
#cmd = 'f2py -c wrapper.f90 *.o -lgfortran -lgomp -liomp5 -lmkl_rt -lpthread -lm --fcompiler=gnu95 --f90flags="-w -O3 -funroll-loops -m64 -Wl,--no-as-needed -openmp -fopenmp -g -fbacktrace"'
os.system(cmd)
os.system('mv dynamic_topmodel_tools.so ../.')
exit()

#Create the signature file
cmd = 'f2py wrapper.f90 -h dynamic_topmodel_tools.pyf -m dynamic_topmodel_tools --overwrite-signature'# only: wrapper :'
#cmd = 'f2py tmp.f90 -h dynamic_topmodel_tools.pyf -m dynamic_topmodel_tools --overwrite-signature only: update :'
os.system(cmd)

#Create driver library
#cmd = 'f2py -c dynamic_topmodel_tools.pyf *.o -lgfortran -lgomp -liomp5 -lmkl_rt -lpthread -lm --fcompiler=gnu95 --f90flags="-w -O3 -funroll-loops -m64 -Wl,--no-as-needed -openmp -fopenmp -g -fbacktrace"'
cmd = 'f2py -c dynamic_topmodel_tools.pyf *.o'
#cmd = 'f2py -c dynamic_topmodel_tools.pyf *.o -lgfortran --fcompiler=gnu95 --f90flags="-fpic"'
#cmd = 'f2py tmp.f90 -lgomp -liomp5 -lmkl_rt -lpthread -lm -c dynamic_topmodel_tools.pyf --fcompiler=gnu95 --f90flags="-w -O3 -funroll-loops -m64 -Wl,--no-as-needed -openmp -fopenmp -g -fbacktrace"'
os.system(cmd)

#Move to the previous directory
os.system('mv dynamic_topmodel_tools.so ../.')

#Remove the temporary file
'''os.system('rm -f tmp.f90')
os.system('rm -f *.mod')
os.system('rm -f *.pyf')
os.system('rm -f *.o')'''
