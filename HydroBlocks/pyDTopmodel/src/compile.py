import os
os.system('rm -f *.mod')
os.system('rm -f *.pyf')
os.system('rm -f *.o')

#User defined info
#mkl_include = '/opt/intel/mkl/include'
#mkl_library = '/opt/intel/mkl/lib'
#mkl_include = '/opt/apps/intel/13/composer_xe_2013_sp1.1.106/mkl/include'
#mkl_library = '/opt/apps/intel/13/composer_xe_2013_sp1.1.106/mkl/lib/intel64'
mkl_include = '/opt/intel/composer_xe_2015.2.164/mkl/include'
mkl_library = '/opt/intel/composer_xe_2015.2.164/mkl/lib/intel64'
#Determine the conda directory
#dir = "/".join(os.popen('which python').read().split('/')[0:-2])
#mkl_include = '%s/include' % dir
#mkl_include = '%s/lib' % dir
#mkl_library = '%s/lib' % dir

#Compile the code to create the module
cmd = 'gfortran -c dynamic_topmodel_tools.f90 -lgomp -liomp5 -lmkl_rt -lpthread -lm -w -O3 -funroll-loops -m64 -Wl,--no-as-needed -openmp -fopenmp -g -fbacktrace -fpic -o dynamic_topmodel_tools.o -I%s -L%s' % (mkl_include,mkl_library)
os.system(cmd)

#Create a version with the include commented out
cmd = 'f2py dynamic_topmodel_tools_wrapper.f90 -h dynamic_topmodel_tools.pyf -m dynamic_topmodel_tools --overwrite-signature'
os.system(cmd)
#cmd = 'f2py -c dynamic_topmodel_tools.pyf wrapper.f90 *.o -lgomp -liomp5 -lmkl_rt -lpthread -lm --fcompiler=gnu95 --f90flags="-w -O3 -funroll-loops -m64 -Wl,--no-as-needed -openmp -fopenmp -g -fbacktrace"'
cmd = 'f2py -c dynamic_topmodel_tools.pyf dynamic_topmodel_tools_wrapper.f90 *.o -lgomp -lmkl_rt -L%s --fcompiler=gnu95' % (mkl_library)
os.system(cmd)

#Move to the previous directory
os.system('mv dynamic_topmodel_tools.so ../.')

#Remove the temporary file
os.system('rm -f tmp.f90')
os.system('rm -f *.mod')
os.system('rm -f *.pyf')
os.system('rm -f *.o')
