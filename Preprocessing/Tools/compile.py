import os
cmd = 'f2py model_tools.f90 -h model_tools.pyf -m model_tools --overwrite-signature'
os.system(cmd)
#Create driver library
cmd = 'f2py model_tools.f90 -lgfortran -lgomp -c model_tools.pyf --fcompiler=gnu95 --f90flags="-w -fopenmp"'
os.system(cmd)
#Move to the previos directory
os.system('mv model_tools.so ../.')
