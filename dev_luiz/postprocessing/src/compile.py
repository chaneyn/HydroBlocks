import os
cmd = 'f2py upscaling_tools.f90 -h upscaling_tools.pyf -m upscaling_fortran --overwrite-signature'
os.system(cmd)
#Create driver library
cmd = 'f2py upscaling_tools.f90 -lgomp -c upscaling_tools.pyf --fcompiler=gnu95 --f90flags="-w -fopenmp"'
os.system(cmd)
#Move to the previos directory
os.system('mv upscaling_fortran.so ../.')
