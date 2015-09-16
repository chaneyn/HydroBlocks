import os
mkl_flag = False
#Preprocessing:Tools
cwd = os.getcwd()
os.chdir('Preprocessing/Tools')
os.system('python compile.py')
os.chdir(cwd)
#HydroBloks:Dynamic TOPMODEL
cwd = os.getcwd()
os.chdir('HydroBloks/pyDTopmodel/src')
os.system('python compile.py')
os.chdir(cwd)
#HydroBloks:Noah-MP
cwd = os.getcwd()
os.chdir('HydroBloks/pyNoahMP')
os.system('python make.py')
os.chdir(cwd)
