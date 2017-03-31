import os
mkl_flag = False
#Preprocessing:Tools
cwd = os.getcwd()
os.chdir('Preprocessing/Tools')
os.system('python compile.py')
os.chdir(cwd)
#HydroBlocks:Dynamic TOPMODEL
cwd = os.getcwd()
os.chdir('HydroBlocks/pyDTopmodel/src')
os.system('python compile.py')
os.chdir(cwd)
#HydroBlocks:Noah-MP
cwd = os.getcwd()
os.chdir('HydroBlocks/pyNoahMP')
os.system('python make.py')
os.chdir(cwd)
