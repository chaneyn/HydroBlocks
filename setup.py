import os
#HydroBlocks:Noah-MP
cwd = os.getcwd()
os.chdir('HydroBlocks/pyNoahMP')
os.system('python make.py')
os.chdir(cwd)
