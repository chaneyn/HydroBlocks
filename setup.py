import os
#model:Noah-MP
cwd = os.getcwd()
os.chdir('model/pyNoahMP')
os.system('python make.py')
os.chdir(cwd)
