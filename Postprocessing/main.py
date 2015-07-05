import os
import upscaling_python
import datetime
import pickle

#Define the number of cores (1,4,16,...)
#ncores = 40**2
#ncores = 16**2
ncores = 5**2

#Define the parameters
mpi_type = 'openmpi'#'ibrun'
startdate = datetime.datetime(2004,1,1,0)
enddate = datetime.datetime(2006,12,31,23)
#enddate = datetime.datetime(2004,3,31,23)
nt_in = 24*((enddate - startdate).days+1)
dt = 1
res = 1#30
#dir = '/home/freeze/nchaney/LittleWashitaRegion'
dir = '/home/freeze/nchaney/HydroBloks/Test'
#dir = '/scratch/02179/chaneyna/conus30m/CONUS_SIMULATIONS/catchments'
#Determine the number of catchments
file = '%s/catchments.pck' % dir
ncatch = len(pickle.load(open(file)))
#dir = '/scratch/02179/chaneyna/catchments'
metadata = {'startdate':startdate,
            'enddate':enddate,
            'dt':dt, #hours
            'ncatch':ncatch,
            'nt_in':nt_in,
            'nt_out':nt_in/dt,
            'dir':dir,
            'res':res, #arcsec
            'vars':['smc1','prcp','lh'],#,'sh','g','qbase','qsurface'],
            'output_dir':'%s/%darcsec' % (dir,res),
            'dir':dir,
	    'nstripes':10,
            'type':'regional',
            'ncores':ncores}
             
python = '/home1/02179/chaneyna/libraries/python/lib/python2.7/site-packages/mpi4py/bin/python-mpi'

#Empty out the directory if it exists
'''os.system('rm -rf %s' % metadata['output_dir'])
os.system('mkdir -p %s/workspace' % metadata['output_dir'])
os.system('lfs setstripe -c %d %s' % (metadata['nstripes'],metadata['output_dir']))'''

#Save the metadata
pickle.dump(metadata,open('%s/workspace/metadata.pck' % metadata['output_dir'],'w'))

#Create the virtual rasters
'''upscaling_python.Create_Virtual_Rasters(metadata)

#Create the main files
upscaling_python.Create_Upscale_Template(metadata)'''

#Initialize the output files
#os.system('ibrun -np %d %s driver_core.py %s %s' % (31,python,metadata['output_dir'],'file_creation'))
#upscaling_python.Initialize_Output_Files(metadata)

#Initialize the mpi processes
#if mpi_type == 'openmpi':
# os.system('mpirun -n %d python driver_core.py %s' % (ncores,metadata['output_dir']))
#elif mpi_type == 'ibrun':
os.system('ibrun -np %d %s driver_core.py %s %s' % (ncores,python,metadata['output_dir'],'upscale'))

#Create the output files
os.system('ibrun -np %d %s driver_core.py %s %s' % (ncores,python,metadata['output_dir'],'file_creation'))
