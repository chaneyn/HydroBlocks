import os
import upscaling_python
import datetime
import pickle

#Define the number of cores (1,4,16,...)
ncores = 5**2

#Define the parameters
startdate = datetime.datetime(2004,1,1,0)
enddate = datetime.datetime(2004,12,31,23)
nt_in = 24*((enddate - startdate).days+1)
dt = 24
res = 3
dir = '/home/freeze/nchaney/LittleWashitaRegion'
metadata = {'startdate':startdate,
            'enddate':enddate,
            'dt':dt, #hours
            'ncatch':9,
            'nt_in':nt_in,
            'nt_out':nt_in/dt,
            'dir':dir,
            'res':res, #arcsec
            'vars':['smc1','prcp','lh'],
            'output_dir':'%s/%darsec' % (dir,res),
            'dir':dir,}
             

#Empty out the directory if it exists
os.system('rm -rf %s' % metadata['output_dir'])
os.system('mkdir -p %s/workspace' % metadata['output_dir'])

#Save the metadata
pickle.dump(metadata,open('%s/workspace/metadata.pck' % metadata['output_dir'],'w'))

#Create the virtual rasters
upscaling_python.Create_Virtual_Rasters(metadata)

#Create the main files
upscaling_python.Create_Upscale_Template(metadata)

#Initialize the output files
upscaling_python.Initialize_Output_Files(metadata)

#Initialize the mpi processes
os.system('mpirun -n %d python driver_core.py %s' % (ncores,metadata['output_dir']))
