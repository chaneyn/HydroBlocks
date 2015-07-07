from mpi4py import MPI
import numpy as np
import json
import os
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
nhrus = np.arange(10,110,40)
workspace = '/scratch/02179/chaneyna/convergence_analysis/workspace'

#Iterate through the hru numbers
for nhru in nhrus[rank::size]:
 print "Rank %d || Running model for %d hrus" % (rank,nhru)
 nhru_nc = 2*nhru/3
 nhru_c = nhru/3
 diff = nhru - nhru_nc - nhru_c
 nhru_nc = nhru_nc + diff
 #Create a directory for this iteration
 dir = '%s/nhru_%d' % (workspace,nhru)
 os.system('mkdir -p %s' % dir)
 #Update the metadata file
 metadata = json.load(open('metadata_template.json')) 
 metadata['nhru_c'] = nhru_c
 metadata['nhru_nc'] = nhru_nc 
 metadata['input_file'] = '%s/input_file.nc' % dir
 metadata['output_file'] = '%s/output_file.nc' % dir
 metadata['soil_file'] = '%s/soil_file.txt' % dir
 metadata['workspace'] = '/scratch/02179/chaneyna/convergence_analysis/LittleWashita/workspace'
 #Save the metadata file
 file = '%s/metadata.json' % (dir,)
 json.dump(metadata,open(file,'w'))
 #Create the input files
 os.system('python ../Preprocessing/Driver.py %s' % file)
 #Run hydrobloks
 os.system('python ../Model/Driver.py %s' % file)
 
 

