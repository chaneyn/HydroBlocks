from mpi4py import MPI
import json
import fiona
import glob
import numpy as np
import os
import sys

def Read_Metadata_File(file):

 import json
 metadata = json.load(open(file))['HydroBlocks']
 mdp = json.load(open(file))['Preprocessing']
 metadata['cid'] = mdp['cid']

 return metadata

def Run_HydroBlocks(metadata_file,edir,cid,rdir):#,hb):

 import datetime
 from dateutil.relativedelta import relativedelta
 import sys
 sys.path.append("model")# % hb)
 import HydroBlocks

 #Read in the metadata file
 metadata = Read_Metadata_File(metadata_file)
 info = metadata
 info['ncores'] = 1
 #info['enddate']['year'] = 2017
 info['startdate']['year'] = 2015
 info['startdate']['month'] = 1#2
 info['startdate']['day'] = 1
 info['enddate']['year'] = 2017
 info['enddate']['month'] = 12
 info['enddate']['day'] = 31
 info['restart']['flag'] = True
 info['segment']['years_per_segment'] = 3
 info['dt'] = 3600#300
 #info['dt'] = 1800#3600#300
 info['cdir'] = '%s/%s' % (edir,cid)
 info['Qobs_file'] = '%s/obs/obs.pck' % rdir
 info['routing_file'] = '%s/%s/octopy.pck' % (edir,cid)
 info['input_file'] = '%s/%s/input_file_routing.nc' % (edir,cid)
 info['output'] = {"dir":"%s/output_data/%s" % (edir,cid),
     "vars":["smc1","trad","sh","lh","runoff","prcp","inundation","smc_root"],
     "routing_vars":["Q","A","Qc","Qf"]}
 info['restart'] = {"flag": True, 
                    "dir":"%s/restart_data/%d" % (edir,info['cid'])}
 info['dz'] = [0.05, 0.1, 0.1, 0.1, 0.2, 0.2, 0.25, 0.25, 0.25, 0.5]

 #Define idate and fdate
 idate = datetime.datetime(metadata['startdate']['year'],metadata['startdate']['month'],metadata['startdate']['day'],0)
 fdate = datetime.datetime(metadata['enddate']['year'],metadata['enddate']['month'],metadata['enddate']['day'],0) + datetime.timedelta(days=1)
 
 #Run the segments for the model
 sidate = idate
 sfdate = idate
 while sidate < fdate:
  sfdate = sidate + relativedelta(years=metadata['segment']['years_per_segment'])
  if sfdate > fdate: sfdate = fdate
  #Set the parameters
  info['idate'] = sidate
  info['fdate'] = sfdate 
  info['MPI'] = MPI
  #Run the model
  #Initialize
  HB = HydroBlocks.initialize(info)
  #Run the model
  HB.run(info)
  #Finalize
  HB.finalize()
  #Update initial time step
  sidate = sfdate

 #Delete HB and library
 del HB
 del HydroBlocks

 return

def run(comm,edir):

 #Get some general info
 rdir = "/".join(edir.split('/')[0:-2])
 size = int(comm.Get_size())
 rank = int(comm.Get_rank())

 #Get cid configuration
 dfile = '%s/input_data/shp/domain.shp' % rdir
 fp = fiona.open(dfile,'r')
 cids = np.array(range(1,len(list(fp))+1))
 fp.close()


 #Iterate per catchment
 for cid in cids[rank::size]:
  print('cid:',cid,flush=True)
  #Create output directory
  os.system('mkdir -p %s/output_data/%d' % (edir,cid))
  #os.chdir('%s' % (edir,))
  mdfile = '%s/%s/metadata.json' % (edir,cid)
  Run_HydroBlocks(mdfile,edir,cid,rdir)
