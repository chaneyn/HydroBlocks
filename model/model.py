from mpi4py import MPI
import json
import fiona
import glob
import numpy as np
import os
import sys
import time
import copy
import model.HydroBlocks as HydroBlocks
from model.pyRouting import routing as HBrouting

def Read_Metadata_File(file):

 metadata = json.load(open(file))

 return metadata

def Run_HydroBlocks(metadata,edir,cids,rdir):
 #print('Run HB',flush=True)
 import datetime
 from dateutil.relativedelta import relativedelta
 #import sys
 #print("%s/model" % hb)
 #sys.path.append("%s/model" % hb)
 comm = MPI.COMM_WORLD
 rank = comm.Get_rank()
 size = comm.Get_size()

 #Initialize HB container
 HBdb = {}

 #Initialize model
 for cid in cids[rank::size]:
  info = metadata
  info['ncores'] = 1
  info['cid'] = cid
  info['cdir'] = '%s/%s' % (edir,cid)
  info['routing_file'] = '%s/%s/octopy.pck' % (edir,cid)
  info['input_file'] = '%s/%s/input_file.nc' % (edir,cid)
  if metadata["routing_module"]["type"] == 'kinematic':
   info['output'] = {"dir":"%s/output_data/%s" % (edir,cid),
      "vars":info['output']['vars'],
      "routing_vars":info['output']['routing_vars']}
  else:
   info['output'] = {"dir":"%s/output_data/%s" % (edir,cid),
      "vars":info['output']['vars']}
  info['restart'] = {"flag":info['restart']['flag'],
                     "dir":"%s/restart_data/%d" % (edir,cid)}
  os.system('mkdir -p %s' % (info['restart']['dir']))

  #Define idate and fdate
  idate = datetime.datetime(metadata['startdate']['year'],metadata['startdate']['month'],metadata['startdate']['day'],0)
  fdate = datetime.datetime(metadata['enddate']['year'],metadata['enddate']['month'],metadata['enddate']['day'],0) + datetime.timedelta(days=1)
  
  sidate = idate
  sfdate = sidate + relativedelta(years=metadata['segment']['years_per_segment'])
  if sfdate > fdate: sfdate = fdate
  #Set the parameters
  info['idate'] = sidate
  info['fdate'] = sfdate
  #Determine cid/rank mapping
  #Initialize model
  print(' Initialize model',flush=True)
  #import model.HydroBlocks as HydroBlocks
  HBdb[cid] = HydroBlocks.initialize(copy.deepcopy(info),MPI)
  #del HydroBlocks

 #Determine cid/rank mapping
 determine_cid_rank_mapping(cids,HBdb,rank,size)

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
  #print(' Initialize model',flush=True)
  #HB = HydroBlocks.initialize(info)
  print(' Run the model',flush=True)
  date = sidate
  HBdb[cid].noahmp.dzwt[:] = 0.0
  i = 0
  tic = time.time()
  while date < sfdate:
   i=i+1
   run_timestep(cids,info,date,tic,HBdb,rank,size)
   #Update time step
   date = date + HBdb[cid].dt_timedelta
  print(' Finalize the model',flush=True)
  for cid in cids[rank::size]:
   HBdb[cid].finalize(info)
  #Update initial time step
  sidate = sfdate

 #Delete HB and library
 for cid in cids[rank::size]:
  del HBdb[cid]

 return

def run_timestep(cids,info,date,tic,HBdb,rank,size):

  #Update input data
  update_input(cids,HBdb,rank,size,date)

  #Calculate initial NOAH water balance
  initialize_water_balance(cids,rank,size,HBdb)

  #Update model
  update_model(cids,rank,size,date,HBdb)

  #Finalize water balance
  finalize_water_balance(cids,rank,size,HBdb)

  #Update the water balance error
  update_water_balance_error(cids,rank,size,HBdb)

  #Update the energy balance error
  update_energy_balance_error(cids,rank,size,HBdb)

  #Update time and date
  update_time_and_date(cids,rank,size,HBdb,info,date,tic)

  #Update output
  update_output(cids,rank,size,HBdb,date)

  #Update time step
  update_timestep(cids,rank,size,HBdb)

  #Output statistics
  update_statistics(cids,rank,size,HBdb,date)

  return

def update_statistics(cids,rank,size,HBdb,date):

  for cid in cids[rank::size]:
   string = '|%s|%s|%s|%s|%s|%s|%s|' % \
        ('Date:%s' % date.strftime("%Y-%m-%d_%H:%M"),\
         'Runtime:%.2f(s)'%(HBdb[cid].runtime),\
         'Acc_ET:%.2f(mm)'%HBdb[cid].acc_et,\
         'Acc_P:%.2f(mm)'%HBdb[cid].acc_prcp,\
         'Acc_Q:%.2f(mm)'%HBdb[cid].acc_q,\
         'Acc_WBE:%.2f(mm)' % HBdb[cid].acc_errwat,\
         'Acc_EBE:%.2f(J/m2)' % HBdb[cid].acc_erreng)
  print(string,flush=True)

  return

def update_timestep(cids,rank,size,HBdb):

  for cid in cids[rank::size]:
   #Update time step
   HBdb[cid].itime = HBdb[cid].itime + 1

  return

def update_output(cids,rank,size,HBdb,date):

  for cid in cids[rank::size]:
   HBdb[cid].update_output(date)

  return

def update_time_and_date(cids,rank,size,HBdb,info,date,tic):

  for cid in cids[rank::size]:
   #Update time and date
   HBdb[cid].date = date
   info['date'] = date
   HBdb[cid].runtime = HBdb[cid].runtime0 + time.time()-tic

  return

def update_energy_balance_error(cids,rank,size,HBdb):

  for cid in cids[rank::size]:
   HBdb[cid].calculate_energy_balance_error()

  return

def update_water_balance_error(cids,rank,size,HBdb):

  for cid in cids[rank::size]:
   HBdb[cid].calculate_water_balance_error()

  return

def finalize_water_balance(cids,rank,size,HBdb):

  for cid in cids[rank::size]:
   HBdb[cid].finalize_water_balance()

  return

def update_model(cids,rank,size,date,HBdb):

  #Update NoahMP
  for cid in cids[rank::size]:
   HBdb[cid].update_noahmp(date)
  
  #Update routing
  update_routing(cids,rank,size,HBdb)

  return

def update_routing(cids,rank,size,HBdb):

  #Update channel source/sink term
  update_channel_source_sink(cids,rank,size,HBdb)

  #Update routing scheme
  update_routing_scheme(cids,rank,size,HBdb)

  #Update hru inundation
  update_hru_inundation(cids,rank,size,HBdb)

  return

def update_hru_inundation(cids,rank,size,HBdb):

  for cid in cids[rank::size]:
   for hru in HBdb[cid].hrus:
    #Only allow fct of the inundated height to infiltrate
    HBdb[cid].routing.hru_inundation[hru] = HBdb[cid].routing.hband_inundation[HBdb[cid].hbands[hru]]

  return

def update_routing_scheme(cids,rank,size,HBdb):

  #Explicit solution
  flag_constant_Kvn = True

  #Initialize variables for the time step (this can be simplified)
  for cid in cids[rank::size]:
   HBdb[cid].routing.A0_org = np.copy(HBdb[cid].routing.A0)
   dt_routing = HBdb[cid].routing.dt_routing
   dt = HBdb[cid].routing.dt
   HBdb[cid].routing.nt = int(dt/dt_routing)
   nt = int(dt/dt_routing)

  for it in range(nt):

   #update routing (kinematic)
   update_routing_kinematic(cids,HBdb,rank,size,flag_constant_Kvn)

   #update routing (particle tracker)
   #update_routing_particle_tracker(cids,HBdb,rank,size)
   #exit()

  #Update local routing inundation
  HBrouting.calculate_routing_inundation(cids,HBdb,rank,size)

  return

def update_routing_particle_tracker(cids,HBdb,rank,size):
 
  #Send/receive velocity fields
  HBrouting.exchange_velocity_fields()

  #Push water down the channel network
  HBrouting.update_particle_tracker_macroscale_polygon()

  #Send/receive volume of water that enters/leaves each channel
  HBrouting.exchange_water_volumes()

  return

def update_routing_kinematic(cids,HBdb,rank,size,flag_constant_Kvn):

  #Exchange boundary conditions
  HBrouting.exchange_bcs_v3(cids,HBdb,rank,size)

  #Update local routing solution
  HBrouting.update_macroscale_polygon_routing(cids,HBdb,rank,size,flag_constant_Kvn)

  return

def update_channel_source_sink(cids,rank,size,HBdb):

  for cid in cids[rank::size]:
   HBdb[cid].update_channel_source_sink()

  return

def initialize_water_balance(cids,rank,size,HBdb):

  for cid in cids[rank::size]:
   HBdb[cid].initialize_water_balance()

  return

def update_input(cids,HBdb,rank,size,date):

  for cid in cids[rank::size]:
   HBdb[cid].update_input(date)

  return

def determine_cid_rank_mapping(cids,HBdb,rank,size):
 
  #Determine what rank has which cid
  for cid in cids[rank::size]:
   HBdb[cid].comm = HBdb[cid].MPI.COMM_WORLD
   HBdb[cid].size = HBdb[cid].comm.Get_size()
   HBdb[cid].rank = HBdb[cid].comm.Get_rank()

  if rank != 0:
   dest = 0
   db_ex = {}
   for cid in cids[rank::size]:
    self = HBdb[cid]
    db_ex[self.cid] = self.rank
   self.comm.send(db_ex,dest=dest,tag=11)
  elif rank == 0:
   db = {}
   for cid in cids[rank::size]:
    self = HBdb[cid]
    db[self.cid] = 0
   for i in range(1,self.size):
    db_ex = self.comm.recv(source=i,tag=11)
    for key in db_ex:
     db[key] = db_ex[key]
  #Wait until completed
  self.comm.Barrier()

  #Send the list to all the cores now
  if self.rank == 0:
   for i in range(1,self.size):
    self.comm.send(db,dest=i,tag=11)
  if self.rank != 0:
   db = self.comm.recv(source=0,tag=11)

  #Memorize links
  for cid in cids[rank::size]:
   self = HBdb[cid]
   self.cid_rank_mapping = db
   self.routing.cid_rank_mapping = db

  return

def run(comm,metadata_file):

 #Get some general info
 metadata = Read_Metadata_File(metadata_file)
 rdir = metadata['rdir']
 edir = '%s/experiments/simulations/%s' % (rdir,metadata['experiment'])
 size = int(comm.Get_size())
 rank = int(comm.Get_rank())

 #Get cid configuration
 dfile = '%s/data/shp/domain.shp' % rdir
 fp = fiona.open(dfile,'r')
 cids = np.array(range(1,len(list(fp))+1))
 fp.close()

 #Iterate per catchment
 #for cid in cids[rank::size]:
 #Create output directory
 #os.system('mkdir -p %s/output_data/%d' % (edir,cid))
 #Run_HydroBlocks(metadata,edir,cid,rdir)
 Run_HydroBlocks(metadata,edir,cids,rdir)
