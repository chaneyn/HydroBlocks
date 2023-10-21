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
 MPdb = HydroBlocks_MacroscalePolygons(comm,cids)

 #Initialize model
 for cid in MPdb.cids:
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
  MPdb.HBdb[cid] = HydroBlocks.initialize(copy.deepcopy(info),MPI)
  #del HydroBlocks

 #Determine cid/rank mapping
 determine_cid_rank_mapping(MPdb)

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
  MPdb.HBdb[cid].noahmp.dzwt[:] = 0.0
  i = 0
  tic = time.time()
  while date < sfdate:
   i=i+1
   run_timestep(MPdb.cids,info,date,tic,MPdb.mpi_rank,MPdb.mpi_size,MPdb.HBdb)
   #Update time step
   date = date + MPdb.HBdb[cid].dt_timedelta
  print(' Finalize the model',flush=True)
  for cid in MPdb.cids:
   MPdb.HBdb[cid].finalize(info)
  #Update initial time step
  sidate = sfdate

 #Delete HB and library
 for cid in MPdb.cids:
  del MPdb.HBdb[cid]

 return

class HydroBlocks_MacroscalePolygons:
    
 def __init__(self,comm,cids):
    
  #Save MPI information
  self.mpi_comm = MPI.COMM_WORLD
  self.mpi_rank = comm.Get_rank()
  self.mpi_size = comm.Get_size()

  #Learn cids for the given container
  ncids = int(np.floor(cids.size/self.mpi_size))
  remainder = cids.size % self.mpi_size
  imax = 0
  for rank in range(self.mpi_size):
   if rank <= remainder-1:ncids = int(np.floor(cids.size/self.mpi_size))+1
   else: ncids = int(np.floor(cids.size/self.mpi_size))
   imin = imax
   imax = imin+ncids
   if rank == self.mpi_rank:
    cids_core = cids[imin:imax]
    break
  self.cids = cids_core

  #Initialize the dictionary for the HB macroscale polygons
  self.HBdb = {}

  return

def run_timestep(cids,info,date,tic,rank,size,HBdb):

  #Update input data
  update_input(cids,HBdb,date)

  #Calculate initial NOAH water balance
  initialize_water_balance(cids,HBdb)

  #Update model
  update_model(cids,rank,size,date,HBdb)

  #Finalize water balance
  finalize_water_balance(cids,HBdb)

  #Update the water balance error
  update_water_balance_error(cids,HBdb)

  #Update the energy balance error
  update_energy_balance_error(cids,HBdb)

  #Update time and date
  update_time_and_date(cids,HBdb,info,date,tic)

  #Update output
  update_output(cids,HBdb,date)

  #Update time step
  update_timestep(cids,HBdb)

  #Output statistics
  update_statistics(cids,HBdb,date)

  return

def update_statistics(cids,HBdb,date):

  for cid in cids:
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

def update_timestep(cids,HBdb):

  for cid in cids:
   #Update time step
   HBdb[cid].itime = HBdb[cid].itime + 1

  return

def update_output(cids,HBdb,date):

  for cid in cids:
   HBdb[cid].update_output(date)

  return

def update_time_and_date(cids,HBdb,info,date,tic):

  for cid in cids:
   #Update time and date
   HBdb[cid].date = date
   info['date'] = date
   HBdb[cid].runtime = HBdb[cid].runtime0 + time.time()-tic

  return

def update_energy_balance_error(cids,HBdb):

  for cid in cids:
   HBdb[cid].calculate_energy_balance_error()

  return

def update_water_balance_error(cids,HBdb):

  for cid in cids:
   HBdb[cid].calculate_water_balance_error()

  return

def finalize_water_balance(cids,HBdb):

  for cid in cids:
   HBdb[cid].finalize_water_balance()

  return

def update_model(cids,rank,size,date,HBdb):

  #Update NoahMP
  for cid in cids:
   HBdb[cid].update_noahmp(date)
  
  #Update routing
  if HBdb[cid].routing_module == 'kinematic':update_routing(cids,rank,size,HBdb)

  return

def update_routing(cids,rank,size,HBdb):

  #Update channel source/sink term
  update_channel_source_sink(cids,HBdb)

  #Update routing scheme
  update_routing_scheme(cids,rank,size,HBdb)

  #Update hru inundation
  update_hru_inundation(cids,HBdb)

  return

def update_hru_inundation(cids,HBdb):

  for cid in cids:
   for hru in HBdb[cid].hrus:
    #Only allow fct of the inundated height to infiltrate
    HBdb[cid].routing.hru_inundation[hru] = HBdb[cid].routing.hband_inundation[HBdb[cid].hbands[hru]]

  return

def update_routing_scheme(cids,rank,size,HBdb):

  #Explicit solution
  flag_constant_Kvn = True

  #Initialize variables for the time step (this can be simplified)
  for cid in cids:
   HBdb[cid].routing.A0_org = np.copy(HBdb[cid].routing.A0)
   dt_routing = HBdb[cid].routing.dt_routing
   dt = HBdb[cid].routing.dt
   HBdb[cid].routing.nt = int(dt/dt_routing)
   nt = int(dt/dt_routing)
   HBdb[cid].routing.Q0sum = np.zeros(HBdb[cid].routing.Q0.size)

  for it in range(nt):

   #update routing (kinematic)
   update_routing_kinematic(cids,HBdb,rank,size,flag_constant_Kvn)

   for cid in cids:
    HBdb[cid].routing.Q0sum += HBdb[cid].routing.Q0

   #update routing (particle tracker)
   #update_routing_particle_tracker(cids,HBdb,rank,size)
   #exit()

  #Compute Q0 as average of all sub Q0
  for cid in cids:
   HBdb[cid].routing.Q0[:] = HBdb[cid].routing.Q0sum/nt

  #Update local routing inundation
  HBrouting.calculate_routing_inundation(cids,HBdb)

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
  HBrouting.update_macroscale_polygon_routing(cids,HBdb,flag_constant_Kvn)

  #Calculate water balance
  #calculate_river_water_balance(cids,HBdb)

  #self.Q1[:] = self.Q0[:]

  return

def calculate_river_water_balance(cids,HBdb):

  for cid in cids:
   dt = HBdb[cid].routing.dt_routing
   ic_outlets = HBdb[cid].routing.outlets_array[:,1]
   Qout_outlets = HBdb[cid].routing.Q0[ic_outlets]
   ic_inlets = HBdb[cid].routing.inlets_array[:,1]
   Qin_inlets = HBdb[cid].routing.bcs[ic_inlets]
   Qin_overland = HBdb[cid].routing.qss*HBdb[cid].routing.c_length
   Vin = dt*np.sum(Qin_inlets) #Volume of water from flow into inlets
   Vout = dt*np.sum(Qout_outlets) #Volume of water from flow out of outlets
   V1 = np.sum(HBdb[cid].routing.A1*HBdb[cid].routing.c_length) #Volume of water in channels at t1
   V0 = np.sum(HBdb[cid].routing.A0*HBdb[cid].routing.c_length) #Volume of water in channels at t0 
   Vin_overland = dt*np.sum(Qin_overland) #Volume of water from overland flow input
   #V1 = V0 + Vin + Vin_overland - Vout
   print(cid,V1-(V0+Vin-Vout+Vin_overland))
   
  return

def update_channel_source_sink(cids,HBdb):

  for cid in cids:
   HBdb[cid].update_channel_source_sink()

  return

def initialize_water_balance(cids,HBdb):

  for cid in cids:
   HBdb[cid].initialize_water_balance()

  return

def update_input(cids,HBdb,date):

  for cid in cids:
   HBdb[cid].update_input(date)

  return

def determine_cid_rank_mapping(MPdb):
 
  cids = MPdb.cids
  HBdb = MPdb.HBdb
  rank = MPdb.mpi_rank

  #Determine what rank has which cid
  for cid in cids:
   HBdb[cid].comm = MPdb.mpi_comm
   HBdb[cid].size = MPdb.mpi_size
   HBdb[cid].rank = MPdb.mpi_rank

  if rank != 0:
   dest = 0
   db_ex = {}
   for cid in cids:
    self = HBdb[cid]
    db_ex[self.cid] = self.rank
   self.comm.send(db_ex,dest=dest,tag=11)
  elif rank == 0:
   db = {}
   for cid in cids:
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
  for cid in cids:
   self = HBdb[cid]
   self.cid_rank_mapping = db
   #self.routing.cid_rank_mapping = db

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
