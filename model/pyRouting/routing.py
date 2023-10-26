import os
import glob
import json
import netCDF4 as nc
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
import scipy.sparse as sparse
import scipy.sparse.linalg
import pickle
import copy
import numba
import time
import sys
import h5py

class particle_tracker:

 def __init__(self,MPI,cid,dt,nhband,nhru,cdir,dt_routing,HB):

  self.itime = 0
  self.comm = MPI.COMM_WORLD
  self.size = self.comm.Get_size()
  self.rank = self.comm.Get_rank()
  self.name = MPI.Get_processor_name()

  #Read in all the relevant data (this should be moved to the preprocessor)
  self.db = pickle.load(open('%s/octopy.pck' % cdir ,'rb'))

  #Define the variables
  self.dt = dt
  self.dt_routing = dt_routing
  self.cid = cid
  self.uhs = self.db['uhs'][:] #unit hydrograph
  self.uh_travel_distance = self.db['uh_travel_time'] #travel distance per hband (travel time is bug)
  self.c_length = self.db['c_length'][:]
  self.c_n = self.db['c_n'][:]
  self.fp_n = self.db['c_n'][:] #NWC - 10/2023 -> set to the same for now, model too "jumpy"
  m =  self.fp_n < self.c_n
  self.fp_n[m] = self.c_n[m]
  self.c_slope = self.db['c_slope'][:]
  self.c_slope[self.c_slope < 10**-3] = 10**-3
  self.c_bankfull = self.db['c_bankfull'][:]
  self.nchannel = self.c_length.size
  self.A0 = 10**-5*np.ones(self.c_length.size)
  self.u0 = 10**-5*np.ones(self.c_length.size)
  self.A1 = np.zeros(self.c_length.size)
  self.A1[:] = self.A0[:]
  self.bcs = np.zeros(self.c_length.size)
  self.qss = np.zeros(self.c_length.size)
  self.qss_remainder = np.zeros(self.c_length.size)
  self.dA = np.zeros(self.c_length.size)
  self.Q0 = np.zeros(self.c_length.size)
  self.Q1 = np.zeros(self.c_length.size)
  self.Vin = np.zeros(self.c_length.size)
  self.Vout = np.zeros(self.c_length.size)
  self.Qin = np.zeros(self.c_length.size)
  self.Qout = np.zeros(self.c_length.size)
  self.Kvn0 = np.zeros(self.c_length.size)
  self.reach2hband = np.asarray(self.db['reach2hband'].todense())
  self.reach2hband_inundation = np.copy(self.reach2hband)
  self.reach2hband_inundation[:] = 0.0
  #self.hdw = self.db['hdw']
  fp = nc.Dataset('%s/input_file.nc' % cdir)
  #Read in downstream channels database
  self.downstream_channels = fp['stream_network']['downstream_channels'][:]
  #self.downstream_channels[:,:,5:] = -9999
  '''if cid == 5:
   ucids = np.unique(self.downstream_channels[:,1,:])
   ucids = ucids[(ucids > 0) & (ucids != cid)]
   for ucid in ucids:
    m = self.downstream_channels[:,1,:] == ucid
    channels = np.unique(self.downstream_channels[:,0,:][m])
    print(ucid,channels.size)
  exit()'''
  #cids to send to 
  tmp = np.unique(fp['stream_network']['outlets'][:,2])
  tmp = tmp[tmp != -9999]
  self.scids_hdw = tmp
  #cids to receive from
  tmp = np.unique(fp['stream_network']['inlets'][:,2:6])
  tmp = tmp[tmp != -9999]
  self.rcids_hdw = tmp
  self.outlets_array = fp['stream_network']['outlets'][:]
  self.inlets_array = fp['stream_network']['inlets'][:]
  #topology
  topology = fp['stream_network']['topology'][:]
  tmp = -9999*np.ones((topology.size,5),dtype=np.int32)
  for i in range(topology.size):
   m = topology == i
   if np.sum(m) >= 1:
    tmp[i,:np.sum(m)] = np.where(m)[0]
  self.topology = tmp[:]
  fp.close()
  self.hdb = self.db['hdb']
  self.hdb['hband'] = self.hdb['hband'].astype(np.int64)
  self.hdb['M'][self.hdb['M'] == 0] = 10**-5
  self.c_width = self.hdb['W'][:,0]#self.db['c_width'][:]
  self.A0_org = np.copy(self.A0)
  self.u0_org = np.copy(self.u0)
  self.Ac0 = np.zeros(self.c_length.size)
  self.Ac1 = np.zeros(self.c_length.size)
  self.dAc0 = np.zeros(self.c_length.size)
  self.dAc1 = np.zeros(self.c_length.size)
  self.Af0 = np.zeros(self.c_length.size)
  self.Af1 = np.zeros(self.c_length.size)
  self.dAf0 = np.zeros(self.c_length.size)
  self.dAf1 = np.zeros(self.c_length.size)
  self.Ac = np.zeros(self.c_length.size)
  self.Af = np.zeros(self.c_length.size)
  self.Qc = np.zeros(self.c_length.size)
  self.Qf = np.zeros(self.c_length.size)
  self.hru_inundation = np.zeros(nhru)
  self.hru_runoff_inundation = np.zeros(nhru)
  self.hband_inundation = np.zeros(nhband)
  self.hband_inundation1 = np.zeros(nhband)
  self.fct_infiltrate = 0.1

  #Define channel hbands per reach
  self.hband_channel = np.zeros(self.c_length.size).astype(np.int32)
  m = (self.hdb['hand'] == 0) & (self.hdb['W'] > 0)
  self.hband_channel[:] = self.hdb['hband'][m]

  return

class kinematic:

 def __init__(self,MPI,cid,dt,nhband,nhru,cdir,dt_routing,HB):

  self.itime = 0
  self.comm = MPI.COMM_WORLD
  self.size = self.comm.Get_size()
  self.rank = self.comm.Get_rank()
  self.name = MPI.Get_processor_name()
  #self.cid_rank_mapping = cid_rank_mapping

  #Read in all the relevant data (this should be moved to the preprocessor)
  self.db = pickle.load(open('%s/octopy.pck' % cdir ,'rb'))

  #Define the variables
  self.dt = dt
  self.dt_routing = dt_routing
  self.cid = cid
  self.uhs = self.db['uhs'][:] #unit hydrographs
  self.uh_travel_distance = self.db['uh_travel_time'] #travel distance per hband (travel time is bug)
  self.c_length = self.db['c_length'][:]
  self.c_n = self.db['c_n'][:]
  #self.fp_n = self.db['fp_n'][:]
  self.fp_n = self.db['c_n'][:] #NWC - 10/2023 -> set to the same for now, model too "jumpy"
  #self.c_n[:] = 0.035#0.05#1#0.05#0.03
  #self.fp_n[:] = 0.15#0.15#25#0.15#0.15
  m =  self.fp_n < self.c_n
  self.fp_n[m] = self.c_n[m]
  self.c_slope = self.db['c_slope'][:]
  self.c_slope[self.c_slope < 10**-3] = 10**-3
  self.c_bankfull = self.db['c_bankfull'][:]
  self.nchannel = self.c_length.size
  self.A0 = 10**-5*np.ones(self.c_length.size)
  self.u0 = 10**-5*np.ones(self.c_length.size)
  self.A1 = np.zeros(self.c_length.size)
  self.bcs = np.zeros(self.c_length.size)
  #self.qin = np.zeros(self.c_length.size)
  #self.qout = np.zeros(self.c_length.size)
  self.qss = np.zeros(self.c_length.size)
  self.qss_remainder = np.zeros(self.c_length.size)
  self.dA = np.zeros(self.c_length.size)
  self.Q0 = np.zeros(self.c_length.size)
  self.Q1 = np.zeros(self.c_length.size)
  self.Kvn0 = np.zeros(self.c_length.size)
  self.reach2hband = np.asarray(self.db['reach2hband'].todense())
  self.reach2hband_inundation = np.copy(self.reach2hband)
  self.reach2hband_inundation[:] = 0.0
  #self.scids_hdw = self.db['scids_hdw']
  #self.rcids_hdw = self.db['rcids_hdw']
  #self.hdw = self.db['hdw']
  fp = nc.Dataset('%s/input_file.nc' % cdir)
  #cids to send to 
  tmp = np.unique(fp['stream_network']['outlets'][:,2])
  tmp = tmp[tmp != -9999]
  self.scids_hdw = tmp
  #cids to receive from
  tmp = np.unique(fp['stream_network']['inlets'][:,2:6])
  tmp = tmp[tmp != -9999]
  self.rcids_hdw = tmp
  self.outlets_array = fp['stream_network']['outlets'][:]
  self.inlets_array = fp['stream_network']['inlets'][:]
  #topology
  topology = fp['stream_network']['topology'][:]
  tmp = -9999*np.ones((topology.size,5),dtype=np.int32)
  for i in range(topology.size):
   m = topology == i
   if np.sum(m) >= 1:
    tmp[i,:np.sum(m)] = np.where(m)[0]
  self.topology = tmp[:]
  fp.close()
  self.hdb = self.db['hdb']
  self.hdb['hband'] = self.hdb['hband'].astype(np.int64)
  self.hdb['M'][self.hdb['M'] == 0] = 10**-5
  self.c_width = self.hdb['W'][:,0]#self.db['c_width'][:]
  self.A0_org = np.copy(self.A0)
  self.u0_org = np.copy(self.u0)
  self.Ac0 = np.zeros(self.c_length.size)
  self.Ac1 = np.zeros(self.c_length.size)
  self.dAc0 = np.zeros(self.c_length.size)
  self.dAc1 = np.zeros(self.c_length.size)
  self.Af0 = np.zeros(self.c_length.size)
  self.Af1 = np.zeros(self.c_length.size)
  self.dAf0 = np.zeros(self.c_length.size)
  self.dAf1 = np.zeros(self.c_length.size)
  self.Ac = np.zeros(self.c_length.size)
  self.Af = np.zeros(self.c_length.size)
  self.Qc = np.zeros(self.c_length.size)
  self.Qf = np.zeros(self.c_length.size)
  self.hru_inundation = np.zeros(nhru)
  self.hru_runoff_inundation = np.zeros(nhru)
  self.hband_inundation = np.zeros(nhband)
  self.hband_inundation1 = np.zeros(nhband)
  self.fct_infiltrate = 0.1

  #Define channel hbands per reach
  self.hband_channel = np.zeros(self.c_length.size).astype(np.int32)
  m = (self.hdb['hand'] == 0) & (self.hdb['W'] > 0)
  self.hband_channel[:] = self.hdb['hband'][m]

  #Pause until all cores have all their data
  #self.comm.Barrier()
	
  return

 def update(self,dt):

  #Explicit solution
  #self.dt_routing = 100 #seconds
  flag_constant_Kvn = False#True
  dt_routing = self.dt_routing
  nt = int(dt/dt_routing)
  A0_org = np.copy(self.A0)
  for it in range(nt):
   #Exchange boundary conditions
   #self.exchange_bcs()
   #self.exchange_bcs_v2()
   #Update solution
   (self.Q0,self.u0,self.A1) = update_solution_explicit(self.c_slope,self.c_n,self.u0,self.A0,self.topology,self.c_length,self.qss,self.bcs,self.dt_routing,self.fp_n,self.hdb['Ac'],self.hdb['Af'],self.hdb['Pc'],self.hdb['Pf'],self.hdb['W'],self.hdb['M'],A0_org,flag_constant_Kvn)
   #self.A0[:] = self.A1[:]
   #self.Q1[:] = self.Q0[:]

  #Zero out qss
  self.qss[:] = 0.0
  #Calculate area-weighted average inundation height per hband
  A = self.hdb['Af'] + self.hdb['Ac']
  #Add channel and flood cross sectional area
  A1 = self.A1[:]
  W = self.hdb['W']
  M = self.hdb['M']
  hand = self.hdb['hand']
  hband = self.hdb['hband']
  (self.hband_inundation[:],self.reach2hband_inundation[:]) = calculate_inundation_height_per_hband(A,A1,W,M,hand,hband,self.reach2hband_inundation,self.reach2hband)
  tmp = np.sum(self.reach2hband*self.reach2hband_inundation,axis=1)/self.c_length
  tmp1 = np.max(np.abs(A1-tmp))
  if tmp1 > 10**-3:
   argmax = np.argmax(np.abs(A1-tmp))
   print(self.itime,np.abs(A1-tmp)[argmax],A1[argmax],tmp[argmax],flush=True)

  #Calculate Qc and Qf
  self.Ac[:] = self.reach2hband[range(self.nchannel),self.hband_channel]*self.reach2hband_inundation[range(self.nchannel),self.hband_channel]/self.c_length
  self.Af[:] = self.A1 - self.Ac
  self.Qc[:] = self.u0*self.Ac
  self.Qf[:] = self.u0*self.Af
  self.Qc[self.Qc < 0.0] = 0.0
  self.Qf[self.Qf < 0.0] = 0.0

  return

 def exchange_bcs_v2(self,):

  #Send headwater data
  crm = self.cid_rank_mapping #Where each cid resides
  for ucid in self.scids_hdw:
   dest = crm[ucid]
   m = self.outlets_array[:,2] == ucid
   channels_ucid = self.outlets_array[m,3].data
   channels_cid = self.outlets_array[m,1].data
   Q0_bcs_ucid = self.Q0[channels_cid]
   db_ex = {'cid':self.cid,'scid_hdw':ucid,
          'channels_ucid':channels_ucid,
          'Q0_bcs_ucid':Q0_bcs_ucid}
   self.comm.send(db_ex,dest=dest,tag=11)

  #Receive headwater data
  recv = {}
  for ucid in self.rcids_hdw:
   if ucid not in recv:recv[ucid] = {}
   source = crm[ucid]
   db_ex = self.comm.recv(source=source,tag=11)
   for var in db_ex:
    recv[ucid][var] = db_ex[var]

  #Update the boundary conditions
  bcs = self.bcs
  self.bcs[:] = 0.0
  for ucid in self.rcids_hdw:
   channels_ucid = recv[ucid]['channels_ucid']
   Q0_bcs_ucid = recv[ucid]['Q0_bcs_ucid']
   self.bcs[channels_ucid] += Q0_bcs_ucid

  #Wait until all are done
  self.comm.Barrier()

  return

@numba.jit(nopython=True,cache=True)
def update_solution_explicit(c_slope,c_n,u0,A0,topology,c_length,qss,bcs,dt_routing,
     fp_n,Ac,Af,Pc,Pf,W,M,A0_org,flag_constant_Kvn,Kvn0):

 #Extract info
 #bcs_c = bcs/c_length
 maxu = 10.0
 minu = 10**-5
 dt = dt_routing

 if flag_constant_Kvn == True:
  #Determine velocity
  #Kvn = calculate_compound_convenyance(Ac,Af,Pc,Pf,W,M,A0_org,c_n,fp_n)
  Kvn = A0_org**(5.0/3.0)/(W[:,0])**(2.0/3.0)/c_n
  #Kvn = (0.2*Kvn + 0.8*Kvn0) #Time smooth conveyance to avoid jumpiness
  u0 = np.zeros(Kvn.size)
  u0[A0_org > 0.0] = Kvn[A0_org > 0.0]*c_slope[A0_org > 0.0]**0.5/A0_org[A0_org > 0.0]
 else:
  #Determine velocity
  #Kvn = calculate_compound_convenyance(Ac,Af,Pc,Pf,W,M,A0,c_n,fp_n)
  Kvn = A0**(5.0/3.0)/(W[:,0])**(2.0/3.0)/c_n
  #Kvn = (0.2*Kvn + 0.8*Kvn0) #Time smooth conveyance to avoid jumpiness
  u0 = np.zeros(Kvn.size)
  u0[A0 > 0.0] = Kvn[A0 > 0.0]*c_slope[A0 > 0.0]**0.5/A0[A0 > 0.0]


 #Constrain velocity
 u0[u0 > maxu] = maxu
 u0[u0 < minu] = minu
 u0[:] = 2.0

 #Compute Q0in
 Q0in = Compute_Q0in(topology,u0,A0)

 #Compute Q0out
 Q0out = A0*u0

 #A1 = A0 + source/sink + boundary conditions - Qout + Qin
 A1 = A0 + dt*qss + dt*bcs/c_length - dt*(u0*A0)/c_length + dt*Q0in/c_length

 #Curate A1 (issues with conservation of mass)
 #A1[A1 < 0] = 0.0

 return (Q0out,u0,A1,Kvn,Q0in)


@numba.jit(nopython=True,cache=True,nogil=True,fastmath=True)
def Compute_Q0in(topology,u0,A0):

 Q0in = np.zeros(A0.size)
 for i in range(Q0in.size):
  for j in range(topology.shape[1]):
   if topology[i,j] == -9999:continue
   else: Q0in[i] += u0[topology[i,j]]*A0[topology[i,j]]

 return Q0in

@numba.jit(nopython=True,cache=True)
def calculate_hydraulic_radius(A,P,W,A1):

 Rh = np.zeros(A.shape[0])
 for i in range(A.shape[0]):
  A0 = 0.0
  P0 = 0.0
  for j in range(A.shape[1]):
   W1 = W[i,j]
   if A[i,j] > A1[i]:break
   if A[i,j+1] == 0.0:break
   A0 = A[i,j]
   P0 = P[i,j]
  #Calculate the height above the segment
  h = (A1[i] - A0)/np.sum(W[i,0:j+1])
  #Calculate the wetted perimeter
  P1 = P0 + 2*h + W1
  #Calculate hydraulic radius
  Rh[i] = A1[i]/P1

 return Rh

@numba.jit(nopython=True,cache=True,nogil=True,fastmath=True)
def calculate_compound_convenyance(Ac,Af,Pc,Pf,W,M,A1,cn,fpn):
 
 Kvn = np.zeros(Ac.shape[0])
 #Determine the level for which we need to back out info
 for i in range(Ac.shape[0]):
  for j in range(Ac.shape[1]-1):
   if (Af[i,j+1]+Ac[i,j+1]) > A1[i]:break
   if (Af[i,j+1]+Ac[i,j+1]) == 0.0:break
  if j == 0:
   A0 = Ac[i,j]
   P0 = Pc[i,j]
   #Calculate the height above the segment
   h = (A1[i] - A0)/W[i,0]
   #Calculate the wetted perimeter
   P1 = P0 + 2*h + W[i,0]
   #Calculate compound conveyance
   Kvn[i] = A1[i]**(5.0/3.0)/P1**(2.0/3.0)/cn[i]
  else:
   #Calculate the height above the segment (quadratic equation)
   c = -(A1[i] - (Af[i,j] + Ac[i,j]))
   b = np.sum(W[i,0:j])
   a = 1.0/M[i,j]
   h = (-b + (b**2.0 - 4.0*a*c)**0.5)/(2.0*a) 
   #h2 = (-b - (b**2.0 - 4.0*a*c)**0.5)/(2.0*a) 
   #print('inside1:',h1,h2,A1[i],(Af[i,j] + Ac[i,j]),a,b,c,M[i,j])
   #Calculate channel cross sectional area
   Ac1 = Ac[i,1] + W[i,0]*h
   #Calculate floodplain cross sectional area
   Af1 = Af[i,1] + np.sum(W[i,1:j])*h + h*h/M[i,j]
   #Calculate channel wetted perimeter
   Pc1 = Pc[i,1]
   #Calculate floodplain wetted perimieter
   Pf1 = Pf[i,j] + 2*(h**2 + (h/M[i,j])**2)**0.5
   #print(h,M[i,j],Pc1,Pf1,cn[i],fpn[i])
   #Calculate conveyances
   Kvnc = (1.0/cn[i])*(Ac1**(5.0/3.0)/Pc1**(2.0/3.0))
   Kvnf = (1.0/fpn[i])*(Af1**(5.0/3.0)/Pf1**(2.0/3.0))
   #Calculate compound conveyance
   Kvn[i] = Kvnc + Kvnf
   if np.isnan(Kvn[i]) == 1:
    print('inside',A1[i],Kvn[i],a,b,c,h,Ac1,Af1,Pc1,Pf1,Kvnc,Kvnf)

 return Kvn

@numba.jit(nopython=True,cache=True)
def calculate_inundation_height_per_hband(A,A1,W,M,hand,hband,reach2hband_inundation,reach2hband):

 #Zero out reach2hband_inundation
 reach2hband_inundation[:] = 0.0
 #Determine the inundation height per hand value per basin
 for i in range(A.shape[0]):
  A0 = 0.0
  for j in range(1,A.shape[1]):
   if A[i,j] > A1[i]:break
   if A[i,j+1] == 0.0:break
   A0 = A[i,j]
  #Calculate the height above the segment
  if j == 1:htop = (A1[i] - A0)/W[i,0]
  else:
   c = -(A1[i] - A0)
   b = np.sum(W[i,0:j-1])
   a = 1.0/M[i,j-1]
   htop = (-b + (b**2.0 - 4.0*a*c)**0.5)/(2.0*a) #trapezoidal
  #Based on htop calculate the water heights (Wrong)
  h = np.zeros(j)
  h[0] = hand[i,j-1] + htop
  if j > 1:h[j-1] = htop*htop/M[i,j-1]/W[i,j-1]
  if j > 2:
   h[1:j-1] = hand[i,j-1] - hand[i,1:j-1] - (hand[i,2:j]-hand[i,1:j-1])/2 + htop
  #Add to hband sum (times area)
  idxs = hband[i,0:j]
  #Place the inundation level
  for k in range(idxs.size):
   idx = idxs[k]
   reach2hband_inundation[i,idx] = h[k]
 areas = np.sum(reach2hband,axis=0)
 hband_inundation = np.zeros(areas.size)	
 hband_inundation[:] = np.sum(reach2hband*reach2hband_inundation,axis=0)/np.sum(reach2hband,axis=0)

 return (hband_inundation,reach2hband_inundation)

@numba.jit(nopython=True,cache=True,nogil=True,fastmath=True)
def compute_qss(reach2hband,crunoff,c_length,qss):
  for i in range(qss.size):
   qss[i] = np.sum(reach2hband[i,:]*crunoff/1000.0)/c_length[i] #m2/s
  #qss[:] += reach2hband.dot(crunoff/1000.0)/c_length #m2/s
  return qss

def exchange_bcs_v3(cids,hbdb,rank,size):

  cids_core = cids
  db_ex_local = {} #exchange between cids on the same process (minimize MPI overhead)

  #Send headwater data
  for cid in cids_core:
   self = hbdb[cid].routing
   crm = hbdb[cid].cid_rank_mapping #Where each cid resides
   db_ex_local[cid] = {}
   for ucid in self.scids_hdw:
    dest = crm[ucid]
    m = self.outlets_array[:,2] == ucid
    channels_ucid = self.outlets_array[m,3].data
    channels_cid = self.outlets_array[m,1].data
    Q0_bcs_ucid = self.Q0[channels_cid]
    db_ex = {'cid':self.cid,'scid_hdw':ucid,
           'channels_ucid':channels_ucid,
           'Q0_bcs_ucid':Q0_bcs_ucid}
    if ucid in cids_core:
     db_ex_local[cid][ucid] = copy.deepcopy(db_ex)
    else:
     tag = int('%s%s' % (str(cid).ljust(4,'0'),str(ucid).ljust(4,'0')))
     self.comm.send(db_ex,dest=dest,tag=tag)

  #Wait until all are done
  #self.comm.Barrier()
  
  #Receive headwater data
  recv = {}
  for cid in cids_core:
   self = hbdb[cid].routing
   crm = hbdb[cid].cid_rank_mapping #Where each cid resides
   recv[cid] = {}
   for ucid in self.rcids_hdw:
    if ucid not in recv[cid]:recv[cid][ucid] = {}
    source = crm[ucid]
    if ucid not in cids_core:
     tag = int('%s%s' % (str(ucid).ljust(4,'0'),str(cid).ljust(4,'0')))
     db_ex = self.comm.recv(source=source,tag=tag)
    else:
     db_ex = db_ex_local[ucid][cid]
    for var in db_ex:
     recv[cid][ucid][var] = db_ex[var]

  #Update the boundary conditions
  for cid in cids_core:
   self = hbdb[cid].routing
   bcs = self.bcs
   self.bcs[:] = 0.0
   for ucid in self.rcids_hdw:
    channels_ucid = recv[cid][ucid]['channels_ucid']
    Q0_bcs_ucid = recv[cid][ucid]['Q0_bcs_ucid']
    for ic in range(channels_ucid.size):
     self.bcs[channels_ucid[ic]] += Q0_bcs_ucid[ic]

  return

def update_macroscale_polygon_routing(cids,HBdb,flag_constant_Kvn):

  for cid in cids:
   self = HBdb[cid].routing
   A0_org = self.A0_org
   self.A0[:] = self.A1[:]
   #Update solution
   (self.Q0,self.u0,self.A1,self.Kvn0,self.Q0in) = update_solution_explicit(self.c_slope,self.c_n,self.u0,self.A0,self.topology,self.c_length,self.qss,self.bcs,self.dt_routing,self.fp_n,self.hdb['Ac'],self.hdb['Af'],self.hdb['Pc'],self.hdb['Pf'],self.hdb['W'],self.hdb['M'],A0_org,flag_constant_Kvn,self.Kvn0)
   #self.A0[:] = self.A1[:]
   #self.Q1[:] = self.Q0[:]

  return

def calculate_routing_inundation(cids,HBdb):

  for cid in cids:
   self = HBdb[cid].routing
   #Zero out qss
   #self.qss[:] = 0.0
   #Calculate area-weighted average inundation height per hband
   A = self.hdb['Af'] + self.hdb['Ac']
   #Add channel and flood cross sectional area
   A1 = self.A0[:]
   W = self.hdb['W']
   M = self.hdb['M']
   hand = self.hdb['hand']
   hband = self.hdb['hband']
   (self.hband_inundation[:],self.reach2hband_inundation[:]) = calculate_inundation_height_per_hband(A,A1,W,M,hand,hband,self.reach2hband_inundation,self.reach2hband)
   tmp = np.sum(self.reach2hband*self.reach2hband_inundation,axis=1)/self.c_length
   tmp1 = np.max(np.abs(A1-tmp))
   if tmp1 > 10**-3:
    argmax = np.argmax(np.abs(A1-tmp))
    print(self.itime,np.abs(A1-tmp)[argmax],A1[argmax],tmp[argmax],flush=True)

   #Calculate Qc and Qf
   self.Ac[:] = self.reach2hband[range(self.nchannel),self.hband_channel]*self.reach2hband_inundation[range(self.nchannel),self.hband_channel]/self.c_length
   self.Af[:] = self.A0 - self.Ac
   self.Qc[:] = self.u0*self.Ac
   self.Qf[:] = self.u0*self.Af
   self.Qc[self.Qc < 0.0] = 0.0
   self.Qf[self.Qf < 0.0] = 0.0

  return

def exchange_velocity_fields(cids,HBdb):

 #Compute velocities on all reaches
 for cid in cids:
  self = HBdb[cid].routing
  self.u0[:] = 2.0 #m/s (Set to 2.0 m/s on all reaches for now)
  #Send velocities
  for ucid in self.particle_tracker_db_send:
   if ucid in cids:continue
   dest = HBdb[cid].cid_rank_mapping[ucid]
   tag = int('%s%s' % (str(cid).ljust(4,'0'),str(ucid).ljust(4,'0')))
   #print('send',cid,ucid,dest,tag,flush=True)
   tmp = self.u0[self.particle_tracker_db_send[ucid]].astype(np.float16)
   self.comm.send(tmp,dest=dest,tag=tag)

 #Receive velocities for the downstream reaches
 for cid in cids:
  self = HBdb[cid].routing
  for ucid in self.particle_tracker_db_receive:
   if ucid in cids:
    self.downstream_u0[self.cid_mapping[ucid-1],self.particle_tracker_db_receive[ucid]] = HBdb[ucid].routing.u0[HBdb[ucid].routing.particle_tracker_db_send[cid]]
   else:
    source = HBdb[cid].cid_rank_mapping[ucid]
    tag = int('%s%s' % (str(ucid).ljust(4,'0'),str(cid).ljust(4,'0')))
    #print('receive',cid,ucid,source,tag,flush=True)
    tmp = self.comm.recv(source=source,tag=tag)
    self.downstream_u0[self.cid_mapping[ucid-1],self.particle_tracker_db_receive[ucid]] = tmp
  
 #Wait
 #self.comm.Barrier()

 return

def update_particle_tracker_macroscale_polygon(cids,HBdb):

 #Push water downstream using u0 fields
 for cid in cids:
  HBr = HBdb[cid].routing
  dt = HBr.dt_routing
  HBr.Vin[:] = 0.0
  #Determine how much water in the channel leaves a given channel
  cl = dt*HBr.u0# HBr.c_length/HBr.u0
  m = cl < HBr.c_length
  f = np.ones(HBr.c_length.size)
  f[m] = cl[m]/HBr.c_length[m]
  Vc = f*(HBr.A0*HBr.c_length)#+dt*HBr.qss*HBr.c_length)
  snake_length = f*HBr.c_length
  HBr.Vout[:] = Vc[:]#0.0
  HBr.downstream_Vin[:] = 0.0
  HBr.downstream_Vout[:] = 0.0
  (HBr.Vin,HBr.Vout,HBr.downstream_Vin,HBr.downstream_Vout) = push_water_downstream(HBr.downstream_u0,HBr.downstream_channels,HBr.dt_routing,HBr.cid_mapping,HBr.u0,cid,HBr.c_length,HBr.downstream_c_length,HBr.Vin,HBr.Vout,HBr.downstream_Vin,HBr.downstream_Vout,Vc,snake_length)

 return

@numba.jit(nopython=True,cache=True,nogil=True,fastmath=True)
def push_water_downstream(downstream_u0,downstream_channels,dt,cid_mapping,u0,cid,c_length,downstream_c_length,Vin,Vout,downstream_Vin,downstream_Vout,Vc,snake_length):

 for ic in range(downstream_channels.shape[0]):
  t = 0.0
  for jc in range(downstream_channels.shape[2]):
   cid1 = downstream_channels[ic,1,jc]
   if ((cid1 == -9999) | (cid1 == -1)):break
   channel = downstream_channels[ic,0,jc]
   if cid1 == cid:
    u = u0[channel]
    cl = c_length[channel]
   else:
    u = downstream_u0[cid_mapping[cid1-1],channel]
    cl = downstream_c_length[cid_mapping[cid1-1],channel]
   t += cl/u
   if t >= dt:
    #Determine fraction that ends up in the last channel
    cc = cl-(t - dt)*u
    f = cc/snake_length[ic]
    if (f > 1) | (jc == 0):
     if cid1 == cid: 
      Vin[channel] += Vc[ic]
     else:
      downstream_Vin[cid_mapping[cid1-1],channel] += Vc[ic]
    else:
     #Split between channels (probably necessary to make possible to split into multiple channels instead; should be ok for now)
     #bottom channel
     if cid1 == cid:
       Vin[channel] += f*Vc[ic]
     else:
      downstream_Vin[cid_mapping[cid1-1],channel] += f*Vc[ic]
     #top channel
     cid2 = downstream_channels[ic,1,jc-1]
     channel2 = downstream_channels[ic,0,jc-1]
     if cid2 == cid:
      Vin[channel2] += (1-f)*Vc[ic]
     else:
      downstream_Vin[cid_mapping[cid2-1],channel2] += (1-f)*Vc[ic]
    break
     
   else:
    if cid1 == cid:
     Vin[channel] += Vc[ic]
     Vout[channel] += Vc[ic]
    else:
     downstream_Vin[cid_mapping[cid1-1],channel] += Vc[ic]
     downstream_Vout[cid_mapping[cid1-1],channel] += Vc[ic]

 return (Vin,Vout,downstream_Vin,downstream_Vout)

def exchange_water_volumes(cids,HBdb):

 #Send computed Vin/Vout to the corresponding cid
 for cid in cids:
  self = HBdb[cid].routing
  #Send Vin/Vout
  for ucid in self.particle_tracker_db_receive:
   if ucid in cids:continue
   dest = HBdb[cid].cid_rank_mapping[ucid]
   tag = int('%s%s' % (str(cid).ljust(4,'0'),str(ucid).ljust(4,'0')))
   ics = self.particle_tracker_db_receive[ucid]
   iucid = self.cid_mapping[ucid-1]
   db = {'Vin':self.downstream_Vin[iucid,ics].astype(np.float16),
          'Vout':self.downstream_Vout[iucid,ics].astype(np.float16)}
   self.comm.send(db,dest=dest,tag=tag)

 #Receive computed Vin/Vout
 for cid in cids:
  self = HBdb[cid].routing
  for ucid in self.particle_tracker_db_send:
   ics = self.particle_tracker_db_send[ucid]
   if ucid in cids:
    iucid2 = HBdb[ucid].routing.cid_mapping[cid-1]
    ics2 = HBdb[ucid].routing.particle_tracker_db_receive[cid]
    HBdb[cid].routing.Vin[ics] += HBdb[ucid].routing.downstream_Vin[iucid2,ics2]
    HBdb[cid].routing.Vout[ics] += HBdb[ucid].routing.downstream_Vout[iucid2,ics2]
   else:
    source = HBdb[cid].cid_rank_mapping[ucid]
    tag = int('%s%s' % (str(ucid).ljust(4,'0'),str(cid).ljust(4,'0')))
    db = self.comm.recv(source=source,tag=tag)
    HBdb[cid].routing.Vin[ics] += db['Vin'][:]
    HBdb[cid].routing.Vout[ics] += db['Vout'][:]

 #Calculate Qin,Qout,A
 for cid in cids:
  HBr = HBdb[cid].routing
  dt = HBr.dt_routing
  HBr.Qin[:] = HBr.Vin/dt
  HBr.Qout[:] = HBr.Vout/dt
  #qss should probable be before pushing water
  HBr.A1[:] = HBr.A0[:] - dt*HBr.Qout/HBr.c_length + dt*HBr.Qin/HBr.c_length + dt*HBr.qss
  #HBr.A1[:] = -dt*HBr.Qout/HBr.c_length + dt*HBr.Qin/HBr.c_length 
  #A1 = A0 + dt*qss + dt*bcs/c_length - dt*(u0*A0)/c_length + dt*Q0in/c_length

 return
