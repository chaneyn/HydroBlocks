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

class kinematic:

 def __init__(self,MPI,cid,cid_rank_mapping,dt,nhband,nhru,cdir,Qobs_file,dt_routing):

  self.itime = 0
  self.comm = MPI.COMM_WORLD
  self.size = self.comm.Get_size()
  self.rank = self.comm.Get_rank()
  self.name = MPI.Get_processor_name()
  self.cid_rank_mapping = cid_rank_mapping

  #Read in all the relevant data (this should be moved to the preprocessor)
  #self.db = pickle.load(open('/stor/soteria/hydro/private/nc153/projects/Octopy/parallel/tmp/%d.pck' ,'rb'))
  self.db = pickle.load(open('%s/octopy.pck' % cdir ,'rb'))

  #Read in the discharge input data
  #self.Qobs = pickle.load(open(Qobs_file,'rb'))
  #self.Qobs = pickle.load(open('/home/nc153/soteria/projects/hydroblocks_inter_catchment/regions/SGP_OK/obs/GSAL_2004-2019.pck','rb'))
  #self.Qobs2 = pickle.load(open('/home/nc153/soteria/projects/hydroblocks_inter_catchment/regions/SGP_OK/obs/chikaskia/chikaskia_2015-2017.pck','rb'))

  #Define the variables
  self.dt = dt
  self.dt_routing = dt_routing
  self.cid = cid
  self.uhs = self.db['uhs'][:] #unit hydrographs
  self.uh_travel_distance = self.db['uh_travel_time'] #travel distance per hband (travel time is bug)
  self.c_length = self.db['c_length'][:]
  self.c_n = self.db['c_n'][:]
  self.fp_n = self.db['fp_n'][:]
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
  self.reach2hband = np.asarray(self.db['reach2hband'].todense())
  self.reach2hband_inundation = np.copy(self.reach2hband)
  self.reach2hband_inundation[:] = 0.0
  self.scids_hdw = self.db['scids_hdw']
  self.rcids_hdw = self.db['rcids_hdw']
  self.hdw = self.db['hdw']
  self.hdb = self.db['hdb']
  self.hdb['hband'] = self.hdb['hband'].astype(np.int64)
  self.hdb['M'][self.hdb['M'] == 0] = 10**-5
  self.c_width = self.hdb['W'][:,0]#self.db['c_width'][:]
  self.LHS = self.db['LHS']
  tmp_topology = self.db['LHS'][:]
  tmp_topology.setdiag(0)
  tmp_topology = tmp_topology.todense()
  topology = np.zeros((tmp_topology.shape[0],3),dtype=np.int32)
  topology[:] = -9999
  for i in range(tmp_topology.shape[0]):
   ids = np.where(tmp_topology[i,:] != 0)[1]
   topology[i,0:ids.size] = ids[:]
  self.topology = topology
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
  self.comm.Barrier()
	
  return

 def update(self,dt):

  #Explicit solution
  #self.dt_routing = 100 #seconds
  dt_routing = self.dt_routing
  nt = int(dt/dt_routing)
  for it in range(nt):
   #Exchange boundary conditions
   self.exchange_bcs()
   #Update solution
   (self.Q0,self.u0,self.A1) = update_solution_explicit(self.c_slope,self.c_n,self.u0,self.A0,self.topology,self.c_length,self.qss,self.bcs,self.dt_routing,self.fp_n,self.hdb['Ac'],self.hdb['Af'],self.hdb['Pc'],self.hdb['Pf'],self.hdb['W'],self.hdb['M'])
   self.A0[:] = self.A1[:]
   self.Q1[:] = self.Q0[:]

  #Zero out qss
  self.qss[:] = 0.0
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


 def exchange_bcs(self,):

  #Send headwater data
  crm = self.cid_rank_mapping #Where each cid resides
  for ucid in self.scids_hdw:
   dest = crm[ucid]
   db_ex = {'cid':self.cid,'scid_hdw':ucid,
          'Q0':self.Q0[self.hdw[ucid][self.cid]['outlet']]}
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
   for ic in range(len(self.hdw[self.cid][ucid]['inlet'])):
    inlet = self.hdw[self.cid][ucid]['inlet'][ic]
    self.bcs[inlet] += recv[ucid]['Q0'][ic]

  #Wait until all are done
  self.comm.Barrier()

  return

 def update_solution_implicit(self,):

  #Extract info
  hdb = self.hdb
  A0 = self.A0[:]
  c_slope = self.c_slope
  c_n = self.c_n
  LHS = self.LHS
  c_length = self.c_length
  A0_org = self.A0_org[:]
  qss = self.qss
  bcs = self.bcs
  maxu = 2.0
  minu = 0.1
  dt = self.dt

  #Determine velocity
  #This effectively linearizes the sub-grid network solver by setting the velocity with the previous time step A (Another way to do this, is converge first locally every change in BCs)
  Kvn = calculate_compound_convenyance(hdb['Ac'],hdb['Af'],hdb['Pc'],hdb['Pf'],hdb['W'],hdb['M'],A0_org,self.c_n,self.fp_n)
  u1 = np.zeros(Kvn.size)
  u1[A0_org > 0.0] = Kvn[A0_org > 0.0]*c_slope[A0_org > 0.0]**0.5/A0_org[A0_org > 0.0]
  #Constrain velocity
  u1[u1 > maxu] = maxu
  u1[u1 < minu] = minu
  #Fill non-diagonals
  tmp = -dt*u1
  LHS = LHS.multiply(tmp)
  #Fill diagonal
  tmp = c_length + dt*u1
  LHS.setdiag(tmp)
  #Set right hand side
  RHS0 = c_length*A0_org
  RHS2 = np.zeros(A0.size)
  RHS2[:] = dt*bcs
  #Iterate qss values to ensure A1 is above 0
  RHS1 = dt*qss*c_length
  RHS = (RHS0 + RHS1 + RHS2)#/(c_length + dt*u1)
  A1 = scipy.sparse.linalg.spsolve(LHS.tocsr(),RHS,use_umfpack=False)
  #Curate A1 (issues with conservation of mass)
  A1[A1 < 0] = 0.0
  #Calculate difference with previous iteration
  dif1 = np.max(np.abs(A0 - A1))
  #Reset A0
  A0[:] = A1[:]
  #Calculate Q1
  Q1 = A0*u1
  #Reset Q0
  self.Q0[:] = Q1[:]
  self.u0[:] = u1[:]
 
  #Update data in db
  self.A0[:] = A0[:]
  self.Q1[:] = Q1[:]
  self.A1[:] = A1[:]
  self.dif0 = dif1

  return

 def update_solution_explicit_original(self,it):

  #Extract info
  hdb = self.hdb
  A0 = self.A0[:]
  c_slope = self.c_slope
  c_slope[c_slope < 10**-3] = 10**-3
  c_n = self.c_n
  LHS = self.LHS
  c_length = self.c_length
  A0_org = A0[:]#elf.A0_org[:]
  qss = self.qss
  bcs = self.bcs/c_length
  maxu = 10.0
  minu = 0.1
  dt = self.dt_routing

  #Determine velocity
  #This effectively linearizes the sub-grid network solver by setting the velocity with the previous time step A (Another way to do this, is converge first locally every change in BCs)
  Kvn = calculate_compound_convenyance(hdb['Ac'],hdb['Af'],hdb['Pc'],hdb['Pf'],hdb['W'],hdb['M'],A0_org,self.c_n,self.fp_n)
  u1 = np.zeros(Kvn.size)
  u1[A0_org > 0.0] = Kvn[A0_org > 0.0]*c_slope[A0_org > 0.0]**0.5/A0_org[A0_org > 0.0]
  #Constrain velocity
  u1[u1 > maxu] = maxu
  u1[u1 < minu] = minu
  #Ammending LHS from the implicit solver for this (HACK)
  LHS.setdiag(0)
  #LHS = np.array(LHS.todense())
  if (it == 0):
   Q0in = LHS*(u1*A0)
   if (self.rank == 0):print(Q0in[0:3])
  exit()
  #A1 = A0 + source/sink + boundary conditions - Qout + Qin
  A1 = A0 + dt*qss + dt*bcs - dt*(u1*A0)/c_length + dt*(LHS*(u1*A0))/c_length
  #Curate A1 (issues with conservation of mass)
  A1[A1 < 0] = 0.0
  #Reset A0
  A0[:] = A1[:]
  #Calculate Q1
  Q1 = A0*u1
  #Update data in db
  self.Q0[:] = Q1[:]
  self.u0[:] = u1[:]
  self.A0[:] = A0[:]
  self.Q1[:] = Q1[:]
  self.A1[:] = A1[:]

  return

@numba.jit(nopython=True,cache=True)
def update_solution_explicit(c_slope,c_n,u0,A0,topology,c_length,qss,bcs,dt_routing,
     fp_n,Ac,Af,Pc,Pf,W,M):

 #Extract info
 bcs_c = bcs/c_length
 maxu = 10.0
 minu = 0.1
 dt = dt_routing

 #Determine velocity
 #Kvn = calculate_compound_convenyance(hdb['Ac'],hdb['Af'],hdb['Pc'],hdb['Pf'],hdb['W'],hdb['M'],A0,c_n,fp_n)
 Kvn = calculate_compound_convenyance(Ac,Af,Pc,Pf,W,M,A0,c_n,fp_n)
 u0 = np.zeros(Kvn.size)
 u0[A0 > 0.0] = Kvn[A0 > 0.0]*c_slope[A0 > 0.0]**0.5/A0[A0 > 0.0]

 #Constrain velocity
 u0[u0 > maxu] = maxu
 u0[u0 < minu] = minu

 #Compute Q0in
 Q0in = Compute_Q0in(topology,u0,A0)

 #A1 = A0 + source/sink + boundary conditions - Qout + Qin
 A1 = A0 + dt*qss + dt*bcs_c - dt*(u0*A0)/c_length + dt*Q0in/c_length

 #Curate A1 (issues with conservation of mass)
 A1[A1 < 0] = 0.0

 #Calculate Q0
 Q0 = A0*u0

 return (Q0,u0,A1)


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
