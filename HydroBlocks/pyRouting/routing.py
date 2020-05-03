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

 def __init__(self,MPI,cid,cid_rank_mapping,dt):

  print("we are running!")
  self.comm = MPI.COMM_WORLD
  self.size = self.comm.Get_size()
  self.rank = self.comm.Get_rank()
  self.name = MPI.Get_processor_name()
  self.cid_rank_mapping = cid_rank_mapping

  #Read in all the relevant data (this should be moved to the preprocessor)
  self.db = pickle.load(open('/stor/soteria/hydro/private/nc153/projects/Octopy/parallel/tmp/%d.pck' % cid,'rb'))

  #Define the variables
  self.dt = dt
  self.cid = cid
  self.c_length = self.db['c_length'][:]
  self.c_n = 1.0#self.db['c_n'][:]
  self.c_slope = self.db['c_slope'][:]
  self.nchannel = self.c_length.size
  self.A0 = 10**-5*np.ones(self.c_length.size)
  self.u0 = 10**-5*np.ones(self.c_length.size)
  self.A1 = np.zeros(self.c_length.size)
  self.bcs = np.zeros(self.c_length.size)
  self.qin = np.zeros(self.c_length.size)
  self.qout = np.zeros(self.c_length.size)
  self.dA = np.zeros(self.c_length.size)
  self.Q0 = np.zeros(self.c_length.size)
  self.Q1 = np.zeros(self.c_length.size)
  self.reach2hru = self.db['reach2hru']
  self.scids_hdw = self.db['scids_hdw']
  self.rcids_hdw = self.db['rcids_hdw']
  self.hdw = self.db['hdw']
  self.hdb = self.db['hdb']
  self.LHS = self.db['LHS']
  self.A0_org = self.A0[:]
  self.u0_org = self.u0[:]
  #dict_keys(['hdb','LHS','tsolve', 'tcount','cdst', 'scids_hdw', 'rcids_hdw', 'hdw', 'runoff', 'reach2hru'])

  #Pause until all cores have all their data
  self.comm.Barrier()
	
  return

 def update(self,dt):

  max_niter = 2
  #Update solution
  for itr in range(max_niter):
   #Exchange boundary conditions
   self.exchange_bcs()
   #Update solution
   self.update_solution(itr,max_niter)

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
   #if (cid == db_ex['scid_hdw']) & (db_ex['cid'] == ucid):
   for var in db_ex:
    recv[ucid][var] = db_ex[var]
   #else:
   # if db_ex['scid_hdw'] not in recv:recv[db_ex['scid_hdw']] = {}
   # if db_ex['cid'] not in recv[db_ex['scid_hdw']]:recv[db_ex['scid_hdw']][db_ex['cid']] = {}
   # for var in db_ex:
   #  recv[db_ex['scid_hdw']][db_ex['cid']][var] = db_ex[var]

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

 def update_solution(self,it,max_niter):

  #Extract info
  hdb = self.hdb
  A0 = self.A0
  c_slope = self.c_slope
  c_slope[c_slope < 5*10**-3] = 5*10**-3
  #c_slope[c_slope < 0.	0i
  c_n = self.c_n
  LHS = self.LHS
  c_length = self.c_length
  A0_org = self.A0_org
  u0_org = self.u0_org
  qin = self.qin
  qout = self.qout
  bcs = self.bcs
  Q0 = self.Q0
  u0 = self.u0
  dA = self.dA
  maxu = 10#0.1
  minu = 10**-5#0.1
  dt = self.dt

  #Determine hydraulic radius
  rh = calculate_hydraulic_radius(hdb['A'],hdb['P'],hdb['W'],A0)
  #Determine velocity
  u1 = rh**(2.0/3.0)*c_slope**0.5/c_n
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
  RHS1 = dt*qin*c_length - dt*qout*c_length
  RHS2 = np.zeros(A0.size)
  RHS2[:] = dt*bcs[:]
  RHS = (RHS0 + RHS1 + RHS2)#/(c_length + dt*u1)
  #Ax = b
  A1 = scipy.sparse.linalg.spsolve(LHS.tocsr(),RHS,use_umfpack=False)
  #QC
  A0[A0 < 0] = 0.0
  dif1 = np.mean(np.abs(A0 - A1))
  if (dif1 < 10**-10) | (it == max_niter-1):
   #Reset A0
   A0[:] = A1[:]
   #Determine hydraulic radius
   rh = calculate_hydraulic_radius(hdb['A'],hdb['P'],hdb['W'],A0)
   #Determine velocity
   u1 = rh**(2.0/3.0)*c_slope**0.5/c_n
   u1[u1 > maxu] = maxu
   u1[u1 < minu] = minu
   #Calculate Q1
   Q1 = A0*u1
   dif0 = -9999
   #Reset Q0
   Q0[:] = Q1[:]
   u0[:] = u1[:]
  else:
   #Reset A0
   A0[:] = A1[:]
   dif0 = dif1
   #Determine hydraulic radius
   rh = calculate_hydraulic_radius(hdb['A'],hdb['P'],hdb['W'],A0)
   #Determine velocity
   u1 = rh**(2.0/3.0)*c_slope**0.5/c_n
   u1[u1 > maxu] = maxu
   u1[u1 < minu] = minu
   #Calculate Q1
   Q1 = A0*u1
   #Reset Q0,u0
   Q0[:] = Q1[:]
   u0[:] = u1[:]
 
  #Update data in db
  self.Q0[:] = Q0[:]
  self.u0[:] = u0[:]
  self.A0[:] = A0[:]
  self.Q1[:] = Q1[:]
  self.A1[:] = A1[:]
  self.dA[:] = dA[:]
  self.bcs[:] = bcs[:]
  self.qin[:] = qin[:]
  self.qout[:] = qout[:]

  return
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
  h = (A1[i] - A0)/W1
  #Calculate the wetted perimeter
  P1 = P0 + 2*h + W1
  #Calculate hydraulic radius
  Rh[i] = A1[i]/P1

 return Rh
