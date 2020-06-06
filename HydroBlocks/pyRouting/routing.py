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

 def __init__(self,MPI,cid,cid_rank_mapping,dt,nhru,hru_area):

  self.itime = 0
  self.comm = MPI.COMM_WORLD
  self.size = self.comm.Get_size()
  self.rank = self.comm.Get_rank()
  self.name = MPI.Get_processor_name()
  self.cid_rank_mapping = cid_rank_mapping

  #Read in all the relevant data (this should be moved to the preprocessor)
  self.db = pickle.load(open('/stor/soteria/hydro/private/nc153/projects/Octopy/parallel/tmp/%d.pck' % cid,'rb'))

  #Define the variables
  self.hru_area = hru_area
  self.dt = dt
  self.cid = cid
  self.uhs = self.db['uhs'][:] #unit hydrographs
  self.uh_travel_distance = self.db['uh_travel_time'] #travel distance per hru (travel time is bug)
  self.c_length = self.db['c_length'][:]
  self.c_n = self.db['c_n'][:]
  self.fp_n = self.db['fp_n'][:]
  self.c_n[:] = 0.035#0.05#1#0.05#0.03
  self.fp_n[:] = 0.15#0.15#25#0.15#0.15
  self.c_slope = self.db['c_slope'][:]
  self.c_bankfull = self.db['c_bankfull'][:]
  self.nchannel = self.c_length.size
  self.A0 = 10**-5*np.ones(self.c_length.size)
  self.u0 = 10**-5*np.ones(self.c_length.size)
  self.A1 = np.zeros(self.c_length.size)
  self.bcs = np.zeros(self.c_length.size)
  #self.qin = np.zeros(self.c_length.size)
  #self.qout = np.zeros(self.c_length.size)
  self.qss = np.zeros(self.c_length.size)
  self.dA = np.zeros(self.c_length.size)
  self.Q0 = np.zeros(self.c_length.size)
  self.Q1 = np.zeros(self.c_length.size)
  self.reach2hru = np.asarray(self.db['reach2hru'].todense())
  self.reach2hru_inundation = np.copy(self.reach2hru)
  self.reach2hru_inundation[:] = 0.0
  self.scids_hdw = self.db['scids_hdw']
  self.rcids_hdw = self.db['rcids_hdw']
  self.hdw = self.db['hdw']
  self.hdb = self.db['hdb']
  self.hdb['M'][self.hdb['M'] == 0] = 10**-5
  self.c_width = self.hdb['W'][:,0]#self.db['c_width'][:]
  self.LHS = self.db['LHS']
  self.A0_org = self.A0[:]
  self.u0_org = self.u0[:]
  self.Ac0 = np.zeros(self.c_length.size)
  self.Ac1 = np.zeros(self.c_length.size)
  self.dAc0 = np.zeros(self.c_length.size)
  self.dAc1 = np.zeros(self.c_length.size)
  self.Af0 = np.zeros(self.c_length.size)
  self.Af1 = np.zeros(self.c_length.size)
  self.dAf0 = np.zeros(self.c_length.size)
  self.dAf1 = np.zeros(self.c_length.size)
  self.hru_inundation = np.zeros(nhru)
  self.hru_inundation1 = np.zeros(nhru)
  self.hru_inundation_available = np.zeros(nhru)

  #Define channel hrus per reach
  self.hru_channel = np.zeros(self.c_length.size).astype(np.int32)
  m = (self.hdb['hand'] == 0) & (self.hdb['W'] > 0)
  self.hru_channel[:] = self.hdb['hru'][m]

  #Pause until all cores have all their data
  self.comm.Barrier()
	
  return

 def update(self,dt):

  #5. Scale HRU changes from NoahMP to individual HRUs per reach
  #6. Balance HRU inundation height per reach
  #7. Extract Ac (and use to compute qin - qout)
  #Determine Ac and Af (compute from inundation per reach)

  #Update data
  self.A0_org[:] = self.A0[:]
  self.u0_org[:] = self.u0[:]

  max_niter = 1
  #Update solution
  for itr in range(max_niter):
   #Exchange boundary conditions
   self.exchange_bcs()
   #Update solution
   self.update_solution(itr,max_niter)

  #Zero out qss
  self.qss[:] = 0.0

  #Calculate area-weighted average inundation height per hru 
  A = self.hdb['Af'] + self.hdb['Ac']
  #Add channel and flood cross sectional area
  #A1 = self.A0
  A1 = self.A0[:]# + self.Af[:]
  W = self.hdb['W']
  M = self.hdb['M']
  hand = self.hdb['hand']
  hru = self.hdb['hru'].astype(np.int64)
  #self.hru_inundation[:] = calculate_inundation_height_per_hru(A,A1,W,hand,hru,self.reach2hru_inundation,self.reach2hru)
  (self.hru_inundation[:],self.reach2hru_inundation[:]) = calculate_inundation_height_per_hru(A,A1,W,M,hand,hru,self.reach2hru_inundation,self.reach2hru)
  tmp = np.sum(self.reach2hru*self.reach2hru_inundation,axis=1)/self.c_length
  tmp1 = np.max(np.abs(A1-tmp))
  if tmp1 > 10**-3:
   argmax = np.argmax(np.abs(A1-tmp))
   print(self.itime,np.abs(A1-tmp)[argmax],A1[argmax],tmp[argmax],flush=True)
  #if self.itime == 89:exit()

  #if (self.rank == 2) & (self.hru_inundation[2] > self.hru_inundation[1]):
  '''if (self.itime == 40) & (self.rank == 2):
   #print(self.c_bankfull)
   #print(self.c_bankfull-self.hru_inundation[self.hru_channel])
   tmp = []
   for ic in range(self.hru_channel.size):
    tmp = self.c_bankfull[ic]-self.reach2hru_inundation[ic,self.hru_channel[ic]]
    if tmp < 0:
     print(self.reach2hru_inundation[ic,self.hru_channel[ic]],self.c_bankfull[ic],self.c_width[ic],A1[ic])
     print(self.hru_channel[ic],self.hru_channel[ic]+2)
     print(self.hru_inundation[self.hru_channel[ic]:self.hru_channel[ic]+3])
     print(self.reach2hru_inundation[ic,self.hru_channel[ic]:self.hru_channel[ic]+3])
     print(hru[ic,:])
  if (self.itime == 40):exit()'''

  #Define change in Ac after inundation
  #Ac0 = []
  #for ic in range(self.hru_channel.size):
  # tmp = self.reach2hru[ic,self.hru_channel[ic]]*self.reach2hru_inundation[ic,self.hru_channel[ic]]/self.c_length[ic]
  # Ac0.append(tmp)
  #Ac0 = np.array(Ac0)
  #Af = A1 - Ac0
  #self.Ac0[:] = Ac0[:]
  #self.dAc0[:] = Ac0 - self.A0
  #self.

  #Add change to source/sink term for next time step
  #self.qss[:] += self.dAc0[:]/self.dt

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

 def update_solution(self,it,max_niter):

  #Extract info
  hdb = self.hdb
  A0 = self.A0
  c_slope = self.c_slope
  c_slope[c_slope < 10**-3] = 10**-3
  c_n = self.c_n
  LHS = self.LHS
  c_length = self.c_length
  A0_org = self.A0_org
  u0_org = self.u0_org
  #qin = self.qin
  #qout = self.qout
  qss = self.qss
  bcs = self.bcs
  Q0 = self.Q0
  u0 = self.u0
  dA = self.dA
  maxu = 2.0#10.0#0.1#0.1#1.0#0.1#1.0#2.0
  minu = 0.1#10.0#1.0#2.0#1.0#0.1#10**-6#0.1
  dt = self.dt

  #Determine hydraulic radius
  #rh = calculate_hydraulic_radius_rect(self.c_bankfull,self.c_width,A0)
  #rh = calculate_hydraulic_radius(hdb['Ac'],hdb['Af'],hdb['Pc'],hdb['Pf'],hdb['W'],A0)
  Kvn = calculate_compound_convenyance(hdb['Ac'],hdb['Af'],hdb['Pc'],hdb['Pf'],hdb['W'],hdb['M'],A0,self.c_n,self.fp_n)
  #Determine velocity
  u1 = np.zeros(Kvn.size)
  u1[A0 > 0.0] = Kvn[A0 > 0.0]*c_slope[A0 > 0.0]**0.5/A0[A0 > 0.0]
  #u1 = Kvn*c_slope**0.5/A0
  #if np.sum(np.isnan(Kvn)) != 0:
  # print('problem!',flush=True)
  # #m = np.isnan(Kvn) == 1
  # #print(u1[m])
  # #print('Kvn',Kvn[m])
  # #print(A0[m])
  # exit()
  #u1 = rh**(2.0/3.0)*c_slope**0.5/c_n
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
  #RHS1 = dt*qin*c_length - dt*qout*c_length
  RHS1 = dt*qss*c_length
  RHS2 = np.zeros(A0.size)
  RHS2[:] = dt*bcs[:]
  RHS = (RHS0 + RHS1 + RHS2)#/(c_length + dt*u1)
  #Ax = b
  A1 = scipy.sparse.linalg.spsolve(LHS.tocsr(),RHS,use_umfpack=False)
  #QC
  #A0[A0 < 0] = 0.0
  A1[A1 < 0] = 0.0
  #A1[A1 <= 0] = 10**-5
  dif1 = np.mean(np.abs(A0 - A1))
  #if it == 0:print(dif1)
  #print(it,dif1)
  if (dif1 < 10**-10) | (it == max_niter-1):
   #Reset A0
   A0[:] = A1[:]
   #Determine hydraulic radius
   #rh = calculate_hydraulic_radius_rect(self.c_bankfull,self.c_width,A0)
   #rh = calculate_hydraulic_radius(hdb['Ac'],hdb['Af'],hdb['Pc'],hdb['Pf'],hdb['W'],A0)
   Kvn = calculate_compound_convenyance(hdb['Ac'],hdb['Af'],hdb['Pc'],hdb['Pf'],hdb['W'],hdb['M'],A0,self.c_n,self.fp_n)
   #Determine velocity
   u1 = np.zeros(Kvn.size)
   u1[A0 > 0.0] = Kvn[A0 > 0.0]*c_slope[A0 > 0.0]**0.5/A0[A0 > 0.0]
   #u1 = Kvn*c_slope**0.5/A0
   #u1 = rh**(2.0/3.0)*c_slope**0.5/c_n
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
   #rh = calculate_hydraulic_radius_rect(self.c_bankfull,self.c_width,A0)
   #rh = calculate_hydraulic_radius(hdb['A'],hdb['P'],hdb['W'],A0)
   #rh = calculate_hydraulic_radius(hdb['Ac'],hdb['Af'],hdb['P'],hdb['W'],A0)
   #rh = calculate_hydraulic_radius(hdb['Ac'],hdb['Af'],hdb['Pc'],hdb['Pf'],hdb['W'],A0)
   Kvn = calculate_compound_convenyance(hdb['Ac'],hdb['Af'],hdb['Pc'],hdb['Pf'],hdb['W'],hdb['M'],A0,self.c_n,self.fp_n)
   #Determine velocity
   u1 = np.zeros(Kvn.size)
   u1[A0 > 0.0] = Kvn[A0 > 0.0]*c_slope[A0 > 0.0]**0.5/A0[A0 > 0.0]
   #u1 = Kvn*c_slope**0.5/A0
   #u1 = rh**(2.0/3.0)*c_slope**0.5/c_n
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
  #self.dA[:] = dA[:]
  self.bcs[:] = bcs[:]
  #self.qin[:] = qin[:]
  #self.qout[:] = qout[:]
  self.qss[:] = qss[:]

  return

def calculate_hydraulic_radius_rect(c_bankfull,c_width,A):

 Arh = np.copy(A)
 #Calculate stage
 h = Arh/c_width
 #Calculate wetted perimeter
 P = 2*h + c_width
 #Restrict to bankfull under flooding conditions
 #m = h > c_bankfull
 #Arh[m] = c_bankfull[m]*c_width[m]
 #P[m] = 2*c_bankfull[m] + c_width[m]
 #Calculate hydraulic radius
 Rh = Arh/P

 return Rh

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

@numba.jit(nopython=True,cache=True)
def calculate_compound_convenyance(Ac,Af,Pc,Pf,W,M,A1,cn,fpn):
 
 Kvn = np.zeros(Ac.shape[0])
 #Determine the level for which we need to back out info
 for i in range(Ac.shape[0]):
  for j in range(Ac.shape[1]):
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
def calculate_inundation_height_per_hru(A,A1,W,M,hand,hru,reach2hru_inundation,reach2hru):

 #hru_sum = np.zeros(reach2hru.shape[1])
 #hru_sum[:] = 0.0
 #Zero out reach2hru_inundation
 reach2hru_inundation[:] = 0.0
 #Determine the inundation height per hand value per basin
 for i in range(A.shape[0]):
  A0 = 0.0
  for j in range(1,A.shape[1]):
   if A[i,j] > A1[i]:break
   if A[i,j+1] == 0.0:break
   A0 = A[i,j]
  #Calculate the height above the segment
  #htop = (A1[i] - A0)/np.sum(W[i,0:j]) (rectangular)
  #print(A1[i],A0,W[i,0])
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
  # print(j,A[i,:])
  # print(j,A0,A1[i],htop,h[0]*W[i,0],h[j-1]*W[i,j-1])
  if j > 2:
   #print(j)
   #print(hand[i,j-1],hand[i,1:j-1],(hand[i,2:j]+hand[i,1:j-1])/2,htop)
   h[1:j-1] = hand[i,j-1] - hand[i,1:j-1] - (hand[i,2:j]-hand[i,1:j-1])/2 + htop
   #h[1:j-1] = hand[i,j-1] - hand[i,1:j-1] + htop
   #h[1:j-1] = hand[i,j-1] - hand[i,1:j-1] - (hand[i,2:j]+hand[i,1:j-1])/2 + htop
   #print(h)
  #h0 = hand[i,j-1] - hand[i,:j]
  #h = htop + h0 #m
  #Add to hru sum (times area)
  idxs = hru[i,0:j]
  #Place the inundation level
  for k in range(idxs.size):
   idx = idxs[k]
   reach2hru_inundation[i,idx] = h[k]
   #a = hru_sum[idx]
   #b = reach2hru[i,idx]*h[k]
   #c = a + b
   #hru_sum[idx] = c
 areas = np.sum(reach2hru,axis=0)
 hru_inundation = np.zeros(areas.size)	
 #for k in range(hru_inundation.size):
 # if areas[k] > 0:hru_inundation[k] = hru_sum[k]/areas[k] #meters
 hru_inundation[:] = np.sum(reach2hru*reach2hru_inundation,axis=0)/np.sum(reach2hru,axis=0)

 #return hru_inundation
 return (hru_inundation,reach2hru_inundation)

