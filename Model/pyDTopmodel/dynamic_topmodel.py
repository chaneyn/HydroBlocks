import numpy as np
import time
import dynamic_topmodel_tools as dtt
import scipy
import scipy.sparse
import scipy.sparse.linalg
import copy

class Dynamic_Topmodel:

 def __init__(self,ngroups,nhru_outlet):

  #Initialize the variables and parameters
  self.itime = 0 #timestep
  self.ngroups = ngroups
  self.m = np.zeros(ngroups,dtype=np.float32)
  self.si = np.zeros(ngroups,dtype=np.float32)
  self.si1 = np.zeros(ngroups,dtype=np.float32)
  self.sdmax = np.zeros(ngroups,dtype=np.float32)
  self.dt = 0.0 #seconds
  self.area = np.zeros(ngroups,dtype=np.float32)
  self.dem = np.zeros(ngroups,dtype=np.float32)
  self.carea = np.zeros(ngroups,dtype=np.float32)
  self.channel = np.zeros(ngroups,dtype=np.float32)
  self.T0 = np.zeros(ngroups,dtype=np.float32)
  self.beta = np.zeros(ngroups,dtype=np.float32)
  self.dx = np.zeros(ngroups,dtype=np.float32)
  self.sideep = np.zeros(ngroups,dtype=np.float32)
  self.pct = np.zeros(ngroups,dtype=np.float32)
  self.dsi = np.zeros(ngroups,dtype=np.float32)
  self.sti = np.zeros(ngroups,dtype=np.float32)

  #Define surface variables
  self.q_surface = np.zeros(ngroups,dtype=np.float32)
  self.storage_surface = np.zeros(ngroups,dtype=np.float32)
  self.qin_surface = np.zeros(ngroups,dtype=np.float32)
  self.qout_surface = np.zeros(ngroups,dtype=np.float32)
  self.storage1_surface = np.zeros(ngroups,dtype=np.float32)
  self.qout1_surface = np.zeros(ngroups,dtype=np.float32)
  self.qout1_surface[:] = -9999.0
  self.qin1_surface = np.zeros(ngroups,dtype=np.float32)
  self.qsurf = np.zeros(ngroups,dtype=np.float32)
  self.recharge_surface = np.zeros(ngroups,dtype=np.float32)
  self.celerity_surface = np.zeros(ngroups,dtype=np.float32)
  self.recharge1_surface = np.zeros(ngroups,dtype=np.float32)
  self.celerity1_surface = np.zeros(ngroups,dtype=np.float32)
  self.surface_velocity = np.zeros(ngroups,dtype=np.float32)
 
  #Current time step
  self.q_subsurface = np.zeros(ngroups,dtype=np.float32)
  self.r = np.zeros(ngroups,dtype=np.float32)
  self.qout = np.zeros(ngroups,dtype=np.float32)
  self.qin = np.zeros(ngroups,dtype=np.float32)
  self.ex = np.zeros(ngroups,dtype=np.float32)
  self.c = np.zeros(ngroups,dtype=np.float32)

  #Previous time step
  self.c1 = np.zeros(ngroups,dtype=np.float32)
  self.r1 = np.zeros(ngroups,dtype=np.float32)
  self.qout1 = np.zeros(ngroups,dtype=np.float32)
  self.qout1[:] = -9999.0
  self.qin1 = np.zeros(ngroups,dtype=np.float32)

  #Outlet information
  self.qin_outlet = np.zeros(nhru_outlet,dtype=np.float32)
  self.qin_outlet_surface = np.zeros(nhru_outlet,dtype=np.float32)
  self.area_outlet = np.zeros(nhru_outlet,dtype=np.float32)
  self.nhru_outlet = nhru_outlet

  #weights
  self.flow_matrix = []

 def update(self,ncores):
  
  #Update the subsurface runoff
  #self.update_subsurface_fortran(ncores)

  #Update the surface runoff
  self.update_surface_fortran(ncores)

  if self.itime % 10 == 0: print self.itime,self.qout_surface

  return

 def update_surface_fortran(self,ncores):

  #Set the recharge to be the sum of surface and excess runoff
  self.recharge1_surface[:] = self.recharge_surface
  self.recharge_surface[:] = self.qsurf + self.ex

  #Initialize q,qout1,qin1,c,and c1
  if self.qout1_surface[0] == -9999.0:
   self.qout1_surface[:] = 0.0
   self.qin1_surface[:] = 0.0
   self.a = 1.67
   self.n = 0.030
   self.b = np.tan(self.beta)**0.5/self.n
   self.celerity_surface[:] = Calculate_Celerity_Surface(self.a,self.b,self.storage_surface)
  
  #Update the celerity
  self.celerity1_surface[:] = self.celerity_surface[:]
  self.celerity_surface[:] = Calculate_Celerity_Surface(self.a,self.b,self.storage_surface)

  #Solve for the given time step
  '''(self.storage_surface,self.storage_surface1,self.qout_surface,self.qout1_surface,
             self.qin_surface,self.qin1_surface,self.celerity_surface,
             self.celerity1_surface) = Update(self.recharge_surface,self.storage_surface,
             self.qout_surface,self.qin_surface,
             self.recharge1_surface,self.storage1_surface,
             self.qout1_surface,self.qin1_surface,
             self.area,self.dx,self.dt,self.celerity_surface,self.celerity1_surface,
             self.flow_matrix,
             self.qin_outlet_surface,self.area_outlet,ncores)'''
  #F = scipy.sparse.csr_matrix(self.flow_matrix.todense().T)
  dtt.update(self.recharge_surface,self.storage_surface,self.qout_surface,self.qin_surface,
             self.recharge1_surface,self.storage1_surface,self.qout1_surface,self.qin1_surface,
            self.area,self.dx,self.dt,self.celerity_surface,self.celerity1_surface,
             #self.flow_matrix.data,self.flow_matrix.indices,self.flow_matrix.indptr,
             F.data,F.indices,F.indptr,
             ncores)

  #Correct the surface storage
  #Determine the amount of water that is "missing"
  #missing = np.sum(self.area[self.storage_surface < 0.0]*self.storage_surface[self.storage_surface < 0.0])
  #self.storage_surface[-1] = self.storage_surface[-1] + missing/self.area[-1]
  self.storage_surface[self.storage_surface < 0] = 0.0

  return

 def update_subsurface_fortran(self,ncores):

  #Initialize q,qout1,qin1,c,and c1
  if self.qout1[0] == -9999.0:
   self.q_subsurface[:] = Calculate_Flux_Subsurface(self.si,self.T0,self.beta,self.m,self.sdmax)
   self.q_subsurface[self.q_subsurface < 0.0] = 0.0
   self.qout1[:] = self.q_subsurface[:]
   self.qin1[:] = 0.0
   self.c[:] = Calculate_Celerity_Subsurface(self.m,self.q_subsurface)
   self.c1[:] = self.c[:]

  #Update the celerity
  self.q_subsurface[:] = Calculate_Flux_Subsurface(self.si,self.T0,self.beta,self.m,self.sdmax)
  self.c1[:] = self.c[:]
  self.c[:] = Calculate_Celerity_Subsurface(self.m,self.q_subsurface)

  #Set deficit in the form that the solver wants
  si = np.copy(-self.si)
  si1 = np.copy(-self.si1)
 
  #Solve for the given time step
  '''dtt.update(self.r,si,self.qout,self.qin,self.q_subsurface,
             self.r1,si1,self.qout1,self.qin1,
             self.area,self.dx,self.dt,self.c,self.c1,
             self.w.data,self.w.indices,self.flow_matrix.indptr,
             self.qin_outlet,self.area_outlet,ncores,0)'''
  (si,si1,self.qout,self.qout1,self.qin,self.qin1,self.c,
             self.c1) = Update(self.r,si,self.qout,self.qin,
             self.r1,si1,self.qout1,self.qin1,
             self.area,self.dx,self.dt,self.c,self.c1,
             self.flow_matrix,
             self.qin_outlet,self.area_outlet,ncores)

  #Revert the deficits to their original form
  self.si[:] = -si
  self.si1[:] = -si1

  #Set the excess runoff
  self.ex[:] = 0
  #print self.si,self.dt
  self.ex[self.si < 0] = -self.si[self.si < 0]/self.dt
  self.si[self.si < 0] = 0.0

  #Memorize the recharge for the time step
  self.r1[:] = self.r

  return

def Update(recharge,storage,qout,qin,recharge1,storage1,qout1,qin1,
                  area,dx,dt,celerity,celerity1,flow_matrix,
                  qin_outlet,area_outlet,nthreads):

 #Determine the appropriate time step
 dt_minimum = np.min(np.abs(dx/celerity))
 #ntt = 1*(int(np.ceil(dt/dt_minimum)) + 1)
 ntt = int(np.ceil(dt/dt_minimum))
 if ntt == 0:ntt = 1
 dtt = dt/ntt

 #Initialize variables to average the sub time steps
 qout_ = np.zeros(storage.size,dtype=np.float32)

 #Define some constatns
 w = 0.5
 I = scipy.sparse.identity(storage.size)
 F = flow_matrix
 scarea = area/dx
 dummy1 = scipy.sparse.dia_matrix(F.shape)
 dummy2 = scipy.sparse.dia_matrix(F.shape)
 dummy2.setdiag(scarea)

 for itime in xrange(ntt):

  #Solve the kinematic wave for this time step
  #Define the constants
  denominator = (1 + dtt*w*celerity/dx)
  numerator1 = qout1 + dtt*w*celerity*recharge + dtt*(1.0 - w)*celerity1*((qin1 - qout1)/dx + recharge1)
  numerator2 = dtt*w*celerity/dx
  p1 = numerator1/denominator
  p2 = numerator2/denominator
  dummy1.setdiag(p2/scarea)
  b = p1
  A = ((F*dummy1).T*dummy2).T

  #Solve for this time step
  tmp = scipy.sparse.linalg.spsolve((I-A).T,b)
  #print tmp
  qout[:] = tmp[:]

  #Set all negative fluxes to 0 
  qout[qout < 0.0] = 0.0

  #Calculate qin
  qin = (scarea*qout*F)/scarea

  #Adjust the storages
  storage = storage + dtt*((qin - qout)/dx + recharge)

  #Set the next time step's info
  qout1[:] = qout[:]
  qin1[:] = qin[:]

 return (storage,storage1,qout,qout1,qin,qin1,celerity,celerity1)

def Calculate_Flux_Subsurface(si,T0,beta,m,sdmax):

 #return T0*np.sin(beta)*(np.exp(-si/m*np.cos(beta)) - np.exp(-sdmax/m*np.cos(beta)))
 return T0*np.tan(beta)*(np.exp(-si/m))

def Calculate_Surface_Velocity(a,b,h):

 return b*h**a

def Calculate_Celerity_Subsurface(m,q):

 return q/m

def Calculate_Flux_Surface(storage_surface,surface_velocity):

 return surface_velocity*storage_surface

def Calculate_Celerity_Surface(a,b,h):

 return a*b*h**(a-1)
