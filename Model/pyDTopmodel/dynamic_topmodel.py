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
  self.w = []#np.empty((ngroups,ngroups),dtype=np.float32,order='F')
  #self.wfull = np.empty((ngroups,ngroups),dtype=np.float32,order='F')

 def update(self,ncores):
  
  #Update the subsurface runoff
  self.update_subsurface_fortran(ncores)

  #Update the surface runoff
  #tic = time.time()
  self.update_surface_fortran(ncores)
  #print time.time() - tic

  return

 def update_surface_fortran(self,ncores):

  #Set the recharge to be the sum of surface and excess runoff
  self.recharge1_surface[:] = self.recharge_surface
  self.recharge_surface[:] = self.qsurf + self.ex
  #self.recharge_surface[:] = 0#self.qsurf + self.ex

  #Initialize q,qout1,qin1,c,and c1
  if self.qout1_surface[0] == -9999.0:
   #self.storage_surface[:] = np.array([0.05,0.05,0.05])
   #self.q_surface[:] = np.array([0.01,0.01,0.01])#Calculate_Flux_Surface(self.storage_surface,self.surface_velocity)
   self.qout1_surface[:] = 0.0#self.q_surface[:]
   self.qin1_surface[:] = 0.0
   self.celerity_surface[:] = 0.0#Calculate_Celerity_Surface(self.surface_velocity)

  #Solve for the given time step
  #dtt.update(self.recharge_surface,self.storage_surface,self.qout_surface,self.qin_surface,self.q_surface,
  #           self.recharge1_surface,self.storage1_surface,self.qout1_surface,self.qin1_surface,
  #           self.area,self.dx,self.dt,self.celerity_surface,self.celerity1_surface,
  #           self.w.data,self.w.indices,self.w.indptr,
  #           self.qin_outlet_surface,self.area_outlet,ncores,1)
  (self.storage_surface,self.storage_surface1,self.qout_surface,self.qout1_surface,
             self.qin_surface,self.qin1_surface,self.celerity_surface,
             self.celerity1_surface) = Update(self.recharge_surface,self.storage_surface,
             self.qout_surface,self.qin_surface,self.q_surface,
             self.recharge1_surface,self.storage1_surface,
             self.qout1_surface,self.qin1_surface,
             self.area,self.dx,self.dt,self.celerity_surface,self.celerity1_surface,
             self.w,
             self.qin_outlet_surface,self.area_outlet,ncores,'surface',
             self.T0,self.beta,self.m,self.sdmax,self.surface_velocity,self.itime)

  #Correct the surface storage
  #print self.storage_surface
  #self.storage_surface[self.storage_surface < 0] = 0.0

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

  #Set deficit in the form that the solver wants
  si = np.copy(-self.si)
  si1 = np.copy(-self.si1)
 
  #Solve for the given time step
  #dx = np.copy(self.dx)
  #dx[:] = 1.0
  '''dtt.update(self.r,si,self.qout,self.qin,self.q_subsurface,
             self.r1,si1,self.qout1,self.qin1,
             self.area,self.dx,self.dt,self.c,self.c1,
             self.w.data,self.w.indices,self.w.indptr,
             self.qin_outlet,self.area_outlet,ncores,0)'''
  #t0 = time.time()
  (si,si1,self.qout,self.qout1,self.qin,self.qin1,self.c,
             self.c1) = Update(self.r,si,self.qout,self.qin,self.q_subsurface,
             self.r1,si1,self.qout1,self.qin1,
             self.area,self.dx,self.dt,self.c,self.c1,
             #self.w.data,self.w.indices,self.w.indptr,
             self.w,
             self.qin_outlet,self.area_outlet,ncores,'subsurface',
             self.T0,self.beta,self.m,self.sdmax,self.surface_velocity,self.itime)
  #print time.time() - t0
  #print self.qout*self.area

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

def Update(recharge,storage,qout,qin,q,recharge1,storage1,qout1,qin1,
                  area,dx,dt,celerity,celerity1,flow_matrix,
                  qin_outlet,area_outlet,nthreads,model_type,
                  T0,beta,m,sdmax,surface_velocity,timestamp):

 #Determine the appropriate time step
 dt_minimum = np.min(np.abs(dx/celerity))
 ntt = 4*(int(np.ceil(dt/dt_minimum)) + 1)
 #ntt = 1*(int(np.ceil(dt/dt_minimum)) + 1)
 ntt = int(np.ceil(dt/dt_minimum))
 if ntt == 0:ntt = 1
 #if (ntt > 100):ntt = 100
 dtt = dt/ntt

 #Initialize variables to average the sub time steps
 #storage_ = np.zeros(storage.size,dtype=np.float32)
 qout_ = np.zeros(storage.size,dtype=np.float32)
 #storage1_ = np.copy(storage)
 #qout1_ = np.copy(qout)
 #qin1_ = np.copy(qin)
 #qin_outlet_ = np.zeros(storage.size,dtype=np.float32)
 #qin_ = np.zeros(storage.size,dtype=np.float32)

 #Define some constatns
 w = 1.0
 I = scipy.sparse.identity(storage.size)
 F = flow_matrix[0:storage.size,0:storage.size]
 scarea = area/dx
 dummy1 = scipy.sparse.dia_matrix(F.shape)
 dummy2 = scipy.sparse.dia_matrix(F.shape)
 dummy2.setdiag(scarea)
 if timestamp % 10 == 0: print timestamp,storage,qout,ntt

 for itime in xrange(ntt):

  #Update the flux and celerity
  if model_type == 'subsurface':
   q[:] = Calculate_Flux_Subsurface(-storage,T0,beta,m,sdmax)
   q[q < 0.0] = 0.0
   celerity[:] = Calculate_Celerity_Subsurface(m,q)
   #qout[:] = q[:]
  if model_type == 'surface':
   a = 1.67
   n = 0.030
   b = np.tan(beta)**0.5/n
   storage[storage < 0] = 0.0
   #surface_velocity = Calculate_Surface_Velocity(a,b,storage)
   #qout[:] = Calculate_Flux_Surface(storage,surface_velocity)
   celerity[:] = Calculate_Celerity_Surface(a,b,storage)
   recharge[:] = -recharge[:]
  #Solve the kinematic wave for this time step
  #Define the constants
  denominator = (1 + dtt*w*celerity/dx)
  numerator1 = qout1 + dtt*w*celerity*recharge + dtt*(1.0 - w)*celerity1*((qin1 - qout1)/dx + recharge1)
  numerator2 = dtt*w*celerity/dx
  p1 = numerator1/denominator
  p2 = numerator2/denominator
  qout[:] = qout1[:]
  '''eps = 10.0
  iter = 0
  while ((eps > 10**-20) & (iter < 1000)):
   qin[:] = (scarea*qout*F)/scarea
   qout[:] = p1 + p2*qin
   eps = np.max(np.abs(qout - qout_))
   qout_[:] = qout[:]
   iter = iter + 1'''
  #print model_type,qout,recharge
  #if model_type == 'surface': exit()
  dummy1.setdiag(p2/scarea)
  b = p1
  #A = F.multiply(p2.T).multiply(area_sp)
  #t0 = time.time()
  #A = ((F*dummy2).T*dummy1).T
  A = ((F*dummy2).T*dummy1).T

  #Solve for this time step
  qout[:] = scipy.sparse.linalg.spsolve((I-A).T,b)
  #qout[:] = scipy.sparse.linalg.bicg(I-A,b,tol=10**-15)[0]
  #if model_type == 'surface':exit()
  #qout = p1 + p2*qin1
  #Set all negative fluxes to 0 
  qout[qout < 0.0] = 0.0
  #qout[qout > (recharge+storage/dtt)/storage] = 0.0

  #Calculate qin
  qin = (scarea*qout*F)/scarea #Changed F

  #Adjust the storages
  storage = storage - dtt*((qout - qin)/dx + recharge)
  #storage[storage < 0] = 0.0

  #Set the next time step's info
  qout1[:] = qout[:]
  qin1[:] = qin[:]
  celerity1[:] = celerity[:]

  #Update celerity?
  #print 'here',area*qout
 #print area*qout

 #Save the variables
 #qout1[:] = qout1_[:]
 #storage1[:] = storage1_[:]
 #qin1[:] = qin1_[:]

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

 maxc = 100
 tmp = a*b*h**(a-1)
 tmp[tmp > maxc/3600.0] = maxc/3600.0
 return tmp
