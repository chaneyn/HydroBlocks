import numpy as np
import time
import dynamic_topmodel_tools as dtt
import scipy
import scipy.sparse
import scipy.sparse.linalg

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
  #self.update_surface_fortran(ncores)
  #print time.time() - tic

  return

 def update_surface_fortran(self,ncores):

  #Set the recharge to be the sum of surface and excess runoff
  self.recharge1_surface[:] = self.recharge_surface
  self.recharge_surface[:] = self.qsurf + self.ex
  #self.recharge_surface[:] = 0#self.qsurf + self.ex

  #Estimate the flux
  self.q_surface = self.Calculate_Flux_Surface(self.storage_surface)
  self.q_surface[self.q_surface < 0] = 0.0

  #Set the celerities
  self.celerity1_surface[:] = self.celerity_surface
  self.celerity_surface[:] = self.Calculate_Celerity_Surface()

  #Solve for the given time step
  dtt.update(self.recharge_surface,self.storage_surface,self.qout_surface,self.qin_surface,self.q_surface,
             self.recharge1_surface,self.storage1_surface,self.qout1_surface,self.qin1_surface,
             self.area,self.dx,self.dt,self.celerity_surface,self.celerity1_surface,
             self.w.data,self.w.indices,self.w.indptr,
             self.qin_outlet_surface,self.area_outlet,ncores,1)
  #dtt.update(self.recharge_surface,self.storage_surface,self.qout_surface,self.qin_surface,
  #	     self.recharge1_surface,self.storage1_surface,self.qout1_surface,self.qin1_surface,
  #	     self.area,self.dx,self.dt,self.celerity_surface,self.celerity1_surface,
  #           self.w.data,self.w.indices,self.w.indptr,ncores)

  #Correct the surface storage
  #print self.storage_surface
  #self.storage_surface[self.storage_surface < 0] = 0.0

  return

 def update_subsurface_fortran(self,ncores):

  #Compute the flux
  #if self.itime == 0:
  self.q_subsurface[:] = self.Calculate_Flux_Subsurface(self.si)
  self.q_subsurface[self.q_subsurface < 0.0] = 0.0
  #print self.q_subsurface[:]

  #Calculate the celerities
  self.c1[:] = self.c
  self.c[:] = self.Calculate_Celerity_Subsurface(self.m,self.q_subsurface)

  #Set deficit in the form that the solver wants
  si = np.copy(-self.si)
  si1 = np.copy(-self.si1)
  #print 'before:',self.si
  #self.r[np.isnan(self.r) == 1] = 0
  #print self.r
 
  #Solve for the given time step
  #dx = np.copy(self.dx)
  #dx[:] = 1.0
  '''dtt.update(self.r,si,self.qout,self.qin,self.q_subsurface,
             self.r1,si1,self.qout1,self.qin1,
             self.area,self.dx,self.dt,self.c,self.c1,
             self.w.data,self.w.indices,self.w.indptr,
             self.qin_outlet,self.area_outlet,ncores,0)'''
  Update(self.r,si,self.qout,self.qin,self.q_subsurface,
             self.r1,si1,self.qout1,self.qin1,
             self.area,self.dx,self.dt,self.c,self.c1,
             #self.w.data,self.w.indices,self.w.indptr,
             self.w,
             self.qin_outlet,self.area_outlet,ncores,0)

  #print 'after:',-si

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

 def Calculate_Flux_Subsurface(self,si):

  return self.T0*np.sin(self.beta)*(np.exp(-si/self.m*np.cos(self.beta)) - np.exp(-self.sdmax/self.m*np.cos(self.beta)))
  #return 10*self.T0#np.sin(self.beta)*(np.exp(-si/self.m*np.cos(self.beta)) - np.exp(-self.sdmax/self.m*np.cos(self.beta)))
  #return self.T0*np.tan(self.beta)*(np.exp(-si/self.m))

 def Calculate_Celerity_Subsurface(self,m,q):
 
  return q/m

 def Calculate_Flux_Surface(self,storage_surface):

  return self.surface_velocity*storage_surface

 def Calculate_Celerity_Surface(self,):

  return self.surface_velocity

def Update(recharge,storage,qout,qin,q,recharge1,storage1,qout1,qin1,
                  area,dx,dt,celerity,celerity1,flow_matrix,
                  qin_outlet,area_outlet,nthreads,model_type):

 #Determine the appropriate time step
 dt_minimum = np.min(dx/celerity)
 ntt = 2*(int(np.ceil(dt/dt_minimum)) + 1)
 dtt = dt/ntt

 #Initialize variables to average the sub time steps
 storage_ = np.zeros(storage.size,dtype=np.float32)
 qout_ = np.zeros(storage.size,dtype=np.float32)
 storage1_ = np.copy(storage)
 qout1_ = np.copy(qout)
 qin1_ = np.copy(qin)
 qin_outlet_ = np.zeros(storage.size,dtype=np.float32)
 qin_ = np.zeros(storage.size,dtype=np.float32)

 #Define some constatns
 w = 1.0#0.5
 I = scipy.sparse.identity(storage.size)
 F = flow_matrix[0:storage.size,0:storage.size]

 #Process the smaller time steps
 for itime in xrange(ntt):

  #Solve the kinematic wave for this time step
  #Define the constants
  denominator = (1 + dtt*w*celerity/dx)
  numerator1 = qout1 + dtt*w*celerity*recharge + dtt*(1.0 - w)*celerity1*((qin1 - qout1)/dx + recharge1)
  numerator2 = dtt*w*celerity/dx
  #F = F.T
  p1 = numerator1/denominator
  p2 = scipy.sparse.csr_matrix(numerator2/denominator)
  #F = scipy.sparse.csr_matrix(np.ones(F.shape))
  b = p1.T
  A = F.multiply(scipy.sparse.csr_matrix(p2/area).T).multiply(scipy.sparse.csr_matrix(area))
  #A = F.multiply(scipy.sparse.csr_matrix(p2))
  #A = (part2.multiply(part3)).multiply(F.T)
  #b = (scipy.sparse.csr_matrix(part1).multiply(part3).T).todense()
  Aprime = I - A

  #Solve for this time step
  #print Aprime.shape
  #print p1
  #print b.shape,p1.size
  qout = scipy.sparse.linalg.bicgstab(Aprime,b,tol=10**-20)[0]#,x0 = qout1)
  #print qout[1]
  #qout = qout[0]
  #qout = scipy.sparse.linalg.spsolve(Aprime,b)
  #print 'in',area*q
  #print 'out',area*qout

  #Set all negative fluxes to 0 
  qout[qout < 0.0] = 0.0

  #Calculate qin
  qin = (area*qout*F)/area

  #Adjust the storages
  storage = storage + dtt*((qin - qout)/dx + recharge)

  #Set the next time step's info
  qout1[:] = qout[:]
  qin1[:] = qin[:]

  #Update celerity?

 #Save the variables
 qout1[:] = qout1_[:]
 storage1[:] = storage1_[:]
 qin1[:] = qin1_[:]

 return
