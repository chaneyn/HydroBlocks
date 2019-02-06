import numpy as np
import scipy.sparse as sparse
import time

class richards:

 def __init__(self,nhru,nsoil):

  self.theta = np.zeros((nhru,nsoil))
  self.thetar = np.zeros(nhru)
  self.thetas = np.zeros(nhru)
  self.b = np.zeros(nhru)
  self.satpsi = np.zeros(nhru)
  self.ksat = np.zeros(nhru)
  self.dem = np.zeros(nhru)
  self.slope = np.zeros(nhru)
  #self.hand = np.zeros(nhru)
  self.area = np.zeros(nhru)
  self.dz = np.zeros((nhru,nsoil))
  self.hdiv = np.zeros((nhru,nsoil))
  self.m = np.zeros(nhru)

  #Initialize the width array
  self.width = []
  self.I = []

  return

 def calculate_soil_moisture_potential(self,il):
  
  eps = 0.01
  theta = self.theta[:,il]
  thetar = self.thetar
  thetas = self.thetas
  b = self.b
  satpsi = self.satpsi #meters
  m = (theta <= (1+eps)*thetar)
  theta[m] = (1+eps)*thetar[m]
  #psi = 1000.0*np.ones(theta.shape) #limit
  #m = theta > thetar
  #psi[m] = satpsi[m]*((theta[m] - thetar[m])/(thetas[m] - thetar[m]))**(-b[m])
  #psi[m] = satpsi[m]*((theta[m] - thetar[m])/(thetas[m] - thetar[m]))**(-b[m])

  with np.errstate(invalid='ignore'):
    psi = satpsi*((theta-thetar)/(thetas-thetar))**-b

  return psi

 def calculate_hydraulic_conductivity(self,psi,il):
 
  af = 1.0 #safe
  #sdz = np.cumsum(self.dz,axis=1)-self.dz/2.
  #df = np.exp(-self.m[:,np.newaxis]/sdz)[:,il]
  #Ksat_x = af*df*self.ksat[:] #lateral saturated hydraulic conductivity (multiply times anisotropy factor) [m/s]
  Ksat_x = af*self.ksat[:] #lateral saturated hydraulic conductivity (multiply times anisotropy factor) [m/s]
  with np.errstate(invalid='ignore'):
   K_x = Ksat_x*(psi/self.satpsi)**(-2-3/self.b)

  return K_x

 def calculate_transmissivity(self,psi,ztop,zbot):
  
  af = 10.0  #safe
  #af = 2.0
  m = np.copy(self.m)
  #m[:] = 1000.0
  Ksat_x = af*self.ksat[:] #lateral saturated hydraulic conductivity (multiply times anisotropy factor) [m/s]
  with np.errstate(invalid='ignore', divide='ignore'):
   K_x = Ksat_x*np.true_divide(psi,self.satpsi)**(-2-np.true_divide(3.,self.b))
   K_x[~np.isfinite(K_x)] = np.nan
   #Calculate transmissivity at top layer (exponential decay)
   Ttop = m*K_x*np.exp(-ztop/m)
   #Calculate transmissivity at bottom of layer (exponential decay)
   Tbot = m*K_x*np.exp(-zbot/m)
   T = Ttop - Tbot
    
  return T

 def calculate_hydraulic_head(self,psi,depth):
  
  h = self.dem - depth - psi

  return h

 #def calculate_divergence_dense(self,h,K_x,dz):
 def calculate_divergence_dense(self,h,T):
 
  dh = h[:,np.newaxis] - h[np.newaxis,:]
  #dx = self.dx #meters (distance between two adjacent grid cells)
  dx = (np.abs(self.dem[:,np.newaxis] - self.dem[np.newaxis,:])**2 + self.dx**2)**0.5
  w = np.array(self.width.todense())
  #w[:] = 0.95*self.area[0]
  #dx[:] = 10.0
  #Calculate distance between centers (assuming hrus are rectangles... very large assumption)
  tmp = np.copy(w)
  tmp = self.area/w
  dx = (tmp + tmp.T)/2
  #dx[:] = 30.0#self.area[0]**0.5
  area = self.area 
  with np.errstate(invalid='ignore',divide='ignore'):
   #Khat = (K_x[:,np.newaxis]*K_x[np.newaxis,:]*(w+w.T))/(K_x[:,np.newaxis]*w.T + K_x[np.newaxis,:]*w)
   That = np.true_divide((2*T[:,np.newaxis]*T[np.newaxis,:]),(T[:,np.newaxis] + T[np.newaxis,:]))
   That[~np.isfinite(That)] = np.nan
   #[mm/s] = [mm/m]*[m/s]*[m]/[m]*[m]*[m]/[m2]
   calc_div = -1000.0*That*np.true_divide(dh,dx)*np.true_divide(w,area) # mm/s
   calc_div[~np.isfinite(calc_div)] = np.nan
  #print(calc_div)
  return calc_div

 #def calculate_divergence_sparse(self,h,K_x,dz):
 def calculate_divergence_sparse(self,h,T):

  #Define the boolean matrix (connections or not?)
  I = self.I
  #Calculate dh
  h1 = (I != 0).multiply(sparse.csr_matrix(h))
  dh = h1.T - h1
  #Calculate dx
  d1 = (I != 0).multiply(sparse.csr_matrix(self.dem))
  dx = d1.T - d1#**2 + self.dx**2)
  dx = dx.power(2)
  dx.data += self.dx**2
  dx = dx.power(0.5)
  #dx = (np.abs(d1.T - d1)**2 + self.dx**2)**0.5
  #dx = (np.abs(d1.T - d1)**2 + self.dx**2)**0.5
  #dx = (np.abs(self.dem[:,np.newaxis] - self.dem[np.newaxis,:])**2 + self.dx**2)**0.5
  #Calculate the effective hydraulic conductivity 
  #k1 = (I != 0).multiply(sparse.csr_matrix(K_x))
  #n = 2*k1.T.multiply(k1)
  #d = k1.T+k1
  #Khat = n.multiply(d.power(-1))
  t1 = (I != 0).multiply(sparse.csr_matrix(T))
  n = 2*t1.T.multiply(t1)
  d = t1.T+t1
  That = n.multiply(d.power(-1)).tocsr()
  print(That.count_nonzero,self.width.count_nonzero)
  #Calculate the flux
  #[m/s] = [m/s]*[m]/[m]*[m]*[m]/[m2]
  print(That.multiply(dh).shape)
  print(That.multiply(dh).multiply(self.width).count_nonzero)
  print(1.0/self.area)
  print(That.multiply(dh).multiply(self.width).multiply(1.0/self.area).count_nonzero)
  print(That.multiply(dh).multiply(self.width).multiply(1.0/self.area).multiply(dx.power(-1)).count_nonzero)
  return -That.multiply(dh).multiply(self.width).multiply(1.0/self.area).multiply(dx.power(-1)).multiply(1000) #mm/s
  #return -Khat.multiply(dh).multiply(self.width).multiply(dz/self.dx/self.area).multiply(1000) #mm/s

 def update(self,type='dense'):

  #Determine if sparse or not
  #if self.nhru <= 1000: type = 'dense'

  #Iterate per layer
  for il in range(self.theta.shape[1]):
   #Calculate soil moisture potential
   psi = self.calculate_soil_moisture_potential(il)
   #Calculate hydraulic conductivity
   #K_x = self.calculate_hydraulic_conductivity(psi)
   #ztop = np.cumsum(self.dz[0:il],axis=1)
   zbot = np.sum(self.dz[:,0:il+1],axis=1)
   ztop = zbot - self.dz[:,il]
   T = self.calculate_transmissivity(psi,ztop,zbot)
   #print('T:',T)
   #Calculate hydraulic head
   #h = self.calculate_hydraulic_head(psi,(ztop+zbot)/2,pw)
   h = self.calculate_hydraulic_head(psi,ztop)
   #Calculate the divergence
   #if type == 'dense':
   #q = self.calculate_divergence_dense(h,K_x)#,self.dz[:,il])
   q = self.calculate_divergence_dense(h,T)
   #print(q)
   #print(np.sum(q,axis=1))
   self.hdiv[:,il] = np.sum(q,axis=0) #mm/s
   #print(3600*self.hdiv)
   #print self.hdiv[:,il]
   #elif type == 'sparse': (Needs more work)
   # q = self.calculate_divergence_sparse(h,T)#,self.dz[:,il])
   # self.hdiv[:,il] = q.sum(axis=0) #mm/s

  return

