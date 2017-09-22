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
  #self.hand = np.zeros(nhru)
  self.area = np.zeros(nhru)
  self.dz = np.zeros((nhru,nsoil))
  self.hdiv = np.zeros((nhru,nsoil))

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
  psi = satpsi*((theta-thetar)/(thetas-thetar))**-b

  return psi

 def calculate_hydraulic_conductivity(self,psi):
  
  Ksat_x = 10.0*self.ksat[:] #lateral saturated hydraulic conductivity (multiply times anisotropy factor) [m/s]
  K_x = Ksat_x*(psi/self.satpsi)**(-2-3/self.b)

  return K_x

 def calculate_hydraulic_head(self,psi):
  
  h = self.dem - psi

  return h

 def calculate_divergence_dense(self,h,K_x,dz):
 
  dh = h[:,np.newaxis] - h[np.newaxis,:]
  dx = self.dx #meters (distance between two adjacent grid cells)
  w = np.array(self.width.todense())
  area = self.area 
  #Khat = (K_x[:,np.newaxis]*K_x[np.newaxis,:]*(w+w.T))/(K_x[:,np.newaxis]*w.T + K_x[np.newaxis,:]*w)
  Khat = (2*K_x[:,np.newaxis]*K_x[np.newaxis,:])/(K_x[:,np.newaxis] + K_x[np.newaxis,:])
  #[mm/s] = [mm/m]*[m/s]*[m]/[m]*[m]*[m]/[m2]
  return -1000.0*Khat*dh/dx*w*dz/area #mm/s

 def calculate_divergence_sparse(self,h,K_x,dz):

  #Define the boolean matrix (connections or not?)
  I = self.I
  #Calculate dh
  h1 = (I != 0).multiply(sparse.csr_matrix(h))
  dh = h1.T - h1
  #Calculate the effective hydraulic conductivity 
  k1 = (I != 0).multiply(sparse.csr_matrix(K_x))
  n = 2*k1.T.multiply(k1)
  d = k1.T+k1
  Khat = n.multiply(d.power(-1))
  #Calculate the flux
  #[m/s] = [m/s]*[m]/[m]*[m]*[m]/[m2]
  return -Khat.multiply(dh).multiply(self.width).multiply(dz/self.dx/self.area).multiply(1000) #mm/s

 def update(self,type='sparse'):

  #Determine if sparse or not
  if self.nhru <= 100: type = 'dense'

  #Iterate per layer
  for il in xrange(self.theta.shape[1]):
   #Calculate soil moisture potential
   psi = self.calculate_soil_moisture_potential(il)
   #Calculate hydraulic conductivity
   K_x = self.calculate_hydraulic_conductivity(psi)
   #Calculate hydraulic head
   h = self.calculate_hydraulic_head(psi)
   #Calculate the divergence
   if type == 'dense':
    q = self.calculate_divergence_dense(h,K_x,self.dz[:,il])
    self.hdiv[:,il] = np.sum(q,axis=0) #mm/s
    #print self.hdiv[:,il]
   elif type == 'sparse':
    q = self.calculate_divergence_sparse(h,K_x,self.dz[:,il])
    self.hdiv[:,il] = q.sum(axis=0) #mm/s

  return

