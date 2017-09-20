import numpy as np

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

  return

 def calculate_soil_moisture_potential(self,il):
  
  theta = self.theta[:,il]
  thetar = self.thetar
  thetas = self.thetas
  b = self.b
  satpsi = self.satpsi #meters
  psi = 1000.0*np.ones(theta.shape) #limit
  m = theta > thetar
  psi[m] = satpsi[m]*((theta[m] - thetar[m])/(thetas[m] - thetar[m]))**(-b[m])

  return psi

 def calculate_hydraulic_conductivity(self,psi):
  
  Ksat_x = 10.0*self.ksat[:] #lateral saturated hydraulic conductivity (multiply times anisotropy factor) [m/s]
  K_x = Ksat_x*(psi/self.satpsi)**(-2-3/self.b)

  return K_x

 def calculate_hydraulic_head(self,psi):
  
  h = self.dem - psi

  return h

 def update(self,):

  #Note:This needs to be updated for a sparse matrix configuration
  #Note:The current setup will be very slow with dense matrices above nhru ~ 100

  #Iterate per layer
  for il in xrange(self.theta.shape[1]):
   #Calculate soil moisture potential
   psi = self.calculate_soil_moisture_potential(il)
   #Calculate hydraulic conductivity
   K_x = self.calculate_hydraulic_conductivity(psi)
   #Calculate hydraulic head
   h = self.calculate_hydraulic_head(psi)
   #Calculate the divergence
   dh = h[:,np.newaxis] - h[np.newaxis,:]
   dx = self.dx #meters (distance between two adjacent grid cells)
   w = np.array(self.width.todense())
   area = self.area
   dz = self.dz[:,il] #meters
   Khat = (K_x[:,np.newaxis]*K_x[np.newaxis,:]*(w+w.T))/(K_x[:,np.newaxis]*w.T + K_x[np.newaxis,:]*w)
   Khat[np.isnan(Khat)] = 0.0
   #[m/s] = [m/s]*[m]/[m]*[m]*[m]/[m2]
   q = -Khat*dh/dx*w*dz/area #m/s
   self.hdiv[:,il] = 1000*np.sum(q,axis=0) #mm/s

  return

