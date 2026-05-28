import numpy as np
import scipy.sparse as sparse
import time
import numba

class richards:

 def __init__(self,nhru,nsoil,vsp_flag):

  # Initialize arrays for soil moisture and hydraulic properties
  self.theta = np.zeros((nhru,nsoil))
  if vsp_flag==True:
   self.thetar = np.zeros((nhru,nsoil)) #laura svp
   self.thetas = np.zeros((nhru,nsoil)) #laura svp
   self.b = np.zeros((nhru,nsoil)) #laura svp
   self.satpsi = np.zeros((nhru,nsoil)) #laura svp
   self.ksat = np.zeros((nhru,nsoil)) #laura svp
  else:
   self.thetar = np.zeros(nhru)
   self.thetas = np.zeros(nhru)
   self.b = np.zeros(nhru)
   self.satpsi = np.zeros(nhru)
   self.ksat = np.zeros(nhru)
  # Initialize arrays for hru properties
  self.dem = np.zeros(nhru)
  self.slope = np.zeros(nhru)
  #self.hand = np.zeros(nhru)
  self.area = np.zeros(nhru)
  self.dz = np.zeros((nhru,nsoil))
  self.m = np.zeros(nhru)
  # Initialize arrays for flow calculations
  self.hdiv = np.zeros((nhru,nsoil))
  self.hdiv_heat = np.zeros((nhru,nsoil))
  #Initialize the width array
  self.width = []
  self.I = []

  return

 def update_numba(self,vsp_flag):
  af = self.af
  theta = self.theta
  dz = self.dz
  hdiv = self.hdiv
  thetar = self.thetar
  thetas = self.thetas
  b = self.b
  satpsi = self.satpsi
  m = self.m
  ksat = self.ksat
  #hand = self.dem
  hand = self.dem1
  w = self.w
  dx = self.dx
  area = self.area
  if vsp_flag==True: #laura svp
   self.hdiv[:] = update_workhorse_vsp(theta,dz,hdiv,thetar,thetas,b,satpsi,m,ksat,hand,w,dx,area,
                                       af, self.flag_sat) #divergence computed with vertical variable soil properties
  else:
   self.hdiv[:] = update_workhorse(theta,dz,hdiv,thetar,thetas,b,satpsi,m,ksat,hand,w,dx,area)

  return

class richards_hbands:

 def __init__(self,nhru,nhband,nsoil,vsp_flag): #laura
  
  # Initialize arrays for soil moisture and hydraulic properties
  self.theta = np.zeros((nhru,nsoil))
  if vsp_flag==True:
   self.thetar = np.zeros((nhru,nsoil)) #laura svp
   self.thetas = np.zeros((nhru,nsoil)) #laura svp
   self.b = np.zeros((nhru,nsoil)) #laura svp
   self.satpsi = np.zeros((nhru,nsoil)) #laura svp
   self.ksat = np.zeros((nhru,nsoil)) #laura svp
  else:
   self.thetar = np.zeros(nhru)
   self.thetas = np.zeros(nhru)
   self.b = np.zeros(nhru)
   self.satpsi = np.zeros(nhru)
   self.ksat = np.zeros(nhru)
  # Initialize arrays for hru properties
  self.dem = np.zeros(nhru)
  self.demhband = np.zeros(nhband) #laura added
  self.dem1hband=np.zeros(nhband) #laura added
  self.slope = np.zeros(nhband)
  #self.hand = np.zeros(nhru)
  self.area = np.zeros(nhband)
  self.dz = np.zeros((nhband,nsoil))
  self.m = np.zeros(nhband)
  # Initialize arrays for flow calculations
  self.hdiv = np.zeros((nhband,nsoil))
  self.hdiv_heat = np.zeros((nhband,nsoil))
  #Initialize the width array
  self.width = {} #laura
  self.I = {}#laura

  return

 def update_numba(self,vsp_flag):
  af = self.af
  theta = self.theta
  dz = self.dz
  hdiv = self.hdiv
  thetar = self.thetar
  thetas = self.thetas
  b = self.b
  satpsi = self.satpsi
  m = self.m
  ksat = self.ksat
  #hand = self.dem
  hand = self.dem1hband
  w = self.w
  dx = self.dx
  area = self.area
  ncsbasins=self.ncsbasins #laura, number of characteristic subbasins
  #Compute divergence independently per characteristic subbasin, laura
  aux=0
  div=np.empty(self.hdiv.shape)
  for bas in range(1,ncsbasins+1):
   w_bas=w['Basin%s' %bas]
   dx_bas=dx['Basin%s' %bas]
   init=aux
   fin=aux+w_bas.shape[0]
   if vsp_flag==False:
    div[init:fin,:]=update_workhorse(theta[init:fin,:],dz[init:fin,:],hdiv[init:fin,:],
                                     thetar[init:fin],thetas[init:fin],b[init:fin],
                                     satpsi[init:fin],m[init:fin],ksat[init:fin],hand[init:fin],
                                     w_bas,dx_bas,area[init:fin],af,self.flag_sat) #divergence computed with uniform soil properties
    aux=fin
   else:
    div[init:fin,:]=update_workhorse_vsp(theta[init:fin,:],dz[init:fin,:],hdiv[init:fin,:],
                                         thetar[init:fin],thetas[init:fin],b[init:fin],
                                         satpsi[init:fin],m[init:fin],ksat[init:fin],hand[init:fin],
                                         w_bas,dx_bas,area[init:fin],af,self.flag_sat)
    aux=fin #laura, added to fix flerchinger
  self.hdiv=div

  return

@numba.jit(nopython=True,cache=True)
def update_workhorse_vsp(theta,dz,hdiv,thetar,thetas,b,satpsi,m,ksat,hand,w,dx,area,af,flag_sat):
 # flag_sat passed as parameter; Dupuit-Forchheimer approximation when True
 #Iterate per layer
 for il in range(theta.shape[1]):
  #Calculate soil moisture potential
  psi = calculate_soil_moisture_potential(il,theta,thetar[:,il],thetas[:,il],b[:,il],satpsi[:,il]) #laura svp
  zbot = np.sum(dz[:,0:il+1],axis=1)
  ztop = zbot - dz[:,il]
  T = calculate_transmissivity(psi,ztop,zbot,m,ksat[:,il],satpsi[:,il],b[:,il],af)#laura svp
  #Calculate hydraulic head
  h = calculate_hydraulic_head(hand,psi,ztop)
  #Calculate the divergence
  q = calculate_divergence(h,T,w,dx,area)
  #Apply sat/unsaturated separation
  if flag_sat:
   eps = 0.01
   for i in range(q.shape[0]):
    if theta[i, il] <= (1 - eps) * thetas[i, il]:
     q[i, :] = 0 #no horizontal flow in unsaturated layers
     q[:, i] = 0
  # q[i, j] stores the divergence contribution for source HRU i toward neighbor j.
  # Sum across each row so the integrated mass closes with area-normalized fluxes.
  hdiv[:,il] = np.sum(q,axis=0) #mm/s - sum over all connections to get divergence at each HRU
 return hdiv
 
@numba.jit(nopython=True,cache=True)
def update_workhorse(theta,dz,hdiv,thetar,thetas,b,satpsi,m,ksat,hand,w,dx,area,af):

 #Iterate per layer
 for il in range(theta.shape[1]):
  #Calculate soil moisture potential
  psi = calculate_soil_moisture_potential(il,theta,thetar,thetas,b,satpsi)
  zbot = np.sum(dz[:,0:il+1],axis=1)
  ztop = zbot - dz[:,il]
  T = calculate_transmissivity(psi,ztop,zbot,m,ksat,satpsi,b,af)
  #Calculate hydraulic head
  h = calculate_hydraulic_head(hand,psi,ztop)
  #Calculate the divergence
  q = calculate_divergence(h,T,w,dx,area)
  hdiv[:,il] = np.sum(q,axis=0) #mm/s

 return hdiv

@numba.jit(nopython=True,cache=True)
def calculate_soil_moisture_potential(il,theta,thetar,thetas,b,satpsi):
  
 eps = 0.01
 theta = theta[:,il]
 m = (theta <= (1+eps)*thetar)
 theta[m] = (1+eps)*thetar[m]
 psi = satpsi*((theta-thetar)/(thetas-thetar))**-b

 return psi

#@numba.jit(nopython=True,cache=True)
#def calculate_transmissivity(psi,ztop,zbot,m,ksat,satpsi,b,af):
  
 #Ksat_x = af*ksat #lateral saturated hydraulic conductivity (multiply times anisotropy factor) [m/s]
 #K_x = Ksat_x*np.true_divide(psi,satpsi)**(-2-np.true_divide(3.,b))
 #Calculate transmissivity at top layer (exponential decay)
 #Ttop = m*K_x*np.exp(-ztop/m)
 #Calculate transmissivity at bottom of layer (exponential decay)
 #Tbot = m*K_x*np.exp(-zbot/m)
 #T = Ttop - Tbot
  
 #return T

@numba.jit(nopython=True,cache=True)
def calculate_transmissivity(psi,ztop,zbot,m,ksat,satpsi,b,af):
  
 Ksat_x = af*ksat #lateral saturated hydraulic conductivity (multiply times anisotropy factor) [m/s]
 K_x = Ksat_x*np.true_divide(psi,satpsi)**(-2-np.true_divide(3.,b))
 #Correct hydraulic conductivity if layer is deeper than 2.0 meters
 depth_threshold = 2.0
 #Calculate transmissivity at top layer (exponential decay)
 Ttop = np.abs(np.where(zbot>depth_threshold, m*K_x*np.exp(-(ztop-depth_threshold)/m),K_x*ztop))
 #Calculate transmissivity at bottom of layer (exponential decay)
 Tbot = np.abs(np.where(zbot>depth_threshold, m*K_x*np.exp(-(zbot-depth_threshold)/m),K_x*zbot))
 # Compute transmissivity
 T = np.abs(Ttop - Tbot)
  
 return T

@numba.jit(nopython=True,cache=True)
def calculate_hydraulic_head(hand,psi,depth):
  
 h = hand - depth - psi

 return h

@numba.jit(nopython=True,cache=True)
def calculate_divergence(h,T,w,dx,area):
 
 #Calculate dh
 dh = calculate_dh(h)
 #Calculate That
 That = calculate_That(T)
 #[mm/s] = [mm/m]*[m/s]*[m]/[m]*[m]*[m]/[m2]
 calc_div = -1000.0*That*dh/dx*w/area # mm/s
 # sign convention: positive dh means flow from i to j, negative dh means flow from j to i. (Darcy's law: q = -K * A * (h[j] - h[i]) / dx[i,j])
 # negative sign because flow is from high to low head, but dh is calculated as h[i] - h[j]
 # The negative sign in the formula accounts for this convention, ensuring that positive divergence corresponds to net outflow from node i.

 return calc_div

@numba.jit(nopython=True,cache=True)
def calculate_dh(h):

 dh = np.zeros((h.size,h.size))
 for i in range(h.size):
  for j in range(h.size):
   dh[i,j] = h[i] - h[j]

 return dh

@numba.jit(nopython=True,cache=True)
def calculate_That(T):

 That = np.zeros((T.size,T.size))
 for i in range(T.size):
  for j in range(T.size):
   That[i,j] = (2*T[i]*T[j])/(T[i] + T[j])

 return That
