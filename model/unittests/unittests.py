import warnings
warnings.filterwarnings("ignore")
import unittest
import sys
sys.path.append('../')
from pyRichards import richards
#import pyDTopmodel.dynamic_topmodel_tools as dtt
import numpy as np
import scipy
import scipy.sparse
import time

class Richards(unittest.TestCase):

  def test_soil_moisture_potential(self):

   eps = 0.01
   nhru = 100
   nsoil = 5
   model = richards.richards(nhru,nsoil)
   #define parameters
   model.theta[:] = np.random.uniform(low=0.05,high=0.5,size=nhru)[:,np.newaxis]
   model.thetar[:] = 0.1
   model.thetas[:] = 0.5 
   model.b[:] = np.random.uniform(low=3,high=8,size=nhru)
   model.satpsi[:] = 0.1
   #choose layer
   psi = []
   bpsi = []
   for il in xrange(nsoil):
    psi.append(model.calculate_soil_moisture_potential(il))
    tmp = model.theta[:,il]
    m = (tmp <= (1+eps)*model.thetar)
    tmp[m] = (1+eps)*model.thetar[m]
    bpsi.append(model.satpsi*((model.theta[:,il]-model.thetar)/(model.thetas-model.thetar))**(-model.b))
   psi = np.array(psi)
   bpsi = np.array(bpsi)
   self.assertTrue(np.allclose(psi,bpsi))

  def test_hydraulic_conductivity(self):

   nhru = 100
   nsoil = 5
   model = richards.richards(nhru,1)
   model.ksat[:] = 10**-3*np.random.uniform(low=10**-6,high=1,size=nhru) #mm/s
   model.satpsi[:] = np.random.uniform(low=0.01,high=1.0,size=nhru)
   model.b[:] = np.random.uniform(low=3,high=8,size=nhru)
   psi = np.random.uniform(low=0.001,high=10**10,size=nhru)
   #model
   kx = model.calculate_hydraulic_conductivity(psi)
   #baseline
   bkx = 10*model.ksat*(psi/model.satpsi)**(-2-3/model.b)
   self.assertTrue(np.allclose(kx,bkx,rtol=1e-10, atol=1e-50))
 
  def test_dense_sparse_comparison(self):

   nhru = 5000
   nsoil = 1
   model = richards.richards(nhru,nsoil)
   #define parameters
   model.nhru = nhru
   model.theta[:] = np.random.uniform(low=0.1,high=0.5,size=(nhru,nsoil))
   model.thetar[:] = 0.1
   model.thetas[:] = 0.5
   model.b[:] = np.random.uniform(low=0.1,high=0.5,size=nhru)
   model.satpsi[:] = np.random.uniform(low=0.001,high=1.0,size=nhru)
   model.ksat[:] = 10**-3*np.random.uniform(low=10**-6,high=1,size=nhru) #mm/s
   model.dz[:] = np.random.uniform(low=0.01,high=100.0,size=nsoil)
   model.dx = 30.0
   model.area[:] = np.random.uniform(low=900,high=90000,size=nhru)
   model.dem[:] = np.random.uniform(low=0.0,high=1000.0,size=nhru)
   tmp = np.random.randint(low=-1000,high=2,size=(nhru,nhru))
   tmp[tmp <= 0] = 0
   tmp[range(nhru),range(nhru)] = 1
   tsymm = 30*(tmp + tmp.T)/2
   model.width = scipy.sparse.csr_matrix(tsymm)
   model.I = model.width.copy()
   model.I[model.I != 0] = 1
   print '%d connections out of %d possible connections' % (np.sum(model.width != 0),nhru*nhru)
   #update (dense)
   tic = time.time()
   model.update(type='dense')
   print 'dense',time.time()-tic
   hdiv = np.copy(model.hdiv)
   #update (sparse)
   tic = time.time()
   model.update(type='sparse')
   print 'sparse',time.time()-tic
   bhdiv = np.copy(model.hdiv)
   #compare
   self.assertTrue(np.allclose(hdiv,bhdiv,rtol=1e-10, atol=1e-50))

suite = unittest.TestLoader().loadTestsFromTestCase(Richards)
unittest.TextTestRunner(verbosity=2).run(suite)
#suite = unittest.TestLoader().loadTestsFromTestCase(DynamicTopmodel)
#unittest.TextTestRunner(verbosity=2).run(suite)
