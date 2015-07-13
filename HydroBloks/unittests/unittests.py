import warnings
warnings.filterwarnings("ignore")
import unittest
import sys
sys.path.append('../')
from pyDTopmodel import dynamic_topmodel
import pyDTopmodel.dynamic_topmodel_tools as dtt
import numpy as np
import scipy
import scipy.sparse

class DynamicTopmodel(unittest.TestCase):

  def test_celerity_subsurface(self):
    m = 1.0
    q = 10.0
    output = dynamic_topmodel.Calculate_Celerity_Subsurface(m,q)
    self.assertEqual(output,10.0)

  def test_flux_subsurface(self):
    si = np.array([1.0,])
    T0 = np.array([0.1,])
    beta = np.array([0.01,])
    m = np.array([0.1,])
    sdmax = np.array([10.0,])
    output = dynamic_topmodel.Calculate_Flux_Subsurface(si,T0,beta,m,sdmax)
    self.assertTrue(np.allclose(output,np.array([4.54218782e-08,]),rtol=1e-05, atol=1e-08))

  def test_kinematic_wave_solution_python_explicit_onlyrecharge(self):

    nhru = 2
    qout = np.zeros(nhru)
    qin = np.zeros(nhru)
    recharge1 = np.zeros(nhru)
    storage1 = np.zeros(nhru)
    qout1 = np.zeros(nhru)
    qin1 = np.zeros(nhru)
    celerity1 = np.zeros(nhru)
    qin_outlet = np.zeros(nhru)
    area_outlet = np.zeros(nhru)
    storage = 0.0*np.ones(nhru)
    recharge = 1.0*np.ones(nhru)
    celerity = 1.0*np.ones(nhru)
    area = 1.0*np.ones(nhru)
    dx = 1.0*np.ones(nhru)
    dt = 1.0
    nthreads = 1
    maxntt = 1
    w = 1.0
    flow_matrix = scipy.sparse.csr_matrix(np.array([[0.,1.],[0.,0.]]))
    (storage,storage1,qout,qout1,qin,qin1,celerity,celerity1) = dynamic_topmodel.Update(recharge,
     storage,qout,qin,recharge1,storage1,qout1,qin1,area,dx,dt,celerity,celerity1,flow_matrix,
     qin_outlet,area_outlet,nthreads,maxntt,w)
    qout_true = np.array([0.5,0.75])
    self.assertTrue(np.allclose(qout,qout_true,rtol=1e-05, atol=1e-08))

  def test_kinematic_wave_solution_python_explicit_norecharge(self):

    nhru = 2
    qout = np.zeros(nhru)
    qin = np.zeros(nhru)
    recharge = np.zeros(nhru)
    recharge1 = np.zeros(nhru)
    storage1 = np.zeros(nhru)
    qin1 = np.zeros(nhru)
    celerity1 = np.zeros(nhru)
    qin_outlet = np.zeros(nhru)
    area_outlet = np.zeros(nhru)
    storage = 1.0*np.ones(nhru)
    celerity = 1.0*np.ones(nhru)
    area = 1.0*np.ones(nhru)
    qout1 = 1.0*np.ones(nhru)
    dx = 1.0*np.ones(nhru)
    dt = 1.0
    nthreads = 1
    maxntt = 1
    w = 1.0
    flow_matrix = scipy.sparse.csr_matrix(np.array([[0.,1.],[0.,0.]]))
    (storage,storage1,qout,qout1,qin,qin1,celerity,celerity1) = dynamic_topmodel.Update(recharge,
     storage,qout,qin,recharge1,storage1,qout1,qin1,area,dx,dt,celerity,celerity1,flow_matrix,
     qin_outlet,area_outlet,nthreads,maxntt,w)
    qout_true = np.array([0.5,0.75])
    self.assertTrue(np.allclose(qout,qout_true,rtol=1e-05, atol=1e-08))

  def test_kinematic_wave_solution_python_explicit(self):

    nhru = 2
    qout = np.zeros(nhru)
    qin = np.zeros(nhru)
    recharge1 = np.zeros(nhru)
    storage1 = np.zeros(nhru)
    qin1 = np.zeros(nhru)
    celerity1 = np.zeros(nhru)
    qin_outlet = np.zeros(nhru)
    area_outlet = np.zeros(nhru)
    recharge = 1.0*np.ones(nhru)
    storage = 1.0*np.ones(nhru)
    celerity = 1.0*np.ones(nhru)
    area = 1.0*np.ones(nhru)
    qout1 = 1.0*np.ones(nhru)
    dx = 1.0*np.ones(nhru)
    dt = 1.0
    nthreads = 1
    maxntt = 1
    w = 1.0
    flow_matrix = scipy.sparse.csr_matrix(np.array([[0.,1.],[0.,0.]]))
    (storage,storage1,qout,qout1,qin,qin1,celerity,celerity1) = dynamic_topmodel.Update(recharge,
     storage,qout,qin,recharge1,storage1,qout1,qin1,area,dx,dt,celerity,celerity1,flow_matrix,
     qin_outlet,area_outlet,nthreads,maxntt,w)
    qout_true = np.array([1.0,1.5])
    self.assertTrue(np.allclose(qout,qout_true,rtol=1e-05, atol=1e-08))

  def test_kinematic_wave_solution_python_implicit_onlyrecharge(self):

    nhru = 2
    qout = np.zeros(nhru)
    qin = np.zeros(nhru)
    qin1 = np.zeros(nhru)
    qout1 = np.zeros(nhru)
    qin_outlet = np.zeros(nhru)
    area_outlet = np.zeros(nhru)
    recharge = 1.0*np.ones(nhru)
    recharge1 = 1.0*np.ones(nhru)
    storage = 1.0*np.ones(nhru)
    storage1 = 1.0*np.ones(nhru)
    celerity = 1.0*np.ones(nhru)
    celerity1 = 1.0*np.ones(nhru)
    area = 1.0*np.ones(nhru)
    dx = 1.0*np.ones(nhru)
    dt = 1.0
    nthreads = 1
    maxntt = 1
    w = 0.5
    flow_matrix = scipy.sparse.csr_matrix(np.array([[0.,1.],[0.,0.]]))
    (storage,storage1,qout,qout1,qin,qin1,celerity,celerity1) = dynamic_topmodel.Update(recharge,
     storage,qout,qin,recharge1,storage1,qout1,qin1,area,dx,dt,celerity,celerity1,flow_matrix,
     qin_outlet,area_outlet,nthreads,maxntt,w)
    qout_true = np.array([0.666666,0.888888])
    self.assertTrue(np.allclose(qout,qout_true,rtol=1e-05, atol=1e-08))

  def test_kinematic_wave_solution_mkl_explicit(self):

    nhru = 2
    storage_mask_subsurface = np.ones(nhru)
    qout = np.zeros(nhru)
    qin = np.zeros(nhru)
    recharge1 = np.zeros(nhru)
    storage1 = np.zeros(nhru)
    qin1 = np.zeros(nhru)
    celerity1 = np.zeros(nhru)
    qin_outlet = np.zeros(nhru)
    area_outlet = np.zeros(nhru)
    storage = 0.0*np.ones(nhru)
    recharge = 1.0*np.ones(nhru)
    celerity = 1.0*np.ones(nhru)
    qout1 = 1.0*np.ones(nhru)
    area = 1.0*np.ones(nhru)
    dx = 1.0*np.ones(nhru)
    dt = 1.0
    ncores = 1
    maxntt = 1
    isw = 1.0
    flow_matrix = scipy.sparse.csr_matrix(np.array([[0.,1.],[0.,0.]]).T)
    flow_matrix.setdiag(flow_matrix.diagonal())
    #Initialize the solver
    dtt.initialize(flow_matrix.indices,flow_matrix.indptr)
    #Solve the system of equations
    dtt.update(recharge,storage,qout,qin,
             recharge1,storage1,qout1,qin1,
             area,dx,dt,celerity,celerity1,storage_mask_subsurface,
             flow_matrix.data,flow_matrix.indices,flow_matrix.indptr,
             ncores,maxntt,isw)
    #Finalize the solver
    dtt.finalize()
    #Compare
    qout_true = np.array([1.0,1.5])
    self.assertTrue(np.allclose(qout,qout_true,rtol=1e-05, atol=1e-08))

  def test_kinematic_wave_solution_mkl_explicit_onlyrecharge(self):

    nhru = 2
    storage_mask_subsurface = np.ones(nhru)
    qout = np.zeros(nhru)
    qin = np.zeros(nhru)
    recharge1 = np.zeros(nhru)
    storage1 = np.zeros(nhru)
    qout1 = np.zeros(nhru)
    qin1 = np.zeros(nhru)
    celerity1 = np.zeros(nhru)
    qin_outlet = np.zeros(nhru)
    area_outlet = np.zeros(nhru)
    storage = 0.0*np.ones(nhru)
    recharge = 1.0*np.ones(nhru)
    celerity = 1.0*np.ones(nhru)
    area = 1.0*np.ones(nhru)
    dx = 1.0*np.ones(nhru)
    dt = 1.0
    ncores = 1
    maxntt = 1
    isw = 1.0
    flow_matrix = scipy.sparse.csr_matrix(np.array([[0.,1.],[0.,0.]]).T)
    flow_matrix.setdiag(flow_matrix.diagonal())
    #Initialize the solver
    dtt.initialize(flow_matrix.indices,flow_matrix.indptr)
    #Solve the system of equations
    dtt.update(recharge,storage,qout,qin,
             recharge1,storage1,qout1,qin1,
             area,dx,dt,celerity,celerity1,storage_mask_subsurface,
             flow_matrix.data,flow_matrix.indices,flow_matrix.indptr,
             ncores,maxntt,isw)
    #Finalize the solver
    dtt.finalize()
    #Compare
    qout_true = np.array([0.5,0.75])
    self.assertTrue(np.allclose(qout,qout_true,rtol=1e-05, atol=1e-08))


  def test_kinematic_wave_solution_mkl_implicit_onlyrecharge(self):

    nhru = 2
    storage_mask_subsurface = np.ones(nhru)
    qout = np.zeros(nhru)
    qin = np.zeros(nhru)
    storage1 = np.zeros(nhru)
    qout1 = np.zeros(nhru)
    qin1 = np.zeros(nhru)
    qin_outlet = np.zeros(nhru)
    area_outlet = np.zeros(nhru)
    storage = 0.0*np.ones(nhru)
    recharge1 = 1.0*np.ones(nhru)
    recharge = 1.0*np.ones(nhru)
    celerity = 1.0*np.ones(nhru)
    celerity1 = np.ones(nhru)
    area = 1.0*np.ones(nhru)
    dx = 1.0*np.ones(nhru)
    dt = 1.0
    ncores = 1
    maxntt = 1
    isw = 0.5
    flow_matrix = scipy.sparse.csr_matrix(np.array([[0.,1.],[0.,0.]]).T)
    flow_matrix.setdiag(flow_matrix.diagonal())
    #Initialize the solver
    dtt.initialize(flow_matrix.indices,flow_matrix.indptr)
    #Solve the system of equations
    dtt.update(recharge,storage,qout,qin,
             recharge1,storage1,qout1,qin1,
             area,dx,dt,celerity,celerity1,storage_mask_subsurface,
             flow_matrix.data,flow_matrix.indices,flow_matrix.indptr,
             ncores,maxntt,isw)
    #Finalize the solver
    dtt.finalize()
    #Compare
    qout_true = np.array([0.666666,0.888888])
    self.assertTrue(np.allclose(qout,qout_true,rtol=1e-05, atol=1e-08))


  def test_kinematic_wave_solution_mkl_implicit(self):

    nhru = 2
    storage_mask_subsurface = np.ones(nhru)
    qout = np.zeros(nhru)
    qin = np.zeros(nhru)
    storage1 = np.zeros(nhru)
    qin1 = np.zeros(nhru)
    qin_outlet = np.zeros(nhru)
    area_outlet = np.zeros(nhru)
    storage = 0.0*np.ones(nhru)
    qout1 = np.ones(nhru)
    recharge1 = 1.0*np.ones(nhru)
    recharge = 1.0*np.ones(nhru)
    celerity = 1.0*np.ones(nhru)
    celerity1 = np.ones(nhru)
    area = 1.0*np.ones(nhru)
    dx = 1.0*np.ones(nhru)
    dt = 1.0
    ncores = 1
    maxntt = 1
    isw = 0.5
    flow_matrix = scipy.sparse.csr_matrix(np.array([[0.,1.],[0.,0.]]).T)
    flow_matrix.setdiag(flow_matrix.diagonal())
    #Initialize the solver
    dtt.initialize(flow_matrix.indices,flow_matrix.indptr)
    #Solve the system of equations
    dtt.update(recharge,storage,qout,qin,
             recharge1,storage1,qout1,qin1,
             area,dx,dt,celerity,celerity1,storage_mask_subsurface,
             flow_matrix.data,flow_matrix.indices,flow_matrix.indptr,
             ncores,maxntt,isw)
    #Finalize the solver
    dtt.finalize()
    #Compare
    qout_true = np.array([1.0,1.33333])
    self.assertTrue(np.allclose(qout,qout_true,rtol=1e-05, atol=1e-08))


  def test_kinematic_wave_solution_mkl_implicit_varydx(self):

    nhru = 2
    storage_mask_subsurface = np.ones(nhru)
    qout = np.zeros(nhru)
    qin = np.zeros(nhru)
    storage1 = np.zeros(nhru)
    qin1 = np.zeros(nhru)
    qin_outlet = np.zeros(nhru)
    area_outlet = np.zeros(nhru)
    storage = 0.0*np.ones(nhru)
    qout1 = np.ones(nhru)
    recharge1 = 1.0*np.ones(nhru)
    recharge = 1.0*np.ones(nhru)
    celerity = 1.0*np.ones(nhru)
    celerity1 = np.ones(nhru)
    area = 1.0*np.ones(nhru)
    dx = 2.0*np.ones(nhru)
    dt = 1.0
    ncores = 1
    maxntt = 1
    isw = 0.5
    flow_matrix = scipy.sparse.csr_matrix(np.array([[0.,1.],[0.,0.]]).T)
    flow_matrix.setdiag(flow_matrix.diagonal())
    #Initialize the solver
    dtt.initialize(flow_matrix.indices,flow_matrix.indptr)
    #Solve the system of equations
    dtt.update(recharge,storage,qout,qin,
             recharge1,storage1,qout1,qin1,
             area,dx,dt,celerity,celerity1,storage_mask_subsurface,
             flow_matrix.data,flow_matrix.indices,flow_matrix.indptr,
             ncores,maxntt,isw)
    #Finalize the solver
    dtt.finalize()
    #Compare
    qout_true = np.array([1.4,1.68])
    self.assertTrue(np.allclose(qout,qout_true,rtol=1e-05, atol=1e-08))


  def test_kinematic_wave_solution_mkl_implicit_varydt(self):

    nhru = 2
    storage_mask_subsurface = np.ones(nhru)
    qout = np.zeros(nhru)
    qin = np.zeros(nhru)
    storage1 = np.zeros(nhru)
    qin1 = np.zeros(nhru)
    qin_outlet = np.zeros(nhru)
    area_outlet = np.zeros(nhru)
    storage = 0.0*np.ones(nhru)
    qout1 = 1.0*np.ones(nhru)
    recharge1 = 1.0*np.ones(nhru)
    recharge = 1.0*np.ones(nhru)
    celerity = 1.0*np.ones(nhru)
    celerity1 = 1.0*np.ones(nhru)
    area = 1.0*np.ones(nhru)
    dx = 1.0*np.ones(nhru)
    dt = 2.0
    ncores = 1
    maxntt = 1
    isw = 0.5
    flow_matrix = scipy.sparse.csr_matrix(np.array([[0.,1.],[0.,0.]]).T)
    flow_matrix.setdiag(flow_matrix.diagonal())
    #Initialize the solver
    dtt.initialize(flow_matrix.indices,flow_matrix.indptr)
    #Solve the system of equations
    dtt.update(recharge,storage,qout,qin,
             recharge1,storage1,qout1,qin1,
             area,dx,dt,celerity,celerity1,storage_mask_subsurface,
             flow_matrix.data,flow_matrix.indices,flow_matrix.indptr,
             ncores,maxntt,isw)
    #Finalize the solver
    dtt.finalize()
    #Compare
    qout_true = np.array([1.0,1.5])
    self.assertTrue(np.allclose(qout,qout_true,rtol=1e-05, atol=1e-08))


  def test_kinematic_wave_solution_mkl_implicit_varyc(self):

    nhru = 2
    storage_mask_subsurface = np.ones(nhru)
    qout = np.zeros(nhru)
    qin = np.zeros(nhru)
    storage1 = np.zeros(nhru)
    qin1 = np.zeros(nhru)
    qin_outlet = np.zeros(nhru)
    area_outlet = np.zeros(nhru)
    storage = 0.0*np.ones(nhru)
    qout1 = np.ones(nhru)
    recharge1 = 1.0*np.ones(nhru)
    recharge = 1.0*np.ones(nhru)
    celerity = 2.0*np.ones(nhru)
    celerity1 = 2.0*np.ones(nhru)
    area = 1.0*np.ones(nhru)
    dx = 1.0*np.ones(nhru)
    dt = 1.0
    ncores = 1
    maxntt = 1
    isw = 0.5
    flow_matrix = scipy.sparse.csr_matrix(np.array([[0.,1.],[0.,0.]]).T)
    flow_matrix.setdiag(flow_matrix.diagonal())
    #Initialize the solver
    dtt.initialize(flow_matrix.indices,flow_matrix.indptr)
    #Solve the system of equations
    dtt.update(recharge,storage,qout,qin,
             recharge1,storage1,qout1,qin1,
             area,dx,dt,celerity,celerity1,storage_mask_subsurface,
             flow_matrix.data,flow_matrix.indices,flow_matrix.indptr,
             ncores,maxntt,isw)
    #Finalize the solver
    dtt.finalize()
    #Compare
    qout_true = np.array([1.0,1.5])
    self.assertTrue(np.allclose(qout,qout_true,rtol=1e-05, atol=1e-08))

  def test_kinematic_wave_solution_mkl_implicit_varyarea(self):

    nhru = 2
    storage_mask_subsurface = np.ones(nhru)
    qout = np.zeros(nhru)
    qin = np.zeros(nhru)
    storage1 = np.zeros(nhru)
    qin1 = np.zeros(nhru)
    qin_outlet = np.zeros(nhru)
    area_outlet = np.zeros(nhru)
    storage = 0.0*np.ones(nhru)
    qout1 = np.ones(nhru)
    recharge1 = 1.0*np.ones(nhru)
    recharge = 1.0*np.ones(nhru)
    celerity = 1.0*np.ones(nhru)
    celerity1 = np.ones(nhru)
    area = np.array([1.0,2.0])#[1.0*np.ones(nhru)
    dx = 1.0*np.ones(nhru)
    dt = 1.0
    ncores = 1
    maxntt = 1
    isw = 0.5
    flow_matrix = scipy.sparse.csr_matrix(np.array([[0.,1.],[0.,0.]]).T)
    flow_matrix.setdiag(flow_matrix.diagonal())
    #Initialize the solver
    dtt.initialize(flow_matrix.indices,flow_matrix.indptr)
    #Solve the system of equations
    dtt.update(recharge,storage,qout,qin,
             recharge1,storage1,qout1,qin1,
             area,dx,dt,celerity,celerity1,storage_mask_subsurface,
             flow_matrix.data,flow_matrix.indices,flow_matrix.indptr,
             ncores,maxntt,isw)
    #Finalize the solver
    dtt.finalize()
    #Compare
    qout_true = np.array([1.0,1.166666])
    self.assertTrue(np.allclose(qout,qout_true,rtol=1e-05, atol=1e-08))

  def test_kinematic_wave_solution_mkl_implicit_3hru(self):

    nhru = 3
    storage_mask_subsurface = np.ones(nhru)
    qout = np.zeros(nhru)
    qin = np.zeros(nhru)
    storage1 = np.zeros(nhru)
    qin1 = np.zeros(nhru)
    qin_outlet = np.zeros(nhru)
    area_outlet = np.zeros(nhru)
    storage = 0.0*np.ones(nhru)
    qout1 = np.ones(nhru)
    recharge1 = 1.0*np.ones(nhru)
    recharge = 1.0*np.ones(nhru)
    celerity = 1.0*np.ones(nhru)
    celerity1 = np.ones(nhru)
    area = 1.0*np.ones(nhru)
    dx = 1.0*np.ones(nhru)
    dt = 2.0
    ncores = 1
    maxntt = 1
    isw = 0.5
    flow_matrix = scipy.sparse.csr_matrix(np.array([[0.,1.,0],[0.,0.,1.],[0.,0.,0.]]).T)
    flow_matrix.setdiag(flow_matrix.diagonal())
    #Initialize the solver
    dtt.initialize(flow_matrix.indices,flow_matrix.indptr)
    #Solve the system of equations
    dtt.update(recharge,storage,qout,qin,
             recharge1,storage1,qout1,qin1,
             area,dx,dt,celerity,celerity1,storage_mask_subsurface,
             flow_matrix.data,flow_matrix.indices,flow_matrix.indptr,
             ncores,maxntt,isw)
    #Finalize the solver
    dtt.finalize()
    #Compare
    qout_true = np.array([1.0,1.5,1.75])
    self.assertTrue(np.allclose(qout,qout_true,rtol=1e-05, atol=1e-08))

  def test_kinematic_wave_solution_mkl_implicit_3hru_returnflow3(self):

    nhru = 3
    storage_mask_subsurface = np.ones(nhru)
    qout = np.zeros(nhru)
    qin = np.zeros(nhru)
    storage1 = np.zeros(nhru)
    qin1 = np.zeros(nhru)
    qin_outlet = np.zeros(nhru)
    area_outlet = np.zeros(nhru)
    storage = 0.0*np.ones(nhru)
    qout1 = np.ones(nhru)
    recharge1 = 1.0*np.ones(nhru)
    recharge = 1.0*np.ones(nhru)
    celerity = 1.0*np.ones(nhru)
    celerity1 = np.ones(nhru)
    area = 1.0*np.ones(nhru)
    dx = 1.0*np.ones(nhru)
    dt = 2.0
    ncores = 1
    maxntt = 1
    isw = 0.5
    flow_matrix = scipy.sparse.csr_matrix(np.array([[0.,1.,0],[0.,0.,1.],[0.,0.,0.5]]).T)
    flow_matrix.setdiag(flow_matrix.diagonal())
    #Initialize the solver
    dtt.initialize(flow_matrix.indices,flow_matrix.indptr)
    #Solve the system of equations
    dtt.update(recharge,storage,qout,qin,
             recharge1,storage1,qout1,qin1,
             area,dx,dt,celerity,celerity1,storage_mask_subsurface,
             flow_matrix.data,flow_matrix.indices,flow_matrix.indptr,
             ncores,maxntt,isw)
    #Finalize the solver
    dtt.finalize()
    #Compare
    qout_true = np.array([1.0,1.5,2.333333])
    self.assertTrue(np.allclose(qout,qout_true,rtol=1e-05, atol=1e-08))

  def test_kinematic_wave_solution_mkl_implicit_3hru_returnflow23(self):

    nhru = 3
    storage_mask_subsurface = np.ones(nhru)
    qout = np.zeros(nhru)
    qin = np.zeros(nhru)
    storage1 = np.zeros(nhru)
    qin1 = np.zeros(nhru)
    qin_outlet = np.zeros(nhru)
    area_outlet = np.zeros(nhru)
    storage = 0.0*np.ones(nhru)
    qout1 = np.ones(nhru)
    recharge1 = 1.0*np.ones(nhru)
    recharge = 1.0*np.ones(nhru)
    celerity = 1.0*np.ones(nhru)
    celerity1 = np.ones(nhru)
    area = 1.0*np.ones(nhru)
    dx = 1.0*np.ones(nhru)
    dt = 2.0
    ncores = 1
    maxntt = 1
    isw = 0.5
    flow_matrix = scipy.sparse.csr_matrix(np.array([[0.,1.,0],[0.,0.5,1.],[0.,0.,0.5]]).T)
    flow_matrix.setdiag(flow_matrix.diagonal())
    #Initialize the solver
    dtt.initialize(flow_matrix.indices,flow_matrix.indptr)
    #Solve the system of equations
    dtt.update(recharge,storage,qout,qin,
             recharge1,storage1,qout1,qin1,
             area,dx,dt,celerity,celerity1,storage_mask_subsurface,
             flow_matrix.data,flow_matrix.indices,flow_matrix.indptr,
             ncores,maxntt,isw)
    #Finalize the solver
    dtt.finalize()
    #Compare
    qout_true = np.array([1.0,2.0,2.6666666])
    self.assertTrue(np.allclose(qout,qout_true,rtol=1e-05, atol=1e-08))

  def test_kinematic_wave_solution_mkl_implicit_3hru_returnflow123(self):

    nhru = 3
    storage_mask_subsurface = np.ones(nhru)
    qout = np.zeros(nhru)
    qin = np.zeros(nhru)
    storage1 = np.zeros(nhru)
    qin1 = np.zeros(nhru)
    qin_outlet = np.zeros(nhru)
    area_outlet = np.zeros(nhru)
    storage = 0.0*np.ones(nhru)
    qout1 = np.ones(nhru)
    recharge1 = 1.0*np.ones(nhru)
    recharge = 1.0*np.ones(nhru)
    celerity = 1.0*np.ones(nhru)
    celerity1 = np.ones(nhru)
    area = 1.0*np.ones(nhru)
    dx = 1.0*np.ones(nhru)
    dt = 2.0
    ncores = 1
    maxntt = 1
    isw = 0.5
    flow_matrix = scipy.sparse.csr_matrix(np.array([[0.5,1.,0],[0.,0.5,1.],[0.,0.,0.5]]).T)
    flow_matrix.setdiag(flow_matrix.diagonal())
    #Initialize the solver
    dtt.initialize(flow_matrix.indices,flow_matrix.indptr)
    #Solve the system of equations
    dtt.update(recharge,storage,qout,qin,
             recharge1,storage1,qout1,qin1,
             area,dx,dt,celerity,celerity1,storage_mask_subsurface,
             flow_matrix.data,flow_matrix.indices,flow_matrix.indptr,
             ncores,maxntt,isw)
    #Finalize the solver
    dtt.finalize()
    #Compare
    qout_true = np.array([1.333333333,2.2222222222,2.81481481481])
    self.assertTrue(np.allclose(qout,qout_true,rtol=1e-05, atol=1e-08))


  def test_kinematic_wave_solution_mkl_implicit_3hru_backflow(self):

    nhru = 3
    storage_mask_subsurface = np.ones(nhru)
    qout = np.zeros(nhru)
    qin = np.zeros(nhru)
    storage1 = np.zeros(nhru)
    qin1 = np.zeros(nhru)
    qin_outlet = np.zeros(nhru)
    area_outlet = np.zeros(nhru)
    storage = 0.0*np.ones(nhru)
    qout1 = np.ones(nhru)
    recharge1 = 1.0*np.ones(nhru)
    recharge = 1.0*np.ones(nhru)
    celerity = 1.0*np.ones(nhru)
    celerity1 = np.ones(nhru)
    area = 1.0*np.ones(nhru)
    dx = 1.0*np.ones(nhru)
    dt = 2.0
    ncores = 1
    maxntt = 1
    isw = 0.5
    flow_matrix = scipy.sparse.csr_matrix(np.array([[0.,1.,0],[0.5,0.,0.5],[0.,0.,0]]).T)
    flow_matrix.setdiag(flow_matrix.diagonal())
    #Initialize the solver
    dtt.initialize(flow_matrix.indices,flow_matrix.indptr)
    #Solve the system of equations
    dtt.update(recharge,storage,qout,qin,
             recharge1,storage1,qout1,qin1,
             area,dx,dt,celerity,celerity1,storage_mask_subsurface,
             flow_matrix.data,flow_matrix.indices,flow_matrix.indptr,
             ncores,maxntt,isw)
    #Finalize the solver
    dtt.finalize()
    #Compare
    qout_true = np.array([1.42857142857,1.71428571429,1.42857142857])
    self.assertTrue(np.allclose(qout,qout_true,rtol=1e-05, atol=1e-08))

suite = unittest.TestLoader().loadTestsFromTestCase(DynamicTopmodel)
unittest.TextTestRunner(verbosity=2).run(suite)
