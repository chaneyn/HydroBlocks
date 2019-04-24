import numpy as np
import matplotlib.pyplot as plt
import time
import netCDF4 as nc
import scipy.sparse as sparse
import scipy.sparse.linalg

def solve_pde(Q0,h0,ql,wc,n,dx,dt,nx,Q_lbc,Q_rbc,max_niter,slope):
 
 C = np.zeros(Q0.size)
 A = np.zeros((nx,nx))
 b = np.zeros((nx))
 Q0_org = np.copy(Q0)
 h0_org = np.copy(h0)
 Q1 = np.copy(Q0)
 h1 = np.copy(h0)
 for it in range(max_niter):
  #Zero out arrays
  A[:] = 0.0
  C[:] = 0.0
  b[:] = 0.0
  #Define the flows
  #Compute conveyance
  dKdh = 5.0/3.0*wc*h0**(2.0/3.0)/n
  #Compute celerity
  C[:] = slope**0.5/wc*dKdh #Use value from previous time step
  #Create A matrix
  #np.fill_diagonal(A[1:,:],-C[1:]/dx)
  np.fill_diagonal(A[1:,:],-C[1:]/dx)
  np.fill_diagonal(A,1/dt + C/dx)
  #Create b array
  #Include boundary conditions 
  b[:] = Q0_org[:]/dt + C*ql[1:-1]
  b[0] = b[0] + C[0]/dx*Q_lbc
  #b[-1] = b[-1] + C[-1]/dx*Q_rbc
  #Compute Q1
  Q1[:] = np.linalg.solve(A,b)
  #Update the water balance
  h1[0] = h0_org[0] + dt*(Q_lbc/dx/wc - Q1[0]/dx/wc)
  h1[1:] = h0_org[1:] + dt*(Q1[0:-1]/dx/wc - Q1[1:]/dx/wc)
  #Reset h0
  h0[:] = h1[:]

 return (Q1,h1)

def numerical_solution_btuw(x,ts,dx,dt,max_niter,wc,n,Q_lbc,Q_rbc,Q0,h0,ql,slope):

 nx = x.size
 output = [np.copy(Q0),]
 output = {'Q':[np.copy(Q0),],'h':[np.copy(h0),]}
 for t in ts[1:]:
  #Reset h0
  if t == ts[1]:
   h0[:] = h0[:] + dt*ql[1:-1]/wc
  #tic = time.time()
  (Q1,h1) = solve_pde(Q0,h0,ql,wc,n,dx,dt,nx,Q_lbc,Q_rbc,max_niter,slope) 
  #Reset Q0
  Q0[:] = Q1[:]
  #Reset h0
  h0[:] = h1[:]# + dt*ql[1:-1]/wc
  ql[:] = 0.0
  #Add output
  output['Q'].append(np.copy(Q0))
  output['h'].append(np.copy(h0))

 for var in output:
  output[var] = np.array(output[var])

 return output

#Read in the stream network information
file = '../../../ReynoldsCreek/input_file.nc'
fp = nc.Dataset(file)
grp = fp['stream_network']
dbc = {}
for var in grp.variables:
 dbc[var] = grp[var][:]

#Assemble the connectivity array (Best to define this reordering in the database creation)
corg = np.arange(dbc['topology'].size)
cdst = dbc['topology'][:]
m = cdst != -1
#cdst[m] = np.abs(cdst[m] - np.max(corg))
#corg = np.abs(corg - np.max(corg))
nc = cdst.size
cmatrix = sparse.coo_matrix((np.ones(cdst[m].size),(corg[m],cdst[m])),shape=(nc,nc),dtype=np.float32)
cmatrix = cmatrix.tocsr().T
'''corg = np.array([0,1,2])
cdst = np.array([-1,0,1])
m = cdst != -1
nc = cdst.size
cmatrix = sparse.coo_matrix((np.ones(cdst[m].size),(corg[m],cdst[m])),shape=(nc,nc),dtype=np.float32)
cmatrix = cmatrix.tocsr().T'''

#Identify headwaters
cup = -1*np.ones(cdst.size)
for id in corg:
 n = np.sum(cdst == id)
 if n > 0:cup[id] = 1

#import copy
c_length = dbc['length'][:]
c_slope = dbc['slope'][:]
c_manning = dbc['manning'][:]
c_width = dbc['width'][:]
c_n = dbc['manning'][:]
c_n[:] = 0.01
c_width[:] = 100.0
c_slope[:] = 0.001
c_length[:] = 1000.0
#c_length = 1000.0*np.ones(cdst.size)
#c_n = 0.1*np.ones(cdst.size)
#c_width = 1000*np.ones(cdst.size)
#c_slope = 0.01*np.ones(cdst.size)
hinit = np.zeros(c_length.size)
#hinit[cup == -1] = 10.0
hinit[:] = 0.01
h0 = np.copy(hinit)
h1 = np.copy(hinit)
#c_width[:] = 1.0

#Filler
tmax = 10**6
dt = 3600 #s
nt = int(tmax/dt)
#Define initial conditions
Qinit = np.zeros(nc)
Qinit[:] = c_width*h0**(5.0/3.0)*c_slope**0.5/c_n
#Qinit[:] = 10**-1
#Qinit[:] = 0.1
Q0 = Qinit[:]
print('Q0',Q0)
out = {'Q':[],'h':[]}
dif0 = -9999
max_niter = 100
for t in range(nt):
 h0_org = np.copy(h0)
 Q0_org = np.copy(Q0)
 for it in range(max_niter):
  print(t,it)
  #Compute conveyance
  dKdh = np.zeros(h0.size)
  m = h0 > 0
  dKdh[m] = 5.0/3.0*c_width[m]*h0[m]**(2.0/3.0)/c_n[m]
  Q0_org[~m] = 0.9*Q0_org[~m]
  #dKdh = 5.0/3.0*c_width*h0**(2.0/3.0)/c_n
  #Q0_org[~m] = 0.01*Q0_org[~m]#h0_org[~m]*c_length[~m]*c_width[~m]
  #Compute celerity
  c = c_slope**0.5/c_width*dKdh #Use value from previous time step
  #c = 0.1 #m/s
  #c = 0.001
  #Fill non-diagonals
  A = cmatrix.multiply(-c/c_length)
  #Fill diagonal
  A.setdiag(1.0/dt + c/c_length)
  #Define Q_nm1
  b = Q0_org/dt
  #Ax = b
  Q1 = scipy.sparse.linalg.spsolve(A,b,use_umfpack=True)
  #Update h
  B = cmatrix.multiply(1)
  B.setdiag(-1)
  #Find Q that lead to negative water heights and correct them
  #tmp = h0_org[:] + dt*np.array(B.multiply(Q1).sum(axis=1))[:,0]/c_length/c_width
  tmp = h0_org + -dt*Q1/c_length/c_width
  m = tmp < 0
  Q1[m] = c_length[m]*c_width[m]/dt*h0_org[m]
  #Compute flux in/out
  dQ = np.array(B.multiply(Q1).sum(axis=1))[:,0]
  h1[:] = h0_org[:] + dt*dQ/c_length/c_width
  #QC
  #print(np.mean(np.abs(Q0 - Q1)),np.min(h1))
  dif1 = np.mean(np.abs(Q0 - Q1))
  if (dif1 < 10**-5) | (dif1 == dif0):
   #Reset h0,Q0
   h0[:] = h1[:]
   Q0[:] = Q1[:]
   dif0 = -9999
   break
  else:
   #Reset h0,Q0
   h0[:] = h1[:]
   Q0[:] = Q1[:]
   dif0 = dif1
  #if(np.sum(h0 < 0) == 0):break
 #Append to output
 out['Q'].append(np.copy(Q1))
 out['h'].append(np.copy(h1))
for var in out:
 out[var] = np.array(out[var])
m = cdst == -1
#print(np.sum(np.trapz(out['Q'][:,m].T,dx=dt)),"m3")
print(dt*np.sum(out['Q'][:,m]))
print(np.sum((hinit-out['h'][-1,:])*c_length*c_width))
plt.subplot(121)
plt.plot(out['Q'])#[:,0])
plt.legend(['C1','C2'])
plt.subplot(122)
plt.plot(out['h'])#[:,0])
plt.show()
exit()

dt = 450.0
dx = 1000
x = np.arange(dx,100000,dx)
t = np.arange(0,dt*10+dt,dt)
wc = 100 #m (channel width)
n = 10**-2 #manning coeffient
Q_lbc = 0.0
Q_rbc = 0.0
Q0 = np.zeros(x.size)
h0 = np.zeros(x.size)
ql = np.zeros(x.size+2)
Q0[:] = 0.0
h0[:] = 0.0
#ql[1] = 10.0/2/dx
#ql[10:15] = 5.00/dx#0.003/2/dx
ql[1:20] = 1.00/dx*np.sin(np.linspace(0,np.pi,19))
print('Input (m3)',dt*np.trapz(ql,dx=dx))
slope = np.zeros(x.size)
slope[:] = 10**-1
it = 10#t.size-1
max_niter = 100
output = numerical_solution_btuw(x,t,dx,dt,max_niter,wc,n,Q_lbc,Q_rbc,Q0,np.copy(h0),ql,slope)
#Compute the water under the curve
for t in range(output['h'].shape[0]):
 print(t,np.trapz(wc*output['h'][t,:],dx=dx))
plt.plot(output['h'][:,:].T)
plt.show()
