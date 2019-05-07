import numpy as np
import matplotlib.pyplot as plt
import time
import netCDF4 as nc
import scipy.sparse as sparse
import scipy.sparse.linalg
import pickle

#Read in the stream network information
file = '../../../ReynoldsCreek/input_file.nc'
fp = nc.Dataset(file)
grp = fp['stream_network']
dbc = {}
for var in grp.variables:
 dbc[var] = grp[var][:]

#Read in the runoff output
file = '/Users/nchaney/Projects/hydroblocks_jan2019/ReynoldsCreek/2015-01-01.nc'
fp = nc.Dataset(file)
runoff = fp['data']['runoff'][:]

#Read in reach/hru area
file = '/Users/nchaney/Projects/hydroblocks_jan2019/ReynoldsCreek/test.pck'
db = pickle.load(open(file,'rb'))

#Create a matrix
reach2hru = np.zeros((dbc['topology'][:].size,runoff.shape[1]))
for reach in db:
 for hru in db[reach]:
  reach2hru[reach-1,hru] = db[reach][hru]
reach2hru = sparse.csr_matrix(reach2hru)
area = np.sum(reach2hru,axis=0)
#pct = area/np.sum(area)
#rf = np.sum(pct*runoff,axis=1)

#Assemble the connectivity array (Best to define this reordering in the database creation)
corg = np.arange(dbc['topology'].size)
cdst = dbc['topology'][:]
m = cdst != -1
#cdst[m] = np.abs(cdst[m] - np.max(corg))
#corg = np.abs(corg - np.max(corg))
nc = cdst.size
cmatrix = sparse.coo_matrix((np.ones(cdst[m].size),(corg[m],cdst[m])),shape=(nc,nc),dtype=np.float32)
cmatrix = cmatrix.tocsr().T

#Identify headwaters
cup = -1*np.ones(cdst.size)
for id in corg:
 n = np.sum(cdst == id)
 if n > 0:cup[id] = 1

#import copy
c_length = dbc['length'][:]
c_slope = dbc['slope'][:]
c_width = dbc['width'][:]
c_n = dbc['manning'][:]
np.random.seed(1)
c_width[:] = np.random.choice(np.linspace(1,1000.0,c_width.size))#100.0
c_n[:] = np.random.choice(np.linspace(0.001,0.1,c_width.size))
c_slope[:] = np.random.choice(np.linspace(0.001,0.1,c_width.size))
c_length[:] = np.random.choice(np.linspace(1.0,1000.0,c_width.size))
hinit = np.zeros(c_length.size)
#hinit[cup == -1] = 0.1
hinit[:] = 0.0
h0 = np.copy(hinit)
h1 = np.copy(hinit)

#Filler
dt = 10800 #s
#tmax = 100*3600*24
tmax = dt*runoff.shape[0]
nt = int(tmax/dt)
#Define initial conditions
Qinit = np.zeros(nc)
Qinit[:] = c_width*h0**(5.0/3.0)*c_slope**0.5/c_n
Q0 = Qinit[:]
qin = np.zeros(c_length.size)
qout = np.zeros(c_length.size)

print('Q0',Q0)
out = {'Q':[],'h':[],'qin':[],'qout':[]}
dif0 = -9999
max_niter = 2
for t in range(nt):
 h0_org = np.copy(h0)
 Q0_org = np.copy(Q0)
 for it in range(max_niter):
  qin[:] = reach2hru.dot(runoff[t,:]/1000.0/dt)/c_length/c_width #m/s
  qin[qin < 0] = 0.0
  print(t,it)
  h = np.zeros(h0.size)
  u = h0**(2.0/3.0)*c_slope**0.5/c_n
  #Fill non-diagonals
  A = cmatrix.multiply(-c_width*dt*u)
  #Fill diagonal
  A.setdiag(c_width*c_length + c_width*dt*u)
  #Set right hand side
  b = c_width*c_length*h0_org + dt*qin*c_width*c_length - dt*qout*c_width*c_length
  #Ax = b
  h1 = scipy.sparse.linalg.spsolve(A,b,use_umfpack=True)
  #QC
  h0[h0 < 0] = 0.0
  dif1 = np.mean(np.abs(h0 - h1))
  if (dif1 < 10**-10) | (it == max_niter-1):
   #Reset h0
   h0[:] = h1[:]
   u = h0**(2.0/3.0)*c_slope**0.5/c_n
   #Calculate Q1
   Q1 = c_width*h0**(5.0/3.0)*c_slope**0.5/c_n
   dif0 = -9999
   break
  else:
   #Reset h0
   h0[:] = h1[:]
   dif0 = dif1 
   Q1 = c_width*h0**(5.0/3.0)*c_slope**0.5/c_n
 #Append to output
 out['Q'].append(np.copy(Q1))
 out['h'].append(np.copy(h1))
 out['qin'].append(np.copy(qin))
 out['qout'].append(np.copy(qout))
for var in out:
 out[var] = np.array(out[var])
m = cdst == -1
#print(np.sum(np.trapz(out['Q'][:,m].T,dx=dt)),"m3")
dVh = -np.sum(c_length*c_width*np.diff(out['h'],axis=0),axis=1)
dVh += np.sum(c_length*c_width*dt*out['qin'],axis=1)[1:]
dVh -= np.sum(c_length*c_width*dt*out['qout'],axis=1)[1:]
dVQ = dt*np.sum(out['Q'][1:,m],axis=1)
print(np.sum(dVh),np.sum(dVQ))
print(np.sum(dVQ)/np.sum(area))
plt.subplot(121)
plt.plot(np.log10(out['Q'][:,0]))
plt.legend(['C1','C2'])
plt.subplot(122)
plt.plot(out['h'][:,0])
plt.show()
exit()
