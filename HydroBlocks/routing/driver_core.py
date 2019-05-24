import numpy as np
import matplotlib.pyplot as plt
import time
import netCDF4 as nc
import scipy.sparse as sparse
import scipy.sparse.linalg
import pickle
import numba
import sys

@numba.jit(nopython=True,cache=True)
def calculate_hydraulic_radius(A,P,W,A1):

 Rh = np.zeros(A.shape[0])
 for i in range(A.shape[0]):
  A0 = 0.0
  P0 = 0.0
  for j in range(A.shape[1]):
   W1 = W[i,j]
   if A[i,j] > A1[i]:break
   if A[i,j+1] == 0.0:break
   A0 = A[i,j]
   P0 = P[i,j]
  #Calculate the height above the segment
  h = (A1[i] - A0)/W1
  #Calculate the wetted perimeter
  P1 = P0 + 2*h + W1
  #Calculate hydraulic radius
  Rh[i] = A1[i]/P1

 return Rh

#cid = 7100454060
cid = int(sys.argv[1])#7100454060
#Read in the stream network information
file = '/home/nc153/soteria/projects/hydroblocks_inter_catchment/run_marion_county/input_data/catchments/%d/input_file_enhanced.nc' % cid
fp = nc.Dataset(file)
grp = fp['stream_network']
dbc = {}
for var in grp.variables:
 dbc[var] = grp[var][:]
ucids = np.unique(dbc['cid'])

#Create basin mapping
basin_mapping = {}
for ucid in ucids:
 basin_mapping[ucid] = {}
 m = dbc['cid'] == ucid
 bs = np.where(m)[0]
 for b in bs:
  ob = dbc['channelid'][b]+1
  basin_mapping[ucid][ob] = b+1

#Read in the runoff output
#file = '/Users/nchaney/Projects/hydroblocks_jan2019/ReynoldsCreek/2015-01-01.nc'
#file = '../../../hydroblocks_inter_catchment/run_marion_county/output_data/%d/output/2002-01-01.nc' % cid
#Read in the runoff output
nhru = 0
hru_mapping = {}
runoff = {}
for ucid in ucids:
 hru_mapping[ucid] = {}
 file = '/home/nc153/soteria/projects/hydroblocks_inter_catchment/run_marion_county/output_data/%s/output/2002-01-01.nc' % ucid
 fp = nc.Dataset(file)
 runoff[ucid] = fp['data']['runoff'][:]
 #Add to hru mapping
 for hru in range(runoff[ucid].shape[1]):
  hru_mapping[ucid][hru] = hru + nhru
 nhru += runoff[ucid].shape[1]
 nt = runoff[ucid].shape[0]
 fp.close()

#Create a new mapping for the hrus
hrus = np.zeros((nhru,2))

#Map the runoff to the new hru mapping and create new hru maping
tmp = np.zeros((nt,nhru))
imin = 0
for ucid in ucids:
 imax = imin + runoff[ucid].shape[1]
 tmp[:,imin:imax] = runoff[ucid]
 hrus[imin:imax,0] = np.arange(runoff[ucid].shape[1])
 hrus[imin:imax,1] = ucid
 imin = imax
runoff = tmp[:]

#Read in reach/hru area
#file = '/Users/nchaney/Projects/hydroblocks_jan2019/ReynoldsCreek/test.pck'
db = {}
for ucid in ucids:
 file = '/home/nc153/soteria/projects/hydroblocks_inter_catchment/run_marion_county/input_data/catchments/%d/routing_info.pck' % ucid
 odb = pickle.load(open(file,'rb'))['reach_hru_area']
 for ob in odb:
  if ob not in basin_mapping[ucid]:continue
  b = basin_mapping[ucid][ob]
  db[b] = {}
  #Map to new basin id
  for ohru in odb[ob]:
   #Map to new hru id
   hru = hru_mapping[ucid][ohru]
   db[b][hru] = odb[ob][ohru]

#Read in the reach/hand database and remap to new topology (THERE COULD BE AN ORDERING PROBLEM)
hdb = {}
for ucid in ucids:
 file = '/home/nc153/soteria/projects/hydroblocks_inter_catchment/run_marion_county/input_data/catchments/%d/routing_info.pck' % ucid
 ohdb = pickle.load(open(file,'rb'))['reach_cross_section']
 for var in ohdb:
  if var not in hdb:hdb[var] = np.zeros((dbc['topology'].size,ohdb[var].shape[1]))
  m = dbc['cid'] == ucid
  hdb[var][m,:] = ohdb[var][dbc['channelid'][m],:]

#Create a matrix
reach2hru = np.zeros((dbc['topology'][:].size,runoff.shape[1]))
for reach in db:
 for hru in db[reach]:
  reach2hru[reach-1,hru] = db[reach][hru]
reach2hru = sparse.csr_matrix(reach2hru)
area = np.sum(reach2hru,axis=0)

#Assemble the connectivity array (Best to define this reordering in the database creation)
corg = np.arange(dbc['topology'].size)
cdst = dbc['topology'][:]
m = cdst != -1
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
#np.random.seed(1)
#c_width[:] = np.random.choice(np.linspace(1,1000.0,c_width.size))#100.0
#c_n[:] = np.random.choice(np.linspace(0.001,0.1,c_width.size))
#c_slope[:] = np.random.choice(np.linspace(0.001,0.1,c_width.size))
#c_length[:] = np.random.choice(np.linspace(1.0,1000.0,c_width.size))
#c_width[:] = 1.0#np.random.choice(np.linspace(1,1000.0,c_width.size))#100.0
#c_n[:] = 0.04#np.random.choice(np.linspace(0.001,0.1,c_width.size))
#c_slope[:] = 0.01#np.random.choice(np.linspace(0.001,0.1,c_width.size))
#c_length[:] = 1000.0#np.random.choice(np.linspace(1.0,1000.0,c_width.size))
Ainit = np.zeros(c_length.size)
#hinit[cup == -1] = 0.1
Ainit[:] = 0.1
#Ainit[cup == -1] = 0.1
A0 = np.copy(Ainit)
A1 = np.copy(Ainit)

#Filler
dt = 10800 #s
#tmax = 100*3600*24
tmax = dt*250#runoff.shape[0]
#maxc = 10.0 #m/s
nt = int(tmax/dt)
#Define initial conditions
Qinit = np.zeros(nc)
Q0 = Qinit[:]
qin = np.zeros(c_length.size)
qout = np.zeros(c_length.size)

out = {'Q':[],'A':[],'qin':[],'qout':[]}
dif0 = -9999
max_niter = 10
for t in range(nt):
 print(cid,t)
 A0_org = np.copy(A0)
 Q0_org = np.copy(Q0)
 for it in range(max_niter):
  #Compute inflows
  #qin[:] = reach2hru.dot(runoff[t,:]/1000.0/dt)/c_length #m/s
  #if t == 0:qin[:] = 10**-6
  #else:qin[:] = 0.0
  qin[qin < 0] = 0.0
  #Determine hydraulic radius
  rh = calculate_hydraulic_radius(hdb['A'],hdb['P'],hdb['W'],A0)
  #rh = A0/c_width
  #Determine velocity
  u = rh**(2.0/3.0)*c_slope**0.5/c_n
  #u[:] = 0.1
  #Define maximum velocity from maximum celerity (5/3u for now;note that this not true)
  #c = 5.0/3.0*u
  #u[c > maxc] = 3.0/5.0*maxc
  #Fill non-diagonals
  LHS = cmatrix.multiply(-dt*u)
  #Fill diagonal
  LHS.setdiag(c_length + dt*u)
  #Set right hand side
  RHS = c_length*A0_org + dt*qin*c_length - dt*qout*c_length
  #Ax = b
  A1 = scipy.sparse.linalg.spsolve(LHS,RHS,use_umfpack=True)
  #QC
  A0[A0 < 0] = 0.0
  dif1 = np.mean(np.abs(A0 - A1))
  if (dif1 < 10**-10) | (it == max_niter-1):
   #Reset A0
   A0[:] = A1[:]
   h = A0/c_length #rectangular channel
   #Determine hydraulic radius
   #rh = A0/c_width
   rh = calculate_hydraulic_radius(hdb['A'],hdb['P'],hdb['W'],A0)
   #Determine velocity
   u = rh**(2.0/3.0)*c_slope**0.5/c_n
   #u[:] = 0.1
   #Define maximum velocity from maximum celerity (5/3u for now;note that this not true)
   #c = 5.0/3.0*u
   #u[c > maxc] = 3.0/5.0*maxc
   #Calculate Q1
   Q1 = A0*u
   dif0 = -9999
   break
  else:
   #Reset A0
   A0[:] = A1[:]
   dif0 = dif1 

 #Append to output
 out['Q'].append(np.copy(Q1))
 out['A'].append(np.copy(A1))
 out['qin'].append(np.copy(qin))
 out['qout'].append(np.copy(qout))
for var in out:
 out[var] = np.array(out[var])
m = cdst == -1
#print(np.sum(np.trapz(out['Q'][:,m].T,dx=dt)),"m3")
dVh = -np.sum(c_length*np.diff(out['A'],axis=0),axis=1)
dVh += np.sum(c_length*dt*out['qin'],axis=1)[1:]
dVh -= np.sum(c_length*dt*out['qout'],axis=1)[1:]
dVQ = dt*np.sum(out['Q'][1:,m],axis=1)
print(np.sum(dVh),np.sum(dVQ))
print(np.sum(dVQ)/np.sum(area))
'''plt.subplot(121)
plt.plot(out['Q'][:,:])
#plt.plot(np.log10(out['Q'][:,0]))
plt.legend(['C1','C2'])
plt.subplot(122)
plt.plot(out['A'][:,:])
#plt.plot(out['A'][:,0])
plt.show()'''
pickle.dump(out,open('workspace/%s.pck' % cid,'wb'))
