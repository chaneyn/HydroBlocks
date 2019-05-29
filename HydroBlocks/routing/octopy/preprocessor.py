import os
import glob
import json
import netCDF4 as nc
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
import scipy.sparse as sparse
import scipy.sparse.linalg
import pickle
import numba
import time
from mpi4py import MPI

#Get some general info
dir = os.getcwd()

#Determine communication information
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()

rdir = '/home/nc153/soteria/projects/hydroblocks_inter_catchment/regions/ohio_basin'
#Define the catchment
tmp = glob.glob('%s/input_data/domain/*' % rdir)
tmp.remove('%s/input_data/domain/domain_database.pck' % rdir)
cids = []
for cid in tmp:
 cids.append(int(cid.split('/')[-1]))
cid = cids[rank]

#Determine catchments that rely on this catchment (ids to send to)
scids = []
for ucid in cids:
 if ucid == cid:continue
 file = '%s/input_data/domain/%d/input_file_enhanced_test.nc' % (rdir,ucid)
 fp = nc.Dataset(file)
 grp = fp['stream_network']
 ucids = np.unique(grp['cid'][:])
 fp.close()
 if cid in ucids:scids.append(ucid)

#Read in the stream network information
file = '%s/input_data/domain/%d/input_file_enhanced_test.nc' % (rdir,cid)
fp = nc.Dataset(file)
grp = fp['stream_network']
dbc = {}
for var in grp.variables:
 dbc[var] = grp[var][:]
ucids = np.unique(dbc['cid'])
fp.close()

#Define ids to receive from
rcids = list(ucids)
rcids.remove(cid)

#Read in the other info
odbc = {}
for ucid in ucids:
 file = '%s/input_data/domain/%d/input_file_enhanced_test.nc' % (rdir,ucid)
 fp = nc.Dataset(file)
 grp = fp['stream_network']
 odbc[ucid] = {}
 for var in ['channelid','cid']:
  odbc[ucid][var] = grp[var][:]

#Construct mapping of other catchment network to that of the current network
mapping_ucid = {}
for ucid in ucids:
 mapping_ucid[ucid] = {'cid':[],'ocid':[]}
 idxs = np.where(dbc['cid'] == ucid)[0]
 for idx in idxs:
  m = (odbc[ucid]['cid'] == ucid) & (odbc[ucid]['channelid'] == dbc['channelid'][idx])
  mapping_ucid[ucid]['cid'].append(idx)
  mapping_ucid[ucid]['ocid'].append(np.where(m)[0][0])
 #Convert to arrays
 for var in mapping_ucid[ucid]:
  mapping_ucid[ucid][var] = np.array(mapping_ucid[ucid][var])

#Assemble database
odb = {'mapping_ucid':mapping_ucid}
pickle.dump(odb,open('input/%s.pck' % cid,'wb'))
