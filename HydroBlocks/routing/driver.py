import glob
import os
#cwd = os.getcwd()
#os.chdir('/home/nc153/soteria/projects/hydroblocks_inter_catchment/run_marion_county/input_data/catchments')
cids = glob.glob('/home/nc153/soteria/projects/hydroblocks_inter_catchment/run_marion_county/input_data/catchments/*')
for cid in cids:
 cid = cid.split('/')[-1]
 if 'catchment' in cid:continue
 print(cid)
 os.system('python driver_core.py %s' % (cid,))
