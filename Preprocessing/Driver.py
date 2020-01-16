import datetime
import Preprocessing
import sys
import time

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

def Read_Metadata_File(file):

 import json
 metadata = json.load(open(file))['Preprocessing']

 return metadata


#Read in the metadata file
metadata_file = sys.argv[1]
metadata = Read_Metadata_File(metadata_file)
info = metadata
info['idate'] = datetime.datetime(metadata['startdate']['year'],
                           metadata['startdate']['month'],
                           metadata['startdate']['day'],0)
info['fdate'] = datetime.datetime(metadata['enddate']['year'],
                           metadata['enddate']['month'],
                           metadata['enddate']['day'],0) + datetime.timedelta(days=1) - datetime.timedelta(seconds=info['dt'])


#Cluster the data
tic = time.time()
Preprocessing.Prepare_Model_Input_Data(info)
print("Elapsed time: ",time.time() - tic)

