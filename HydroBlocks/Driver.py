import datetime
import HydroBlocks# as HB
import sys
import cPickle as pickle

def Read_Metadata_File(file):

 import json
 metadata = json.load(open(file))

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
                           metadata['enddate']['day'],23)

#Initialize
HB = HydroBlocks.initialize(info)

#Run the model
HB.run(info)

#Finalize
HB.finalize()
