import datetime
from dateutil.relativedelta import relativedelta
import HydroBlocks# as HB
import sys
import pickle

def Read_Metadata_File(file):

 import json
 metadata = json.load(open(file))['HydroBlocks']

 return metadata

#Read in the metadata file
metadata_file = sys.argv[1]
metadata = Read_Metadata_File(metadata_file)
info = metadata
info["dz"] = [0.05, 0.1, 0.1, 0.1, 0.2, 0.2, 0.25, 0.25, 0.25, 0.5,1.0,2.0]
#info['enddate']['year'] = 2002
#info['enddate']['month'] = 2
#info['enddate']['day'] = 1
#info['subsurface_module'] = 'none'
#info['segment']['years_per_segment'] = 2017-2002+1

#Define idate and fdate
idate = datetime.datetime(metadata['startdate']['year'],metadata['startdate']['month'],metadata['startdate']['day'],0)
fdate = datetime.datetime(metadata['enddate']['year'],metadata['enddate']['month'],metadata['enddate']['day'],0) + datetime.timedelta(days=1)

#Run the segments for the model
sidate = idate
sfdate = idate
while sidate < fdate:
 sfdate = sidate + relativedelta(years=metadata['segment']['years_per_segment'])
 if sfdate > fdate: sfdate = fdate
 #Set the parameters
 info['idate'] = sidate
 info['fdate'] = sfdate
 #Run the model
 #Initialize
 HB = HydroBlocks.initialize(info)
 #Run the model
 HB.run(info)
 #Finalize
 HB.finalize()
 #Update initial time step
 sidate = sfdate



