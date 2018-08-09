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

#Define idate and fdate
idate = datetime.datetime(metadata['startdate']['year'],metadata['startdate']['month'],metadata['startdate']['day'],0)
fdate = datetime.datetime(metadata['enddate']['year'],metadata['enddate']['month'],metadata['enddate']['day'],0) + datetime.timedelta(days=1)

#Run the segments for the model
sidate = idate
sfdate = idate
while sidate < fdate:
 sfdate = sidate + relativedelta(years=metadata['segment']['years_per_segment'])
 if sfdate > fdate:sfdate = fdate
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



