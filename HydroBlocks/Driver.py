import datetime
from dateutil.relativedelta import relativedelta
import HydroBlocks# as HB
import sys
import cPickle as pickle

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

'''
#Define the info
info = {
    'icatch':metadata['catchment_id'],
    'input_file':metadata['input_file'],
    'output_file':metadata['output_file'],
    #'soil_file':metadata['soil_file'],
    'workspace':metadata['workspace'],
    'surface_flow_flag':metadata['surface_flow_flag'],
    'subsurface_module':metadata['subsurface_module'],
    "hwu_flag":metadata['hwu_flag'],
    "hwu_sf_flag":metadata['hwu_sf_flag'],
    "hwu_gw_flag":metadata['hwu_gw_flag'],
    "hwu_agric_flag":metadata['hwu_agric_flag'],
    "hwu_domest_flag":metadata['hwu_domest_flag'],
    "hwu_indust_flag":metadata['hwu_indust_flag'],
    "hwu_lstock_flag":metadata['hwu_lstock_flag'],
    'dt':metadata['dt'],#seconds
    'dx':metadata['dx'],#meters
    'nsoil':metadata['nsoil'],
    'ncores':metadata['parallel_ncores'],
    'idate':idate,
    'fdate':fdate,
    'model_type':metadata['model_type'],
    'create_mask_flag':metadata['create_mask_flag'],
    'mkl_flag':metadata['mkl_flag']
    }
'''



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



