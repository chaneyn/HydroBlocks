import datetime
import Preprocessing
import sys

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



'''
#Define the info
hydrobloks_info = {
    'icatch':metadata['catchment_id'],
	'input_file':metadata['input_file'],
	'output_file':metadata['output_file'],
    'workspace':metadata['workspace'],
	'surface_flow_flag':metadata['surface_flow_flag'],
	'subsurface_module':metadata['subsurface_module'],
	'dt':metadata['dt'],#seconds
	'dx':metadata['dx'],#meters
	'nsoil':metadata['nsoil'],
	'ncores':metadata['parallel_ncores'],
	'idate':idate,
	'fdate':fdate,
	'nclusters':metadata['nhru'],
	'model_type':metadata['model_type'],
    'create_mask_flag':metadata['create_mask_flag'],
    'covariates':metadata['covariates'],
    'hwu_flag':metadata['hwu_flag'],
    'hwu_sf_flag':metadata['hwu_sf_flag'],
    'hwu_gw_flag':metadata['hwu_gw_flag'],
    'hwu_agric_flag':metadata['hwu_agric_flag'],
    'hwu_domest_flag':metadata['hwu_domest_flag'],
    'hwu_indust_flag':metadata['hwu_indust_flag'],
    'hwu_lstock_flag':metadata['hwu_lstock_flag'],
    }
'''


#Cluster the data
Preprocessing.Prepare_Model_Input_Data(info)

