import datetime
import HydroBlocks# as HB
import sys

def Read_Metadata_File(file):

 import json
 metadata = json.load(open(file))

 return metadata

#Read in the metadata file
metadata_file = sys.argv[1]
metadata = Read_Metadata_File(metadata_file)
ncores = metadata['parallel_ncores']

#Define the dates
idate = datetime.datetime(metadata['startdate']['year'],
			   metadata['startdate']['month'],
			   metadata['startdate']['day'],0)
fdate = datetime.datetime(metadata['enddate']['year'],
			   metadata['enddate']['month'],
			   metadata['enddate']['day'],23)

#Define the info
info = {
        'icatch':metadata['catchment_id'],
	'input_file':metadata['input_file'],
	'output_file':metadata['output_file'],
        #'soil_file':metadata['soil_file'],
        'workspace':metadata['workspace'],
	'surface_flow_flag':metadata['surface_flow_flag'],
	'subsurface_flow_flag':metadata['subsurface_flow_flag'],
        'hwu_flag':metadata['hwu_flag'],
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

#Initialize
HB = HydroBlocks.initialize(info)

#Run the model
HB.run(info)
 
#Finalize
HB.finalize()
