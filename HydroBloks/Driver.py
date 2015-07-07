#import sys
#sys.path.append('Preprocessing')
#sys.path.append('Model')
import datetime
import HydroBloks as HB
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
hydrobloks_info = {
        'icatch':metadata['catchment_id'],
	'input_file':metadata['input_file'],
	'output_file':metadata['output_file'],
        'soil_file':metadata['soil_file'],
        'workspace':metadata['workspace'],
	'surface_flow_flag':metadata['surface_flow_flag'],
	'subsurface_flow_flag':metadata['subsurface_flow_flag'],
	'dt':metadata['dt'],#seconds
	'dtt':metadata['dtt'],#seconds
	'dx':metadata['dx'],#meters
	'nsoil':metadata['nsoil'],
	'ncores':metadata['parallel_ncores'],
	'idate':idate,
	'fdate':fdate,
	'nclusters_nc':metadata['nhru_nc'],
	'nclusters_c':metadata['nhru_c'],
	'nclusters':metadata['nhru_nc'] + metadata['nhru_c'],
	'model_type':metadata['model_type'],
	}

#Run the model
HB.Run_Model(hydrobloks_info)

#Setup output redirection
#fout,ferr = os.open('/dev/null',os.O_RDWR|os.O_CREAT),os.open('/dev/null',os.O_RDWR|os.O_CREAT)
#so,se = os.dup(1),os.dup(2)

#Flush out the output
#sys.stdout.flush()
#sys.stderr.flush()

#Redirect the output back to the terminal 
#os.dup2(so, 1),os.dup2(se,2)
#os.dup2(so,1)
