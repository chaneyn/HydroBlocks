import warnings
import collections
warnings.filterwarnings('ignore')
#import gdal
import os
import pickle
import numpy as np
#import matplotlib.pyplot as plt
import sys
import glob
sys.stdout.flush()
import json

import geospatialtools.pedotransfer as pedotransfer
import geospatialtools.gdal_tools as gdal_tools
import geospatialtools.terrain_tools as terrain_tools
#import geospatialtools.netcdf_tools as nc_io
import time
import datetime
from random import shuffle
import netCDF4 as nc
from scipy.io import netcdf as scipy_nc
import time
from dateutil.relativedelta import relativedelta
#from rpy2.robjects import r,FloatVector
from osgeo import ogr,osr
import gc
#import sparse
from scipy.sparse import csr_matrix, csc_matrix, find, hstack
import psutil
import rasterio
import fiona
mb = 1024*1024

def driver(comm,metadata_file):

 size = comm.Get_size()
 rank = comm.Get_rank()
 #Read in the metadata
 metadata = json.load(open(metadata_file))
 #Define files
 metadata['output_data'] = '%s/data' % metadata['dir']
 metadata['dem'] = '%s/%s' % (metadata['elevation']['dir'],metadata['elevation']['vars']['dem'])
 metadata['acc'] = '%s/%s' % (metadata['elevation']['dir'],metadata['elevation']['vars']['acc'])
 metadata['fdir'] = '%s/%s' % (metadata['elevation']['dir'],metadata['elevation']['vars']['fdir'])
 metadata['svp'] = True
 #Create the domain decomposition
 if rank == 0:
  domain_decomposition(metadata)
 comm.Barrier()
 #Correct and finalize the domain decomposition
 correct_domain_decomposition(comm,metadata)
 comm.Barrier()

 #Read in the catchment summary database
 pck_file = '%s/cids/domain_database.pck' % metadata['output_data']
 cdb = pickle.load(open(pck_file,'rb'))
 crange = range(len(cdb))
 
 for ic in crange[rank::size]:

  print("Rank:%d, Catchment:%s - Initializing" % (rank,ic),flush=True)

  cid = cdb[ic]['cid']
  #cid = 1

  #Define the catchment directory
  cdir = '%s/cids/%d' % (metadata['output_data'],cid)

  #Prepare the workspace for the sub-domain
  os.system('rm -rf %s' % cdir)
  os.system('mkdir -p %s' % cdir)
  prepare_input_data(cdir,cdb[ic],metadata,rank,ic)

 return

def Create_Other_Soil_Properties_Soilgrids_Texture(cdb,workspace,metadata,icatch,log,properties):

 #Create soil texture maps
 S = properties['sand'] #%w
 C = properties['clay'] #%w
 ST = properties['silt'] #%w
 OM = properties['om'] #%w
 badvals = (S == -9999.0) | (C == -9999.0) | (ST == -9999.0) | (OM == -9999.0) | np.isnan(S)
 texture = np.ones(C.shape)*(-9999.0)
 texture[~badvals] = Texture_Class(C[~badvals],ST[~badvals],S[~badvals],OM[~badvals])

 #Look-up table
 file = '/home/nc153/soteria/projects/agu2019/data/noah_soils_lookup.pck'
 lutable = pickle.load(open(file,'rb')) 
 
 #Define unique classes
 utextures = np.unique(texture)
 utextures = utextures[utextures != -9999]
 output  = {'texture_class':texture[:]}

 #Theta saturated
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures: 
  tmp[texture == txt] = lutable['maxsmc'][int(txt)-1]
 output['thetas'] = tmp[:]

 #Theta residual
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['drysmc'][int(txt)-1]
 output['thetar'] = tmp[:]

 #Ksat
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['satdk'][int(txt)-1]
 output['ksat'] = 3600.0*1000.0*tmp[:]

 #Theta33  (1kPa ~ 10.197 cm H20)
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['refsmc'][int(txt)-1]
 output['theta33'] = tmp[:]

 #Theta1500
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['wltsmc'][int(txt)-1]
 output['theta1500'] = tmp[:]

 #Bb
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['bb'][int(txt)-1]
 output['bb'] = tmp[:]
   
 #Bubble pressure
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['satpsi'][int(txt)-1]
 output['psisat'] = tmp[:]
 
 #SATDW
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['satdw'][int(txt)-1]
 output['dsat'] = 3600.0*1000.0*tmp[:]

 #F11
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['f11'][int(txt)-1]
 output['f11'] = tmp[:]

 #QTZ (qtz)
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['qtz'][int(txt)-1]
 output['qtz'] = tmp[:]
  
 mask = gdal_tools.read_raster('%s/mask_latlon.tif' % workspace)

 #Output data
 md = gdal_tools.retrieve_metadata('%s/sand_latlon.tif' % workspace)
 md['nodata'] = -9999.0
 for var in output:
  file = '%s/%s_latlon.tif' % (workspace,var)
  gdal_tools.write_raster(file,md,output[var])

 del output
 gc.collect()
 return

def Create_Other_Soil_Properties_Polaris_Texture(cdb,workspace,metadata,icatch,log,properties):

 #Create soil texture maps
 S = properties['sand'] #%w
 C = properties['clay'] #%w
 ST = properties['silt'] #%w
 OM = properties['om'] #%w
 badvals = (S == -9999.0) | (C == -9999.0) | (ST == -9999.0) | (OM == -9999.0) | np.isnan(S)
 texture = np.ones(C.shape)*(-9999.0)
 texture[~badvals] = Texture_Class(C[~badvals],ST[~badvals],S[~badvals],OM[~badvals])

 #Look-up table
 file = '/home/nc153/soteria/projects/agu2019/data/noah_soils_lookup.pck'
 lutable = pickle.load(open(file,'rb')) 
 
 #Define unique classes
 utextures = np.unique(texture)
 utextures = utextures[utextures != -9999]
 output  = {'texture_class':texture[:]}

 #Theta saturated
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures: 
  tmp[texture == txt] = lutable['maxsmc'][int(txt)-1]
 output['thetas'] = tmp[:]

 #Theta residual
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['drysmc'][int(txt)-1]
 output['thetar'] = tmp[:]

 #Ksat
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['satdk'][int(txt)-1]
 output['ksat'] = 3600.0*1000.0*tmp[:]

 #Theta33  (1kPa ~ 10.197 cm H20)
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['refsmc'][int(txt)-1]
 output['theta33'] = tmp[:]

 #Theta1500
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['wltsmc'][int(txt)-1]
 output['theta1500'] = tmp[:]

 #Bb
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['bb'][int(txt)-1]
 output['bb'] = tmp[:]
   
 #Bubble pressure
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['satpsi'][int(txt)-1]
 output['psisat'] = tmp[:]
 
 #SATDW
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['satdw'][int(txt)-1]
 output['dsat'] = 3600.0*1000.0*tmp[:]

 #F11
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['f11'][int(txt)-1]
 output['f11'] = tmp[:]

 #QTZ (qtz)
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['qtz'][int(txt)-1]
 output['qtz'] = tmp[:]
  
 mask = gdal_tools.read_raster('%s/mask_latlon.tif' % workspace)

 #Output data
 md = gdal_tools.retrieve_metadata('%s/sand_latlon.tif' % workspace)
 md['nodata'] = -9999.0
 for var in output:
  file = '%s/%s_latlon.tif' % (workspace,var)
  gdal_tools.write_raster(file,md,output[var])

 del output
 gc.collect()
 return

def Create_Other_Soil_Properties_Conus_Soil(cdb,workspace,metadata,icatch,log,properties):

 #Look-up table
 file = '/home/nc153/soteria/projects/agu2019/data/noah_soils_lookup.pck'
 lutable = pickle.load(open(file,'rb')) 
 
 #Set the input properties
 texture = properties['texture_class'] #%w
 badvals = (texture == -9999.0) | np.isnan(texture)
 
 #Define unique classes
 utextures = np.unique(texture)
 output  = {}

 #Theta saturated
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures: 
  tmp[texture == txt] = lutable['maxsmc'][int(txt)-1]
 output['thetas'] = tmp[:]

 #Theta residual
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['drysmc'][int(txt)-1]
 output['thetar'] = tmp[:]

 #Ksat
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['satdk'][int(txt)-1]
 output['ksat'] = 3600.0*1000.0*tmp[:]

 #Theta33  (1kPa ~ 10.197 cm H20)
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['refsmc'][int(txt)-1]
 output['theta33'] = tmp[:]

 #Theta1500
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['wltsmc'][int(txt)-1]
 output['theta1500'] = tmp[:]

 #Bb
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['bb'][int(txt)-1]
 output['bb'] = tmp[:]
   
 #Bubble pressure
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['satpsi'][int(txt)-1]
 output['psisat'] = tmp[:]
 
 #SATDW
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['satdw'][int(txt)-1]
 output['dsat'] = 3600.0*1000.0*tmp[:]

 #F11
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['f11'][int(txt)-1]
 output['f11'] = tmp[:]

 #QTZ (qtz)
 tmp = -9999*np.ones(texture.shape)
 for txt in utextures:
  tmp[texture == txt] = lutable['qtz'][int(txt)-1]
 output['qtz'] = tmp[:]
  
 mask = gdal_tools.read_raster('%s/mask_latlon.tif' % workspace)

 #Output data
 md = gdal_tools.retrieve_metadata('%s/sand_latlon.tif' % workspace)
 md['nodata'] = -9999.0
 for var in output:
  file = '%s/%s_latlon.tif' % (workspace,var)
  gdal_tools.write_raster(file,md,output[var])

 del output
 gc.collect()
 return

def Create_Other_Soil_Properties(cdb,workspace,metadata,icatch,log,properties):
 
 #Set the input properties
 S = properties['sand'] #%w
 C = properties['clay'] #%w
 ST = properties['silt'] #%w
 OM = properties['om'] #%w  
 badvals = (S == -9999.0) | (C == -9999.0) | (ST == -9999.0) | (OM == -9999.0) | np.isnan(S)
 

 # Saxton pedotransfers limited to OM < 8%w and clay <60%
 #OM[OM > 8.0] = 8.0
 #C[C > 60.0] = 60.0
 
 output  = {}

 if metadata['soil_database'] == 'soilgrids':
   #Theta saturated
   output['thetas'] = pedotransfer.ThetaS_Saxton2006(S/100,C/100,OM/100)

   #Theta residual
   output['thetar'] = pedotransfer.Residual_Water_Content_Maidment92(output['thetas'],C,S)

   #Ksat
   output['ksat'] = pedotransfer.Ksat_Saxton2006(S/100,C/100,OM/100) # mm/h

   #Bb -- Brooks and Correy parameter lamba to the distribution of pore sizes  b = 1/lambda
   output['bb'] = 1.0/pedotransfer.Lambda_Maidment92(output['thetas'],C,S)

   #Bubble pressure
   output['psisat'] = pedotransfer.Bubbling_Pressure_Maidment92(output['thetas'],C,S)/100.0 # cm -> m

   #Theta33
   output['theta33'] = pedotransfer.Theta_33_Saxton2006(S/100,C/100,OM/100)

   #Theta1500
   output['theta1500'] = pedotransfer.Theta_1500_Saxton2006(S/100,C/100,OM/100)


 if metadata['soil_database'] == 'polaris':
   #Theta saturated
   output['thetas'] = properties['thetas']

   #Theta residual
   output['thetar'] = properties['thetar']

   #Ksat
   output['ksat'] = properties['ksat'] * 10  # cm/h --> mm/h

   #Theta33  (1kPa ~ 10.197 cm H20)
   output['theta33'] = np.power(properties['hb']/33,properties['lambda'])*(output['thetas']-output['thetar'])+output['thetar']

   #Theta1500
   output['theta1500'] = np.power(properties['hb']/1500,properties['lambda'])*(output['thetas']-output['thetar'])+output['thetar']

   #Convert BC parameters to Campbell
   lamda_campbell = (np.log(output['theta33'])-np.log(output['theta1500']))/(np.log(15000)-np.log(330))
   psisat_campbell = (output['theta1500']/output['thetas'])**(1/lamda_campbell)*15000

   #Bb -- Brooks and Corey parameter lamba to the distribution of pore sizes  b = 1/lambda
   #output['bb'] = 1.0/properties['lambda']
   output['bb'] = 1.0/lamda_campbell
   
   #Bubble pressure
   #output['psisat'] = properties['hb']/100.0 # cm --> m
   output['psisat'] = psisat_campbell/100.0 # cm --> m
 
 #SATDW (satdw = bb*satdk*satpsi/maxsmc)
 output['dsat'] = output['bb']*output['ksat']*output['psisat']/output['thetas'] # mm2/h

 #F11 (f11 = log10(satpsi) + bb*log10(maxsmc) + 2.0)
 #output['f11'] = np.log10(output['psisat']) + output['bb']*np.log10(output['thetas']) + 2.0
 output['f11'] = np.copy(output['psisat'])
 output['f11'][:] = -9999.0

 #QTZ (qtz)
 #output['qtz'] = np.copy(output['f11'])
 #output['qtz'][:] = 0.4


 S = properties['sand'] #%w
 C = properties['clay'] #%w
 ST = properties['silt'] #%w
 OM = properties['om'] #%w
 

 #Soil Texture
 #print C, ST, S, OM
 output['texture_class'] = np.ones(C.shape)*(-9999.0)
 output['texture_class'][~badvals] = Texture_Class(C[~badvals],ST[~badvals],S[~badvals],OM[~badvals])

 #QTZ (qtz)
 #From NOAH-LSM look-up table
 output['qtz'] = np.copy(output['psisat'])
 output['qtz'][output['texture_class'] == 1] = 0.92
 output['qtz'][output['texture_class'] == 2] = 0.82
 output['qtz'][output['texture_class'] == 3] = 0.60
 output['qtz'][output['texture_class'] == 4] = 0.25
 output['qtz'][output['texture_class'] == 5] = 0.10
 output['qtz'][output['texture_class'] == 6] = 0.40
 output['qtz'][output['texture_class'] == 7] = 0.60
 output['qtz'][output['texture_class'] == 8] = 0.10
 output['qtz'][output['texture_class'] == 9] = 0.35
 output['qtz'][output['texture_class'] == 10] = 0.52
 output['qtz'][output['texture_class'] == 11] = 0.10
 output['qtz'][output['texture_class'] == 12] = 0.25
  
 mask = gdal_tools.read_raster('%s/mask_latlon.tif' % workspace)

 badvals2 = (output['texture_class'] < 0) | (np.isnan(output['f11'])) | badvals
 for var in output:
  output[var][badvals2]=-9999.0
  m2 = ( mask >= 0 ) & np.invert(badvals2)
  missing_ratio = 1.0 - np.sum(m2)/float(np.sum(mask >= 0))
  if missing_ratio > 0.95:
   import sys
   sys.stderr.write('Error_preprocessing: %s_full_of_nans %s\n' % (var,icatch))
   return
   
 #Output data
 if metadata['svp']==False:
  md = gdal_tools.retrieve_metadata('%s/sand/sand_latlon.tif' % workspace)#laura, svp
  md['nodata'] = -9999.0
  for var in output:
   file = '%s/%s_latlon.tif' % (workspace,var)
   gdal_tools.write_raster(file,md,output[var])
 else:
  md = gdal_tools.retrieve_metadata('%s/sand/sand_latlon_%scm.tif' % (workspace,layer))
  md['nodata'] = -9999.0
  for var in output:
   if var=='f11':
    file = '%s/%s_latlon.tif' % (workspace,var)
    gdal_tools.write_raster(file,md,output[var])
   else:
    if var in ['bb','dsat','qtz','theta1500','theta33','texture_class','psisat','ksat','thetar','thetas']:
     if os.path.isdir('%s/%s' %(workspace,var))==False:
      os.system('mkdir %s/%s' %(workspace,var))
     file_out='%s/%s/%s_latlon_%scm.tif' % (workspace,var,var,layer)
     gdal_tools.write_raster(file_out,md,output[var][layer])

 del output
 gc.collect()
 return

def Create_Other_Soil_Properties_svp(cdb,workspace,metadata,icatch,log,properties):
 
 for layer in properties['sand']:
  #Set the input properties
  S = properties['sand'][layer] #%w
  C = properties['clay'][layer] #%w
  ST = properties['silt'][layer] #%w
  OM = properties['om'][layer] #%w 
  from numpy import inf 
  badvals = (S == -9999.0) | (C == -9999.0) | (ST == -9999.0) | (OM == -9999.0) | np.isnan(S) 

 # Saxton pedotransfers limited to OM < 8%w and clay <60%
 #OM[OM > 8.0] = 8.0
 #C[C > 60.0] = 60.0
 
  output  = {}

  if metadata['soil']['dataset'] == 'soilgrids':
   #Theta saturated
   output['thetas'] = pedotransfer.ThetaS_Saxton2006(S/100,C/100,OM/100)

   #Theta residual
   output['thetar'] = pedotransfer.Residual_Water_Content_Maidment92(output['thetas'],C,S)

   #Ksat
   output['ksat'] = pedotransfer.Ksat_Saxton2006(S/100,C/100,OM/100) # mm/h

   #Bb -- Brooks and Correy parameter lamba to the distribution of pore sizes  b = 1/lambda
   output['bb'] = 1.0/pedotransfer.Lambda_Maidment92(output['thetas'],C,S)

   #Bubble pressure
   output['psisat'] = pedotransfer.Bubbling_Pressure_Maidment92(output['thetas'],C,S)/100.0 # cm -> m

   #Theta33
   output['theta33'] = pedotransfer.Theta_33_Saxton2006(S/100,C/100,OM/100)

   #Theta1500
   output['theta1500'] = pedotransfer.Theta_1500_Saxton2006(S/100,C/100,OM/100)

  if metadata['soil']['dataset'] == 'polaris':
    #Theta saturated
    output['thetas']={}
    output['thetas'][layer] = properties['thetas'][layer]

    #Theta residual
    output['thetar']={}
    output['thetar'][layer] = properties['thetar'][layer]

    #Ksat
    output['ksat']={}
    output['ksat'][layer] = properties['ksat'][layer] * 10  # cm/h --> mm/h

    #Theta33  (1kPa ~ 10.197 cm H20)
    output['theta33']={}
    output['theta33'][layer] = np.float32(np.power(properties['hb'][layer]/33,properties['lambda'][layer])*(output['thetas'][layer]-output['thetar'][layer])+output['thetar'][layer])

    #Theta1500
    output['theta1500']={}
    output['theta1500'][layer] = np.power(properties['hb'][layer]/1500,properties['lambda'][layer])*(output['thetas'][layer]-output['thetar'][layer])+output['thetar'][layer]

    #Convert BC parameters to Campbell
    lamda_campbell = (np.log(output['theta33'][layer])-np.log(output['theta1500'][layer]))/(np.log(15000)-np.log(330))
    psisat_campbell = (output['theta1500'][layer]/output['thetas'][layer])**(1/lamda_campbell)*15000

    #Bb -- Brooks and Corey parameter lamba to the distribution of pore sizes  b = 1/lambda
    #output['bb'] = 1.0/properties['lambda']
    output['bb']={}
    output['bb'][layer] = 1.0/lamda_campbell
   
    #Bubble pressure
    #output['psisat'] = properties['hb']/100.0 # cm --> m
    output['psisat']={}
    output['psisat'][layer] = psisat_campbell/100.0 # cm --> m
 
  #SATDW (satdw = bb*satdk*satpsi/maxsmc)
  output['dsat']={}
  output['dsat'][layer] = output['bb'][layer]*output['ksat'][layer]*output['psisat'][layer]/output['thetas'][layer] # mm2/h

 #F11 (f11 = log10(satpsi) + bb*log10(maxsmc) + 2.0)
 #output['f11'] = np.log10(output['psisat']) + output['bb']*np.log10(output['thetas']) + 2.0
  output['f11'] = np.copy(output['psisat'][layer])
  output['f11'][:] = -9999.0

 #QTZ (qtz)
 #output['qtz'] = np.copy(output['f11'])
 #output['qtz'][:] = 0.4

  S = properties['sand'][layer] #%w
  C = properties['clay'][layer] #%w
  ST = properties['silt'][layer] #%w
  OM = properties['om'][layer] #%w
 
  #Soil Texture
  #print C, ST, S, OM
  output['texture_class']={}
  output['texture_class'][layer] = np.ones(C.shape)*(-9999.0)
  output['texture_class'][layer][~badvals] = Texture_Class(C[~badvals],ST[~badvals],S[~badvals],OM[~badvals])

  #QTZ (qtz)
  #From NOAH-LSM look-up table
  output['qtz']={}
  output['qtz'][layer] = np.copy(output['psisat'][layer])
  output['qtz'][layer][output['texture_class'] == 1] = 0.92
  output['qtz'][layer][output['texture_class'] == 2] = 0.82
  output['qtz'][layer][output['texture_class'] == 3] = 0.60
  output['qtz'][layer][output['texture_class'] == 4] = 0.25
  output['qtz'][layer][output['texture_class'] == 5] = 0.10
  output['qtz'][layer][output['texture_class'] == 6] = 0.40
  output['qtz'][layer][output['texture_class'] == 7] = 0.60
  output['qtz'][layer][output['texture_class'] == 8] = 0.10
  output['qtz'][layer][output['texture_class'] == 9] = 0.35
  output['qtz'][layer][output['texture_class'] == 10] = 0.52
  output['qtz'][layer][output['texture_class'] == 11] = 0.10
  output['qtz'][layer][output['texture_class'] == 12] = 0.25
  
  mask = gdal_tools.read_raster('%s/mask_latlon.tif' % workspace)

  badvals2 = (output['texture_class'][layer] < 0) | (np.isnan(output['f11'])) | badvals
  for var in output:
   if var=='f11':
    output[var][badvals2]=-9999.0
    m2 = ( mask >= 0 ) & np.invert(badvals2)
    missing_ratio = 1.0 - np.sum(m2)/float(np.sum(mask >= 0))
    if missing_ratio > 0.95:
     import sys
     sys.stderr.write('Error_preprocessing: %s_full_of_nans %s\n' % (var,icatch))
     return
   else:
    output[var][layer][badvals2]=-9999.0
    m2= (mask >= 0 ) & np.invert(badvals2)
    missing_ratio = 1.0 - np.sum(m2)/float(np.sum(mask >= 0))
    if missing_ratio > 0.95:
     import sys
     sys.stderr.write('Error_preprocessing: %s_full_of_nans %s\n' % (var,icatch))
     return
   
 #Output data
  md = gdal_tools.retrieve_metadata('%s/sand/sand_latlon_%scm.tif' % (workspace,layer))
  md['nodata'] = -9999.0
  for var in output:
   if var=='f11':
    file = '%s/%s_latlon.tif' % (workspace,var)
    gdal_tools.write_raster(file,md,output[var])
   else:
    if var in ['bb','dsat','qtz','theta1500','theta33','texture_class','psisat','ksat','thetar','thetas']:
     if os.path.isdir('%s/%s' %(workspace,var))==False:
      os.system('mkdir %s/%s' %(workspace,var))
     file_out='%s/%s/%s_latlon_%scm.tif' % (workspace,var,var,layer)
     gdal_tools.write_raster(file_out,md,output[var][layer])

 del output
 gc.collect()
 return


def Texture_Class(clay,silt,sand,om):
 
 #r('library("soiltexture")')

 #Ensure sand/silt/clay add up to 100
 m = (sand != -9999) | (clay != -9999) | (silt != -9999) | (om != -9999)
 sum = sand[m]+silt[m]+clay[m]
 sand[m] = 100*sand[m]/sum
 silt[m] = 100*silt[m]/sum
 clay[m] = 100*clay[m]/sum
 
 oc = 0.58*om*(1000/100.) # OC=0.58*OM [g/kg]
 shape = sand.shape

 '''#Pass to R
 r.assign('clay',FloatVector(clay[m]))
 r.assign('sand',FloatVector(sand[m]))
 r.assign('silt',FloatVector(silt[m]))
 #r.assign('om',FloatVector(om))
 r.assign('oc',FloatVector(oc[m]))
 del sum, sand, silt, clay, om ; gc.collect()
 #Compute the texture class
 # clay, sand, silt in %, OC=0.58*OM [g/kg]
 r('my.text <- data.frame("CLAY" = clay,"SILT" = silt,"SAND" = sand,"OC" = oc)')
 r('rm(clay,sand,silt,om); gc()')
 # Classify according to the USDA classification
 r('output <- TT.points.in.classes(tri.data = my.text,class.sys = "USDA.TT")')
 r('rm(my.text); gc()')
 #output = csr_matrix(np.array(r("output")))
 output = np.array(r("output"))
 r('rm(list=ls()); gc()')
 #print(output)

 #novals = csr_matrix(np.zeros((output.shape[0],1)))
 novals = np.zeros((output.shape[0],1))
 #output = hstack([output,novals]).tocsr()
 output = np.hstack([output,novals])

 tmp = output.sum(axis=1)
 mask1 = find(tmp == 0)[0]
 output[mask1,12] = 1

 mask2 = find(output>1)
 maskl = mask2[0]
 maskc = mask2[1]
 maskl,index = np.unique(maskl,return_index=True)
 output[maskl,:]=0
 #output[maskl,maskc[index]]=1
 output[maskl,12]=1

  
 #Convert to what NOAH wants
 #Cl SiCl SaCl ClLo SiClLo SaClLo Lo SiLo SaLo Si LoSa Sa -9999.0
 mapping = np.array([12,11,10,9,8,7,6,4,3,5,2,1,-9999.0])
 mask = find(output > 0)[1]
 soil_texture  = mapping[mask]
 texture = np.ones(shape)
 texture[:] = -9999.0
 texture[m] = soil_texture'''
 texture = pedotransfer.compute_soil_texture_class(sand,clay)
 #Set -9999 to Loam
 texture[texture == -9999] = 6

 #del output, mask, soil_texture; gc.collect()
 return texture

def domain_decomposition(md):
 
 print("Assembling domain shapefile")
 #Make new lat/lon domain
 if md['domain_decomposition']['type'] == 'latlon':
  #Use to create shapefile of domain
  create_domain_shapefile(md)
 if md['domain_decomposition']['type'] == 'gfdl_gridspec': 
  #Create shapefile from gfdl grid spec file
  create_domain_shapefile_gfdl_gridspec(md)
 #Create summary database
 summarize_domain_decompisition(md)

 return

def create_domain_shapefile_gfdl_gridspec(metadata):
 
 #Extract metadata info
 gridspec_dir = metadata['domain_decomposition']['gridspec_dir']
 tgs = '%s/%s' % (gridspec_dir,metadata['domain_decomposition']['gridspec_name'])
 tlm = '%s/%s' % (gridspec_dir,metadata['domain_decomposition']['landmask_name'])
 ntiles = metadata['domain_decomposition']['ntiles']
 sdir = '%s/data/shp' % metadata['dir']
 os.system('mkdir -p %s' % sdir)
 #Create the shapefile
 driver = ogr.GetDriverByName("ESRI Shapefile")
 ds = driver.CreateDataSource("%s/domain.shp" % sdir)
 srs = osr.SpatialReference()
 srs.ImportFromEPSG(4326)
 layer = ds.CreateLayer("grid", srs, ogr.wkbPolygon)
 layer.CreateField(ogr.FieldDefn("ID", ogr.OFTInteger))
 layer.CreateField(ogr.FieldDefn("TILE", ogr.OFTInteger))
 layer.CreateField(ogr.FieldDefn("X", ogr.OFTInteger))
 layer.CreateField(ogr.FieldDefn("Y", ogr.OFTInteger))
 layer.CreateField(ogr.FieldDefn("LFN", ogr.OFTReal))
 layer.CreateField(ogr.FieldDefn("LFO", ogr.OFTReal))

 count = 0
 for tile in range(1,ntiles+1):
  #Extract xs,ys,and lms
  gs_tile = tgs.replace('$tid',str(tile))
  lm_tile = tlm.replace('$tid',str(tile))
  fp  = nc.Dataset(gs_tile)
  fplm = nc.Dataset(lm_tile)
  #Extract lats/lons
  xs = fp.variables['x'][:]
  ys = fp.variables['y'][:]
  xs[xs > 180.0] = xs[xs > 180.0] - 360.0
  #Extract mask
  lm = fplm.variables['mask'][:]
  fp.close()
  fplm.close()

 #Iterate through each cell
 ncells_x = int((xs.shape[1]-1)/2)
 ncells_y = int((ys.shape[0]-1)/2)
 for icell in range(ncells_x):
  imin = icell*2
  imax = (icell+1)*2
  for jcell in range(ncells_y):
   jmin = jcell*2
   jmax = (jcell+1)*2
   #Skip 87-90N and 87-90S
   minlat = np.min(ys[jmin:jmax+1,imin:imax+1])
   maxlat = np.max(ys[jmin:jmax+1,imin:imax+1])
   minlon = np.min(xs[jmin:jmax+1,imin:imax+1])
   maxlon = np.max(xs[jmin:jmax+1,imin:imax+1])
   if ((minlat >= 86.0) | (maxlat <= -86.0)):continue #Do not deal with poles
   if ((minlon < -90) & (maxlon > 90)):continue #Do not go over dateline
   #Update the count
   count += 1
   #if (count-1) % size != rank:continue
   print(tile,minlat,maxlat,minlon,maxlon)
   #Compute the info
   compute_info(xs,ys,metadata,count,sdir,imin,imax,jmin,jmax,lm,icell,jcell,layer,tile)#,log)

 # Destroy the data source to free resources
 ds.Destroy()

 #Memorize domain decomposition file
 metadata['domain_decomposition']['file'] = "%s/domain.shp" % sdir

 return metadata

def compute_info(xs,ys,metadata,count,cdir,imin,imax,jmin,jmax,lm,icell,jcell,layer,tile):#,log):

 output = {}
 #Extract regional boundaries
 xs_subset = xs[jmin:jmax+1,imin:imax+1]
 #if (np.max(xs_subset) > 100) and (np.min(xs_subset) < -100):
 # xs_subset[xs_subset < 0] = xs_subset[xs_subset < 0] + 360.0

 #Construct list of points
 mpoint = ogr.Geometry(ogr.wkbMultiPoint)
 point = ogr.Geometry(ogr.wkbPoint)
 for i in range(imin,imax+1):
  for j in range(jmin,jmax+1):
   point.AddPoint(xs[j,i],ys[j,i])
   #point.AddPoint(ys[j,i],xs[j,i])
   mpoint.AddGeometry(point)
 poly = mpoint.ConvexHull()

 #Calculate land fraction
 #Retrieve coordinates of envelope
 bbox = poly.GetEnvelope()
 metadata['bbox'] = bbox
 #Calculate the original land fraction
 lfrac_org = lm[jcell,icell]
 lfrac = 1.0 #HACK
 #Calculate the land fraction
 #if (lfrac_org == 0.0):
 # return
 #else:
 # lfrac = calculate_land_fraction(metadata,count,cdir,poly)#log,poly)

 #Add to the info
 # create the feature
 feature = ogr.Feature(layer.GetLayerDefn())
 #Set the calculate land fraction
 feature.SetField("LFN",lfrac)
 #Set the land fraction from the grid spec
 feature.SetField("LFO",lfrac_org)
 #Set the ID
 feature.SetField("ID",count)
 #Set the tile
 feature.SetField("TILE",tile)
 #Set the i id
 feature.SetField("X",icell+1)
 #Set the j id
 feature.SetField("Y",jcell+1)
 # Set the feature geometry using the point
 feature.SetGeometry(poly)
 # Create the feature in the layer (shapefile)
 layer.CreateFeature(feature)
 # Destroy the feature to free resources
 feature.Destroy()

 return

def create_domain_shapefile(md):

 #Extract parameters
 minlat = md['domain_decomposition']['boundaries']['minlat']
 maxlat = md['domain_decomposition']['boundaries']['maxlat']
 minlon = md['domain_decomposition']['boundaries']['minlon']
 maxlon = md['domain_decomposition']['boundaries']['maxlon']
 res = md['domain_decomposition']['res']
 sdir = '%s/data/shp' % md['dir']
 os.system('mkdir -p %s' % sdir)

 #Create the shapefile
 driver = ogr.GetDriverByName("ESRI Shapefile")
 ds = driver.CreateDataSource("%s/domain.shp" % sdir)
 srs = osr.SpatialReference()
 srs.ImportFromEPSG(4326)
 layer = ds.CreateLayer("grid", srs, ogr.wkbPolygon)
 layer.CreateField(ogr.FieldDefn("ID", ogr.OFTInteger))
 layer.CreateField(ogr.FieldDefn("X", ogr.OFTInteger))
 layer.CreateField(ogr.FieldDefn("Y", ogr.OFTInteger))

 #Define lats,lons
 #print(int(np.round((maxlon-minlon)/res)))
 #print(int(np.round((maxlat-minlat)/res)))
 nlat = int(np.round((maxlat-minlat)/res)) + 1
 nlon = int(np.round((maxlon-minlon)/res)) + 1
 lats = np.linspace(minlat,maxlat,nlat)
 lons = np.linspace(minlon,maxlon,nlon)
 #lats = np.arange(minlat,maxlat+res,res)
 #lons = np.arange(minlon,maxlon+res,res)

 #Iterate through each cell
 cid = 1
 for ilat in range(lats.size-1):
  for ilon in range(lons.size-1):
   #Construct list of points
   mpoint = ogr.Geometry(ogr.wkbMultiPoint)
   point = ogr.Geometry(ogr.wkbPoint)
   iss = [0,0,1,1,0]
   jss = [0,1,1,0,0]
   for k in range(len(iss)):
    i = iss[k]
    j = jss[k]
    point.AddPoint(lons[ilon+j],lats[ilat+i])
    mpoint.AddGeometry(point)
   poly = mpoint.ConvexHull()
   #Add to the info
   # create the feature
   feature = ogr.Feature(layer.GetLayerDefn())
   #Set the ID
   feature.SetField("ID",cid)
   #Set the i id
   feature.SetField("X",ilat)
   #Set the j id
   feature.SetField("Y",ilon)
   #Set the feature geometry using the point
   feature.SetGeometry(poly)
   #Create the feature in the layer (shapefile)
   layer.CreateFeature(feature)
   #Destroy the feature to free resources
   feature.Destroy()
   #Update cid
   cid += 1

 #Destroy the data source to free resources
 ds.Destroy()

 #Memorize domain decomposition file
 md['domain_decomposition']['file'] = "%s/domain.shp" % sdir

 return md

def summarize_domain_decompisition(md):

 #Get info for the basin from the WBD dataset
 file_in = md['domain_decomposition']['file']
 cid = str(md['domain_decomposition']['id'])
 file_out = '%s/shp' % md['output_data']
 #minlat = md['boundaries']['minlat']
 #minlon = md['boundaries']['minlon']
 #maxlat = md['boundaries']['maxlat']
 #maxlon = md['boundaries']['maxlon']

 #Prepare the domain directory
 cdir = '%s/cids' % md['output_data']
 #os.system('rm -rf %s' % cdir)
 os.system('mkdir -p %s' % cdir)

 #Extract the region of interest
 #os.system('rm -rf %s' % file_out)
 #os.system('ogr2ogr -spat %.16f %.16f %.16f %.16f %s %s' % (minlon,minlat,maxlon,maxlat,file_out,file_in))
 #os.system('cp %s %s' % (file_in,file_out))
 file_in = file_out

 #Open access to the database
 from osgeo import ogr
 driver = ogr.GetDriverByName("ESRI Shapefile")
 ds = driver.Open(file_out,0)

 #Iterate through each feature getting the necessary info
 layer = ds.GetLayer()
 output = []
 for feature in layer:
   info = {}
   bbox = feature.GetGeometryRef().GetEnvelope()
   info['cid'] = feature.GetField(cid)
   info['bbox'] =  {'minlat':bbox[2],'minlon':bbox[0],'maxlat':bbox[3],'maxlon':bbox[1]}
   output.append(info)

 #Close the shapefile file
 del ds, ogr, layer, bbox, info
 gc.collect()

 #Pickle the database
 pickle.dump(output,open('%s/cids/domain_database.pck' % md['output_data'],'wb'),pickle.HIGHEST_PROTOCOL)
 del output
 gc.collect()

 return

def Create_Mask(cdb,workspace,metadata,icatch,log):

 #Define parameters
 shp_in = '%s/shp' % metadata['output_data']
 shp_out = '%s/shp' % workspace
 cistr = metadata['domain_decomposition']['id']
 lstr = metadata['domain_decomposition']['layer']
 ci = cdb['cid']
 bbox =  cdb['bbox']
 res_latlon = metadata['res_latlon']
 
 #Define the files
 mask_latlon_file = '%s/mask_latlon.tif' % workspace
 tmp_file = '%s/tmp.tif' % workspace
 
 #Rasterize the area
 buff = 0.1
 
 print(' buffer size:',buff,' icatch:',ci,flush=True) 
 minx = bbox['minlon']-buff
 miny = bbox['minlat']-buff
 maxx = bbox['maxlon']+buff
 maxy = bbox['maxlat']+buff
 cache = int(psutil.virtual_memory().available*0.7/mb)
 print(minx,maxx,miny,maxy)

 #Correct coordinates to avoid reprojections
 #Fix coordinates to the dem vrt to avoid inconsistiences
 md = gdal_tools.retrieve_metadata(metadata['dem'])
 res_latlon = np.abs(md['resx'])
 minx = md['minx'] + np.floor((minx-md['minx'])/res_latlon)*res_latlon
 miny = md['miny'] + np.floor((miny-md['miny'])/res_latlon)*res_latlon#np.floor(miny/res_latlon)*res_latlon
 maxx = md['minx'] + np.ceil((maxx-md['minx'])/res_latlon)*res_latlon
 maxy = md['miny'] + np.ceil((maxy-md['miny'])/res_latlon)*res_latlon

 #Rasterize
 os.system("gdal_rasterize -at -ot Float64 --config GDAL_CACHEMAX %i -a_nodata -9999 -init -9999 -tr %.16f %.16f -te %.16f %.16f %.16f %.16f -a %s -l %s %s %s >> %s 2>&1" % (cache, res_latlon,res_latlon,minx,miny,maxx,maxy,cistr,lstr,shp_in,mask_latlon_file,log))

 #del data
 gc.collect()

 return

def Correct_Mask(cdb,workspace,metadata,icatch,log):

 bbox =  cdb['bbox']
 res_latlon = metadata['res_latlon']

 #0. Read in mask
 mask = gdal_tools.read_data('%s/mask_latlon.tif' % workspace).data
 mask_org = np.copy(mask)

 #0.25 Read admin boundaries
 admin = gdal_tools.read_data('%s/admin_latlon.tif' % workspace).data

 #0.5 Copy original mask
 os.system('cp %s %s' % ('%s/mask_latlon.tif' % workspace,'%s/mask_org_latlon.tif' % workspace))

 #1a. Read in elevation data
 dem = gdal_tools.read_data('%s/dem_latlon.tif' % workspace).data

 #1b. Read in arcgis flow direction data
 fdir_arcgis = gdal_tools.read_data('%s/fdir_latlon.tif' % workspace).data

 #1c. Read in accumulation area
 acc = 10**6*gdal_tools.read_data('%s/acc_latlon.tif' % workspace).data #km2->m2

 #2. Update the mask
 m2 = np.copy(mask).astype(np.bool)
 m2[:] = 0
 m2[dem != -9999] = 1
 m2[admin <= 0] = 0
 demns = dem

 #eares = 90 #meter(hack)
 #Calculate slope and aspect
 #res_array = np.copy(demns)
 #res_array[:] = eares
 #(slope,aspect) = terrain_tools.ttf.calculate_slope_and_aspect(np.flipud(demns),res_array,res_array)
 #slope = np.flipud(slope)
 #aspect = np.flipud(aspect)
 #Calculate accumulation area
 fdir = terrain_tools.transform_arcgis_fdir(fdir_arcgis)
 #area = terrain_tools.ttf.calculate_d8_acc_pfdir(demns,m2,eares,fdir)
 ##(area,fdir) = terrain_tools.ttf.calculate_d8_acc(demns,m2,eares)
 #Calculate channel initiation points (2 parameters)
 #C = area/eares*slope**2
 #ipoints = ((C > 100) & (area > 10**5)).astype(np.int32)
 cthrs = metadata['channel_initiation']["athrs"]#Laura #10**6
 #ipoints = ((area > cthrs)).astype(np.int32)
 #ipoints[ipoints == 0] = -9999
 #Create area for channel delineation
 #ac = terrain_tools.ttf.calculate_d8_acc_wipoints_pfdir(demns,m2,ipoints,eares,fdir)
 fdc = fdir
 ac = acc
 ac[m2 == 0] = -9999
 #ac[ac != 0] = area[ac != 0]
 
 #Compute the channels
 (channels,channels_wob,channel_topology,shreve_order) = terrain_tools.ttf.calculate_channels_wocean_wprop(ac,cthrs,cthrs,fdc,m2)
 basins = terrain_tools.ttf.delineate_basins(channels_wob,m2,fdc)
 basins[basins == 0] = -9999
 tmp = gdal_tools.read_data('%s/mask_latlon.tif' % workspace)
 tmp.data = basins.astype(np.float32)
 tmp.write_data('%s/basins_latlon.tif' % workspace)
 #Determine the basins that have the majority in the given area
 #0.Match up with dem
 #mask[(mask == -9999) & (basins != -9999)] = np.max(mask) + 1
 mask_alt = np.copy(mask)
 #1.Compute the area of each basin in each "coarse grid cell"
 ubasins = np.unique(basins)
 ubasins = ubasins[ubasins != -9999]
 ucatchs = np.unique(mask)
 ucatchs = ucatchs[ucatchs != -9999]
 i = 0
 for uc in ucatchs:
  mask_alt[mask == uc] = i
  i += 1 
 count = np.zeros((ucatchs.size,ubasins.size))
 channel_count = np.zeros((ucatchs.size,ubasins.size))
 basin_external = np.ones(ubasins.size)
 for i in range(basins.shape[0]):
  for j in range(basins.shape[1]):
   b = int(basins[i,j])
   r = int(mask_alt[i,j])
   if b != -9999:
    if(mask_alt[i,j] < 0):basin_external[b-1] = 0
    if r != -9999:
     if(channels_wob[i,j] > 0):channel_count[r,b-1] += 1
   if (b == -9999) | (r == -9999):continue
   count[r,b-1] += 1
 #Find the basins that have the largest representation in the given cell or catchment
 argmax = np.argmax(count,axis=0)
 #If none have a basin then set argmax to nan
 scount = np.sum(count,axis=0)
 argmax[scount == 0] = -9999
 mask_v2 = np.copy(mask)
 mask_v2[:] = -9999
 #Remove all basins that are not in the final list
 for i in range(basins.shape[0]):
  for j in range(basins.shape[1]):
   b = basins[i,j]
   if b == -9999 :continue
   if argmax[b-1] == -9999:continue
   mask_v2[i,j] = ucatchs[argmax[b-1]]
 #Set the external to -9999
 #mask_v2[mask_v2 == np.max(mask_v2)] = -9999
 #Output the data
 mask = gdal_tools.read_data('%s/mask_latlon.tif' % workspace)
 mask.data = mask_v2
 mask.write_data('%s/mask_latlon.tif' % workspace)

 gc.collect()
 
 return

def Terrain_Analysis(cdb,workspace,metadata,icatch,log):

 #0. Get the parameters
 md = gdal_tools.retrieve_metadata('%s/mask_latlon.tif' % workspace)
 minx = md['minx']
 miny = md['miny']
 maxx = md['maxx']
 maxy = md['maxy']
 res = abs(md['resx'])
 dem_region = metadata['dem']
 acc_region = metadata['acc']
 fdir_region = metadata['fdir']
 #acc_region = metadata['acc']
 lproj = md['proj4']

 #1. Cutout the region of interest
 dem_latlon_file = '%s/dem_latlon.tif' % workspace
 acc_latlon_file = '%s/acc_latlon.tif' % workspace
 fdir_latlon_file = '%s/fdir_latlon.tif' % workspace
 cache = int(psutil.virtual_memory().available*0.7/mb)
 #os.system('gdalwarp -t_srs \'%s\' -dstnodata -9999.0 -r bilinear -tr %.16f %.16f -te %.16f %.16f %.16f %.16f --config GDAL_CACHEMAX %i  %s %s >> %s 2>&1' % (lproj,res,res,minx,miny,maxx,maxy,cache,dem_region,dem_latlon_file,log))
 os.system('gdalwarp -t_srs \'%s\' -dstnodata -9999.0 -te %.16f %.16f %.16f %.16f --config GDAL_CACHEMAX %i  %s %s >> %s 2>&1' % (lproj,minx,miny,maxx,maxy,cache,dem_region,dem_latlon_file,log))
 os.system('gdalwarp -t_srs \'%s\' -dstnodata -9999.0 -te %.16f %.16f %.16f %.16f --config GDAL_CACHEMAX %i  %s %s >> %s 2>&1' % (lproj,minx,miny,maxx,maxy,cache,acc_region,acc_latlon_file,log))
 os.system('gdalwarp -t_srs \'%s\' -dstnodata -9999.0 -te %.16f %.16f %.16f %.16f --config GDAL_CACHEMAX %i  %s %s >> %s 2>&1' % (lproj,minx,miny,maxx,maxy,cache,fdir_region,fdir_latlon_file,log))

 #1. Cutout the region of interest
 #acc_latlon_file = '%s/acc_latlon.tif' % workspace
 #os.system('gdalwarp -t_srs \'%s\' -dstnodata -9999.0 -te %.16f %.16f %.16f %.16f --config GDAL_CACHEMAX %i  %s %s >> %s 2>&1' % (lproj,minx,miny,maxx,maxy,cache,acc_region,acc_latlon_file,log))

 data = gdal_tools.read_raster(dem_latlon_file)
 metadata = gdal_tools.retrieve_metadata(dem_latlon_file)
 metadata['nodata'] = -9999.0
 data[data == -32768] = -9999.0
 gdal_tools.write_raster(dem_latlon_file,metadata,data)

 data = gdal_tools.read_raster(acc_latlon_file)
 metadata = gdal_tools.retrieve_metadata(acc_latlon_file)
 metadata['nodata'] = -9999.0
 #data[data == -32768] = -9999.0
 gdal_tools.write_raster(acc_latlon_file,metadata,data)

 data = gdal_tools.read_raster(fdir_latlon_file)
 metadata = gdal_tools.retrieve_metadata(fdir_latlon_file)
 metadata['nodata'] = -9999.0
 #data[data == -32768] = -9999.0
 gdal_tools.write_raster(fdir_latlon_file,metadata,data)

 mask = gdal_tools.read_raster('%s/mask_latlon.tif' % workspace)
 m2 = ( mask >= 0 ) & ( data != -9999.0)
 missing_ratio = 1. - np.sum(m2)/float(np.sum(mask >= 0))
 if missing_ratio > 0.95 :
  import sys
  sys.stderr.write('Error_preprocessing: dem_full_of_nans %s\n' % (icatch))
  return

 del data
 gc.collect()
 
 return

def Extract_Dbedrock(cdb,workspace,metadata,icatch,log):
 #0. Get the parameters
 md = gdal_tools.retrieve_metadata('%s/mask_latlon.tif' % workspace)
 minx = md['minx']
 miny = md['miny']
 maxx = md['maxx']
 maxy = md['maxy']
 res = abs( md['resx'])
 lc_region = metadata['dbedrock']
 lproj = md['proj4']

 latlon_file = '%s/dbedrock_latlon.tif' % (workspace)
 cache = int(psutil.virtual_memory().available*0.7/mb)
 cmd = 'gdalwarp -overwrite -t_srs \'%s\' -te %.16f %.16f %.16f %.16f -tr %.16f %.16f --config GDAL_CACHEMAX %i  %s %s >> %s 2>&1' % (lproj,minx,miny,maxx,maxy,res,res,cache,lc_region,latlon_file,log)
 os.system(cmd)

 data = gdal_tools.read_raster(latlon_file)
 m1 = (data == -1) | (data == 255) | (data < 0)
 data[m1] = -9999.0
 m = ( data >= 0 ) & (data < 1.0) 
 data[m] = 1.0 # set the min depth
 md['nodata'] = -9999.0
 gdal_tools.write_raster(latlon_file,md,data)

 mask = gdal_tools.read_raster('%s/mask_latlon.tif' % workspace)
 m2 = ( mask >= 0 ) & ( data != -9999.0)
 missing_ratio = 1. - np.sum(m2)/float(np.sum(mask >= 0))
 if missing_ratio > 0.95 :
  import sys
  sys.stderr.write('Error_preprocessing: dbedrock_full_of_nans %s\n' % (icatch))
  return

 del data
 gc.collect()

 return

def Extract_Land_Cover(cdb,workspace,metadata,icatch,log):

 #0. Get the parameters
 md = gdal_tools.retrieve_metadata('%s/mask_latlon.tif' % workspace)
 minx = md['minx']
 miny = md['miny']
 maxx = md['maxx']
 maxy = md['maxy']
 res = abs(md['resx'])
 lc_region = metadata['landcover']['file']
 lproj = md['proj4']
 #mapping = eval(metadata['landcover']['mapping'])
 mapping = {0:-9999.0,11:17,12:15,21:10,22:13,23:13,24:13,31:16,41:4,42:2,43:5,45:5,46:5,51:7,52:6,71:10,72:19,73:19,74:19,81:12,82:12,90:11,95:11,-9999.0:-9999.0}

 #1. Prepare Land cover data
 lc_latlon_file = '%s/lc_latlon.tif' % workspace
 cache = int(psutil.virtual_memory().available*0.7/mb)
 os.system('gdalwarp -t_srs \'%s\' -dstnodata -9999 -r near -tr %.16f %.16f -te %.16f %.16f %.16f %.16f --config GDAL_CACHEMAX %i %s %s >> %s 2>&1' % (lproj,res,res,minx,miny,maxx,maxy,cache,lc_region,lc_latlon_file,log))
 data = gdal_tools.read_raster(lc_latlon_file)

 #0.5. Map to final 
 lcs = np.unique(data)[:]
 tmp = np.copy(data)
 for lc in lcs:
  tmp[data == lc] = mapping[lc]
 data = tmp

 #2. Export to tif
 md = gdal_tools.retrieve_metadata(lc_latlon_file)
 md['nodata'] = -9999.0
 gdal_tools.write_raster(lc_latlon_file,md,data)
 tmp = '%s/lc_ea2.tif' % workspace
 cache = int(psutil.virtual_memory().available*0.7/mb)
 os.system('gdal_translate --config GDAL_CACHEMAX %i %s %s >> %s 2>&1' % (cache,lc_latlon_file,tmp,log))
 os.system('mv %s %s && rm -rf %s >> %s 2>&1' % (tmp,lc_latlon_file,tmp,log))

 del data, tmp
 gc.collect()

 return

def Extract_Soils(cdb,workspace,metadata,icatch,log):

 #0. Get the parameters
 md = gdal_tools.retrieve_metadata('%s/mask_latlon.tif' % workspace)
 minx = md['minx']
 miny = md['miny']
 maxx = md['maxx']
 maxy = md['maxy']
 res = abs(md['resx'])
 lproj = md['proj4']

 # POLARIS
 if (metadata['soil']['dataset'] == 'polaris'):
  vars = ['clay','sand','silt','om','hb','thetar','thetas','ksat','lambda']

  if metadata['svp']==True:  #laura, svp
   properties = {}
   for var in vars:
    #dir_in=metadata['soil'][var]
    dir_in = metadata['soil']['directory']
    if var in ['thetas','thetar']:
     layers=glob.glob(dir_in+'theta_%s_mean*' %var.split('a')[1])
    else:
     layers=glob.glob(dir_in+'%s_mean*' %var)
    i=0
    for lay in layers:
     #if '30_60' not in lay:continue
     file_in=lay
     avrgd=(float(lay.split('mean_')[1].split('.vrt')[0].split('_')[0])+float(lay.split('mean_')[1].split('.vrt')[0].split('_')[1]))/2
     if i==0:
      os.system('mkdir %s/%s' %(workspace,var))
      properties[var]={}
     i=i+1
     file_out='%s/%s/%s_latlon_%scm.tif' % (workspace,var,var,avrgd)
     cache=int(psutil.virtual_memory().available*0.7/mb)
     os.system('gdalwarp -t_srs \'%s\' -dstnodata -9999 -r bilinear -tr %.16f %.16f -te %.16f %.16f %.16f %.16f --config GDAL_CACHEMAX %i %s %s >> %s 2>&1' % (lproj,res,res,minx,miny,maxx,maxy,cache,file_in,file_out,log))
     properties[var][str(avrgd)] = gdal_tools.read_raster(file_out)
     if var in ['om','ksat','hb']:
      badvals = (properties[var][str(avrgd)] == -9999.0)
      properties[var][str(avrgd)] = np.power(10.,properties[var][str(avrgd)])
      properties[var][str(avrgd)][badvals] = -9999.0
      #end svp block here, laura

 # Write out the data
 mask = gdal_tools.read_raster('%s/mask_latlon.tif' % workspace)
 if metadata['soil']['dataset'] != 'conus-soil':
  if metadata['svp']==True:
   for layer in properties['clay']: #laura svp
    badvals = ((properties['clay'][layer]==-9999.0) | (properties['sand'][layer]==-9999.0)) | ((properties['silt'][layer]==-9999.0) | (properties['om'][layer]==-9999.0))
    for var in ['clay','sand','silt','om']:
     file_out = '%s/%s/%s_latlon_%scm.tif' % (workspace,var,var,layer) #laura svp
     properties[var][layer][badvals] = -9999.0
     
     m2 = ( mask >= 0 ) & np.invert(badvals)
     missing_ratio = 1.0 -np.sum(m2)/float(np.sum(mask >= 0))
     if missing_ratio > 0.95 :
      #os.system('rm -rf %s' % file_out)
      import sys
      sys.stderr.write('Warning_preprocessing: %s_layer_%s_full_of_nans %s\n' % (var,layer,icatch)) #laura svp
      continue
      #return
     if var not in ['hb','lambda','thetas','thetar','ksat']:
      md = gdal_tools.retrieve_metadata(file_out)
      md['nodata'] = -9999.0
      gdal_tools.write_raster(file_out,md,properties[var][layer]) 
      #end block svp, laura

 else:
  badvals = ((properties['clay']==-9999.0) | (properties['sand']==-9999.0)) | ((properties['silt']==-9999.0))
  for var in ['clay','sand','silt']:
   file_out = '%s/%s_latlon.tif' % (workspace,var)
   properties[var][badvals] = -9999.0
   m2 = ( mask >= 0 ) & np.invert(badvals)
   missing_ratio = 1.0 -np.sum(m2)/float(np.sum(mask >= 0))
   if missing_ratio > 0.95 :
    os.system('rm -rf %s' % file_out)
    import sys
    sys.stderr.write('Error_preprocessing: %s_full_of_nans %s\n' % (var,icatch))
    return
   md = gdal_tools.retrieve_metadata(file_out)
   md['nodata'] = -9999.0
   gdal_tools.write_raster(file_out,md,properties[var])

 #Create the missing properties
 Create_Other_Soil_Properties_svp(cdb,workspace,metadata,icatch,log,properties) #modified function svp, laura

 for var in ['hb','lambda']:
   os.system('rm -rf %s/%s_latlon.tif' % (workspace,var))
 del properties
 gc.collect()

 return

def Extract_Meteorology(cdb,workspace,metadata,icatch,log):

 #Get the parameters
 md = gdal_tools.retrieve_metadata('%s/mask_latlon.tif' % workspace)
 cminlon = md['minx']
 cminlat = md['miny']
 cmaxlon = md['maxx']
 cmaxlat = md['maxy']

 #Define time info
 startdate = datetime.datetime.strptime(metadata['meteo']['startdate'],'%d%b%Y')
 enddate = datetime.datetime.strptime(metadata['meteo']['enddate'],'%d%b%Y')
  
 vars = ['tair','spfh','psurf','wind','swdown','lwdown','precip']
 shuffle(vars)
 
 #Process each variable
 for var in vars:

  #Get parameters
  var_name = metadata['meteo']['vars'][var]['name']

  #Get metadata
  file = metadata['meteo']['vars'][var]['file'] 
  file = file.replace('$YEAR',str(startdate.year))
  file = file.replace('$MTH','%02d' % startdate.month)
  fp = nc.Dataset(file,'r')
  lats = fp.variables['lat'][:]
  lons = fp.variables['lon'][:]
  undef = fp.variables[var_name]._FillValue
  fp.close()
  del fp
  gc.collect()
  
  #Set up domain (with buffer)
  iminlat = np.argmin(np.abs(lats - cminlat)) - 1
  imaxlat = np.argmin(np.abs(lats - cmaxlat)) + 1
  iminlon = np.argmin(np.abs(lons - cminlon)) - 1
  imaxlon = np.argmin(np.abs(lons - cmaxlon)) + 1
  if iminlat < 0:iminlat = 0
  if imaxlat >= lats.size:imaxlat = lats.size-1
  if iminlon < 0:iminlon = 0
  if imaxlon >= lons.size:imaxlon = lons.size-1
  minlat = lats[iminlat]
  maxlat = lats[imaxlat]
  minlon = lons[iminlon]
  maxlon = lons[imaxlon]
  resy = (lats[-1]-lats[0])/len(lats)
  resx = (lons[-1]-lons[0])/len(lons)
  res = (resx+resy)/2.

  #Determine the box size
  nlon = int(np.round((maxlon - minlon)/res + 1))
  nlat = int(np.round((maxlat - minlat)/res + 1))

  #Set the metadata
  md = {'nlat':nlat,'nlon':nlon,'minlat':minlat,'minlon':minlon,'maxlat':maxlat,'maxlon':maxlon,'res':res}
  md['undef'] = -9999.0

  #Read in data and create local copy
  date = startdate
  dt = relativedelta(months=1)
  tstep = metadata['meteo']['tstep']
  nts = int(((enddate-startdate).days+1)*(24/int(tstep.split('h')[0])))
  data = np.zeros((nts,nlat,nlon))
  while date <= enddate:

   #Open access to file
   file = metadata['meteo']['vars'][var]['file']
   #print var, date.year
   file = file.replace('$YEAR',str(date.year))
   file = file.replace('$MTH','%02d' % date.month)  
   fp = nc.Dataset(file,'r')
   tic = time.time()
   #Extract the data
   if date == startdate:
    tmp = fp.variables[var_name][:,iminlat:imaxlat+1,iminlon:imaxlon+1]
    it = 0
    ft = it+tmp.shape[0]
    data[it:ft,:,:]=tmp
   else:
    tmp = fp.variables[var_name][:,iminlat:imaxlat+1,iminlon:imaxlon+1]
    it = ft
    ft = it+tmp.shape[0]
    data[it:ft,:,:]=tmp
   print(icatch,var,date,tmp.shape,time.time()-tic,flush=True)

   fp.close()
   #Update the time step
   date = date + dt

  del fp, tmp
  gc.collect()
  #Correct for undefined values
  correction = {'precip':0.0,'tair':273.0,'wind':2.0,'spfh':0.01,'lwdown':200.0,'swdown':0.0,'psurf':90000}
  data[data == undef] = correction[var]
  data[data < -1000] = correction[var]
  #Open the output file
  #Update the units using the conversion factor
  data = metadata['meteo']['vars'][var]['factor']*data
  ncfile = '%s/%s.nc' % (workspace,var)
  md['file']=ncfile
  md['nt']=data.shape[0]
  md['tinitial'] = datetime.datetime(startdate.year,startdate.month,1,0)
  md['tinitial_all'] = md['tinitial']
  md['tstep'] = tstep
  md['vars'] = [var]
  fp = Create_NETCDF_File(md)
  
  #Write the data
  fp.variables[var][:] = data
  #Close the file
  fp.close()
  del fp, data
  gc.collect()

  #Create a sample grid using the mask
  mask_latlon_file = '%s/mask_latlon.tif' % (workspace)
  file_coarse = '%s/%s_latlon_coarse.tif' % (workspace,var)
  cache = int(psutil.virtual_memory().available*0.7/mb)
  os.system('gdalwarp -tr %.16f %.16f -te %.16f %.16f %.16f %.16f --config GDAL_CACHEMAX %i %s %s >> %s 2>&1' % (res,res,minlon-res/2,minlat-res/2,maxlon+res/2,maxlat+res/2,cache,mask_latlon_file,file_coarse,log))

  #Define the coarse and fine scale mapping
  maskij = gdal_tools.read_raster(file_coarse)
  metadata_maskij = gdal_tools.retrieve_metadata(file_coarse)
  for i in np.arange(maskij.shape[0]):
   maskij[i,:] = np.arange(i*maskij.shape[1],(i+1)*maskij.shape[1])
  metadata_maskij['nodata'] = -9999.0
  gdal_tools.write_raster(file_coarse,metadata_maskij,np.flipud(maskij))
  del maskij
  gc.collect()

  #Get the parameters
  md = gdal_tools.retrieve_metadata('%s/mask_latlon.tif' % workspace)
  minx = md['minx']
  miny = md['miny']
  maxx = md['maxx']
  maxy = md['maxy']
  res = abs(md['resx'])
  lproj = md['proj4']
  
  #Regrid and downscale
  file_in = file_coarse
  file_out = '%s/%s_latlon_fine.tif' % (workspace,var)
  cache = int(psutil.virtual_memory().available*0.7/mb)
  os.system('gdalwarp -t_srs \'%s\' -dstnodata -9999 -tr %.16f %.16f -te %.16f %.16f %.16f %.16f --config GDAL_CACHEMAX %i %s %s >> %s 2>&1' % (lproj,res,res,minx,miny,maxx,maxy,cache,file_in,file_out,log))

 return

def Extract_Meteorology_Daily(cdb,workspace,metadata,icatch,log):

 #Get the parameters
 md = gdal_tools.retrieve_metadata('%s/mask_latlon.tif' % workspace)
 cminlon = md['minx']
 cminlat = md['miny']
 cmaxlon = md['maxx']
 cmaxlat = md['maxy']

 #Define time info
 startdate = datetime.datetime.strptime(metadata['meteo']['startdate'],'%d%b%Y')
 enddate = datetime.datetime.strptime(metadata['meteo']['enddate'],'%d%b%Y')
  
 vars = ['tair','spfh','psurf','wind','swdown','lwdown','precip']
 shuffle(vars)
 
 #Process each variable
 for var in vars:

  #Get parameters
  #var_name = metadata['meteo']['vars'][var]['name']
  var_name = var

  #Get metadata
  #file = metadata['meteo']['vars'][var]['file'] 
  file = '%s/%s' % (metadata['meteo']['dir'],metadata['meteo']['file_template'])
  file = file.replace('$YEAR',str(startdate.year))
  file = file.replace('$MTH','%02d' % startdate.month)
  file = file.replace('$DAY','%02d' % startdate.day)
  fp = nc.Dataset(file,'r')
  lats = fp.variables['lat'][:]
  lons = fp.variables['lon'][:]
  undef = fp.variables[var_name]._FillValue
  fp.close()
  del fp
  gc.collect()
  
  #Set up domain (with buffer)
  iminlat = np.argmin(np.abs(lats - cminlat)) - 1
  imaxlat = np.argmin(np.abs(lats - cmaxlat)) + 1
  iminlon = np.argmin(np.abs(lons - cminlon)) - 1
  imaxlon = np.argmin(np.abs(lons - cmaxlon)) + 1
  if iminlat < 0:iminlat = 0
  if imaxlat >= lats.size:imaxlat = lats.size-1
  if iminlon < 0:iminlon = 0
  if imaxlon >= lons.size:imaxlon = lons.size-1
  minlat = lats[iminlat]
  maxlat = lats[imaxlat]
  minlon = lons[iminlon]
  maxlon = lons[imaxlon]
  resy = (lats[-1]-lats[0])/len(lats)
  resx = (lons[-1]-lons[0])/len(lons)
  res = (resx+resy)/2.

  #Determine the box size
  nlon = int(np.round((maxlon - minlon)/res + 1))
  nlat = int(np.round((maxlat - minlat)/res + 1))
  
  #Set the metadata
  md = {'nlat':nlat,'nlon':nlon,'minlat':minlat,'minlon':minlon,'maxlat':maxlat,'maxlon':maxlon,'res':res}
  md['undef'] = -9999.0

  #Read in data and create local copy
  date = startdate
  dt = relativedelta(days=1)
  tstep = metadata['meteo']['tstep']
  nts = int(((enddate-startdate).days+1)*(24/int(tstep.split('h')[0])))
  data = np.zeros((nts,nlat,nlon))
  while date <= enddate:

   #Open access to file
   file = '%s/%s' % (metadata['meteo']['dir'],metadata['meteo']['file_template'])
   #file = metadata['meteo']['vars'][var]['file']
   #print var, date.year
   file = file.replace('$YEAR',str(date.year))
   file = file.replace('$MTH','%02d' % date.month)  
   file = file.replace('$DAY','%02d' % date.day)  
   fp = nc.Dataset(file,'r')
   #print file 
   tic = time.time()
   #Extract the data
   if date == startdate:
    tmp = fp.variables[var_name][:,iminlat:imaxlat+1,iminlon:imaxlon+1]
    it = 0
    ft = it+tmp.shape[0]
    data[it:ft,:,:]=tmp
   else:
    tmp = fp.variables[var_name][:,iminlat:imaxlat+1,iminlon:imaxlon+1]
    it = ft
    ft = it+tmp.shape[0]
    data[it:ft,:,:]=tmp
   print(icatch,var,date,tmp.shape,time.time()-tic,flush=True)

   fp.close()
   #Update the time step
   date = date + dt

  del fp, tmp
  gc.collect()
  #Correct for undefined values
  correction = {'precip':0.0,'tair':273.0,'wind':2.0,'spfh':0.01,'lwdown':200.0,'swdown':0.0,'psurf':90000}
  data[data == undef] = correction[var]
  data[data < -1000] = correction[var]
  #Open the output file
  #Update the units using the conversion factor
  #data = metadata['meteo']['vars'][var]['factor']*data
  data = 1.0*data
  ncfile = '%s/%s.nc' % (workspace,var)
  md['file']=ncfile
  md['nt']=data.shape[0]
  md['tinitial'] = datetime.datetime(startdate.year,startdate.month,1,0)
  md['tinitial_all'] = md['tinitial']
  md['tstep'] = tstep
  md['vars'] = [var]
  fp = Create_NETCDF_File(md)
  
  #Write the data
  fp.variables[var][:] = data
  #Close the file
  fp.close()
  del fp, data
  gc.collect()

  #Create a sample grid using the mask
  mask_latlon_file = '%s/mask_latlon.tif' % (workspace)
  file_coarse = '%s/%s_latlon_coarse.tif' % (workspace,var)
  cache = int(psutil.virtual_memory().available*0.7/mb)
  os.system('gdalwarp -tr %.16f %.16f -te %.16f %.16f %.16f %.16f --config GDAL_CACHEMAX %i %s %s >> %s 2>&1' % (res,res,minlon-res/2,minlat-res/2,maxlon+res/2,maxlat+res/2,cache,mask_latlon_file,file_coarse,log))

  #Define the coarse and fine scale mapping
  maskij = gdal_tools.read_raster(file_coarse)
  metadata_maskij = gdal_tools.retrieve_metadata(file_coarse)
  for i in np.arange(maskij.shape[0]):
   maskij[i,:] = np.arange(i*maskij.shape[1],(i+1)*maskij.shape[1])
  metadata_maskij['nodata'] = -9999.0
  gdal_tools.write_raster(file_coarse,metadata_maskij,np.flipud(maskij))
  del maskij
  gc.collect()

  #Get the parameters
  md = gdal_tools.retrieve_metadata('%s/mask_latlon.tif' % workspace)
  minx = md['minx']
  miny = md['miny']
  maxx = md['maxx']
  maxy = md['maxy']
  res = abs(md['resx'])
  lproj = md['proj4']
  
  #Regrid and downscale
  file_in = file_coarse
  file_out = '%s/%s_latlon_fine.tif' % (workspace,var)
  cache = int(psutil.virtual_memory().available*0.7/mb)
  os.system('gdalwarp -t_srs \'%s\' -dstnodata -9999 -tr %.16f %.16f -te %.16f %.16f %.16f %.16f --config GDAL_CACHEMAX %i %s %s >> %s 2>&1' % (lproj,res,res,minx,miny,maxx,maxy,cache,file_in,file_out,log))

 return

def prepare_input_data(cdir,cdb,metadata,rank,icatch):
 
 #Create the workspace
 #workspace = '%s/workspace' % cdir
 workspace = cdir
 #os.system('mkdir -p %s' % workspace)

 #Define the log
 log = '%s/log.txt' % workspace

 #Create the mask
 print(rank,'Preparing catchment mask',time.ctime(),icatch,flush=True)
 Create_Mask(cdb,workspace,metadata,icatch,log)

 #Create the administrative boundaries map
 print(rank,'Preparing administrative boundaries mask',time.ctime(),icatch,flush=True)
 Create_Administrative_Boundaries(cdb,workspace,metadata,icatch,log)

 #Terrain analysis
 print(rank,'Preparing dem products',time.ctime(),icatch,flush=True)
 Terrain_Analysis(cdb,workspace,metadata,icatch,log)
 if os.path.isfile('%s/dem_latlon.tif' % workspace) == False: return

 #Correct the mask
 print(rank,'Correcting catchment mask',time.ctime(),icatch,flush=True)
 Correct_Mask(cdb,workspace,metadata,icatch,log)
   
 #Create land cover product
 print(rank,'Preparing land cover data',time.ctime(),icatch,flush=True)
 Extract_Land_Cover(cdb,workspace,metadata,icatch,log)
 
 #Create soil product
 print(rank,'Preparing the soil data',time.ctime(),icatch,flush=True)
 Extract_Soils(cdb,workspace,metadata,icatch,log)

 #Create depth to bedrock product
 print(rank,'Preparing the depth to the bedrock',time.ctime(),icatch,flush=True)
 Extract_Dbedrock(cdb,workspace,metadata,icatch,log)

 #Create meteorology product
 print(rank,'Preparing the meteorological data',time.ctime(),icatch,flush=True)
 #Extract_Meteorology(cdb,workspace,metadata,icatch,log)
 Extract_Meteorology_Daily(cdb,workspace,metadata,icatch,log)

 return

def Create_NETCDF_File(md):

 nlat = md['nlat']
 nlon = md['nlon']
 res = md['res']
 minlon = md['minlon'] + res/2
 minlat = md['minlat'] + res/2
 undef = md['undef']
 nt = md['nt']
 tstep = md['tstep']
 tinitial = md['tinitial']
 tinitial_all = md['tinitial_all']
 vars = md['vars']
 if 'vars_info' in md:vars_info = md['vars_info']
 else:vars_info = vars
 file = md['file']

 #Determine the initial it
 it = int((tinitial-tinitial_all).total_seconds()/3600.0)
 if nt > 0:t = np.arange(it,nt+it)

 #Prepare the netcdf file
 #Create file
 f = nc.Dataset(file, 'w')

 #Define dimensions
 f.createDimension('lon',nlon)
 f.createDimension('lat',nlat)
 if nt > 0:f.createDimension('t',None)#len(t))

 #Longitude
 f.createVariable('lon','d',('lon',))
 f.variables['lon'][:] = np.linspace(minlon,minlon+res*(nlon-1),nlon)
 f.variables['lon'].units = 'degrees_east'
 f.variables['lon'].long_name = 'Longitude'
 f.variables['lon'].res = res

 #Latitude
 f.createVariable('lat','d',('lat',))
 f.variables['lat'][:] = np.linspace(minlat,minlat+res*(nlat-1),nlat)
 f.variables['lat'].units = 'degrees_north'
 f.variables['lat'].long_name = 'Latitude'
 f.variables['lat'].res = res

 #Time
 if nt > 0:
  times = f.createVariable('t','d',('t',))
  f.variables['t'][:] = t
  f.variables['t'].units = '%s since %04d-%02d-%02d %02d:00:00.0' % (tstep,tinitial_all.year,tinitial_all.month,tinitial_all.day,tinitial_all.hour)
  f.variables['t'].long_name = 'Time'

 #Data
 i = 0
 for var in vars:
  if nt > 0:f.createVariable(var,'f',('t','lat','lon'),fill_value=undef)#,zlib=True)
  else: f.createVariable(var,'f',('lat','lon'),fill_value=undef)#,zlib=True)
  f.variables[var].long_name = vars_info[i]
  i = i + 1

 return f

def flip(m, axis):
    if not hasattr(m, 'ndim'):
        m = asarray(m)
    indexer = [slice(None)] * m.ndim
    try:
        indexer[axis] = slice(None, None, -1)
    except IndexError:
        raise ValueError("axis=%i is invalid for the %i-dimensional input array"
                         % (axis, m.ndim))
    return m[tuple(indexer)]

def correct_domain_decomposition(comm,metadata):

 size = comm.Get_size()
 rank = comm.Get_rank()

 #Read in the catchment summary database
 pck_file = '%s/cids/domain_database.pck' % metadata['output_data']
 cdb = pickle.load(open(pck_file,'rb'))
 crange = range(len(cdb))
 odb = {}
 for ic in crange[rank::size]:

  print("Rank:%d, Catchment:%s - Initializing" % (rank,ic),flush=True)

  cid = cdb[ic]['cid']

  #Define the catchment directory
  cdir = '%s/cids/%d' % (metadata['output_data'],cid)
  workspace = cdir

  #Define the log
  log = '%s/log.txt' % cdir

  #Prepare the workspace for the sub-domain
  os.system('rm -rf %s' % cdir)
  os.system('mkdir -p %s' % cdir)
  #Create the mask
  print(rank,'Preparing catchment mask',time.ctime(),cid,flush=True)
  Create_Mask(cdb[ic],cdir,metadata,cid,log)

  print(rank,'Preparing administrative boundaries mask',time.ctime(),cid,flush=True)
  Create_Administrative_Boundaries(cdb[ic],workspace,metadata,cid,log)

  #Terrain analysis
  print(rank,'Preparing dem products',time.ctime(),cid,flush=True)
  Terrain_Analysis(cdb[ic],cdir,metadata,cid,log)

  #Correct the mask
  print(rank,'Correcting catchment mask',time.ctime(),cid,flush=True)
  Correct_Mask(cdb[ic],cdir,metadata,cid,log)

  #Create land cover product
  print(rank,'Preparing land cover data',time.ctime(),cid,flush=True)
  Extract_Land_Cover(cdb[ic],cdir,metadata,cid,log)

  #Create soil product
  print(rank,'Preparing the soil data',time.ctime(),cid,flush=True)
  Extract_Soils(cdb[ic],cdir,metadata,cid,log)

  #Create meteo file
  #file_in = '/gpfs/f5/gfdl_b/proj-shared/Nathaniel.Chaney/datasets/PCF/1hr/tair.tif'
  file_in = '%s/tair.tif' % (metadata['meteo']['dir'],)
  md = gdal_tools.retrieve_metadata('%s/mask_latlon.tif' % cdir)
  minx = md['minx']
  miny = md['miny']
  maxx = md['maxx']
  maxy = md['maxy']
  res = abs(md['resx'])
  lproj = md['proj4']
  file_out = '%s/meteo_latlon.tif' % cdir
  cache = int(psutil.virtual_memory().available*0.7/mb)
  os.system('gdalwarp -t_srs \'%s\' -dstnodata -9999 -r bilinear -tr %.16f %.16f -te %.16f %.16f %.16f %.16f --config GDAL_CACHEMAX %i %s %s >> %s 2>&1' % (lproj,res,res,minx,miny,maxx,maxy,cache,file_in,file_out,log))

  #Determine number of pixels
  file = '%s/mask_latlon.tif' % workspace
  mask = rasterio.open(file).read(1)
  file = '%s/sand/sand_latlon_2.5cm.tif' % workspace
  sand = rasterio.open(file).read(1)
  file = '%s/meteo_latlon.tif' % workspace
  meteo = rasterio.open(file).read(1)
  file = '%s/lc_latlon.tif' % workspace
  lc = rasterio.open(file).read(1)
  npx_mask = np.sum(mask == cid)
  npx_sand = np.sum(sand != -9999)
  npx_meteo = np.sum(meteo != -9999)
  npx_lc = np.sum(lc != -9999)
  odb[cid] = min(npx_mask,npx_sand,npx_meteo,npx_lc)
  #odb[cid] = min(npx_mask,npx_meteo,npx_lc)
  print(npx_mask,npx_sand,npx_meteo,npx_lc,flush=True)
  #print(npx_mask,npx_meteo,npx_lc,flush=True)

 #Broadcast and collect
 if rank == 0:
   print(rank,'Cleaning up cid map',time.ctime(),flush=True)
   for i in range(1,size):
    odb2 = comm.recv(source=i, tag=11)
    for key in odb2:odb[key] = odb2[key]
 else:
    comm.send(odb,dest=0, tag=11)

 #Create new shapefile and summary
 if rank == 0:
  skeys = np.array(sorted(odb.keys()))
  count = 1
  odb2 = collections.OrderedDict()
  for key in skeys:
   if odb[key] > 250:
    odb2[key] = count
    count += 1
   
  '''odb2 = {}
  count = 1
  for key in odb:
   #if odb[key] > 0:
   if odb[key] > 250:#Minimum number of valid pixels in the cid of the domain to count
    odb2[key] = count
    count += 1'''
    
  odb = odb2
  dfile = '%s/shp/domain.shp' % metadata['output_data']
  dfile2 = '%s/shp/domain2.shp' % metadata['output_data']
  fp = fiona.open(dfile,'r')
  fp2 = fiona.open(dfile2,'w',crs=fp.crs,driver='ESRI Shapefile',schema=fp.schema)
  for poly in fp.values():
   ID = poly['properties']['ID']
   if ID in odb:
    poly['id'] = odb[ID]-1
    poly['properties']['ID'] = odb[ID]
    fp2.write(poly)
  fp.close()
  fp2.close()
  mdir = metadata['output_data']
  gdir = mdir
  #gdir = '%s/general' % metadata['dir']
  os.system('mv %s/shp/domain2.shp %s/shp/domain.shp' % (mdir,gdir))
  os.system('mv %s/shp/domain2.shx %s/shp/domain.shx' % (mdir,gdir))
  os.system('mv %s/shp/domain2.prj %s/shp/domain.prj' % (mdir,gdir))
  os.system('mv %s/shp/domain2.cpg %s/shp/domain.cpg' % (mdir,gdir))
  os.system('mv %s/shp/domain2.dbf %s/shp/domain.dbf' % (mdir,gdir))
  #perform new summary
  summarize_domain_decompisition(metadata)

 #Wait until rank 0 completes
 comm.Barrier()

 #Remove all cdir
 for ic in crange[rank::size]:

  cid = cdb[ic]['cid']

  #Define the catchment directory
  cdir = '%s/cids/%d' % (metadata['output_data'],cid)
  os.system('rm -rf %s' % cdir)

 comm.Barrier()

 return

def Create_Administrative_Boundaries(cdb,workspace,metadata,icatch,log):

 #Define parameters
 shp_in = metadata['administrative_boundaries']
 shp_out = '%s/admin_shp' % workspace
 ci = cdb['cid']
 bbox =  cdb['bbox']
 res_latlon = metadata['res_latlon']
 
 #Define the files
 admin_latlon_file = '%s/admin_latlon.tif' % workspace
 tmp_file = '%s/tmp.tif' % workspace
 
 #Rasterize the area
 buff = 0.1
 
 print(' buffer size:',buff,' icatch:',ci,flush=True) 
 minx = bbox['minlon']-buff
 miny = bbox['minlat']-buff
 maxx = bbox['maxlon']+buff
 maxy = bbox['maxlat']+buff
 cache = int(psutil.virtual_memory().available*0.7/mb)

 #Correct coordinates to avoid reprojections
 #Fix coordinates to the dem vrt to avoid inconsistiences
 md = gdal_tools.retrieve_metadata(metadata['dem'])
 res_latlon = np.abs(md['resx'])
 minx = md['minx'] + np.floor((minx-md['minx'])/res_latlon)*res_latlon
 miny = md['miny'] + np.floor((miny-md['miny'])/res_latlon)*res_latlon#np.floor(miny/res_latlon)*res_latlon
 maxx = md['minx'] + np.ceil((maxx-md['minx'])/res_latlon)*res_latlon
 maxy = md['miny'] + np.ceil((maxy-md['miny'])/res_latlon)*res_latlon

 #Rasterize
 os.system("gdal_rasterize -at -ot Float64 --config GDAL_CACHEMAX %i -a_nodata -9999 -init -9999 -tr %.16f %.16f -te %.16f %.16f %.16f %.16f -burn 1 %s %s >> %s 2>&1" % (cache, res_latlon,res_latlon,minx,miny,maxx,maxy,shp_in,admin_latlon_file,log))

 #del data
 gc.collect()

 return


