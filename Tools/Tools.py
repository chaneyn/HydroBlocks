import os
import numpy as np
import cPickle as pickle
import ordereddict
import sys
sys.path.append('Model')
import HydroBloks as HB
import datetime
import time
import netCDF4 as nc
import gdal_tools

def Create_LHS_parameter_set(nsamples):

 from SALib.sample import latin_hypercube
 from SALib.util import scale_samples, read_param_file
 import random as rd

 # Set random seed (does not affect quasi-random Sobol sampling)
 seed = 1
 np.random.seed(seed)
 rd.seed(seed)

 #Define parameters and ranges
 parameters = ordereddict.OrderedDict()
 parameters['log10m'] = [np.log10(0.001),np.log10(0.1)]
 parameters['lnTe'] = [np.log(np.exp(-8.0)/3600.0),np.log(np.exp(8.0)/3600.0)]
 parameters['log10soil'] = [np.log10(1.0),np.log10(2.00)]
 parameters['sdmax'] = [0.1,2.0] #dtopmodel

 #Make directory
 if os.path.exists('LHS') == False:
  os.mkdir('LHS')

 #Prepare file with parameter range
 fp = open('LHS/parameters.txt','w')
 vars = []
 for var in parameters:
  vars.append(var)
  fp.write('%s %f %f\n' % (var,parameters[var][0],parameters[var][1]))
 fp.close()

 #Read the parameter range file and generate samples
 param_file = 'LHS/parameters.txt'
 pf = read_param_file(param_file)

 #Generate samples (choose method here)
 param_values = latin_hypercube.sample(nsamples, pf['num_vars'])

 #Samples are given in range [0, 1] by default. Rescale them to your parameter bounds.
 scale_samples(param_values, pf['bounds'])

 #Save parameters to file         
 np.savetxt('LHS/LHS_sampling.txt', param_values, delimiter=' ',header=" ".join(vars))

 return

def Initialize_Output():

 variables = {
           #Fluxes
           'qout_subsurface':[],
           'qout_surface':[],
           #State variables
           'smc1':[], #soil moisture content (layer 1)
           'smc2':[], #soil moisture content (...)
           'smc3':[], #soil moisture content
           'smc4':[], #soil moisture content
           #'smc5':[], #soil moisture content
           'swd':[], #water deficit
           'sstorage':[], #storage surface
           'sndpth':[], #snow depth
	   'swe':[], #snow water equivalent (mm)
           #Water balance
 	   'et':[], #evapotranspiration
           'qexcess':[],
           'qsurface':[],
           'prcp':[], 
           #Energy fluxes
           'g':[], #ground heat flux (w/m2)
           'sh':[], #sensible heat flux (w/m2)
           'lh':[], #latent heat flux (w/m2)
 	   'swnet':[], #absorved shortwave radiation (w/m2)
           'lwnet':[], #net upward longwave radiation (w/m2)
           }

 misc = {
      'sldpth':[], #soil depth
      'dates':[], #dates
      }
           
 return {'variables':variables,'misc':misc}

def update_output(output,NOAH,TOPMODEL,date):
 
 #Miscellanous
 output['misc']['sldpth'] = NOAH.sldpth
 output['misc']['dates'].append(date)
 output['misc']['pct'] = TOPMODEL.pct
 output['misc']['dx'] = TOPMODEL.dx
 output['misc']['area'] = TOPMODEL.area
 output['misc']['outlet_hsu'] = TOPMODEL.outlet_hsu
 metadata = {
             'g':{'description':'Ground heat flux','units':'W/m2'},
             'sh':{'description':'Sensible heat flux','units':'W/m2'},
             'lh':{'description':'Latent heat flux','units':'W/m2'},
             'lwnet':{'description':'Net longwave radiation','units':'W/m2'},
             'swnet':{'description':'Absorbed shortwave radiation','units':'W/m2'},
             'et':{'description':'Evapotranspiration','units':'mm/s'},
             'qexcess':{'description':'Excess runoff','units':'mm/s'},
             'qsurface':{'description':'Surface runoff','units':'mm/s'},
             'prcp':{'description':'Precipitation','units':'mm/s'},
	     'smc1':{'description':'Soil water content','units':'m3/m3'},
             'smc2':{'description':'Soil water content','units':'m3/m3'},
             'smc3':{'description':'Soil water content','units':'m3/m3'},
             'smc4':{'description':'Soil water content','units':'m3/m3'},
             #'smc5':{'description':'Soil water content','units':'m3/m3'},
             'swd':{'description':'Soil water deficit','units':'mm'},
             'sstorage':{'description':'Surface storage','units':'mm'},
             'sndpth':{'description':'Snow depth','units':'mm'},
             'swe':{'description':'Snow water equivalent','units':'mm'},
             'qout_subsurface':{'description':'Subsurface flux','units':'m2/s'},
             'qout_surface':{'description':'Surface flux','units':'m2/s'},
             }
 output['misc']['metadata'] = metadata
 #Energy fluxes
 output['variables']['g'].append(np.copy(NOAH.ssoil)) #W/m2
 output['variables']['sh'].append(np.copy(NOAH.fsh)) #W/m2
 output['variables']['lh'].append(np.copy(NOAH.fcev + NOAH.fgev + NOAH.fctr)) #W/m2
 output['variables']['swnet'].append(np.copy(NOAH.fsa)) #W/m2
 output['variables']['lwnet'].append(np.copy(NOAH.fira)) #W/m2
 #Water balance
 output['variables']['et'].append(np.copy(NOAH.ecan + NOAH.etran + NOAH.esoil)) #mm/s
 imax = TOPMODEL.outlet_hsu
 qexcess = np.copy(NOAH.runsb)
 qexcess[imax] = qexcess[imax] + 10**3*(TOPMODEL.dx[imax]*(TOPMODEL.qout[imax])/np.sum(TOPMODEL.area))
 output['variables']['qexcess'].append(np.copy(qexcess)) #mm/s
 output['variables']['qsurface'].append(np.copy(NOAH.runsf)) #mm/s
 output['variables']['prcp'].append(np.copy(NOAH.prcp)) #mm/s
 #State variables
 output['variables']['smc1'].append(np.copy(NOAH.smc[:,0])) #m3/m3
 output['variables']['smc2'].append(np.copy(NOAH.smc[:,1])) #m3/m3
 output['variables']['smc3'].append(np.copy(NOAH.smc[:,2])) #m3/m3
 output['variables']['smc4'].append(np.copy(NOAH.smc[:,3])) #m3/m3
 output['variables']['swd'].append(np.copy(10**3*TOPMODEL.si)) #mm
 output['variables']['sstorage'].append(np.copy(10**3*TOPMODEL.storage_surface)) #mm
 output['variables']['sndpth'].append(np.copy(10**3*NOAH.sndpth)) #mm
 output['variables']['swe'].append(np.copy(NOAH.swe)) #mm
 #Fluxes
 output['variables']['qout_subsurface'].append(np.copy(TOPMODEL.qout)) #m2/s
 output['variables']['qout_surface'].append(np.copy(TOPMODEL.qout_surface)) #m2/s

 return output

def run_model(parameters,cdir,idate,fdate,ncores,wbd):

 #Read in the pickled data
 data = pickle.load(open('%s/input/data.pck' % cdir))
 ncells = len(data['hsu'].keys())
 dt = 3600.
 nsoil = 20
 dt_timedelta = datetime.timedelta(seconds=dt)
  
 #Initialize the output
 output = Initialize_Output()
 
 #Define the parameter files
 info = {}
 info['VEGPARM'] = 'Model/pyNoahMP/data/VEGPARM.TBL'#'data/VEGPARM.TBL'
 info['GENPARM'] = 'Model/pyNoahMP/data/GENPARM.TBL'#'data/GENPARM.TBL'
 info['MPTABLE'] = 'Model/pyNoahMP/data/MPTABLE.TBL'#'pyNoahMP/data/MPTABLE.TBL'
 info['SOILPARM'] = 'Model/pyNoahMP/data/SOILPARM.TBL'#'data/SOILPARM.TBL'

 #Initialize the model
 NOAH = HB.Initialize_Model(ncells,dt,nsoil,data,parameters,info,wbd)
 output['misc']['zsoil'] = NOAH.zsoil

 #Initialize Topmodel
 TOPMODEL = HB.Initialize_DTopmodel(ncells,dt,data,parameters) 

 #Run the model
 meteorology = data['meteorology']
 date = idate
 i = 0
 dE = 0
 tic = time.time()
 #errwat = 0
 r = 0
 dr = 0
 et,etran,esoil,ecan = 0,0,0,0
 prcp = 0
 q = 0
 errwat = 0
 NOAH.dzwt[:] = 0.0
 TOPMODEL.ex[:] = 0.0
 while date <= fdate:
  if (date.hour == 0) and (date.day == 1):
   print date,np.sum(TOPMODEL.pct*TOPMODEL.si),time.time() - tic,et,prcp,'q:%f'%q,np.sum(TOPMODEL.pct*NOAH.smcwtd),np.sum(TOPMODEL.pct*NOAH.zwt),'WB ERR:%f' % np.sum(TOPMODEL.pct*errwat),'etran:%f'%etran,'ecan:%f'%ecan,'esoil:%f'%esoil
  if (date.month == 1) and (date.day == 1) and (date.hour == 0):
   et = 0
   prcp = 0
   q = 0
   ecan,etran,esoil = 0,0,0
 
  #Update input data
  Update_Input(NOAH,TOPMODEL,date,meteorology,i)

  #Calculate initial NOAH water balance
  smw = 1000.0*((NOAH.smcmax*np.abs(NOAH.zwt) - TOPMODEL.si) + (NOAH.smcmax*(100-np.abs(NOAH.zwt))))
  beg_wb = np.copy(NOAH.canliq + NOAH.canice + NOAH.swe + NOAH.wa + smw)#TOPMODEL.si*1000)
  dzwt0 = np.copy(NOAH.dzwt)# - NOAH.deeprech)

  #Update model
  (NOAH,TOPMODEL) = HB.Update_Model(NOAH,TOPMODEL,ncores)

  #Calculate final water balance
  smw = 1000.0*((NOAH.smcmax*np.abs(NOAH.zwt) - TOPMODEL.si) + (NOAH.smcmax*(100.0-np.abs(NOAH.zwt))))
  end_wb = np.copy(NOAH.canliq + NOAH.canice + NOAH.swe + NOAH.wa + smw)# + TOPMODEL.si*1000)

  #Update the water balance error
  tmp = np.copy(end_wb - beg_wb - NOAH.dt*(NOAH.prcp-NOAH.ecan-NOAH.etran-NOAH.esoil-NOAH.runsf-NOAH.runsb) - 1000*dzwt0)
  errwat += tmp# - tmp
  q = q + 3600.*np.sum(TOPMODEL.pct*NOAH.runsb)
  et = et + 3600*np.sum(TOPMODEL.pct*(NOAH.ecan + NOAH.etran + NOAH.esoil))
  etran += 3600*np.sum(TOPMODEL.pct*NOAH.etran)
  ecan += 3600*np.sum(TOPMODEL.pct*NOAH.ecan)
  esoil += 3600*np.sum(TOPMODEL.pct*NOAH.esoil)
  prcp = prcp + 3600*np.sum(TOPMODEL.pct*NOAH.prcp)

  #Update output
  output = update_output(output,NOAH,TOPMODEL,date)

  #Update time step
  date = date + dt_timedelta
  i = i + 1

 #Write the output
 for var in output['variables']:
  output['variables'][var] = np.array(output['variables'][var])
 pickle.dump(output,open('%s/output/output.pck' % cdir,'wb'),pickle.HIGHEST_PROTOCOL)

 #Finalize the model
 HB.Finalize_Model(NOAH,TOPMODEL)

 return

def Update_Input(NOAH,TOPMODEL,date,meteorology,i):

  NOAH.itime = i
  TOPMODEL.itime = i
  NOAH.nowdate[:] = date.strftime('%Y-%m-%d_%H:%M:%S')
  NOAH.julian = (date - datetime.datetime(date.year,1,1,0)).days
  NOAH.yearlen = (datetime.datetime(date.year+1,1,1,0) - datetime.datetime(date.year,1,1,1,0)).days + 1

  #Update meteorology
  NOAH.lwdn[:] = meteorology['nldas_dlwrf'][i,:] #W/m2
  NOAH.swdn[:] = meteorology['nldas_dswrf'][i,:] #W/m2
  NOAH.psfc[:] = meteorology['nldas_pres'][i,:] #Pa
  NOAH.p_ml[:] = meteorology['nldas_pres'][i,:] #Pa
  NOAH.u_ml[:] = (meteorology['nldas_wind'][i,:]**2/2)**0.5 #m/s
  NOAH.v_ml[:] = (meteorology['nldas_wind'][i,:]**2/2)**0.5 #m/s
  NOAH.t_ml[:] = 273.15+meteorology['nldas_tair'][i,:] #K
  estar = 611.0*np.exp(17.27*((NOAH.t_ml - 273.15)/(NOAH.t_ml[:] - 36.0)))
  e = meteorology['nldas_rh'][i,:]*estar/100
  q = 0.622*e/NOAH.psfc
  NOAH.q_ml[:] = q #Kg/Kg
  NOAH.qsfc1d[:] = q #Kg/Kg
  NOAH.prcp[:] = meteorology['nldas_prec'][i,:]/3600.0 #mm/s

  return (NOAH,TOPMODEL)

def create_netcdf_file(file_netcdf,output,nens,input,cdir):

 #Create the file
 fp = nc.Dataset(file_netcdf, 'w', format='NETCDF4')

 #Create the dimensions
 ntime = output['variables']['lh'].shape[0]/24
 nhsu = output['variables']['lh'].shape[1]
 nsoil = 4
 fp.createDimension('hsu',nhsu)
 fp.createDimension('time',ntime)
 fp.createDimension('ensemble',nens)
 fp.createDimension('soil',nsoil)

 #Create the summary group and its variables
 grp = fp.createGroup('summary')
 for var in output['variables']:
  ncvar = grp.createVariable(var,'f4',('time','ensemble'))
  ncvar.description = output['misc']['metadata'][var]['description']
  ncvar.units = output['misc']['metadata'][var]['units']

 #Create the output
 grp = fp.createGroup('catchment')
 for var in output['variables']:
  grp.createVariable(var,'f4',('time','hsu','ensemble'))
  ncvar.description = output['misc']['metadata'][var]['description']
  ncvar.units = output['misc']['metadata'][var]['units']

 #Create the metadata
 grp = fp.createGroup('metadata')
 #dates
 times = grp.createVariable('dates','f8',('time',))
 dates = []
 for date in output['misc']['dates']:
  if date.hour == 0:dates.append(date)
 times.units = 'days since 1900-01-01'
 times.calendar = 'standard'
 times[:] = nc.date2num(np.array(dates),units=times.units,calendar=times.calendar)
 #Soil depth
 soil = grp.createVariable('soil','f4',('soil',))
 soil[:] = np.array(output['misc']['sldpth'])[0,0:nsoil]
 soil.units = 'meters'
 soil.description = 'soil depth'
 #HSU percentage coverage
 pcts = grp.createVariable('pct','f4',('hsu',))
 pcts[:] = 100*np.array(output['misc']['pct'])
 pcts.description = 'hsu percentage coverage'
 pcts.units = '%'
 #HSU Spatial resolution
 dx = grp.createVariable('dx','f4',('hsu',))
 dx[:] = np.array(output['misc']['dx'])
 dx.units = 'meters'
 #HSU area
 area = grp.createVariable('area','f4',('hsu',))
 area[:] = np.array(output['misc']['area'])
 area.units = 'meters squared'
 #HSU ids
 hsu = grp.createVariable('hsu','i4',('hsu',))
 hsus =[]
 for value in xrange(len(output['misc']['pct'])):hsus.append(value)
 hsu[:] = np.array(hsus)
 hsu.description = 'hsu ids'
 #Add wbd metadata
 grp.HUC = input['wbd']['HUC10']
 #Define outlet hsu
 grp.outlet_hsu = int(output['misc']['outlet_hsu'])
 #Define catchment name
 grp.catchment_name = input['wbd']['Name']
 #Define catchment area
 grp.AreaSqKm = input['wbd']['AreaSqKm']
 #Add HSU transition probability matrix
 tp = grp.createVariable('tpm','f4',('hsu','hsu'))
 tp.description = 'Transition probability matrix between hsus'
 tp[:] =  input['tp']

 #Retrieve the conus_albers metadata
 metadata = gdal_tools.retrieve_metadata(input['wbd']['files']['ti']) 
 metadata['nodata'] = -9999.0
 #Save the conus_albers metadata
 fp.createDimension('nx',metadata['nx'])
 fp.createDimension('ny',metadata['ny'])
 hmca = grp.createVariable('hmca','f4',('ny','nx')) 
 hmca.gt = metadata['gt']
 hmca.projection = metadata['projection']
 hmca.description = 'HSU mapping (conus albers)'
 hmca.nodata = metadata['nodata']
 #Save the conus albers mapping
 hsu_map = np.copy(input['hsu_map'])
 hsu_map[np.isnan(hsu_map) == 1] = metadata['nodata']
 hmca[:] = hsu_map
 #Write out the mapping
 file_ca = '%s/workspace/hsu_mapping_conus_albers.tif' % cdir
 gdal_tools.write_raster(file_ca,metadata,hsu_map)

 #Map the mapping to regular lat/lon
 file_ll = '%s/workspace/hsu_mapping_latlon.tif' % cdir
 os.system('rm -f %s' % file_ll)
 res = input['wbd']['bbox']['res']
 minlat = input['wbd']['bbox']['minlat']
 minlon = input['wbd']['bbox']['minlon']
 maxlat = input['wbd']['bbox']['maxlat']
 maxlon = input['wbd']['bbox']['maxlon']
 log = '%s/workspace/log.txt' % cdir
 os.system('gdalwarp -tr %.16f %.16f -dstnodata %.16f -t_srs EPSG:4326 -s_srs EPSG:102039 -te %.16f %.16f %.16f %.16f %s %s >> %s 2>&1' % (res,res,metadata['nodata'],minlon,minlat,maxlon,maxlat,file_ca,file_ll,log))
 #Add the lat/lon mapping
 #Retrieve the lat/lon metadata
 metadata = gdal_tools.retrieve_metadata(file_ll)
 metadata['nodata'] = -9999.0
 #Save the lat/lon metadata
 fp.createDimension('nlon',metadata['nx'])
 fp.createDimension('nlat',metadata['ny'])
 hmll = grp.createVariable('hmll','f4',('nlat','nlon'))
 hmll.gt = metadata['gt']
 hmll.projection = metadata['projection']
 hmll.description = 'HSU mapping (regular lat/lon)'
 hmll.nodata = metadata['nodata']
 #Save the lat/lon mapping
 hsu_map = np.copy(gdal_tools.read_raster(file_ll))
 hsu_map[np.isnan(hsu_map) == 1] = metadata['nodata']
 hmll[:] = hsu_map
 
 #Close the file 
 fp.close()

 return

def update_netcdf(cdir,iens,nens,parameters):
 
 #Define the netcdf file
 file_netcdf = '%s/output/output.nc' % cdir

 #Read in the ensemble output
 file_pickle = '%s/output/output.pck' % cdir
 output = pickle.load(open(file_pickle))

 #Read in the model input data
 file_pickle = '%s/input/data.pck' % cdir
 input = pickle.load(open(file_pickle))

 #Create the netcdf file if necessary
 if iens == 0: create_netcdf_file(file_netcdf,output,nens,input,cdir)

 #Open the file
 fp = nc.Dataset(file_netcdf, 'a')
 
 #Output the data
 for var in output['variables']:
  
  #Compute the daily average
  data = []
  for itime in xrange(output['variables'][var].shape[0]/24):
   data.append(np.mean(output['variables'][var][24*itime:24*(itime+1)],axis=0))
  data = np.array(data)

  #Write the catchment info
  fp.groups['catchment'].variables[var][:,:,iens] = data

  #Compute the catchment summary 
  if var not in ['qout_surface','qout_subsurface']:
   data = np.sum(output['misc']['pct']*data,axis=1)
  else:
   imax = output['misc']['outlet_hsu']
   data = data[:,imax]

  #Write the catchment summary
  fp.groups['summary'].variables[var][:,iens] = data

 #Close the file
 fp.close()

 return
