import os
import numpy as np
import netCDF4 as nc
import datetime
import time
import sys
sys.path.append('Model/pyNoahMP')
sys.path.append('Model/pyDTopmodel')
import scipy.sparse as sparse
import pickle

def Finalize_Model(NOAH,TOPMODEL):

 #Deallocate the NOAH memory
 NOAH.finalize()

 #Delete the objects
 del NOAH
 del TOPMODEL

 return

def Initialize_Model(ncells,dt,nsoil,parameters,info,wbd):

 from NoahMP import model

 model.ncells = ncells
 model.nsoil = nsoil
 model.dt = dt
 model.nsnow = 3
 model.llanduse[:] = 'MODIFIED_IGBP_MODIS_NOAH'
 #model.lsoil[:] = 'STAS'
 model.lsoil[:] = 'CUST'
 model.vegparm_file[:] = info['VEGPARM']#'data/VEGPARM.TBL'
 model.genparm_file[:] = info['GENPARM']#'data/GENPARM.TBL'
 model.soilparm_file[:] = info['SOILPARM']
 model.mptable_file[:] = info['MPTABLE']#'pyNoahMP/data/MPTABLE.TBL'
 #Read in the soil parameter file
 #fp = open(info['SOILPARM'])
 #iline = 0
 #soils_data = {'MAXSMC':[],'DRYSMC':[],'REFSMC':[]}
 #for line in fp:
 # if (iline > 2) & (iline < 15):
 #  tmp = line.split(',')
 #  soils_data['MAXSMC'].append(float(tmp[4]))
 #  soils_data['DRYSMC'].append(float(tmp[2]))
 #  soils_data['REFSMC'].append(float(tmp[5]))
 # iline = iline + 1
 #fp.close()
 #Define the options
 model.idveg = 3#3#4 # dynamic vegetation (1 -> off ; 2 -> on)
 model.iopt_crs = 2#2 # canopy stomatal resistance (1-> Ball-Berry; 2->Jarvis)
 model.iopt_btr = 1#2 # soil moisture factor for stomatal resistance (1-> Noah; 2-> CLM; 3-> SSiB)
 model.iopt_run = 5#5#2#5##2#2 # runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS)
 model.iopt_sfc = 2#2 # surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
 model.iopt_frz = 2#1 # supercooled liquid water (1-> NY06; 2->Koren99)
 model.iopt_inf = 2#1#1#2 # frozen soil permeability (1-> NY06; 2->Koren99)
 model.iopt_rad = 2 # radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
 model.iopt_alb = 1 # snow surface albedo (1->BATS; 2->CLASS)
 model.iopt_snf = 3 # rainfall & snowfall (1-Jordan91; 2->BATS; 3->Noah)]
 model.iopt_tbot = 1#1#1 # lower boundary of soil temperature (1->zero-flux; 2->Noah) 
 model.iopt_stc = 1#1#2 # snow/soil temperature time scheme (only layer 1) 1 -> semi-implicit; 2 -> full implicit (original Noah)
 #Initialize the model
 model.initialize()
 #Set initial info
 model.iz0tlnd = 0
 model.z_ml[:] = 2.0
 #psoil = 10**parameters['log10soil']
 #psoil = parameters['psoil']
 tmp = 0.1*np.ones(nsoil)
 model.sldpth[:] = tmp#psoil*tmp
 model.zsoil[:] = -np.cumsum(model.sldpth[:],axis=1)
 model.zsnso[:] = 0.0
 model.zsnso[:,3::] = model.zsoil[:]

 #Set all to land (can include lake in the future...)
 model.ist[:] = 1
 #Set soil color
 model.isc[:] = 4#4
 model.ice[:] = 0
 #Initialize snow info
 model.isnow[:] = 0
 #Others
 model.dx = 30.0
 model.ice[:] = 0
 model.foln[:] = 1.0
 model.ficeold[:] = 0.0
 model.albold[:] = 0.65
 model.sneqvo[:] = 0.0
 model.ch[:] = 0.0001
 model.cm[:] = 0.0001
 model.canliq[:] =  0.0
 model.canice[:] = 0.0
 model.sndpth[:] = 0.0
 model.swe[:] = 0.0
 model.snice[:] = 0.0
 model.snliq[:] = 0.0
 model.wa[:] = 0.0
 model.wt[:] = model.wa[:]
 model.zwt[:] = 0.0
 model.wslake[:] = 0.0
 model.lfmass[:] = 9.0
 model.rtmass[:] = 500.0
 model.stmass[:] = 3.33
 model.wood[:] = 500.0
 model.stblcp[:] = 1000.0
 model.fastcp[:] = 1000.0
 model.plai[:] = 0.5
 model.tauss[:] = 0.0
 model.smcwtd[:] = model.sh2o[:,1]
 model.deeprech[:] = 0.0
 model.rech[:] = 0.0
 model.eah[:] = 1000.0
 model.fwet[:] = 0.0
 model.psai[:] = 0.1
 model.stc[:] = 285.0
 model.slopetyp[:] = 3
 model.albold[:] = 0.5
 #Define the data
 model.vegtyp[:] = info['input_fp'].groups['parameters'].variables['land_cover'][:]
 model.soiltyp[:] = np.arange(1,ncells+1)
 #info['input_fp'].groups['parameters'].variables['soil_texture_class'][:]
 model.smcmax[:] = info['input_fp'].groups['parameters'].variables['MAXSMC'][:]
 model.smcref[:] = info['input_fp'].groups['parameters'].variables['REFSMC'][:]
 model.smcdry[:] = info['input_fp'].groups['parameters'].variables['DRYSMC'][:]
 for ilayer in xrange(model.sh2o.shape[1]):
  model.sh2o[:,ilayer] = info['input_fp'].groups['parameters'].variables['MAXSMC'][:]
 '''for hsu in data['hsu']:
  ihsu = data['hsu'].keys().index(hsu)
  model.vegtyp[ihsu] = data['hsu'][hsu]['land_cover'] #HERE
  model.soiltyp[ihsu] = ihsu + 1 #data['hsu'][hsu]['soil_texture_class']#1#ihsu+1 #HERE
  model.smcmax[ihsu] = data['hsu'][hsu]['soil_parameters']['MAXSMC'] #soils_data['MAXSMC'][model.soiltyp[ihsu]-1]
  model.sh2o[ihsu] = data['hsu'][hsu]['soil_parameters']['MAXSMC'] #soils_data['MAXSMC'][model.soiltyp[ihsu]-1]
  model.smcref[ihsu] = data['hsu'][hsu]['soil_parameters']['REFSMC'] #soils_data['REFSMC'][model.soiltyp[ihsu]-1]
  model.smcdry[ihsu] = data['hsu'][hsu]['soil_parameters']['DRYSMC'] #soils_data['DRYSMC'][model.soiltyp[ihsu]-1]
 '''
 model.smc[:] = model.sh2o[:]
 model.smcwtd[:] = model.sh2o[:,0]
 #Set lat/lon (declination calculation)
 model.lat[:] = 0.0174532925*(wbd['bbox']['minlat'] + wbd['bbox']['maxlat'])/2
 model.lon[:] = 0.0174532925*(360.0+(wbd['bbox']['minlon'] + wbd['bbox']['maxlon'])/2)

 #Initialize output
 model.tg[:] = 285.0
 model.tv[:] = 285.0
 for ilayer in xrange(model.nsnow-1,model.stc.shape[1]):
  model.stc[:,ilayer] = 285.0
 model.tah[:] = 285.0
 model.t2mv[:] = 285.0
 model.t2mb[:] = 285.0
 model.runsf[:] = 0.0
 model.runsb[:] = 0.0
 #forcing
 model.fveg[:] = 1.0
 model.fvgmax[:] = 1.0
 model.tbot[:] = 285.0

 return model

def Initialize_DTopmodel(ncells,dt,parameters,info):

 #sys.path.append("topmodel")
 import dynamic_topmodel

 #Define some metadata
 nhru_outlet = info['input_fp'].groups['outlet'].groups['summary'].variables['hru_dst'].size

 #Initialize Dynamic Topmodel
 model = dynamic_topmodel.Dynamic_Topmodel(ncells,nhru_outlet)
 dx = 30.0 #meters

 #Set model parameters
 model.dt = dt #seconds
 model.area[:] = 0.0 #meters^2
 model.dx[:] = dx #meters
 model.m[:] = 10**parameters['log10m'] #HERE
 model.sdmax[:] = parameters['sdmax'] #HERE
 #model.sdmax[:] = 1000.0#parameters['sdmax'] #HERE
 dem = []
 #Set cluster information
 model.pct[:] = info['input_fp'].groups['parameters'].variables['area_pct'][:]/100
 model.area[:] = info['input_fp'].groups['parameters'].variables['area'][:]
 model.T0[:] = info['input_fp'].groups['parameters'].variables['SATDK'][:]*(10**parameters['log10m'])
 model.sti[:] = info['input_fp'].groups['parameters'].variables['ti'][:]
 model.beta[:] = info['input_fp'].groups['parameters'].variables['slope'][:]
 model.carea[:] = info['input_fp'].groups['parameters'].variables['carea'][:]
 model.channel[:] = info['input_fp'].groups['parameters'].variables['channel'][:]
 model.surface_velocity[:] = 1000.0/dt #m/s
 model.dem[:] = info['input_fp'].groups['parameters'].variables['dem'][:] 
 #Set outlet information
 model.area_outlet[:] = dx**2*info['input_fp'].groups['outlet'].groups['summary'].variables['counts'][:]
 #model.sdmax[:] = 0.1#np.array(dem) - np.min(dem)
 model.pct = model.pct/np.sum(model.pct)
 ti_mean = np.sum(model.pct*model.sti[:])
 #lnTe = parameters['lnTe']
 
 #Calculate the sti
 #model.T0[:] = (model.dem - np.min(model.dem))*model.T0
 lnT0 = np.log(model.T0)
 #lnT0 = lnTe*lnT0/np.sum(model.pct*lnT0)
 #model.T0[:] = np.exp(lnT0)
 lnTe = np.sum(model.pct*lnT0)
 model.sti = model.sti - (lnT0 - lnTe)

 #Set weight matrix
 #model.outlet_hsu = int(info['input_fp'].groups['metadata'].outlet_hsu)#int(data['outlet']['hsu'])
 #hsu = model.outlet_hsu
 flow_matrix = sparse.csr_matrix((info['input_fp'].groups['flow_matrix'].variables['data'][:],
			          info['input_fp'].groups['flow_matrix'].variables['indices'][:],
			          info['input_fp'].groups['flow_matrix'].variables['indptr'][:]),
 				  dtype=np.float32)
 #flow_matrix.data = info['input_fp'].groups['flow_matrix'].variables['data'][:]
 #flow_matrix.indices = info['input_fp'].groups['flow_matrix'].variables['indices'][:]
 #flow_matrix.indptr = info['input_fp'].groups['flow_matrix'].variables['indptr'][:]
 #data['tp'][hsu,hsu] = data['tp'][model.outlet_hsu,hsu] - model.dx[hsu]**2/model.area[hsu]
 #flow_matrix[hsu,hsu] = flow_matrix[hsu,hsu] - model.dx[hsu]**2/model.area[hsu]
 #model.w = sparse.csr_matrix(data['tp'].T,dtype=np.float32)
 model.w = flow_matrix
 #print model.w.shape
 #exit()
 #model.wfull[:] = data['tp'].T
 #Initialize the soil moisture deficit values
 model.si[:] = 0.0
 
 return model

def Update_Model(NOAH,TOPMODEL,ncores):
 
 #Set the partial pressure of CO2 and O2
 NOAH.co2air[:] = 355.E-6*NOAH.psfc[:]# ! Partial pressure of CO2 (Pa) ! From NOAH-MP-WRF
 NOAH.o2air[:] = 0.209*NOAH.psfc[:]# ! Partial pressure of O2 (Pa)  ! From NOAH-MP-WRF

 #Update NOAH
 NOAH.minzwt[:] = -100.0#-0.1*((TOPMODEL.dem - np.min(TOPMODEL.dem))+0.5)
 ntt = 1#4
 #NOAH.dzwt[:] = NOAH.dzwt[:]/ntt
 NOAH.dzwt[:] = NOAH.dzwt[:]/ntt
 dt = np.copy(NOAH.dt)
 NOAH.dt = dt/ntt
 #print "Update Noah",ncores
 for i in xrange(ntt):
  NOAH.run_model(ncores)
 NOAH.dt = np.copy(dt)

 #Calculate the updated soil moisture deficit
 si0 = np.copy(NOAH.si0)
 si1 = np.copy(NOAH.si1)

 #Calculate the change in deficit
 TOPMODEL.si[:] = si1[:]
 #TOPMODEL.si[:] = NOAH.si1[:]
 TOPMODEL.dsi[:] = np.copy(si1 - si0)
 #TOPMODEL.dsi[:] = NOAH.si1[:] - NOAH.si0[:]
 TOPMODEL.r[:] = -TOPMODEL.dsi[:]/TOPMODEL.dt

 #Add the surface runoff
 TOPMODEL.qsurf[:] = NOAH.runsf[:]/1000.0

 #Update dynamic topmodel
 #print "Update TOPMODEL"
 TOPMODEL.update(ncores)

 #Calculate the change in deficit
 TOPMODEL.sideep[:] = TOPMODEL.sideep.astype(np.float32)
 TOPMODEL.si[:] = TOPMODEL.si.astype(np.float32)
 dsi = np.copy(si1 - TOPMODEL.si)

 #Update the soil moisture values
 NOAH.dzwt[:] = np.copy(dsi+TOPMODEL.dt*TOPMODEL.ex-TOPMODEL.dt*TOPMODEL.r)
 TOPMODEL.ex[:] = 0.0

 return (NOAH,TOPMODEL)

def run_model(info):

 #data = pickle.load(open(info['input']))
 #Define the metadata
 dt = info['dt']
 nsoil = info['nsoil']
 wbd = info['wbd']
 ncores = info['ncores']
 idate = info['idate']
 fdate = info['fdate']
 parameters = info['parameters']
 output_type = info['output_type']
 soil_file = info['soil_file']
 input_file = info['input']
 output_file = info['output']
 dt_timedelta = datetime.timedelta(seconds=dt)

 #Open access to the input netcdf file
 info = {}
 info['input_fp'] = nc.Dataset(input_file)

 #Open access to the output netcdf file
 info['output_fp'] = nc.Dataset(output_file,'w',format='NETCDF4')
  
 #Define the parameter files
 info['VEGPARM'] = 'Model/pyNoahMP/data/VEGPARM.TBL'#'data/VEGPARM.TBL'
 info['GENPARM'] = 'Model/pyNoahMP/data/GENPARM.TBL'#'data/GENPARM.TBL'
 info['MPTABLE'] = 'Model/pyNoahMP/data/MPTABLE.TBL'#'pyNoahMP/data/MPTABLE.TBL'
 info['SOILPARM'] = soil_file
 
 #Set the number of cells (hsus)
 ncells = len(info['input_fp'].dimensions['hsu'])

 #Initialize the output
 output = Initialize_Output(output_type)

 #Initialize the model
 print "Initializing Noah"
 NOAH = Initialize_Model(ncells,dt,nsoil,parameters,info,wbd)
 output['misc']['zsoil'] = NOAH.zsoil

 #Initialize Topmodel
 print "Initializing TOPMODEL"
 TOPMODEL = Initialize_DTopmodel(ncells,dt,parameters,info) 

 #Run the model
 #meteorology = {}
 #for var in info['input_fp'].groups['meteorology'].variables:
 # meteorology[var] = info['input_fp'].groups['meteorology'].variables[var][:]
 #meteorology = data['meteorology']
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
   print date
   testvar = 1
   #print date,np.sum(TOPMODEL.pct*TOPMODEL.si),time.time() - tic,et,prcp,'q:%f'%q,np.sum(TOPMODEL.pct*NOAH.smcwtd),np.sum(TOPMODEL.pct*NOAH.zwt),'WB ERR:%f' % np.sum(TOPMODEL.pct*errwat),'etran:%f'%etran,'ecan:%f'%ecan,'esoil:%f'%esoil
  if (date.month == 1) and (date.day == 1) and (date.hour == 0):
   et = 0
   prcp = 0
   q = 0
   ecan,etran,esoil = 0,0,0
 
  #Update input data
  #Update_Input(NOAH,TOPMODEL,date,meteorology,i)
  Update_Input(NOAH,TOPMODEL,date,info,i)

  #Calculate initial NOAH water balance
  #smw = 1000.0*((NOAH.smcmax*np.abs(NOAH.zwt) - TOPMODEL.si) + (NOAH.smcmax*(100-np.abs(NOAH.zwt))))
  #beg_wb = np.copy(NOAH.canliq + NOAH.canice + NOAH.swe + NOAH.wa + smw)#TOPMODEL.si*1000)
  #dzwt0 = np.copy(NOAH.dzwt)# - NOAH.deeprech)

  #Update model
  (NOAH,TOPMODEL) = Update_Model(NOAH,TOPMODEL,ncores)

  #Calculate final water balance
  #smw = 1000.0*((NOAH.smcmax*np.abs(NOAH.zwt) - TOPMODEL.si) + (NOAH.smcmax*(100.0-np.abs(NOAH.zwt))))
  #end_wb = np.copy(NOAH.canliq + NOAH.canice + NOAH.swe + NOAH.wa + smw)# + TOPMODEL.si*1000)

  #Update the water balance error
  #tmp = np.copy(end_wb - beg_wb - NOAH.dt*(NOAH.prcp-NOAH.ecan-NOAH.etran-NOAH.esoil-NOAH.runsf-NOAH.runsb) - 1000*dzwt0)
  #errwat += tmp# - tmp
  #q = q + dt*np.sum(TOPMODEL.pct*NOAH.runsb)
  #et = et + dt*np.sum(TOPMODEL.pct*(NOAH.ecan + NOAH.etran + NOAH.esoil))
  #etran += dt*np.sum(TOPMODEL.pct*NOAH.etran)
  #ecan += dt*np.sum(TOPMODEL.pct*NOAH.ecan)
  #esoil += dt*np.sum(TOPMODEL.pct*NOAH.esoil)
  #prcp = prcp + dt*np.sum(TOPMODEL.pct*NOAH.prcp)

  #Update output
  #Write to netcdf
  #file_netcdf = hydrobloks_info['output'] 
  ##Update the netcdf file
  #vars = output['variables'].keys()
  #update_netcdf(hydrobloks_info['dir'],0,1,parameters,file_netcdf,output,0)
  Update_Output(info,i,NOAH,TOPMODEL)
  #output = Update_Output(output,NOAH,TOPMODEL,date,output_type)

  #Update time step
  date = date + dt_timedelta
  i = i + 1

 #Finalize the model
 Finalize_Model(NOAH,TOPMODEL)

 #Close the input and output files
 info['input_fp'].close()
 info['output_fp'].close()

 return

def Update_Input(NOAH,TOPMODEL,date,info,i):

  NOAH.itime = i
  dt = NOAH.dt
  TOPMODEL.itime = i
  NOAH.nowdate[:] = date.strftime('%Y-%m-%d_%H:%M:%S')
  NOAH.julian = (date - datetime.datetime(date.year,1,1,0)).days
  NOAH.yearlen = (datetime.datetime(date.year+1,1,1,0) - datetime.datetime(date.year,1,1,1,0)).days + 1

  #Update meteorology
  meteorology = info['input_fp'].groups['meteorology']
  NOAH.lwdn[:] = meteorology.variables['dlwrf'][i,:] #W/m2
  NOAH.swdn[:] = meteorology.variables['dswrf'][i,:] #W/m2
  NOAH.psfc[:] = meteorology.variables['pres'][i,:] #Pa
  NOAH.p_ml[:] = meteorology.variables['pres'][i,:] #Pa
  NOAH.u_ml[:] = (meteorology.variables['wind'][i,:]**2/2)**0.5 #m/s
  NOAH.v_ml[:] = (meteorology.variables['wind'][i,:]**2/2)**0.5 #m/s
  NOAH.t_ml[:] = 273.15+meteorology.variables['tair'][i,:] #K
  estar = 611.0*np.exp(17.27*((NOAH.t_ml[:] - 273.15)/(NOAH.t_ml[:] - 36.0)))
  e = meteorology.variables['rh'][i,:]*estar/100
  q = 0.622*e/NOAH.psfc[:]
  NOAH.q_ml[:] = q #Kg/Kg
  NOAH.qsfc1d[:] = q #Kg/Kg
  NOAH.prcp[:] = meteorology.variables['prec'][i,:]/dt #mm/s ONLY NLDAS
  #if np.mean(meteorology.variables['apcpsfc'][i,:]) >= 0.0:
  # NOAH.prcp[:] = meteorology.variables['apcpsfc'][i,:]/dt #mm/s
  #else:
  # NOAH.prcp[:] = meteorology.variables['prec'][i,:]/dt #mm/s

  return (NOAH,TOPMODEL)

def Initialize_Output(output_type):

 if output_type == 'Full':
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
 elif output_type == 'Summary':
  variables = {
           #State variables
           'smc1':{'mean':[],'std':[]}, #soil moisture content (layer 1)
           'swe':{'mean':[],'std':[]}, #soil moisture content (layer 1)
           'lh':{'mean':[],'std':[]}, #soil moisture content (layer 1)
           'sh':{'mean':[],'std':[]}, #soil moisture content (layer 1)
           'qexcess':{'mean':[],'std':[]}, #soil moisture content (layer 1)
           'qsurface':{'mean':[],'std':[]}, #soil moisture content (layer 1)
           'prcp':{'mean':[],'std':[]}, #soil moisture content (layer 1)
           }

 misc = {
      'sldpth':[], #soil depth
      'dates':[], #dates
      }

 return {'variables':variables,'misc':misc}

def Create_Netcdf_File(info):

 #Extract the pointers to both input and output files
 fp_out = info['output_fp']
 fp_in = info['input_fp']

 #Define the metadata
 metadata = {
             'g':{'description':'Ground heat flux','units':'W/m2'},
             'sh':{'description':'Sensible heat flux','units':'W/m2'},
             'lh':{'description':'Latent heat flux','units':'W/m2'},
             #'lwnet':{'description':'Net longwave radiation','units':'W/m2'},
             #'swnet':{'description':'Absorbed shortwave radiation','units':'W/m2'},
             #'et':{'description':'Evapotranspiration','units':'mm/s'},
             #'qexcess':{'description':'Excess runoff','units':'mm/s'},
             #'qsurface':{'description':'Surface runoff','units':'mm/s'},
             #'prcp':{'description':'Precipitation','units':'mm/s'},
             'smc1':{'description':'Soil water content','units':'m3/m3'},
             #'smc2':{'description':'Soil water content','units':'m3/m3'},
             #'smc3':{'description':'Soil water content','units':'m3/m3'},
             #'smc4':{'description':'Soil water content','units':'m3/m3'},
             'swd':{'description':'Soil water deficit','units':'mm'},
             #'sstorage':{'description':'Surface storage','units':'mm'},
             #'sndpth':{'description':'Snow depth','units':'mm'},
             #'swe':{'description':'Snow water equivalent','units':'mm'},
             'qout_subsurface':{'description':'Subsurface flux','units':'m2/s'},
             #'qout_surface':{'description':'Surface flux','units':'m2/s'},
             }

 #Create the dimensions
 print 'Creating the dimensions'
 ntime = len(fp_in.dimensions['time'])
 nhsu = len(fp_in.dimensions['hsu'])
 nsoil = 4
 fp_out.createDimension('hsu',nhsu)
 fp_out.createDimension('time',ntime)
 fp_out.createDimension('soil',nsoil)

 #Create the output
 print 'Creating the catchment group'
 grp = fp_out.createGroup('catchment')
 for var in metadata:
  ncvar = grp.createVariable(var,'f4',('time','hsu',))
  ncvar.description = metadata[var]['description']
  ncvar.units = metadata[var]['units']

 #Create the outlet group 
 print 'Creating the outlet group'
 grp = fp_out.createGroup('outlet')
 grp.createDimension('hru_outlet',len(fp_in.groups['outlet'].groups['summary'].dimensions['hru']))
 grp.createVariable('hru_org','i4',('hru_outlet',))
 grp.createVariable('hru_dst','i4',('hru_outlet',))
 grp.createVariable('counts','i4',('hru_outlet',))
 for var in ['qin_subsurface_outlet',]:
  ncvar = grp.createVariable(var,'f4',('time','hru_outlet',))

 #Create the metadata
 print 'Creating the metadata group'
 grp = fp_out.createGroup('metadata')
 #dates
 '''times = grp.createVariable('dates','f8',('time',))
 dates = []
 for date in output['misc']['dates']:
  #if date.hour == 0:dates.append(date)
  dates.append(date)
 times.units = 'days since 1900-01-01'
 times.calendar = 'standard'
 times[:] = nc.date2num(np.array(dates),units=times.units,calendar=times.calendar)'''
 #HSU percentage coverage
 print 'Setting the HRU percentage coverage'
 pcts = grp.createVariable('pct','f4',('hsu',))
 pcts[:] = fp_in.groups['parameters'].variables['area_pct'][:]
 pcts.description = 'hsu percentage coverage'
 pcts.units = '%'
 '''#HSU Spatial resolution
 print 'Setting the spatial resolution'
 dx = grp.createVariable('dx','f4',('hsu',))
 dx[:] = np.array(output['misc']['dx'])
 dx.units = 'meters'''
 #HSU area
 print 'Setting the HRU areal coverage'
 area = grp.createVariable('area','f4',('hsu',))
 area[:] = fp_in.groups['parameters'].variables['area'][:]
 area.units = 'meters squared'
 print 'Defining the HRU ids'
 hsu = grp.createVariable('hsu','i4',('hsu',))
 hsus =[]
 for value in xrange(nhsu):hsus.append(value)
 hsu[:] = np.array(hsus)
 hsu.description = 'hsu ids'

 #Create the mapping
 grp = fp_out.createGroup('latlon_mapping')
 grp.createDimension('nlon',len(fp_in.groups['latlon_mapping'].dimensions['nlon']))
 grp.createDimension('nlat',len(fp_in.groups['latlon_mapping'].dimensions['nlat']))
 hmll = grp.createVariable('hmll','f4',('nlat','nlon'))
 hmll.gt = fp_in.groups['latlon_mapping'].variables['hmll'].gt
 hmll.projection = fp_in.groups['latlon_mapping'].variables['hmll'].projection
 hmll.description = 'HSU mapping (regular lat/lon)'
 hmll.nodata = fp_in.groups['latlon_mapping'].variables['hmll'].nodata
 #Save the lat/lon mapping
 hmll[:] = fp_in.groups['latlon_mapping'].variables['hmll'][:]

 return

def Update_Output(info,itime,NOAH,TOPMODEL):

 #Create the netcdf file
 if itime == 0: Create_Netcdf_File(info)

 #Update the variables (catchment)
 grp =info['output_fp'].groups['catchment']
 #NoahMP
 grp.variables['smc1'][itime,:] = np.copy(NOAH.smc[:,0]) #m3/m3
 grp.variables['g'][itime,:] = np.copy(NOAH.ssoil) #W/m2
 grp.variables['sh'][itime,:] = np.copy(NOAH.fsh) #W/m2
 grp.variables['lh'][itime,:] = np.copy(NOAH.fcev + NOAH.fgev + NOAH.fctr) #W/m2
 #TOPMODEL
 grp.variables['swd'][itime,:] = np.copy(10**3*TOPMODEL.si) #mm
 grp.variables['qout_subsurface'][itime,:] = np.copy(TOPMODEL.qout) #m2/s

 #Update the variables (outlet)
 grp = info['output_fp'].groups['outlet']
 grp.variables['qin_subsurface_outlet'][itime,:] = np.copy(TOPMODEL.qin_outlet) #m2/s
 
 return

def update_output(output,NOAH,TOPMODEL,date,output_type):
 
 #Miscellanous
 output['misc']['sldpth'] = np.copy(NOAH.sldpth)
 output['misc']['dates'].append(date)
 output['misc']['pct'] = TOPMODEL.pct
 output['misc']['dx'] = TOPMODEL.dx
 output['misc']['area'] = TOPMODEL.area
 #output['misc']['outlet_hsu'] = TOPMODEL.outlet_hsu
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

 if output_type == 'Full':
  #Energy fluxes
  output['variables']['g'].append(np.copy(NOAH.ssoil)) #W/m2
  output['variables']['sh'].append(np.copy(NOAH.fsh)) #W/m2
  output['variables']['lh'].append(np.copy(NOAH.fcev + NOAH.fgev + NOAH.fctr)) #W/m2
  output['variables']['swnet'].append(np.copy(NOAH.fsa)) #W/m2
  output['variables']['lwnet'].append(np.copy(NOAH.fira)) #W/m2
  #Water balance
  output['variables']['et'].append(np.copy(NOAH.ecan + NOAH.etran + NOAH.esoil)) #mm/s
  #imax = TOPMODEL.outlet_hsu
  #qexcess = np.copy(10**3*TOPMODEL.ex)
  qexcess = np.copy(NOAH.runsb)
  #qexcess[imax] = qexcess[imax] + 10**3*(TOPMODEL.dx[imax]*(TOPMODEL.qout[imax])/np.sum(TOPMODEL.area))
  output['variables']['qexcess'].append(np.copy(qexcess)) #mm/s
  output['variables']['qsurface'].append(np.copy(NOAH.runsf)) #mm/s
  output['variables']['prcp'].append(np.copy(NOAH.prcp)) #mm/s
  #State variables
  output['variables']['smc1'].append(np.copy(NOAH.smc[:,0])) #m3/m3
  output['variables']['smc2'].append(np.copy(NOAH.smc[:,1])) #m3/m3
  output['variables']['smc3'].append(np.copy(NOAH.smc[:,2])) #m3/m3
  output['variables']['smc4'].append(np.copy(NOAH.smc[:,3])) #m3/m3
  #output['variables']['smc5'].append(np.copy(NOAH.smc[:,4])) #m3/m3
  output['variables']['swd'].append(np.copy(10**3*TOPMODEL.si)) #mm
  output['variables']['sstorage'].append(np.copy(10**3*TOPMODEL.storage_surface)) #mm
  output['variables']['sndpth'].append(np.copy(10**3*NOAH.sndpth)) #mm
  output['variables']['swe'].append(np.copy(NOAH.swe)) #mm
  #Fluxes
  output['variables']['qout_subsurface'].append(np.copy(TOPMODEL.qout)) #m2/s
  output['variables']['qout_surface'].append(np.copy(TOPMODEL.qout_surface)) #m2/s

 elif output_type == 'Summary':

  #Compute the mean and standard deviation
  pcts = TOPMODEL.pct
  #smc1
  mean = np.sum(pcts*NOAH.smc[:,0])
  std = np.sum(pcts*(NOAH.smc[:,0] - mean)**2)**0.5
  output['variables']['smc1']['mean'].append(mean)
  output['variables']['smc1']['std'].append(std)
  #swe
  swe = NOAH.swe
  mean = np.sum(pcts*swe)
  std = np.sum(pcts*(swe - mean)**2)**0.5
  output['variables']['swe']['mean'].append(mean)
  output['variables']['swe']['std'].append(std)
  #lh
  lh = NOAH.fcev + NOAH.fgev + NOAH.fctr
  mean = np.sum(pcts*lh)
  std = np.sum(pcts*(lh - mean)**2)**0.5
  output['variables']['lh']['mean'].append(mean)
  output['variables']['lh']['std'].append(std)
  #sh
  sh = NOAH.fsh
  mean = np.sum(pcts*sh)
  std = np.sum(pcts*(sh - mean)**2)**0.5
  output['variables']['sh']['mean'].append(mean)
  output['variables']['sh']['std'].append(std)
  #qexcess
  #imax = TOPMODEL.outlet_hsu
  qexcess = np.copy(NOAH.runsb)
  #qexcess[imax] = qexcess[imax] + 10**3*(TOPMODEL.dx[imax]*(TOPMODEL.qout[imax])/np.sum(TOPMODEL.area))
  mean = np.sum(pcts*qexcess)
  std = np.sum(pcts*(qexcess - mean)**2)**0.5
  output['variables']['qexcess']['mean'].append(mean)
  output['variables']['qexcess']['std'].append(std)
  #qsurface
  qsurf = NOAH.runsf
  mean = np.sum(pcts*qsurf)
  std = np.sum(pcts*(qsurf - mean)**2)**0.5
  output['variables']['qsurface']['mean'].append(mean)
  output['variables']['qsurface']['std'].append(std)
  #prcp
  prcp = NOAH.prcp
  mean = np.sum(pcts*prcp)
  std = np.sum(pcts*(prcp - mean)**2)**0.5
  output['variables']['prcp']['mean'].append(mean)
  output['variables']['prcp']['std'].append(std)

 return output

