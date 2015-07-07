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

 #Close the TOPMODEL solver
 TOPMODEL.dtt.finalize()

 #Delete the objects
 del NOAH
 del TOPMODEL

 return

class HydroBloks:

 def __init__(self,nhru):

  #NOAH-MP
  self.smc1 = np.zeros(nhru,dtype=np.float32)
  self.g = np.zeros(nhru,dtype=np.float32)
  self.sh = np.zeros(nhru,dtype=np.float32)
  self.lh = np.zeros(nhru,dtype=np.float32)
  self.qbase = np.zeros(nhru,dtype=np.float32)
  self.qsurface = np.zeros(nhru,dtype=np.float32)
  self.prcp = np.zeros(nhru,dtype=np.float32)

  #TOPMODEL
  self.swd = np.zeros(nhru,dtype=np.float32)
  self.qout_subsurface = np.zeros(nhru,dtype=np.float32)
  self.qout_surface = np.zeros(nhru,dtype=np.float32)
  self.qout_sstorage = np.zeros(nhru,dtype=np.float32)

  #Summary
  self.dE = 0.0
  self.r = 0.0
  self.dr = 0.0
  self.et = 0.0
  self.etran = 0.0
  self.esoil = 0.0
  self.ecan = 0.0
  self.prcp = 0.0
  self.q = 0.0
  self.errwat = 0.0#np.zeros(nhru,dtype=np.float32)
  self.dzwt0 = np.zeros(nhru,dtype=np.float32)
  self.beg_wb = np.zeros(nhru,dtype=np.float32)
  self.end_wb = np.zeros(nhru,dtype=np.float32)

  return

 def Initialize_Water_Balance(self,NOAH,TOPMODEL):

  smw = 1000.0*((NOAH.smcmax*np.abs(NOAH.zwt) - TOPMODEL.si) + (NOAH.smcmax*(100-np.abs(NOAH.zwt))))
  self.beg_wb = np.copy(NOAH.canliq + NOAH.canice + NOAH.swe + NOAH.wa + smw)#TOPMODEL.si*1000)
  self.dzwt0 = np.copy(NOAH.dzwt)# - NOAH.deeprech)

  return 

 def Finalize_Water_Balance(self,NOAH,TOPMODEL):

  smw = 1000.0*((NOAH.smcmax*np.abs(NOAH.zwt) - TOPMODEL.si) + (NOAH.smcmax*(100.0-np.abs(NOAH.zwt))))
  self.end_wb = np.copy(NOAH.canliq + NOAH.canice + NOAH.swe + NOAH.wa + smw)# + TOPMODEL.si*1000)

  return

 def Calculate_Water_Balance_Error(self,NOAH,TOPMODEL,dt):

  tmp = np.copy(self.end_wb - self.beg_wb - NOAH.dt*(NOAH.prcp-NOAH.ecan-
        NOAH.etran-NOAH.esoil-NOAH.runsf-NOAH.runsb) - 1000*self.dzwt0)
  self.errwat += np.sum(TOPMODEL.pct*tmp)
  self.q = self.q + dt*np.sum(TOPMODEL.pct*NOAH.runsb) + dt*np.sum(TOPMODEL.pct*NOAH.runsf)
  self.et = self.et + dt*np.sum(TOPMODEL.pct*(NOAH.ecan + NOAH.etran + NOAH.esoil))
  self.etran += dt*np.sum(TOPMODEL.pct*NOAH.etran)
  self.ecan += dt*np.sum(TOPMODEL.pct*NOAH.ecan)
  self.esoil += dt*np.sum(TOPMODEL.pct*NOAH.esoil)
  self.prcp = self.prcp + dt*np.sum(TOPMODEL.pct*NOAH.prcp)

  return

def Initialize_Model(ncells,dt,nsoil,info):

 from pyNoahMP.NoahMP import model
 #from NoahMP import model

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
 tmp = 0.1*np.ones(nsoil)
 model.sldpth[:] = tmp
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
 model.dx = info['dx']
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
 model.lat[:] = 0.0174532925*info['input_fp'].groups['metadata'].latitude
 #(wbd['bbox']['minlat'] + wbd['bbox']['maxlat'])/2
 model.lon[:] = 0.0174532925*info['input_fp'].groups['metadata'].longitude
 #(360.0+(wbd['bbox']['minlon'] + wbd['bbox']['maxlon'])/2)

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

def Initialize_DTopmodel(ncells,dt,info):

 #sys.path.append("topmodel")
 from pyDTopmodel import dynamic_topmodel

 #Define some metadata
 nhru_outlet = info['input_fp'].groups['outlet'].groups['summary'].variables['hru_dst'].size

 #Initialize Dynamic Topmodel
 model = dynamic_topmodel.Dynamic_Topmodel(ncells,nhru_outlet)
 dx = info['dx'] #meters

 #Set flags
 model.surface_flow_flag = info['surface_flow_flag']
 model.subsurface_flow_flag = info['subsurface_flow_flag']

 #Set model parameters
 model.dt = dt #seconds
 model.area[:] = 0.0 #meters^2
 model.dx[:] = dx #meters
 model.m[:] = info['input_fp'].groups['parameters'].variables['m'][:]
 model.sdmax[:] = info['input_fp'].groups['parameters'].variables['sdmax'][:]
 dem = []
 #Set cluster information
 model.pct[:] = info['input_fp'].groups['parameters'].variables['area_pct'][:]/100
 model.area[:] = info['input_fp'].groups['parameters'].variables['area'][:]
 model.T0[:] = info['input_fp'].groups['parameters'].variables['SATDK'][:]*model.m
 model.sti[:] = info['input_fp'].groups['parameters'].variables['ti'][:]
 model.beta[:] = info['input_fp'].groups['parameters'].variables['slope'][:]
 model.carea[:] = info['input_fp'].groups['parameters'].variables['carea'][:]
 model.channel[:] = info['input_fp'].groups['parameters'].variables['channel'][:]
 model.dem[:] = info['input_fp'].groups['parameters'].variables['dem'][:] 
 model.mannings[:] = info['input_fp'].groups['parameters'].variables['mannings'][:]
 #Set outlet informationA
 model.area_outlet[:] = dx**2*info['input_fp'].groups['outlet'].groups['summary'].variables['counts'][:]
 model.pct = model.pct/np.sum(model.pct)
 ti_mean = np.sum(model.pct*model.sti[:])
 
 #Calculate the sti
 lnT0 = np.log(model.T0)
 lnTe = np.sum(model.pct*lnT0)
 model.sti = model.sti - (lnT0 - lnTe)

 #Set weight matrix
 model.flow_matrix = sparse.csr_matrix((info['input_fp'].groups['flow_matrix'].variables['data'][:],
			          info['input_fp'].groups['flow_matrix'].variables['indices'][:],
			          info['input_fp'].groups['flow_matrix'].variables['indptr'][:]),
 				  dtype=np.float64)[0:ncells,0:ncells]
 model.flow_matrix.setdiag(model.flow_matrix.diagonal()) #Ensure the zeros are not sparse (for kinematic wave solution).
 model.flow_matrix_T = sparse.csr_matrix(model.flow_matrix.T) #transposed
 model.flow_matrix_T.setdiag(model.flow_matrix_T.diagonal()) #Ensure the zeros are not sparse  (for kinematic wave solution).

 #Initialize the solver
 model.dtt.initialize(model.flow_matrix_T.indices,model.flow_matrix_T.indptr)

 #Initialize the soil moisture deficit values
 model.si[:] = 0.0
 
 return model

def Update_Model(NOAH,TOPMODEL,ncores):
 
 #Set the partial pressure of CO2 and O2
 NOAH.co2air[:] = 355.E-6*NOAH.psfc[:]# ! Partial pressure of CO2 (Pa) ! From NOAH-MP-WRF
 NOAH.o2air[:] = 0.209*NOAH.psfc[:]# ! Partial pressure of O2 (Pa)  ! From NOAH-MP-WRF

 #Update NOAH
 NOAH.minzwt[:] = -100.0#-0.1*((TOPMODEL.dem - np.min(TOPMODEL.dem))+0.5)
 NOAH.run_model(ncores)
 #print "Update Noah",time.time() - t0

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
 TOPMODEL.qsurf[:] = NOAH.runsf[:]/1000.0# + NOAH.runsb[:]/1000.0

 #Update dynamic topmodel
 #print "Update TOPMODEL"
 t0 = time.time()
 TOPMODEL.update(ncores)
 #print"Update TOPMODEL",time.time() - t0

 #Calculate the change in deficit
 TOPMODEL.sideep[:] = TOPMODEL.sideep.astype(np.float32)
 TOPMODEL.si[:] = TOPMODEL.si.astype(np.float32)
 dsi = np.copy(si1 - TOPMODEL.si)

 #Update the soil moisture values
 NOAH.dzwt[:] = np.copy(dsi+TOPMODEL.dt*TOPMODEL.ex)#-TOPMODEL.dt*TOPMODEL.r)
 #NOAH.dzwt[:] = np.copy(dsi)#*TOPMODEL.ex)#-TOPMODEL.dt*TOPMODEL.r)
 #NOAH.runsf[:] -= 1000.0*TOPMODEL.ex
 #print 1000.0*TOPMODEL.ex
 #TOPMODEL.ex[:] = 0.0

 return (NOAH,TOPMODEL)

def Run_Model(info):

 #data = pickle.load(open(info['input']))
 #Define the metadata
 dt = info['dt']
 dtt = info['dtt']
 nsoil = info['nsoil']
 ncores = info['ncores']
 idate = info['idate']
 fdate = info['fdate']
 soil_file = info['soil_file']
 input_file = info['input_file']
 output_file = info['output_file']
 dt_timedelta = datetime.timedelta(seconds=dt)

 #Open access to the input netcdf file
 #info = {'dt':info['dt'],'dx':info['dx']}
 info['input_fp'] = nc.Dataset(input_file)

 #Open access to the output netcdf file
 info['output_fp'] = nc.Dataset(output_file,'w',format='NETCDF4')
  
 #Define the parameter files
 dir = os.path.dirname(os.path.abspath(__file__))
 info['VEGPARM'] = '%s/pyNoahMP/data/VEGPARM.TBL' % dir#'data/VEGPARM.TBL'
 info['GENPARM'] = '%s/pyNoahMP/data/GENPARM.TBL' % dir#'data/GENPARM.TBL'
 info['MPTABLE'] = '%s/pyNoahMP/data/MPTABLE.TBL' % dir#'pyNoahMP/data/MPTABLE.TBL'
 info['SOILPARM'] = soil_file
 
 #Set the number of cells (hsus)
 ncells = len(info['input_fp'].dimensions['hsu'])

 #Initialize the model
 print "Initializing Noah"
 NOAH = Initialize_Model(ncells,dt,nsoil,info)

 #Initialize Topmodel
 print "Initializing TOPMODEL"
 TOPMODEL = Initialize_DTopmodel(ncells,dt,info) 

 #Initialize the coupler
 HB = HydroBloks(ncells)

 #Run the model
 date = idate
 i = 0
 tic = time.time()
 NOAH.dzwt[:] = 0.0
 TOPMODEL.ex[:] = 0.0
 while date <= fdate:

  #Memorize time stamp
  time0 = time.time()

  if (date.hour == 0) and (date.day == 1):
   testvar = 1
   print date,time.time() - tic,'et:%f'%HB.et,'prcp:%f'%HB.prcp,'q:%f'%HB.q,'WB ERR:%f' % HB.errwat
 
  #Update input data
  Update_Input(NOAH,TOPMODEL,date,info,i)

  #Run sub-timesteps
  ntt = int(dt/dtt)
  NOAH.dt = dtt
  TOPMODEL.dt = dtt

  for itt in xrange(ntt):

   #Calculate initial NOAH water balance
   HB.Initialize_Water_Balance(NOAH,TOPMODEL)

   #Update model
   (NOAH,TOPMODEL) = Update_Model(NOAH,TOPMODEL,ncores)

   #Calculate final water balance
   HB.Finalize_Water_Balance(NOAH,TOPMODEL)

   #Update the water balance error
   HB.Calculate_Water_Balance_Error(NOAH,TOPMODEL,NOAH.dt)

   #Update output

  #Redefine other info
  NOAH.dt = dt
  TOPMODEL.dt = dt

  #Update time and date
  info['time'] = time.time() - time0#seconds
  info['date'] = date

  #Update output
  Update_Output(info,i,NOAH,TOPMODEL)

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
  #NOAH.prcp[:] = meteorology.variables['prec'][i,:]/dt #mm/s ONLY NLDAS
  if np.mean(meteorology.variables['apcpsfc'][i,:]) >= 0.0:
   NOAH.prcp[:] = meteorology.variables['apcpsfc'][i,:]/dt #mm/s
  else:
   NOAH.prcp[:] = meteorology.variables['prec'][i,:]/dt #mm/s

  return (NOAH,TOPMODEL)

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
             'qbase':{'description':'Excess runoff','units':'mm/s'},
             'qsurface':{'description':'Surface runoff','units':'mm/s'},
             'prcp':{'description':'Precipitation','units':'mm/s'},
             'smc1':{'description':'Soil water content','units':'m3/m3'},
             #'smc2':{'description':'Soil water content','units':'m3/m3'},
             #'smc3':{'description':'Soil water content','units':'m3/m3'},
             #'smc4':{'description':'Soil water content','units':'m3/m3'},
             'swd':{'description':'Soil water deficit','units':'mm'},
             'sstorage':{'description':'Surface storage','units':'mm'},
             #'sndpth':{'description':'Snow depth','units':'mm'},
             #'swe':{'description':'Snow water equivalent','units':'mm'},
             'qout_subsurface':{'description':'Subsurface flux','units':'m2/s'},
             'qout_surface':{'description':'Surface flux','units':'m2/s'},
             }

 #Create the dimensions
 print 'Creating the dimensions'
 ntime = len(fp_in.dimensions['time'])
 nhsu = len(fp_in.dimensions['hsu'])
 #nsoil = 4
 fp_out.createDimension('hsu',nhsu)
 fp_out.createDimension('time',ntime)
 #fp_out.createDimension('soil',nsoil)

 #Create the output
 print 'Creating the catchment group'
 grp = fp_out.createGroup('catchment')
 for var in metadata:
  ncvar = grp.createVariable(var,'f4',('time','hsu',))
  ncvar.description = metadata[var]['description']
  ncvar.units = metadata[var]['units']

 #Create the outlet group 
 #print 'Creating the outlet group'
 #grp = fp_out.createGroup('outlet')
 #grp.createDimension('hru_outlet',len(fp_in.groups['outlet'].groups['summary'].dimensions['hru']))
 #grp.createVariable('hru_org','i4',('hru_outlet',))
 #grp.createVariable('hru_dst','i4',('hru_outlet',))
 #grp.createVariable('counts','i4',('hru_outlet',))
 #grp.createVariable('area_outlet','i4',('hru_outlet',))
 #for var in ['qin_subsurface_outlet','qin_surface_outlet']:
 # ncvar = grp.createVariable(var,'f4',('time','hru_outlet',))
 #dx = info['dx']
 #grp.variables['area_outlet'][:] = dx**2*info['input_fp'].groups['outlet'].groups['summary'].variables['counts'][:]
 #grp.variables['hru_org'][:] = info['input_fp'].groups['outlet'].groups['summary'].variables['hru_org'][:]
 #grp.variables['hru_dst'][:] = info['input_fp'].groups['outlet'].groups['summary'].variables['hru_dst'][:]

 #Create the metadata
 print 'Creating the metadata group'
 grp = fp_out.createGroup('metadata')
 #Time
 grp.createVariable('time','f8',('time',))
 dates = grp.createVariable('date','f8',('time',))
 dates.units = 'hours since 1900-01-01'
 dates.calendar = 'standard'
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

 #General info
 grp = info['output_fp'].groups['metadata']
 grp.variables['time'][itime] = info['time'] #Seconds
 dates = grp.variables['date']
 dates[itime] = nc.date2num(info['date'],units=dates.units,calendar=dates.calendar)

 #Update the variables (catchment)
 grp = info['output_fp'].groups['catchment']

 #NoahMP
 grp.variables['smc1'][itime,:] = np.copy(NOAH.smc[:,0]) #m3/m3
 grp.variables['g'][itime,:] = np.copy(NOAH.ssoil) #W/m2
 grp.variables['sh'][itime,:] = np.copy(NOAH.fsh) #W/m2
 grp.variables['lh'][itime,:] = np.copy(NOAH.fcev + NOAH.fgev + NOAH.fctr) #W/m2
 grp.variables['qbase'][itime,:] = NOAH.dt*np.copy(NOAH.runsb) #mm
 grp.variables['qsurface'][itime,:] = NOAH.dt*np.copy(NOAH.runsf) #mm
 grp.variables['prcp'][itime,:] = NOAH.dt*np.copy(NOAH.prcp) #mm

 #TOPMODEL
 grp.variables['swd'][itime,:] = np.copy(10**3*TOPMODEL.si) #mm
 grp.variables['qout_subsurface'][itime,:] = np.copy(TOPMODEL.qout) #m2/s
 grp.variables['qout_surface'][itime,:] = np.copy(TOPMODEL.qout_surface) #m2/s
 grp.variables['sstorage'][itime,:] = np.copy(TOPMODEL.storage_surface)

 #Update the variables (outlet)
 #grp = info['output_fp'].groups['outlet']
 #grp.variables['qin_subsurface_outlet'][itime,:] = np.copy(TOPMODEL.qin_outlet) #m2/s
 #grp.variables['qin_surface_outlet'][itime,:] = np.copy(TOPMODEL.qin_outlet_surface) #m2/s
 
 return
