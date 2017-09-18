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

def initialize(info):

 HB = HydroBlocks(info)

 return HB

class HydroBlocks:

 def __init__(self,info):

  #Read in general information
  self.general_information(info)

  #Initialize Noah-MP
  print "Initializing Noah-MP"
  self.initialize_noahmp()

  #Initialize Dynamic TOPMODEL
  print "Initializing Dynamic TOPMODEL"
  self.initialize_dtopmodel()

  #Initialize human water use module
  print "Initializing Human Water Management"
  self.initialize_hwu()

  #Other metrics
  self.dE = 0.0
  self.r = 0.0
  self.dr = 0.0
  self.et = 0.0
  self.etran = 0.0
  self.esoil = 0.0
  self.ecan = 0.0
  self.prcp = 0.0
  self.q = 0.0
  self.errwat = 0.0
  self.erreng = 0.0
  self.dzwt0 = np.zeros(self.nhru,dtype=np.float32)
  self.beg_wb = np.zeros(self.nhru,dtype=np.float32)
  self.end_wb = np.zeros(self.nhru,dtype=np.float32)

  return
 
 def general_information(self,info):

  #Define the metadata
  self.dx = info['dx']
  self.dt = info['dt']
  self.dtt = self.dt#info['dtt']
  self.nsoil = info['nsoil']
  self.ncores = info['ncores']
  self.idate = info['idate']
  self.fdate = info['fdate']
  self.mkl_flag = info['mkl_flag']
  self.dt_timedelta = datetime.timedelta(seconds=self.dt)
  self.input_fp = nc.Dataset(info['input_file'])
  self.output_fp = nc.Dataset(info['output_file'],'w',format='NETCDF4')
  self.nhru = len(self.input_fp.dimensions['hsu'])
  self.surface_flow_flag = info['surface_flow_flag']
  self.subsurface_flow_flag = info['subsurface_flow_flag']
  self.hwu_flag = info['hwu_flag']
  self.create_mask_flag = info['create_mask_flag']

  return

 def initialize_noahmp(self,):

  #Initialize noahmp
  from pyNoahMP.NoahMP import model as noahmp
  self.noahmp = noahmp

  #Initialize parameters
  self.noahmp.ncells = self.nhru
  self.noahmp.nsoil = self.nsoil
  self.noahmp.dt = self.dt
  self.noahmp.nsnow = 3
  self.noahmp.llanduse[:] = 'MODIFIED_IGBP_MODIS_NOAH'
  self.noahmp.lsoil[:] = 'CUST'
  dir = os.path.dirname(os.path.abspath(__file__))
  VEGPARM = '%s/pyNoahMP/data/VEGPARM.TBL' % dir
  GENPARM = '%s/pyNoahMP/data/GENPARM.TBL' % dir
  MPTABLE = '%s/pyNoahMP/data/MPTABLE.TBL' % dir
  self.noahmp.vegparm_file[0:len(VEGPARM)] = VEGPARM
  self.noahmp.genparm_file[0:len(GENPARM)] = GENPARM
  self.noahmp.mptable_file[0:len(MPTABLE)] = MPTABLE
  #Define the options
  self.noahmp.idveg = 3#3#4 # dynamic vegetation (1 -> off ; 2 -> on)
  self.noahmp.iopt_crs = 2#2 # canopy stomatal resistance (1-> Ball-Berry; 2->Jarvis)
  self.noahmp.iopt_btr = 1#2 # soil moisture factor for stomatal resistance (1-> Noah; 2-> CLM; 3-> SSiB)
  self.noahmp.iopt_run = 5#5#2#5##2#2 # runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS)
  self.noahmp.iopt_sfc = 2#2 # surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
  self.noahmp.iopt_frz = 2#1 # supercooled liquid water (1-> NY06; 2->Koren99)
  self.noahmp.iopt_inf = 2#1#1#2 # frozen soil permeability (1-> NY06; 2->Koren99)
  self.noahmp.iopt_rad = 2 # radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
  self.noahmp.iopt_alb = 1 # snow surface albedo (1->BATS; 2->CLASS)
  self.noahmp.iopt_snf = 3 # rainfall & snowfall (1-Jordan91; 2->BATS; 3->Noah)]
  self.noahmp.iopt_tbot = 1#1#1 # lower boundary of soil temperature (1->zero-flux; 2->Noah) 
  self.noahmp.iopt_stc = 1#1#2 # snow/soil temperature time scheme (only layer 1) 1 -> semi-implicit; 2 -> full implicit (original Noah)

  #Allocate memory
  self.noahmp.initialize_general()

  #Set info
  self.noahmp.iz0tlnd = 0
  self.noahmp.z_ml[:] = 10.0#2.0 -- Noemi
  tmp = np.array([0.1,0.3,0.6,1.0,2.0,2.0,2.0,2.0]) #Let's bring this out to the metadata
  self.noahmp.sldpth[:] = tmp
  self.noahmp.zsoil[:] = -np.cumsum(self.noahmp.sldpth[:],axis=1)
  self.noahmp.zsnso[:] = 0.0
  self.noahmp.zsnso[:,3::] = self.noahmp.zsoil[:]

  #Set all to land (can include lake in the future...)
  self.noahmp.ist[:] = 1
  #Set soil color
  self.noahmp.isc[:] = 4#4
  self.noahmp.ice[:] = 0
  #Initialize snow info
  self.noahmp.isnow[:] = 0
  #Others
  self.noahmp.dx = self.dx
  self.noahmp.ice[:] = 0
  self.noahmp.foln[:] = 1.0
  self.noahmp.ficeold[:] = 0.0
  self.noahmp.albold[:] = 0.65
  self.noahmp.sneqvo[:] = 0.0
  self.noahmp.ch[:] = 0.0001
  self.noahmp.cm[:] = 0.0001
  self.noahmp.canliq[:] =  0.0
  self.noahmp.canice[:] = 0.0
  self.noahmp.sndpth[:] = 0.0
  self.noahmp.swe[:] = 0.0
  self.noahmp.snice[:] = 0.0
  self.noahmp.snliq[:] = 0.0
  self.noahmp.wa[:] = 0.0
  self.noahmp.wt[:] = self.noahmp.wa[:]
  self.noahmp.zwt[:] = 0.0
  self.noahmp.wslake[:] = 0.0
  self.noahmp.lfmass[:] = 9.0
  self.noahmp.rtmass[:] = 500.0
  self.noahmp.stmass[:] = 3.33
  self.noahmp.wood[:] = 500.0
  self.noahmp.stblcp[:] = 1000.0
  self.noahmp.fastcp[:] = 1000.0
  self.noahmp.plai[:] = 0.5
  self.noahmp.tauss[:] = 0.0
  self.noahmp.smcwtd[:] = self.noahmp.sh2o[:,1]
  self.noahmp.deeprech[:] = 0.0
  self.noahmp.rech[:] = 0.0
  self.noahmp.eah[:] = 1000.0
  self.noahmp.fwet[:] = 0.0
  self.noahmp.psai[:] = 0.1
  self.noahmp.stc[:] = 285.0
  self.noahmp.slopetyp[:] = 3
  self.noahmp.albold[:] = 0.5
  #Define the data
  self.noahmp.vegtyp[:] = self.input_fp.groups['parameters'].variables['land_cover'][:]
  self.noahmp.soiltyp[:] = np.arange(1,self.noahmp.ncells+1)
  self.noahmp.smcmax[:] = self.input_fp.groups['parameters'].variables['MAXSMC'][:]
  self.noahmp.smcref[:] = self.input_fp.groups['parameters'].variables['REFSMC'][:]
  self.noahmp.smcdry[:] = self.input_fp.groups['parameters'].variables['DRYSMC'][:]
  for ilayer in xrange(self.noahmp.sh2o.shape[1]):
   self.noahmp.sh2o[:,ilayer] = self.input_fp.groups['parameters'].variables['MAXSMC'][:]
  self.noahmp.smc[:] = self.noahmp.sh2o[:]
  self.noahmp.smcwtd[:] = self.noahmp.sh2o[:,0]
  #Initialize the soil parameters
  self.noahmp.bb0[:] = self.input_fp.groups['parameters'].variables['BB'][:]
  self.noahmp.drysmc0[:] = self.input_fp.groups['parameters'].variables['DRYSMC'][:]
  self.noahmp.f110[:] = self.input_fp.groups['parameters'].variables['F11'][:]
  self.noahmp.maxsmc0[:] = self.input_fp.groups['parameters'].variables['MAXSMC'][:]
  self.noahmp.refsmc0[:] = self.input_fp.groups['parameters'].variables['REFSMC'][:]
  self.noahmp.satpsi0[:] = self.input_fp.groups['parameters'].variables['SATPSI'][:]
  self.noahmp.satdk0[:] = self.input_fp.groups['parameters'].variables['SATDK'][:]
  self.noahmp.satdw0[:] = self.input_fp.groups['parameters'].variables['SATDW'][:]
  self.noahmp.wltsmc0[:] = self.input_fp.groups['parameters'].variables['WLTSMC'][:]
  self.noahmp.qtz0[:] = self.input_fp.groups['parameters'].variables['QTZ'][:]

  #Set lat/lon (declination calculation)
  self.noahmp.lat[:] = 0.0174532925*self.input_fp.groups['metadata'].latitude
  self.noahmp.lon[:] = 0.0174532925*self.input_fp.groups['metadata'].longitude

  #Initialize output
  self.noahmp.tg[:] = 285.0
  self.noahmp.tv[:] = 285.0
  for ilayer in xrange(self.noahmp.nsnow-1,self.noahmp.stc.shape[1]):
   self.noahmp.stc[:,ilayer] = 285.0
  self.noahmp.tah[:] = 285.0
  self.noahmp.t2mv[:] = 285.0
  self.noahmp.t2mb[:] = 285.0
  self.noahmp.runsf[:] = 0.0
  self.noahmp.runsb[:] = 0.0
  #forcing
  self.noahmp.fveg[:] = 1.0
  self.noahmp.fvgmax[:] = 1.0
  self.noahmp.tbot[:] = 285.0

  #Define the parameters
  self.noahmp.initialize_parameters()

  return

 def initialize_dtopmodel(self,):

  from pyDTopmodel import dynamic_topmodel as dtopmodel

  #Define some metadata
  nhru_outlet = self.input_fp.groups['outlet'].groups['summary'].variables['hru_dst'].size

  #Initialize Dynamic Topmodel
  self.dtopmodel = dtopmodel.Dynamic_Topmodel(self.nhru,nhru_outlet,self.mkl_flag)

  #Set flags
  self.dtopmodel.surface_flow_flag = self.surface_flow_flag
  self.dtopmodel.subsurface_flow_flag = self.subsurface_flow_flag

  #Set self.dtopmodel parameters
  self.dtopmodel.dt = self.dt #seconds
  self.dtopmodel.area[:] = 0.0 #meters^2
  self.dtopmodel.dx[:] = self.dx #meters
  self.dtopmodel.m[:] = self.input_fp.groups['parameters'].variables['m'][:]
  self.dtopmodel.sdmax[:] = self.input_fp.groups['parameters'].variables['sdmax'][:]
  dem = []
  #Set cluster information
  self.dtopmodel.pct[:] = self.input_fp.groups['parameters'].variables['area_pct'][:]/100
  self.dtopmodel.area[:] = self.input_fp.groups['parameters'].variables['area'][:]
  self.dtopmodel.T0[:] = self.input_fp.groups['parameters'].variables['SATDK'][:]*self.dtopmodel.m
  self.dtopmodel.sti[:] = self.input_fp.groups['parameters'].variables['ti'][:]
  self.dtopmodel.beta[:] = self.input_fp.groups['parameters'].variables['slope'][:]
  self.dtopmodel.carea[:] = self.input_fp.groups['parameters'].variables['carea'][:]
  self.dtopmodel.channel[:] = self.input_fp.groups['parameters'].variables['channel'][:]
  self.dtopmodel.dem[:] = self.input_fp.groups['parameters'].variables['dem'][:] 
  self.dtopmodel.mannings[:] = self.input_fp.groups['parameters'].variables['mannings'][:]
  #Set outlet information
  self.dtopmodel.area_outlet[:] = self.dx**2*self.input_fp.groups['outlet'].groups['summary'].variables['counts'][:]
  self.dtopmodel.pct = self.dtopmodel.pct/np.sum(self.dtopmodel.pct)
  ti_mean = np.sum(self.dtopmodel.pct*self.dtopmodel.sti[:])
 
  #Calculate the sti
  lnT0 = np.log(self.dtopmodel.T0)
  lnTe = np.sum(self.dtopmodel.pct*lnT0)
  self.dtopmodel.sti = self.dtopmodel.sti - (lnT0 - lnTe)

  #Set weight matrix
  self.dtopmodel.flow_matrix = sparse.csr_matrix((self.input_fp.groups['flow_matrix'].variables['data'][:],
			          self.input_fp.groups['flow_matrix'].variables['indices'][:],
			          self.input_fp.groups['flow_matrix'].variables['indptr'][:]),
 				  dtype=np.float64)[0:self.nhru,0:self.nhru]
  self.dtopmodel.flow_matrix.setdiag(self.dtopmodel.flow_matrix.diagonal()) #Ensure the zeros are not sparse (for kinematic wave solution).
  self.dtopmodel.flow_matrix_T = sparse.csr_matrix(self.dtopmodel.flow_matrix.T) #transposed
  self.dtopmodel.flow_matrix_T.setdiag(self.dtopmodel.flow_matrix_T.diagonal()) #Ensure the zeros are not sparse  (for kinematic wave solution).

  #Initialize the solver
  if self.mkl_flag:self.dtopmodel.dtt.initialize(self.dtopmodel.flow_matrix_T.indices,self.dtopmodel.flow_matrix_T.indptr)

  #Initialize the soil moisture deficit values
  self.dtopmodel.si[:] = 0.0
 
  return

 def initialize_hwu(self,):

  from Human_Water_Use import Human_Water_Use as hwu
  self.hwu = hwu(self.noahmp,self.nhru)
  self.hwu.hwu_flag = self.hwu_flag
  self.hwu.irrig = np.zeros(self.nhru,dtype=np.float32)

  return

 def run(self,info):

  info['input_fp'] = self.input_fp
  info['output_fp'] = self.output_fp

  #Run the model
  date = self.idate
  i = 0
  tic = time.time()
  self.noahmp.dzwt[:] = 0.0
  self.dtopmodel.ex[:] = 0.0
  while date <= self.fdate:

   #Memorize time stamp
   time0 = time.time()

   if (date.hour == 0) and (date.day == 1):
    testvar = 1
    print date,time.time() - tic,'et:%f'%self.et,'prcp:%f'%self.prcp,'q:%f'%self.q,'WB ERR:%f' % self.errwat,'ENG ERR:%f' % self.erreng
 
   #Update input data
   self.update_input(date,i)#self.noahmp,self.dtopmodel,date,info,i)

   #Save the original precip
   precip = np.copy(self.noahmp.prcp)
 
   #Calculate initial NOAH water balance
   self.initialize_water_balance()

   #Update model
   self.update()

   #Return precip to original value
   self.noahmp.prcp[:] = precip[:] 

   #Calculate final water balance
   self.finalize_water_balance()

   #Update the water balance error
   self.calculate_water_balance_error()

   #Update the energy balance error
   self.calculate_energy_balance_error()
 
   #Update time and date
   self.date = date
   info['date'] = date

   #Update output
   self.update_output(date,i)

   #Update time step
   date = date + self.dt_timedelta
   i = i + 1

  return

 def update_input(self,date,i):

  self.noahmp.itime = i
  dt = self.dt
  self.dtopmodel.itime = i
  self.noahmp.nowdate[:] = date.strftime('%Y-%m-%d_%H:%M:%S')
  self.noahmp.julian = (date - datetime.datetime(date.year,1,1,0)).days
  self.noahmp.yearlen = (datetime.datetime(date.year+1,1,1,0) - datetime.datetime(date.year,1,1,1,0)).days + 1

  #Update meteorology
  meteorology = self.input_fp.groups['meteorology']
  self.noahmp.lwdn[:] = meteorology.variables['lwdown'][i,:] #W/m2
  self.noahmp.swdn[:] = meteorology.variables['swdown'][i,:] #W/m2
  self.noahmp.psfc[:] = meteorology.variables['psurf'][i,:] #Pa
  self.noahmp.p_ml[:] = meteorology.variables['psurf'][i,:] #Pa
  self.noahmp.u_ml[:] = (meteorology.variables['wind'][i,:]**2/2)**0.5 #m/s
  self.noahmp.v_ml[:] = (meteorology.variables['wind'][i,:]**2/2)**0.5 #m/s
  self.noahmp.t_ml[:] = meteorology.variables['tair'][i,:] #K
  self.noahmp.q_ml[:] = meteorology.variables['spfh'][i,:] #Kg/Kg
  self.noahmp.qsfc1d[:] = meteorology.variables['spfh'][i,:] #Kg/Kg
  self.noahmp.prcp[:] = meteorology.variables['precip'][i,:] #mm/s

  return

 def update(self,):

  #Abstract Water from Human Use 
  if self.hwu.hwu_flag == True:
   if self.hwu.itime > 0: self.hwu.irrigation(self.noahmp,self)
   self.hwu.itime = self.hwu.itime +1

  #Set the partial pressure of CO2 and O2
  self.noahmp.co2air[:] = 355.E-6*self.noahmp.psfc[:]# ! Partial pressure of CO2 (Pa) ! From NOAH-MP-WRF
  self.noahmp.o2air[:] = 0.209*self.noahmp.psfc[:]# ! Partial pressure of O2 (Pa)  ! From NOAH-MP-WRF

  #Update NOAH
  self.noahmp.run_model(self.ncores)

  #Reinitialize dzwt
  self.noahmp.dzwt[:] = 0.0

  if self.subsurface_flow_flag == True:

   #Calculate the updated soil moisture deficit
   si0 = np.copy(NOAH.si0)
   si1 = np.copy(NOAH.si1)

   #Calculate the change in deficit
   self.dtopmodel.si[:] = si1[:]
   self.dtopmodel.dsi[:] = np.copy(si1 - si0)
   self.dtopmodel.r[:] = -self.dtopmodel.dsi[:]/self.dtopmodel.dt

   #Add the surface runoff
   self.dtopmodel.qsurf[:] = self.noahmp.runsf[:]/1000.0

   #Update dynamic topmodel
   self.dtopmodel.update(self.ncores)

   #Calculate the change in deficit
   self.dtopmodel.sideep[:] = self.dtopmodel.sideep.astype(np.float32)
   self.dtopmodel.si[:] = self.dtopmodel.si.astype(np.float32)
   dsi = np.copy(si1 - self.dtopmodel.si)

   #Update the soil moisture values
   self.noahmp.dzwt[:] = dsi+self.dtopmodel.dt*self.dtopmodel.ex-self.dtopmodel.dt*self.dtopmodel.r

  return

 def initialize_water_balance(self,):
 
  smw = np.sum(1000*self.noahmp.sldpth*self.noahmp.smc,axis=1)
  self.beg_wb = np.copy(self.noahmp.canliq + self.noahmp.canice + self.noahmp.swe + self.noahmp.wa + smw)
  self.dzwt0 = np.copy(self.noahmp.dzwt)

  return 

 def finalize_water_balance(self,):

  NOAH = self.noahmp
  smw = np.sum(1000*NOAH.sldpth*NOAH.smc,axis=1) 
  self.end_wb = np.copy(NOAH.canliq + NOAH.canice + NOAH.swe + NOAH.wa + smw)# + TOPMODEL.si*1000)
  
  return

 def calculate_water_balance_error(self,):

  NOAH = self.noahmp
  TOPMODEL = self.dtopmodel
  HWU = self.hwu
  dt = self.dt
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

 def calculate_energy_balance_error(self,):

  NOAH = self.noahmp
  TOPMODEL = self.dtopmodel
  tmp = np.copy(NOAH.sav+NOAH.sag-NOAH.fira-NOAH.fsh-NOAH.fcev-NOAH.fgev-NOAH.fctr-NOAH.ssoil)
  self.erreng += np.sum(TOPMODEL.pct*tmp)

  return

 def update_output(self,date,itime):

  NOAH = self.noahmp
  TOPMODEL = self.dtopmodel
  HWU = self.hwu
  HB = self

  #Create the netcdf file
  if itime == 0: self.create_netcdf_file()

  #General info
  grp = self.output_fp.groups['metadata']
  dates = grp.variables['date']
  dates[itime] = nc.date2num(date,units=dates.units,calendar=dates.calendar)

  #Update the variables (macroscale)
  grp = self.output_fp.groups['macroscale']

  #NoahMP
  cs = np.cumsum(NOAH.sldpth[0,:])
  mask = cs <= 0.1
  pct = NOAH.sldpth[0,mask]/np.sum(NOAH.sldpth[0,mask])
  grp.variables['smc1'][itime,:] = np.sum(pct*NOAH.smc[:,mask],axis=1) #m3/m3
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

  #New
  grp.variables['wtd'][itime,:] = np.copy(NOAH.zwt)
  grp.variables['errwat'][itime,:] = np.copy(HB.errwat)
  grp.variables['totsmc'][itime,:] = smw = np.sum(1000*NOAH.sldpth*NOAH.smc,axis=1)

  return

 def create_netcdf_file(self,):

  #Extract the pointers to both input and output files
  fp_out = self.output_fp
  fp_in = self.input_fp

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
             'wtd':{'description':'WTD','units':'m'},
             'errwat':{'description':'errwat','units':'mm'},
             'totsmc':{'description':'totsmc','units':'??'},
             }

  #Create the dimensions
  print 'Creating the dimensions'
  ntime = len(fp_in.dimensions['time'])
  nhru = len(fp_in.dimensions['hsu'])
  fp_out.createDimension('hru',nhru)
  fp_out.createDimension('time',ntime)

  #Create the output
  print 'Creating the macroscale group'
  grp = fp_out.createGroup('macroscale')
  for var in metadata:
   ncvar = grp.createVariable(var,'f4',('time','hru',))
   ncvar.description = metadata[var]['description']
   ncvar.units = metadata[var]['units']

  #Create the metadata
  print 'Creating the metadata group'
  grp = fp_out.createGroup('metadata')
  #Time
  grp.createVariable('time','f8',('time',))
  dates = grp.createVariable('date','f8',('time',))
  dates.units = 'hours since 1900-01-01'
  dates.calendar = 'standard'

  #HRU percentage coverage
  print 'Setting the HRU percentage coverage'
  pcts = grp.createVariable('pct','f4',('hru',))
  pcts[:] = fp_in.groups['parameters'].variables['area_pct'][:]
  pcts.description = 'hru percentage coverage'
  pcts.units = '%'
  #HRU area
  print 'Setting the HRU areal coverage'
  area = grp.createVariable('area','f4',('hru',))
  area[:] = fp_in.groups['parameters'].variables['area'][:]
  area.units = 'meters squared'
  print 'Defining the HRU ids'
  hru = grp.createVariable('hru','i4',('hru',))
  hrus =[]
  for value in xrange(nhru):hrus.append(value)
  hru[:] = np.array(hrus)
  hru.description = 'hru ids'

  #Create the mapping
  if self.create_mask_flag == True:
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

 def finalize(self,):

  #Deallocate the NOAH memory
  self.noahmp.finalize()

  #Close the TOPMODEL solver
  if self.mkl_flag: self.dtopmodel.dtt.finalize()

  #Delete the objects
  del self.noahmp
  del self.dtopmodel

  #Close the files
  self.input_fp.close()
  self.output_fp.close()

  return
