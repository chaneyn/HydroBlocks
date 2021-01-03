import os
import numpy as np
import netCDF4 as nc
import h5py
import datetime
import time
import sys
import scipy.sparse as sparse
import pickle

def assign_string(nelem,pstring):

 tmp = np.chararray(shape=(nelem,))
 tmp[:] = ' '
 tmp2 = []
 for letter in pstring:
  tmp2.append(letter)
 tmp[0:len(tmp2)] = tmp2

 return tmp

def initialize(info):

 HB = HydroBlocks(info)

 return HB

class HydroBlocks:

 def __init__(self,info):

  #Read in general information
  self.general_information(info)

  #Initialize Noah-MP
  #print("Initializing Noah-MP",flush=True)
  self.initialize_noahmp()
  self.tsno=np.zeros((self.noahmp.ncells,self.noahmp.nsnow),order='F').astype(np.float32) #Initial condition TSNO Laura
  #Initialize subsurface module
  #print("Initializing subsurface module",flush=True)
  self.initialize_subsurface()

  #Initialize human water use module
  #print("Initializing Human Water Management",flush=True)
  #self.initialize_hwu(info)

  #Initialize routing module
  #print("Initializing the routing module",flush=True)
  self.initialize_routing()
  #Other metrics
  self.dE = 0.0
  self.r = 0.0
  self.dr = 0.0
  self.acc_et = 0.0
  self.acc_etran = 0.0
  self.acc_edir = 0.0
  self.acc_ecan = 0.0
  self.acc_prcp = 0.0
  self.acc_q = 0.0
  self.acc_errwat = 0.0
  self.acc_erreng = 0.0
  self.runtime = 0.0
  self.runtime0 = 0.0
  self.dzwt0 = np.zeros(self.nhru,dtype=np.float32)
  self.beg_wb = np.zeros(self.nhru,dtype=np.float32)
  self.end_wb = np.zeros(self.nhru,dtype=np.float32)
  self.itime = 0

  #Restart from initial conditions?
  self.restart()

  return
  
 def restart(self,):

  file_restart = '%s/%s.h5' % (self.metadata['restart']['dir'],self.idate.strftime('%Y-%m-%d'))
  if (os.path.exists(file_restart) == False 
      or self.metadata['restart']['flag'] == False):
    #print("Cold startup",flush=True)
    return

  #Read in the restart information
  fp = h5py.File(file_restart,'r')
  #diagnostics
  self.acc_et = fp['acc_et'][()]
  self.acc_etran = fp['acc_etran'][()]
  self.acc_edir = fp['acc_edir'][()]
  self.acc_ecan = fp['acc_ecan'][()]
  self.acc_prcp = fp['acc_prcp'][()]
  self.acc_q = fp['acc_q'][()]
  self.acc_errwat = fp['acc_errwat'][()]
  self.acc_erreng = fp['acc_erreng'][()]
  self.runtime0 = fp['runtime'][()]
  #noahmp
  self.noahmp.smceq[:] = fp['smceq'][:]
  self.noahmp.albold[:] = fp['albold'][:]
  self.noahmp.albedo[:] = fp['albedo'][:]
  self.noahmp.sneqvo[:] = fp['sneqvo'][:]
  self.noahmp.tslb[:] = fp['tslb'][:]
  self.noahmp.sh2o[:] = fp['sh2o'][:]
  self.noahmp.smois[:] = fp['smois'][:]
  self.noahmp.tah[:] = fp['tah'][:]
  self.noahmp.eah[:] = fp['eah'][:]
  self.noahmp.fwet[:] = fp['fwet'][:]
  self.noahmp.canliq[:] = fp['canliq'][:]
  self.noahmp.canice[:] = fp['canice'][:]
  self.noahmp.tv[:] = fp['tv'][:]
  self.noahmp.tg[:] = fp['tg'][:]
  self.noahmp.qsfc1d[:] = fp['qsfc1d'][:]
  self.noahmp.qsnow[:] = fp['qsnow'][:]
  self.noahmp.isnow[:] = fp['isnow'][:]
  self.noahmp.zsnso[:] = fp['zsnso'][:]
  self.noahmp.snowh[:] = fp['snowh'][:]
  self.noahmp.snow[:] = fp['snow'][:]
  self.noahmp.snice[:] = fp['snice'][:]
  self.noahmp.snliq[:] = fp['snliq'][:]
  self.noahmp.zwt[:] = fp['zwt'][:]
  self.noahmp.wa[:] = fp['wa'][:]
  self.noahmp.wt[:] = fp['wt'][:]
  self.noahmp.wslake[:] = fp['wslake'][:]
  self.noahmp.lfmass[:] = fp['lfmass'][:]
  self.noahmp.rtmass[:] = fp['rtmass'][:]
  self.noahmp.stmass[:] = fp['stmass'][:]
  self.noahmp.wood[:] = fp['wood'][:]
  self.noahmp.stblcp[:] = fp['stblcp'][:]
  self.noahmp.fastcp[:] = fp['fastcp'][:]
  self.noahmp.lai[:] = fp['lai'][:]
  self.noahmp.sai[:] = fp['sai'][:]
  self.noahmp.cm[:] = fp['cm'][:]
  self.noahmp.ch[:] = fp['ch'][:]
  self.noahmp.tauss[:] = fp['tauss'][:]
  self.noahmp.smcwtd[:] = fp['smcwtd'][:]
  self.noahmp.deeprech[:] = fp['deeprech'][:]
  self.noahmp.rech[:] = fp['rech'][:]
  self.noahmp.tsno[:] = fp['tsno'][:]
  self.noahmp.grain[:] = fp['grain'][:]
  self.noahmp.gdd[:] = fp['gdd'][:]
  self.noahmp.fveg[:] = fp['fveg'][:]
  self.noahmp.fvgmax[:] = fp['fvgmax'][:]
  self.noahmp.acsnom[:] = fp['acsnom'][:]
  self.noahmp.acsnow[:] = fp['acsnow'][:]
  self.noahmp.qsfc[:] = fp['qsfc'][:]
  self.noahmp.sfcrunoff[:] = fp['sfcrunoff'][:]
  self.noahmp.udrunoff[:] = fp['udrunoff'][:]
  #routing
  if self.routing_module == 'kinematic':
   self.routing.Q0[:] = fp['Q0'][:]
   self.routing.u0[:] = fp['u0'][:]
   self.routing.A0[:] = fp['A0'][:]
   self.routing.Q1[:] = fp['Q1'][:]
   self.routing.A1[:] = fp['A1'][:]
   self.routing.dA[:] = fp['dA'][:]
   self.routing.bcs[:] = fp['bcs'][:]
   self.routing.qss[:] = fp['qss'][:]
   self.routing.hband_inundation[:] = fp['hband_inundation'][:]
   self.routing.reach2hband_inundation[:] = fp['reach2hband_inundation'][:]
  fp.close()

  return
 
 def general_information(self,info):

  #Parallel information
  self.MPI  = info['MPI']

  #Define the metadata
  self.cid = info['cid']
  self.dt = info['dt']
  self.dtt = self.dt#info['dtt']
  self.nsoil = len(info['dz'])#['nsoil']
  self.idate = info['idate']
  self.fdate = info['fdate']
  self.cdir = info['cdir']
  self.Qobs_file = info['Qobs_file']
  self.dt_timedelta = datetime.timedelta(seconds=self.dt)
  self.input_fp = nc.Dataset(info['input_file'])
  self.dx = self.input_fp.groups['metadata'].dx
  self.nhru = len(self.input_fp.dimensions['hru'])
  self.subsurface_module = info['subsurface_module']
  self.routing_module = info['routing_module']['type']
  self.routing_surface_coupling = info['routing_module']['surface_coupling']
  #self.hwu_flag = info['water_management']['hwu_flag']
  self.area = self.input_fp.groups['parameters'].variables['area'][:]
  self.pct = self.input_fp.groups['parameters'].variables['area_pct'][:]/100
  self.pct = self.pct/np.sum(self.pct)
  self.metadata = info
  self.m = self.input_fp.groups['parameters'].variables['m'][:]  #Noemi
  self.m[:] = 10.0 #m
  self.input_fp_meteo_time = self.input_fp.groups['meteorology'].variables['time']

  #Create a list of all the dates
  dates = []
  date = self.idate
  while date < self.fdate:
   dates.append(date)
   date = date + self.dt_timedelta
  self.dates = np.array(dates)

  #Read in all the meteorology for the given time period
  #determine the first time step for the meteorology
  var = self.input_fp['meteorology'].variables['time']
  ndates = var[:]
  #convert current date to num
  ndate = nc.date2num(self.idate,units=var.units,calendar=var.calendar)
  minitial_time = np.where(ndates == ndate)[0][0]
  #define period to extract
  idate = nc.date2num(self.idate,units=var.units,calendar=var.calendar)
  fdate = nc.date2num(self.fdate,units=var.units,calendar=var.calendar)
  idx  = np.where((ndates >= idate) & (ndates < fdate))[0]
  self.meteodata = {}
  for var in ['lwdown','swdown','psurf','wind','tair','spfh','precip']:
   self.meteodata[var] = self.input_fp['meteorology'][var][idx,:]
  #Define subtimestep
  dt = self.dt
  dt_meteo = (ndates[1] - ndates[0])*3600.0
  self.nt_timestep = int(dt_meteo/dt)

  return

 def initialize_noahmp(self,):

  #Initialize noahmp
  from model.pyNoahMP import NoahMP
  #self.noahmp = pyNoahMP.NoahMP
  self.noahmp = NoahMP

  #Initialize parameters
  self.noahmp.ncells = self.nhru
  self.noahmp.nsoil = self.nsoil
  self.noahmp.nsnow = self.metadata['noahmp_options']['nsnow']
  self.noahmp.dt = self.dt

  #Allocate memory
  #1d,integer
  vars = ['isnow','ice','isltyp','ivgtyp','cropcat','pgs']
  for var in vars:
   exec('self.noahmp.%s = np.zeros(self.nhru).astype(np.int32)' % var)
   exec('self.noahmp.%s[:] = -2147483647' % var)
  #1d,real
  vars = ['z_ml','lwdn','swdn','p_ml','psfc','prcp','t_ml','q_ml','u_ml','v_ml','fsh','ssoil','salb',\
          'fsno','swe','emissi','qsfc1d','tv','tg','canice','canliq','eah','tah','cm','ch',\
          'fwet','sneqvo','albold','qsnow','wslake','zwt','dzwt','wa','wt','smcwtd','deeprech','rech',\
          'lfmass','rtmass','stmass','wood','stblcp','fastcp','tauss','t2mv','t2mb','q2mv',\
          'q2mb','trad','nee','gpp','npp','fvegmp','runsf','runsb','ecan','etran','esoil','fsa','fira',\
          'apar','psn','sav','sag','rssun','rssha','bgap','wgap','tgv','tgb','chv','chb','irc','irg','shc',\
          'shg','evg','ghv','irb','shb','evb','ghb','tr','evc','chleaf','chuc','chv2','chb2','cosz','lat',\
          'lon','fveg','fvgmax','fpice','fcev','fgev','fctr','qsnbot',\
          'foln','smcmax','smcdry','smcref','errwat','zwt0','minzwt',\
          'mozb','fvb','mozv','fvv','zpd','zpdg','tauxb','tauyb','tauxv','tauyv','sfcheadrt',\
          'snow','snowh','canwat','xlat','xlon','tsk','tmn','xice','sai','lai','grain',\
          'gdd','chstar','area','msftx','msfty','qrfs','qsprings','qslat','fdepth',\
          'ht','riverbed','eqzwt','rivercond','pexp','rechclim',\
          'planting','harvest','season_gdd',\
          'hfx','qfx','lh','grdflx','smstav','smstot',\
          'sfcrunoff','udrunoff','albedo','snowc',\
          'acsnom','acsnow',\
          'emiss','qsfc',\
          'z0','znt',\
          'edir','rs',\
          'pah']
  for var in vars:
   exec('self.noahmp.%s = np.zeros(self.nhru,order=\'F\').astype(np.float32)' % var)
   exec('self.noahmp.%s[:] = 9999999999.0' % var)
  #2d,real
  vars = ['sh2o','smc','smceq','zsoil','sldpth','hdiv','smois','tslb','smoiseq',\
          'bexp','smcdry','smcwlt','smcref','smcmax','dksat','dwsat','psisat','quartz']
  for var in vars:
   exec('self.noahmp.%s = np.zeros((self.nhru,self.nsoil),order=\'F\').astype(np.float32)' % var)
   exec('self.noahmp.%s[:] = 9999999999.0' % var)
  self.noahmp.zsnso = np.zeros((self.nhru,self.noahmp.nsoil+self.noahmp.nsnow),order='F').astype(np.float32)
  self.noahmp.stc = np.zeros((self.nhru,self.noahmp.nsoil+self.noahmp.nsnow),order='F').astype(np.float32)
  self.noahmp.snice = np.zeros((self.nhru,self.noahmp.nsnow),order='F').astype(np.float32)
  self.noahmp.snliq = np.zeros((self.nhru,self.noahmp.nsnow),order='F').astype(np.float32)
  self.noahmp.dzs = np.zeros((self.noahmp.nsoil),order='F').astype(np.float32)
  self.noahmp.croptype = np.zeros((self.nhru,5),order='F').astype(np.float32)
  self.noahmp.fndsoilw = False
  self.noahmp.fndsnowh = True
  self.noahmp.tsno = np.zeros((self.nhru,self.noahmp.nsnow),order='F').astype(np.float32)
  self.noahmp.restart = False
  self.noahmp.allowed_to_read = True
  self.noahmp.dx = 0.0
  self.noahmp.dy = 0.0
  self.noahmp.wtddt = 0.0
  self.noahmp.stepwtd = 0
  self.noahmp.gecros_state = np.zeros((self.nhru,60),order='F').astype(np.float32)
  self.noahmp.iz0tlnd = 0#-16497279 #0
  self.noahmp.planting[:] = 126.0
  self.noahmp.harvest[:] = 290.0
  self.noahmp.season_gdd[:] = 1605.0
  self.noahmp.fvgmax[:] = 96
  #Initialize the rest
  self.noahmp.hdiv[:] = 0.0
  self.noahmp.acsnom[:] = 0.0
  self.noahmp.acsnow[:] = 0.0
  self.noahmp.sfcrunoff[:] = 0.0
  self.noahmp.udrunoff[:] = 0.0
  self.noahmp.deeprech[:] = 0.0
  self.noahmp.rech[:] = 0.0
  self.noahmp.rechclim[:] = 0.0
  #Define new ics
  self.noahmp.xice[:] = 0.0
  self.noahmp.snow[:] = 1.0
  self.noahmp.snowh[:] = 0.01
  self.noahmp.isltyp[:] = self.input_fp.groups['parameters'].variables['soil_texture_class'][:]
  self.noahmp.ivgtyp[:] = self.input_fp.groups['parameters'].variables['land_cover'][:]
  self.noahmp.tslb[:] = np.array([266.1,274.0,276.9,279.9])
  self.noahmp.smois[:] = np.array([0.298,0.294,0.271,0.307])
  self.noahmp.dzs[:] = np.array(self.metadata['dz'])
  self.noahmp.tsk[:] = 263.70
  self.noahmp.tmn[:] = 285.0
  self.noahmp.ht[:] = 218.0
  self.noahmp.lai[:] = 2.0
  self.noahmp.xlat[:] = self.input_fp.groups['metadata'].latitude
  self.noahmp.xlon[:] = self.input_fp.groups['metadata'].longitude
  self.noahmp.wtddt = 30.0
  self.noahmp.stepwtd = np.max([np.round(self.noahmp.wtddt*60.0/self.noahmp.dt),1])
  self.noahmp.runsf[:] = 0.0
  self.noahmp.runsb[:] = 0.0
  self.noahmp.canliq[:] = 0.0
  self.noahmp.canice[:] = 0.0
  self.noahmp.canliq[:] = 0.0
  self.noahmp.swe[:] = 0.0
  self.noahmp.wa[:] = 0.0
  #End new ics

  #Define data
  self.noahmp.llanduse = 'MODIFIED_IGBP_MODIS_NOAH'
  self.noahmp.lsoil = 'CUST'
  fdir = os.path.dirname(os.path.abspath(__file__))
  SOILPARM = '%s/pyNoahMP/data/SOILPARM.TBL' % fdir
  GENPARM = '%s/pyNoahMP/data/GENPARM.TBL' % fdir
  MPTABLE = '%s/pyNoahMP/data/MPTABLE.TBL' % fdir
  self.noahmp.soilparm_file = SOILPARM
  self.noahmp.genparm_file = GENPARM
  self.noahmp.mptable_file = MPTABLE
  #Define the noahmp options
  self.noahmp.idveg = self.metadata['noahmp_options']['idveg'] # dynamic vegetation [Look at noahmp .F files for opts]
  self.noahmp.iopt_crs = self.metadata['noahmp_options']['iopt_crs'] #canopy stomatal resistance (1-> Ball-Berry; 2->Jarvis)
  self.noahmp.iopt_btr = self.metadata['noahmp_options']['iopt_btr'] # sm factor for st. resistance (1-> Noah; 2-> CLM; 3-> SSiB)
  self.noahmp.iopt_run = self.metadata['noahmp_options']['iopt_run'] # runoff and gw (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS)
  self.noahmp.iopt_sfc = self.metadata['noahmp_options']['iopt_sfc'] # surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
  self.noahmp.iopt_frz = self.metadata['noahmp_options']['iopt_frz'] # supercooled liquid water (1-> NY06; 2->Koren99)
  self.noahmp.iopt_inf = self.metadata['noahmp_options']['iopt_inf'] # frozen soil permeability (1-> NY06; 2->Koren99)
  self.noahmp.iopt_rad = self.metadata['noahmp_options']['iopt_rad'] # radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
  self.noahmp.iopt_alb = self.metadata['noahmp_options']['iopt_alb'] # snow surface albedo (1->BATS; 2->CLASS)
  self.noahmp.iopt_snf = self.metadata['noahmp_options']['iopt_snf'] # rainfall & snowfall (1-Jordan91; 2->BATS; 3->Noah)]
  self.noahmp.iopt_tbot = self.metadata['noahmp_options']['iopt_tbot'] # lower boundary of soil temperature (1->zero-flux; 2->Noah) 
  self.noahmp.iopt_stc = self.metadata['noahmp_options']['iopt_stc'] # snow/soil temp. time scheme 1 -> semi-implicit; 2 -> implicit
  self.noahmp.iopt_gla = self.metadata['noahmp_options']['iopt_gla']
  self.noahmp.iopt_rsf = self.metadata['noahmp_options']['iopt_rsf']
  self.noahmp.iopt_crop = self.metadata['noahmp_options']['iopt_crop']
  self.noahmp.sf_urban_physics = self.metadata['noahmp_options']['sf_urban_physics']
  self.noahmp.sldpth[:] = np.array(self.metadata['dz'])
  self.noahmp.z_ml[:] = self.metadata['noahmp_options']['z_ml']
  self.noahmp.zsoil = -np.cumsum(self.noahmp.sldpth[:],axis=1)
  self.noahmp.zsnso[:] = 9999999999.0
  self.noahmp.zsnso[:,self.noahmp.nsnow:] = self.noahmp.zsoil[:]
  self.noahmp.dx = self.dx
  self.noahmp.tauss[:] = 0.0
  #Define the data
  hand = self.input_fp.groups['parameters'].variables['hand'][:]
  self.hrus = self.input_fp.groups['parameters'].variables['hru'][:].astype(np.int32)
  self.hbands = self.input_fp.groups['parameters'].variables['hband'][:].astype(np.int32)
  self.nhband = np.unique(self.hbands).size
  self.noahmp.clay_pct  = self.input_fp.groups['parameters'].variables['clay'][:] # Noemi
  self.noahmp.smcwtd[:] = self.noahmp.sh2o[:,0]
  #Initialize the soil parameters
  for isoil in range(self.noahmp.nsoil):
   self.noahmp.bexp[:,isoil] = self.input_fp.groups['parameters'].variables['BB'][:]
   self.noahmp.smcdry[:,isoil] = self.input_fp.groups['parameters'].variables['DRYSMC'][:]
   self.noahmp.smcwlt[:,isoil] = self.input_fp.groups['parameters'].variables['WLTSMC'][:]
   self.noahmp.smcref[:,isoil] = self.input_fp.groups['parameters'].variables['REFSMC'][:]
   self.noahmp.smcmax[:,isoil] = self.input_fp.groups['parameters'].variables['MAXSMC'][:]
   self.noahmp.dksat[:,isoil] = self.input_fp.groups['parameters'].variables['SATDK'][:]
   self.noahmp.dwsat[:,isoil] = self.input_fp.groups['parameters'].variables['SATDW'][:]
   self.noahmp.psisat[:,isoil] = self.input_fp.groups['parameters'].variables['SATPSI'][:]
   self.noahmp.quartz[:,isoil] = self.input_fp.groups['parameters'].variables['QTZ'][:]

  #Set lat/lon (declination calculation)
  self.noahmp.lat[:] = 0.0174532925*self.input_fp.groups['metadata'].latitude
  self.noahmp.lon[:] = 0.0174532925*self.input_fp.groups['metadata'].longitude

  #Define the parameters
  noah = self.noahmp
  self.noahmp.initialize(noah.llanduse,noah.snow,noah.snowh,noah.canwat,\
          noah.isltyp,noah.ivgtyp,noah.xlat,noah.tslb,noah.smois,noah.sh2o,noah.dzs,\
          noah.fndsoilw,noah.fndsnowh,noah.tsk,noah.isnow,noah.tv,noah.tg,noah.canice,\
          noah.tmn,noah.xice,noah.canliq,noah.eah,noah.tah,noah.cm,noah.ch,\
          noah.fwet,noah.sneqvo,noah.albold,noah.qsnow,noah.wslake,noah.zwt,noah.wa,\
          noah.wt,noah.tsno,noah.zsnso,noah.snice,noah.snliq,noah.lfmass,noah.rtmass,\
          noah.stmass,noah.wood,noah.stblcp,noah.fastcp,noah.sai,noah.lai,noah.grain,\
          noah.gdd,noah.croptype,noah.cropcat,noah.t2mv,noah.t2mb,noah.chstar,\
          noah.restart,noah.allowed_to_read,noah.iopt_run,noah.iopt_crop,\
          noah.sf_urban_physics,noah.smoiseq,noah.smcwtd,noah.rech,\
          noah.deeprech,noah.area,noah.dx,noah.dy,noah.msftx,noah.msfty,noah.wtddt,\
          noah.stepwtd,noah.dt,noah.qrfs,noah.qsprings,noah.qslat,noah.fdepth,\
          noah.ht,noah.riverbed,noah.eqzwt,noah.rivercond,noah.pexp,noah.rechclim,\
          noah.gecros_state,\
          noah.soilparm_file,noah.genparm_file,noah.mptable_file)

  return

 def initialize_routing(self,):

  #Determine what rank has which cid
  self.comm = self.MPI.COMM_WORLD
  self.size = self.comm.Get_size()
  self.rank = self.comm.Get_rank()
  if self.rank != 0:
   dest = 0
   db_ex = {self.cid:self.rank}
   self.comm.send(db_ex,dest=dest,tag=11)
  elif self.rank == 0:
   db = {}
   db[self.cid] = 0
   for i in range(1,self.size):
    db_ex = self.comm.recv(source=i,tag=11)
    db[list(db_ex.keys())[0]] = db_ex[list(db_ex.keys())[0]]
  #Wait until completed
  self.comm.Barrier()
  #Send the list to all the cores now
  if self.rank == 0:
   for i in range(1,self.size):
    self.comm.send(db,dest=i,tag=11)
  if self.rank != 0:
   db = self.comm.recv(source=0,tag=11)
  #Memorize links
  self.cid_rank_mapping = db

  if self.routing_module == 'kinematic':self.initialize_kinematic()

  return

 def initialize_kinematic(self,):

  from model.pyRouting import routing

  #Initialize kinematic wave routing
  self.routing = routing.kinematic(self.MPI,self.cid,self.cid_rank_mapping,self.dt,
                                   self.nhband,self.nhru,self.cdir,self.Qobs_file)
  self.routing.calculate_inundation_height_per_hband = routing.calculate_inundation_height_per_hband

  #Initialize hru to reach IRF
  ntmax = 100
  self.routing.IRF = {}
  uhs = self.routing.uhs[:]
  uhs = uhs/np.sum(uhs,axis=1)[:,np.newaxis]
  self.routing.IRF['uh'] = uhs[:]
  #Setup qfuture
  self.routing.IRF['qfuture'] = np.zeros((self.nhband,ntmax-1))

  return

 def update_routing(self,):

  if self.routing_module == 'kinematic':

   #Clean up runoff
   self.noahmp.runsf[self.noahmp.runsf < 0] = 0.0
   self.noahmp.runsb[self.noahmp.runsb < 0] = 0.0

   #Compute total runoff for the HRU
   runoff = self.noahmp.runsf+self.noahmp.runsb

   #Back out the actual runoff for this time step
   runoff_true = runoff - self.routing.hru_inundation*1000/self.dt #mm/s
   runoff_true[runoff_true < 0] = 0.0
   r = np.copy(runoff_true)
   r[:] = 1.0
   r[runoff > 0] = runoff_true[runoff > 0]/runoff[runoff > 0]
   self.noahmp.runsf[:] = r*self.noahmp.runsf[:]
   self.noahmp.runsb[:] = r*self.noahmp.runsb[:]
   self.routing.hru_runoff_inundation[:] = runoff - self.noahmp.runsf - self.noahmp.runsb
   self.routing.hru_runoff_inundation[self.routing.hru_runoff_inundation<0] = 0.0

   #Calculate runoff per hband (This needs to be optimized; Numba?)
   runoff_hband = np.zeros(self.nhband)
   area_hband = np.zeros(self.nhband)
   for hru in self.hrus:
    runoff_hband[self.hbands[hru]] += self.area[hru]*runoff[hru]
    area_hband[self.hbands[hru]] += self.area[hru]
   runoff_hband = runoff_hband/area_hband

   if self.routing_surface_coupling == True:
    #Zero out hband_inundation1
    self.routing.hband_inundation1[:] = 0.0

    #Over inundated hbands set runoff to be inundation and zero out corresponding runoff HERE
    m = self.routing.hband_inundation > 0
    self.routing.hband_inundation1[m] = self.dt*runoff_hband[m]/1000.0 #PROBLEMATIC
    runoff_hband[m] = 0.0

    #Determine the change per hband per reach
    scaling = self.routing.hband_inundation1[m]/(self.routing.fct_infiltrate*self.routing.hband_inundation[m])
    self.routing.reach2hband_inundation[:,m] = (1-self.routing.fct_infiltrate)*self.routing.reach2hband_inundation[:,m] + scaling*self.routing.fct_infiltrate*self.routing.reach2hband_inundation[:,m]
    #self.routing.reach2hband_inundation[:,m] = scaling*self.routing.reach2hband_inundation[:,m]

    #Determine the change per hband per reach (option 2)
    '''for ihband in range(self.routing.reach2hband_inundation.shape[1]):
     #Calculate total volume of water that was removed from the height band
     m = self.routing.reach2hband[:,ihband] > 0
     clength_hband = self.routing.c_length[m]
     Vloss = np.sum(clength_hband*(self.routing.hband_inundation1[ihband] - self.routing.hband_inundation[ihband]))
     #Split up the total volume of water among the reaches
     for ireach in range(self.routing.reach2hband_inundation.shape[0]):
      if m[ireach] == False:continue
      self.routing.reach2hband_inundation[ireach,ihband] = self.routing.reach2hband_inundation[ireach,ihband] - self.routing.c_length[ireach]/np.sum(clength_hband)*Vloss
      #Vloss += 
     print(ihband,np.sum(m),Vloss)'''

    #dif = self.routing.hband_inundation1[m] - self.routing.fct_infiltrate*self.routing.hband_inundation[m]
    #self.routing.reach2hband_inundation[:,m] = 0.5*self.routing.reach2hband_inundation[:,m]
    #self.routing.reach2hband_inundation[:,m] = self.routing.reach2hband_inundation[:,m] + dif
    #self.routing.reach2hband_inundation[(self.routing.reach2hband > 0) & (self.routing.reach2hband_inundation < 0.0)] = 0.0
    #tmp = self.routing.reach2hband_inundation[range(self.routing.reach2hband_inundation.shape[0]),self.routing.hband_channel]
    #tmp[tmp < 0.1] = 0.1
    #self.routing.reach2hband_inundation[range(self.routing.reach2hband_inundation.shape[0]),self.routing.hband_channel] = tmp[:]
    #print(np.unique(self.routing.reach2hband_inundation))
    

   #Recompute cross sectional area
   A = np.sum(self.routing.reach2hband*self.routing.reach2hband_inundation,axis=1)/self.routing.c_length
   Adb = self.routing.hdb['Ac'] + self.routing.hdb['Af']
   (self.routing.hband_inundation[:],self.routing.reach2hband_inundation[:]) = self.routing.calculate_inundation_height_per_hband(Adb,A,self.routing.hdb['W'],self.routing.hdb['M'],self.routing.hdb['hand'],self.routing.hdb['hband'].astype(np.int64),self.routing.reach2hband_inundation,self.routing.reach2hband)

   #Apply convolution to the runoff of each hband
   qr = runoff_hband[:,np.newaxis]*self.routing.IRF['uh']
   #QC to ensure we are conserving mass
   m = np.sum(qr,axis=1) > 0
   if np.sum(m) > 0:
    qr[m,:] = runoff_hband[m,np.newaxis]*qr[m,:]/np.sum(qr[m,:],axis=1)[:,np.newaxis]
   #Add corresponding qfuture to this time step
   crunoff = self.routing.IRF['qfuture'][:,0]
   #Update qfuture
   self.routing.IRF['qfuture'][:,0:-1] = self.routing.IRF['qfuture'][:,1:]
   self.routing.IRF['qfuture'][:,-1] = 0.0
   #Add current time step
   crunoff = crunoff + qr[:,0]
   #Add the rest to qfuture
   self.routing.IRF['qfuture'] = self.routing.IRF['qfuture'] + qr[:,1:]

   #Aggregate the hband runoff at the reaches
   self.routing.qss[:] += self.routing.reach2hband.dot(crunoff/1000.0)/self.routing.c_length #m2/s
 
   if self.routing_surface_coupling == True:
    #Add changes between Ac1 and Ac0 to qss term
    self.routing.qss[:] += (A - self.routing.A0[:])/self.dt
    #self.routing.A0[:] = A[:]

   #Update routing module
   self.routing.itime = self.itime
   self.routing.update(self.dt)

   #Update the hru inundation values
   for hru in self.hrus:
    #Only allow fct of the inundated height to infiltrate
    self.routing.hru_inundation[hru] = self.routing.hband_inundation[self.hbands[hru]]
   
  return

 def initialize_subsurface(self,):

  if self.subsurface_module == 'richards':self.initialize_richards()

  return

 def initialize_richards(self,):
   
  from model.pyRichards import richards
  
  #Initialize richards
  self.richards = richards.richards(self.nhru,self.nsoil)

  #Set other parameters
  self.richards.dx = self.dx
  self.richards.nhru = self.nhru
  self.richards.m[:] = self.input_fp.groups['parameters'].variables['m'][:]
  #self.richards.dem[:] = self.input_fp.groups['parameters'].variables['dem'][:]
  self.richards.dem[:] = self.input_fp.groups['parameters'].variables['hand'][:]
  self.richards.slope[:] = self.input_fp.groups['parameters'].variables['slope'][:]
  #self.richards.hand[:] = self.input_fp.groups['parameters'].variables['hand'][:]
  self.richards.area[:] = self.input_fp.groups['parameters'].variables['area'][:]
  self.richards.width = sparse.csr_matrix((self.input_fp.groups['wmatrix'].variables['data'][:],
                                  self.input_fp.groups['wmatrix'].variables['indices'][:],
                                  self.input_fp.groups['wmatrix'].variables['indptr'][:]),
                                  shape=(self.nhru,self.nhru),dtype=np.float64)
  #self.richards.width_dense = np.array(self.richards.width.todense())
  #with np.errstate(invalid='ignore',divide='ignore'):tmp = self.richards.width_dense/self.richards.area
  #self.richards.dx = (tmp + tmp.T)/2
  self.richards.I = self.richards.width.copy()
  self.richards.I[self.richards.I != 0] = 1
  self.richards.w = np.array(self.richards.width.todense())
  with np.errstate(invalid='ignore', divide='ignore'):
   tmp = self.richards.area/self.richards.w
  dx = (tmp + tmp.T)/2
  self.richards.dx = dx
  
  return

 def initialize_hwu(self,info):

  #from pyHWU.Human_Water_Use import Human_Water_Use as hwu
  #self.hwu = hwu(self,info)
  #self.hwu.area = self.input_fp.groups['parameters'].variables['area'][:]

  #if self.hwu.hwu_flag == True:
  # print("Initializing Human Water Management")
  # self.hwu.initialize_allocation(self,)

  return

 def run(self,info):

  #Run the model
  date = self.idate
  tic = time.time()
  self.noahmp.dzwt[:] = 0.0
  i=0

  while date < self.fdate:
   i=i+1 
   #Update input data
   tic0 = time.time()
   self.noahmp.stc[:,0:self.noahmp.nsnow]=self.tsno #Laura
   self.update_input(date)
   if i==1: #Laura. Initial conditions for tah and eah as defined in NOAH
    self.noahmp.tah=self.noahmp.t_ml
    self.noahmp.eah=(self.noahmp.p_ml*self.noahmp.q_ml)/(0.622+self.noahmp.q_ml)
    self.noahmp.cm[:] = 0.1
    self.noahmp.ch[:] = 0.1
   #Save the original precip
   precip = np.copy(self.noahmp.prcp)

   #Calculate initial NOAH water balance
   self.initialize_water_balance()

   #Update model
   tic0 = time.time()
   self.update(date)
   #if i>1: self.update(date)

   #print('update model',time.time() - tic0,flush=True)
   
   #Return precip to original value
   #self.noahmp.prcp[:] = precip[:] 

   #Calculate final water balance
   self.finalize_water_balance()

   #Update the water balance error
   self.calculate_water_balance_error()

   #Update the energy balance error
   self.calculate_energy_balance_error()
 
   #Update time and date
   self.date = date
   info['date'] = date
   self.runtime = self.runtime0 + time.time()-tic

   #Update output
   tic0 = time.time()
   self.update_output(date)
   #print('update output',time.time() - tic0,flush=True)

   #Update time step
   date = date + self.dt_timedelta
   self.itime = self.itime + 1

   #Output some statistics
   #if (date.minute ==0):# and (date.hour == 0) and (date.day == 1):
   string = '|%s|%s|%s|%s|%s|%s|%s|' % \
        ('Date:%s' % date.strftime("%Y-%m-%d_%H:%M"),\
         'Runtime:%.2f(s)'%(self.runtime),\
         'Acc_ET:%.2f(mm)'%self.acc_et,\
         'Acc_P:%.2f(mm)'%self.acc_prcp,\
         'Acc_Q:%.2f(mm)'%self.acc_q,\
         'Acc_WBE:%.2f(mm)' % self.acc_errwat,\
         'Acc_EBE:%.2f(J/m2)' % self.acc_erreng)
   print(string,flush=True)
   #print(' '*len(string),flush=True)

  return

 def update_input(self,date):

  self.noahmp.itime = self.itime
  dt = self.dt
  self.noahmp.nowdate = assign_string(19,date.strftime('%Y-%m-%d_%H:%M:%S'))
  self.noahmp.julian = (date - datetime.datetime(date.year,1,1,0)).days
  self.noahmp.year = date.year

  #Update meteorology
  i = int(np.floor(self.itime/self.nt_timestep))
  self.noahmp.lwdn[:] = self.meteodata['lwdown'][i,:] #W/m2
  self.noahmp.swdn[:] = self.meteodata['swdown'][i,:] #W/m2
  self.noahmp.psfc[:] = self.meteodata['psurf'][i,:] #Pa
  self.noahmp.p_ml[:] = self.noahmp.psfc[:] #Pa
  self.noahmp.u_ml[:] = self.meteodata['wind'][i,:] #m/s
  self.noahmp.v_ml[:] = 0.0 #m/s
  self.noahmp.t_ml[:] = self.meteodata['tair'][i,:] #K
  self.noahmp.q_ml[:] = self.meteodata['spfh'][i,:] #Kg/Kg
  self.noahmp.qsfc1d[:] = self.noahmp.q_ml[:] #Kg/Kg
  self.noahmp.prcp[:] = self.meteodata['precip'][i,:] #mm/s

  # Update water demands
  '''if self.hwu.hwu_flag == True:
   if (date.hour*3600)%self.hwu.dta == 0:
    water_use = self.input_fp.groups['water_use']
    if self.hwu.hwu_indust_flag == True:
     self.hwu.demand_indust[:]  = water_use.variables['industrial'][i,:] #m/s
     self.hwu.deficit_indust[:] = np.copy(self.hwu.demand_indust[:])
    if self.hwu.hwu_domest_flag == True:
     self.hwu.demand_domest[:]  = water_use.variables['domestic'][i,:] #m/s
     self.hwu.deficit_domest[:] = np.copy(self.hwu.demand_domest[:])
    if self.hwu.hwu_lstock_flag == True:
     self.hwu.demand_lstock[:]  = water_use.variables['livestock'][i,:] #m/s
     self.i    hwu.deficit_lstock[:] = np.copy(self.hwu.demand_lstock[:])'''


  return

 def update(self,date):

  if self.routing_surface_coupling == True:
   self.noahmp.sfcheadrt[:] = self.routing.fct_infiltrate*self.routing.hru_inundation[:]*1000 #mm
  else:
   self.noahmp.sfcheadrt[:] = 0.0

  # Update subsurface
  self.update_subsurface()
  
  # Update NOAH
  n = self.noahmp
  n.update(n.z_ml,n.dt,n.lwdn,n.swdn,n.u_ml,n.v_ml,n.q_ml,n.t_ml,n.prcp,n.psfc,\
           n.nowdate,n.xlat,n.xlon,n.cosz,n.julian,\
           n.itime,n.year,\
           n.dzs,n.dx,\
           n.ivgtyp,n.isltyp,n.fvegmp,n.fvgmax,n.tmn,\
           n.cropcat,\
           n.planting,n.harvest,n.season_gdd,\
           n.idveg,n.iopt_crs,n.iopt_btr,n.iopt_run,n.iopt_sfc,n.iopt_frz,\
           n.iopt_inf,n.iopt_rad,n.iopt_alb,n.iopt_snf,n.iopt_tbot,n.iopt_stc,\
           n.iopt_gla,n.iopt_rsf,n.iopt_crop,\
           n.iz0tlnd,n.sf_urban_physics,\
           n.tsk,n.hfx,n.qfx,n.lh,n.grdflx,n.smstav,\
           n.smstot,n.sfcrunoff,n.udrunoff,n.albedo,n.snowc,n.smois,\
           n.sh2o,n.tslb,n.snow,n.snowh,n.canwat,n.acsnom,\
           n.acsnow,n.emiss,n.qsfc,\
           n.z0,n.znt,\
           n.isnow,n.tv,n.tg,n.canice,n.canliq,n.eah,\
           n.tah,n.cm,n.ch,n.fwet,n.sneqvo,n.albold,\
           n.qsnow,n.wslake,n.zwt,n.wa,n.wt,n.tsno,\
           n.zsnso,n.snice,n.snliq,n.lfmass,n.rtmass,n.stmass,\
           n.wood,n.stblcp,n.fastcp,n.lai,n.sai,n.tauss,\
           n.smoiseq,n.smcwtd,n.deeprech,n.rech,n.grain,n.gdd,n.pgs,\
           n.gecros_state,\
           n.t2mv,n.t2mb,n.q2mv,n.q2mb,\
           n.trad,n.nee,n.gpp,n.npp,n.fvegmp,n.runsf,\
           n.runsb,n.ecan,n.edir,n.etran,n.fsa,n.fira,\
           n.apar,n.psn,n.sav,n.sag,n.rssun,n.rssha,\
           n.bgap,n.wgap,n.tgv,n.tgb,n.chv,n.chb,\
           n.shg,n.shc,n.shb,n.evg,n.evb,n.ghv,\
           n.ghb,n.irg,n.irc,n.irb,n.tr,n.evc,\
           n.chleaf,n.chuc,n.chv2,n.chb2,n.rs,\
           n.pah,\
           n.bexp,n.smcdry,n.smcwlt,n.smcref,n.smcmax,\
           n.dksat,n.dwsat,n.psisat,n.quartz,\
           n.hdiv,\
           n.sfcheadrt\
          )

  self.tsno=self.noahmp.stc[:,0:self.noahmp.nsnow] #Laura

  # Calculate water demands and supplies, and allocate volumes
  #self.hwu.Calc_Human_Water_Demand_Supply(self,date)

  # Abstract Surface Water and Groundwater
  #self.hwu.Water_Supply_Abstraction(self,date)

  # Update routing
  self.update_routing()

  return

 def update_subsurface(self,):

  self.noahmp.dzwt[:] = 0.0

  if self.subsurface_module == 'richards':

   #Assign noahmp variables to subsurface module
   self.richards.theta[:] = self.noahmp.smois[:]
   self.richards.thetar[:] = self.noahmp.smcdry[:,0]
   self.richards.thetas[:] = self.noahmp.smcmax[:,0]
   self.richards.b[:] = self.noahmp.bexp[:,0]
   self.richards.satpsi[:] = self.noahmp.psisat[:,0]
   self.richards.ksat[:] = self.noahmp.dksat[:,0]
   self.richards.dz[:] = self.noahmp.sldpth[:]

   #Update subsurface module
   #0.Update hand value to account for hru inundation (This is a hack to facilitate a non-flooding stream to influence its surrounding hrus)
   if self.routing_module == 'kinematic':
     self.richards.dem1 = self.richards.dem+self.routing.hru_inundation
     m = self.richards.dem1[0:-1] > self.richards.dem1[1:]
     self.richards.dem1[0:-1][m] = self.richards.dem1[1:][m]
   else:
     self.richards.dem1 = self.richards.dem

   #self.richards.update()
   self.richards.update_numba()

   #Assign subsurface module variables to noahmp
   self.noahmp.hdiv[:] = self.richards.hdiv[:]
   
  return

 def initialize_water_balance(self,):
 
  smw = np.sum(1000*self.noahmp.sldpth*self.noahmp.smois,axis=1)
  self.beg_wb = np.copy(self.noahmp.canliq + self.noahmp.canice + self.noahmp.swe + self.noahmp.wa + smw)
  self.dzwt0 = np.copy(self.noahmp.dzwt)

  return 

 def finalize_water_balance(self,):

  NOAH = self.noahmp
  smw = np.sum(1000*NOAH.sldpth*NOAH.smois,axis=1) 
  self.end_wb = np.copy(NOAH.canliq + NOAH.canice + NOAH.swe + NOAH.wa + smw)
  
  return

 def calculate_water_balance_error(self,):

  NOAH = self.noahmp
  #HWU = self.hwu
  dt = self.dt
  if self.subsurface_module == 'richards':
   tmp = np.copy(self.end_wb - self.beg_wb - NOAH.dt*(NOAH.prcp-NOAH.ecan-
         NOAH.etran-NOAH.edir-NOAH.runsf-NOAH.runsb-np.sum(NOAH.hdiv,axis=1)))
  else:
   tmp = np.copy(self.end_wb - self.beg_wb - NOAH.dt*(NOAH.prcp-NOAH.ecan-
         NOAH.etran-NOAH.edir-NOAH.runsf-NOAH.runsb))
  if self.routing_module == 'kinematic':
   tmp = tmp - np.copy(NOAH.sfcheadrt) + np.copy(self.routing.hru_runoff_inundation)*dt
  self.acc_errwat += np.sum(self.pct*tmp)
  self.acc_q += dt*np.sum(self.pct*NOAH.runsb) + dt*np.sum(self.pct*NOAH.runsf)
  self.acc_et += dt*np.sum(self.pct*(NOAH.ecan + NOAH.etran + NOAH.edir))
  self.acc_etran += dt*np.sum(self.pct*NOAH.etran)
  self.acc_ecan += dt*np.sum(self.pct*NOAH.ecan)
  self.acc_edir += dt*np.sum(self.pct*NOAH.edir)
  self.acc_prcp += dt*np.sum(self.pct*NOAH.prcp)

  return

 def calculate_energy_balance_error(self,):

  NOAH = self.noahmp
  tmp = np.copy(NOAH.sav+NOAH.sag-NOAH.fira-NOAH.hfx-NOAH.lh-NOAH.grdflx+NOAH.pah)
  self.acc_erreng += np.sum(self.dt*self.pct*tmp)

  return

 def update_output(self,date):

  NOAH = self.noahmp
  #HWU = self.hwu
  HB = self
  itime = self.itime

  #Create the netcdf file
  if date == self.idate: self.create_netcdf_file()

  #General info
  grp = self.output_fp.groups['metadata']
  dates = grp.variables['date']
  dates[itime] = nc.date2num(date,units=dates.units,calendar=dates.calendar)
  
  #Update the variables
  grp = self.output_fp.groups['data']

  tmp = {}
  #NoahMP
  tmp['lai'] = np.copy(NOAH.lai) #m2/m2
  tmp['sai'] = np.copy(NOAH.sai) #m2/m2
  tmp['smc'] = np.copy(NOAH.smois) #m3/m3
  tmp['g'] = np.copy(NOAH.grdflx) #W/m2
  tmp['sh'] = np.copy(NOAH.hfx) #W/m2
  tmp['lh'] = np.copy(NOAH.lh) #W/m2
  tmp['lwnet'] = np.copy(NOAH.fira) #W/m2
  tmp['swnet'] = np.copy(NOAH.sav + NOAH.sag) #W/m2
  tmp['shb'] = np.copy(NOAH.shb) #W/m2
  tmp['shc'] = np.copy(NOAH.shc) #W/m2
  tmp['shg'] = np.copy(NOAH.shg) #W/m2
  tmp['evb'] = np.copy(NOAH.evb) #W/m2
  tmp['evc'] = np.copy(NOAH.evc) #W/m2
  tmp['evg'] = np.copy(NOAH.evg) #W/m2
  tmp['ecan'] = np.copy(NOAH.ecan) #W/m2
  tmp['etran'] = np.copy(NOAH.etran) #W/m2
  tmp['edir'] = np.copy(NOAH.edir) #W/m2
  tmp['swdn'] = np.copy(NOAH.swdn) #W/m2
  tmp['lwdn'] = np.copy(NOAH.lwdn) #W/m2
  tmp['t2mv'] = np.copy(NOAH.t2mv) #W/m2
  tmp['t2mb'] = np.copy(NOAH.t2mb) #W/m2
  tmp['qbase'] = NOAH.dt*np.copy(NOAH.runsb) #mm
  tmp['qsurface'] = NOAH.dt*np.copy(NOAH.runsf) #mm
  tmp['udrunoff'] = np.copy(NOAH.udrunoff)
  tmp['sfcrunoff'] = np.copy(NOAH.sfcrunoff)
  tmp['runoff'] = NOAH.dt*(np.copy(NOAH.runsf)+np.copy(NOAH.runsb)) #mm 
  tmp['prcp'] = NOAH.dt*np.copy(NOAH.prcp) #mm
  tmp['trad'] = np.copy(NOAH.trad) #K
  tmp['stc'] = np.copy(NOAH.stc[:,3]) #K
  tmp['tv'] = np.copy(NOAH.tv) #K
  tmp['salb'] = np.copy(NOAH.albedo) 
  tmp['fsa'] = np.copy(NOAH.fsa) #W/m2
  tmp['canliq'] = np.copy(NOAH.canliq) #W/m2
  tmp['canice'] = np.copy(NOAH.canice) #W/m2
  tmp['wa'] = np.copy(NOAH.wa) #W/m2
  tmp['wt'] = np.copy(NOAH.wt) #W/m2
  tmp['sav'] = np.copy(NOAH.sav) #W/m2
  tmp['tr'] = np.copy(NOAH.tr) #W/m2
  tmp['irc'] = np.copy(NOAH.irc) #W/m2
  tmp['irg'] = np.copy(NOAH.irg) #W/m2
  tmp['ghv'] = np.copy(NOAH.ghv) #W/m2
  tmp['sag'] = np.copy(NOAH.sag) #W/m2
  tmp['irb'] = np.copy(NOAH.irb) #W/m2
  tmp['ghb'] = np.copy(NOAH.ghb) #W/m2
  tmp['tg'] = np.copy(NOAH.tg) #K
  tmp['tah'] = np.copy(NOAH.tah) #K
  tmp['tgv'] = np.copy(NOAH.tgv) #K
  tmp['tgb'] = np.copy(NOAH.tgb) #K
  tmp['q2mv'] = np.copy(NOAH.q2mv) #K
  tmp['q2mb'] = np.copy(NOAH.q2mb) #K
  tmp['eah'] = np.copy(NOAH.eah) #K
  tmp['fwet'] = np.copy(NOAH.fwet) #K
  tmp['zsnso_sn'] = np.copy(NOAH.zsnso[:,0:self.noahmp.nsnow]) #m
  tmp['snice'] = np.copy(NOAH.snice) #mm
  tmp['snliq'] = np.copy(NOAH.snliq) #mm
  tmp['soil_t'] = np.copy(NOAH.tslb) #K
  tmp['soil_m'] = np.copy(NOAH.smois) #m3/m3
  tmp['soil_w'] = np.copy(NOAH.sh2o) #m3/m3
  tmp['snow_t'] = np.copy(NOAH.tsno) #K
  tmp['snowh'] = np.copy(NOAH.snowh) #mm
  tmp['qsnow'] = np.copy(NOAH.qsnow) #mm/s
  tmp['isnow'] = np.copy(NOAH.isnow) #count
  tmp['fsno'] = np.copy(NOAH.snowc) #fraction
  #tmp['acsnow'] = np.copy(NOAH.acsnow) #mm (after update)
  tmp['chv'] = np.copy(NOAH.chv) #
  tmp['chb'] = np.copy(NOAH.chb) #
  tmp['chleaf'] = np.copy(NOAH.chleaf) #
  tmp['chuc'] = np.copy(NOAH.chuc) #
  tmp['chv2'] = np.copy(NOAH.chv2) #
  tmp['chb2'] = np.copy(NOAH.chb2) #
  tmp['lfmass'] = np.copy(NOAH.lfmass) #
  tmp['rtmass'] = np.copy(NOAH.rtmass) #
  tmp['stmass'] = np.copy(NOAH.stmass) #
  tmp['wood'] = np.copy(NOAH.wood) #
  #tmp['grain'] = np.copy(NOAH.wood) # (after update)
  #tmp['gdd'] = np.copy(NOAH.wood) # (after update)
  tmp['stblcp'] = np.copy(NOAH.stblcp) #
  tmp['fastcp'] = np.copy(NOAH.fastcp) #
  tmp['nee'] = np.copy(NOAH.nee) #
  tmp['gpp'] = np.copy(NOAH.gpp) #
  tmp['npp'] = np.copy(NOAH.npp) #
  tmp['psn'] = np.copy(NOAH.psn) #
  tmp['apar'] = np.copy(NOAH.apar) #
  tmp['emissi'] = np.copy(NOAH.emiss) 
  tmp['cosz'] = np.copy(NOAH.cosz) 
  tmp['zwt'] = np.copy(NOAH.zwt) #m
  tmp['swe'] = np.copy(NOAH.snow) #m
  tmp['totsmc'] = np.sum(NOAH.sldpth*NOAH.smois,axis=1)/np.sum(NOAH.sldpth[0]) #m3/m3
  tmp['hdiv'] = np.copy(NOAH.hdiv)
  tmp['mozb'] = np.copy(NOAH.mozb)
  tmp['mozv'] = np.copy(NOAH.mozv)
  tmp['fvb'] = np.copy(NOAH.fvb)
  tmp['fvv'] = np.copy(NOAH.fvv)
  tmp['fveg'] = np.copy(NOAH.fvegmp)
  tmp['zpd'] = np.copy(NOAH.zpd)
  tmp['zpdg'] = np.copy(NOAH.zpdg)
  tmp['tauxb'] = np.copy(NOAH.tauxb)
  tmp['tauyb'] = np.copy(NOAH.tauyb)
  tmp['tauxv'] = np.copy(NOAH.tauxv)
  tmp['tauyv'] = np.copy(NOAH.tauyv)
  tmp['cm'] = np.copy(NOAH.cm)
  tmp['ch'] = np.copy(NOAH.ch)
  #routing
  if self.routing_module == 'kinematic':
   tmp['sfcheadrt'] = np.copy(NOAH.sfcheadrt)
   tmp['inundation'] = np.copy(self.routing.hru_inundation)


  # root zone
  cs = np.cumsum(NOAH.sldpth[0,:])
  mask = cs <= 0.5
  pct = NOAH.sldpth[0,mask]/np.sum(NOAH.sldpth[0,mask])
  tmp['smc_root'] = np.sum(pct*NOAH.smois[:,mask],axis=1) #m3/m3
  # top soil layer
  tmp['smc1'] = np.copy(NOAH.smois[:,0])
  
  #General
  tmp['errwat'] = np.copy(HB.acc_errwat)

  # Water Management
  '''if self.hwu.hwu_agric_flag == True:
   tmp['demand_agric'] = np.copy(HWU.demand_agric)*NOAH.dt       #m
   tmp['deficit_agric'] = np.copy(HWU.deficit_agric)*NOAH.dt     #m
   tmp['irrig_agric'] = np.copy(HWU.irrigation)*(NOAH.dt/1000.0)   #m

  if self.hwu.hwu_flag == True:
   if self.hwu.hwu_indust_flag == True:
     tmp['demand_indust'] = np.copy(HWU.demand_indust)*NOAH.dt #m
     tmp['deficit_indust'] = np.copy(HWU.deficit_indust)*NOAH.dt #m
     tmp['alloc_indust'] = np.copy(HWU.demand_indust-HWU.deficit_indust)*NOAH.dt #m
   if self.hwu.hwu_domest_flag == True:
     tmp['demand_domest'] = np.copy(HWU.demand_domest)*NOAH.dt #m
     tmp['deficit_domest'] = np.copy(HWU.deficit_domest)*NOAH.dt #m
     tmp['alloc_domest'] = np.copy(HWU.demand_domest-HWU.deficit_domest)*NOAH.dt #m
   if self.hwu.hwu_lstock_flag == True:
     tmp['demand_lstock'] = np.copy(HWU.demand_lstock)*NOAH.dt #m
     tmp['deficit_lstock'] = np.copy(HWU.deficit_lstock)*NOAH.dt #m
     tmp['alloc_lstock'] = np.copy(HWU.demand_lstock-HWU.deficit_lstock)*NOAH.dt #m
   if self.hwu.hwu_sf_flag == True:
     tmp['alloc_sf'] = np.copy(HWU.alloc_sf) #m
   if self.hwu.hwu_gw_flag == True:                                                                                              tmp['alloc_gw'] = np.copy(HWU.alloc_gw) #m'''

  #Output the variables
  sep = 100
  if itime == 0:
   self.output = {}
   for var in self.metadata['output']['vars']:
    shp = grp.variables[var].shape
    if len(shp) == 2:self.output[var] = np.zeros((sep,shp[1]))
    if len(shp) == 3:self.output[var] = np.zeros((sep,shp[1],shp[2]))
  #Fill the data (MISSING!)
  val = itime % sep
  for var in self.metadata['output']['vars']:
    self.output[var][val,:] = tmp[var]
  # self.output[itime] =
  if (itime+1) % sep == 0:
   for var in self.metadata['output']['vars']:
    grp.variables[var][itime-sep+1:itime+1,:] = self.output[var][:]
  if (itime+1) == self.ntime:
   val = int(sep*np.ceil((itime - sep)/sep))
   for var in self.metadata['output']['vars']:
    grp.variables[var][val:itime+1,:] = self.output[var][0:itime-val+1,:]

  #Output routing variables
  if self.routing_module == 'kinematic':
   grp = self.output_fp.groups['data_routing']
   tmp = {}
   tmp['A'] = self.routing.A1[:]
   tmp['Q'] = self.routing.Q1[:]
   tmp['Qc'] = self.routing.Qc[:]
   tmp['Qf'] = self.routing.Qf[:]
   tmp['reach_inundation'] = self.routing.reach2hband_inundation[:]
   sep = 100
   if itime == 0:
    self.output_routing = {}
    for var in self.metadata['output']['routing_vars']:
     shp = grp.variables[var].shape
     self.output_routing[var] = np.zeros((sep,shp[1]))
   #Fill the data (MISSING!)
   val = itime % sep
   for var in self.metadata['output']['routing_vars']:
     self.output_routing[var][val,:] = tmp[var]
   # self.output[itime] =
   if (itime+1) % sep == 0:
    for var in self.metadata['output']['routing_vars']:
     grp.variables[var][itime-sep+1:itime+1,:] = self.output_routing[var][:]
   if (itime+1) == self.ntime:
    val = int(sep*np.ceil((itime - sep)/sep))
    for var in self.metadata['output']['routing_vars']:
     grp.variables[var][val:itime+1,:] = self.output_routing[var][0:itime-val+1,:]

  return

 def create_netcdf_file(self,):

  #Create the output directory if necessary
  os.system('mkdir -p %s' % self.metadata['output']['dir'])

  #Extract the pointers to both input and output files
  ofile = '%s/%s.nc' % (self.metadata['output']['dir'],self.idate.strftime('%Y-%m-%d'))
  self.output_fp = nc.Dataset(ofile,'w',format='NETCDF4')
  fp_out = self.output_fp
  fp_in = self.input_fp

  #Define the metadata
  metadata = {
             "ecan":{'description':'Vegetated canopy latent heat flux','units':'mm/s','dims':('time','hru',),'precision':8},
             "ecan":{'description':'Vegetated canopy latent heat flux','units':'mm/s','dims':('time','hru',),'precision':8},
             "etran":{'description':'Vegetated canopy latent heat flux','units':'mm/s','dims':('time','hru',),'precision':8},
             "edir":{'description':'Vegetated canopy latent heat flux','units':'mm/s','dims':('time','hru',),'precision':8},
             'lai':{'description':'Leaf area index','units':'m2/m2','dims':('time','hru',),'precision':8},
             'sai':{'description':'Steam area index','units':'m2/m2','dims':('time','hru',),'precision':8},
             'g':{'description':'Ground heat flux','units':'W/m2','dims':('time','hru',),'precision':4},
             'sh':{'description':'Sensible heat flux','units':'W/m2','dims':('time','hru',),'precision':4},
             'lh':{'description':'Latent heat flux','units':'W/m2','dims':('time','hru',),'precision':4},
             'lwnet':{'description':'Net longwave radiation','units':'W/m2','dims':('time','hru',),'precision':4},
             'swnet':{'description':'Absorbed shortwave radiation','units':'W/m2','dims':('time','hru',),'precision':4},
             "t2mv":{'description':'Vegetated air temperature','units':'K','dims':('time','hru',),'precision':4},
	     "t2mb":{'description':'Bare air temperature','units':'K','dims':('time','hru',),'precision':4},
	     "fveg":{'description':'Vegetated fraction','units':'','dims':('time','hru',),'precision':8},
	     "mozb":{'description':'Bare stability parameter','units':'','dims':('time','hru',),'precision':4},
	     "mozv":{'description':'Vegetated stability parameter','units':'','dims':('time','hru',),'precision':4},
	     "zpd":{'description':'Vegetated zero plane displacement','units':'','dims':('time','hru',),'precision':4},
	     "zpdg":{'description':'Ground zero plane displacement','units':'','dims':('time','hru',),'precision':4},
             "cm":{'description':'Momentum drag coefficient','units':'','dims':('time','hru',),'precision':4},
             "ch":{'description':'Sensible heat exchange coefficient','units':'','dims':('time','hru',),'precision':4},
	     "tauxb":{'description':'Bare wind stress (e-w)','units':'','dims':('time','hru',),'precision':4},
	     "tauyb":{'description':'Bare wind stress (n-s)','units':'','dims':('time','hru',),'precision':4},
	     "tauxv":{'description':'Bare wind stress (e-w)','units':'','dims':('time','hru',),'precision':4},
	     "tauyv":{'description':'Bare wind stress (n-s)','units':'','dims':('time','hru',),'precision':4},
	     "fvb":{'description':'Bare friction velocity','units':'m/s','dims':('time','hru',),'precision':4},
	     "fvv":{'description':'Vegetated friction velocity','units':'m/s','dims':('time','hru',),'precision':4},
             "shb":{'description':'Bare sensible heat flux','units':'W/m2','dims':('time','hru',),'precision':4},
             "shc":{'description':'Vegetated canopy sensible heat flux','units':'W/m2','dims':('time','hru',),'precision':4},
             "shg":{'description':'Vegetated ground sensible heat flux','units':'W/m2','dims':('time','hru',),'precision':4},
             "evb":{'description':'Bare latent heat flux','units':'W/m2','dims':('time','hru',),'precision':4},
             "evc":{'description':'Vegetated canopy latent heat flux','units':'W/m2','dims':('time','hru',),'precision':4},
             "evg":{'description':'Vegetated ground latent heat flux','units':'W/m2','dims':('time','hru',),'precision':4},
             "swdn":{'description':'Shortwave down','units':'W/m2','dims':('time','hru',),'precision':4},
             "lwdn":{'description':'Longwave down','units':'W/m2','dims':('time','hru',),'precision':4},
             'trad':{'description':'Land surface skin temperature','units':'K','dims':('time','hru',),'precision':2},
             'stc':{'description': 'Snow/soil temperature','units':'K','dims':('time','hru',),'precision':2},
             'tv':{'description': 'Canopy temperature','units':'K','dims':('time','hru',),'precision':2},
             'salb':{'description':'Land surface albedo','units':' ','dims':('time','hru',),'precision':8},
             'fsa':{'description':'Total absorbed solar radiation','units':'W/m2','dims':('time','hru',),'precision':4},
             'canliq':{'description':'Canopy liquid water content','units':'mm','dims':('time','hru',),'precision':4},
             'canice':{'description':'Canopy ice water content','units':'mm','dims':('time','hru',),'precision':4},
             'wa':{'description':'Water in aquifer','units':'kg/m2','dims':('time','hru',),'precision':4},
             'wt':{'description':'Water in aquifer and saturated soil','units':'kg/m2','dims':('time','hru',),'precision':4},
             'sav':{'description':'Solar radiative heat flux absorbed by veg','units':'W/m2','dims':('time','hru',),'precision':4},
             'tr':{'description':'Transpiration heat','units':'W/m2','dims':('time','hru',),'precision':4},
             'irc':{'description':'Canopy net LW rad','units':'W/m2','dims':('time','hru',),'precision':4},
             'irg':{'description':'Ground net LW rad','units':'W/m2','dims':('time','hru',),'precision':4},
             'ghv':{'description':'Ground heat flux + to soil vegetated','units':'W/m2','dims':('time','hru',),'precision':4},
             'sag':{'description':'Solar radiative heat flux absorbed by ground','units':'W/m2','dims':('time','hru',),'precision':4},
             'irb':{'description':'Net LW rad to atm bare','units':'W/m2','dims':('time','hru',),'precision':4},
             'ghb':{'description':'Ground heat flux + to soil bare"','units':'W/m2','dims':('time','hru',),'precision':4},
             'tg':{'description':'Ground temperature','units':'K','dims':('time','hru',),'precision':4},
             'tah':{'description':'Canopy air temperature','units':'K','dims':('time','hru',),'precision':4},
             'tgv':{'description':'Ground surface Temp vegetated','units':'K','dims':('time','hru',),'precision':4},
             'tgb':{'description':'Ground surface Temp bare','units':'K','dims':('time','hru',),'precision':4},
             'q2mv':{'description':'2m mixing ratio vegetated','units':'kg/kg','dims':('time','hru',),'precision':4},
             'q2mb':{'description':'2m mixing ratio bare','units':'kg/kg','dims':('time','hru',),'precision':4},
             'eah':{'description':'Canopy air vapor pressure','units':'Pa','dims':('time','hru',),'precision':4},
             'fwet':{'description':'Wetted or snowed fraction of canopy','units':'fraction','dims':('time','hru'),'precision':4},
             'zsnso_sn':{'description':'Snow layer depths from snow surface','units':'m','dims':('time','hru','snow'),'precision':4},
             'snice':{'description':'Snow layer ice','units':'mm','dims':('time','hru','snow'),'precision':4},
             'snliq':{'description':'Snow layer liquid water','units':'mm','dims':('time','hru','snow'),'precision':4},
             'soil_t':{'description':'soil temperature','units':'K','dims':('time','hru','soil'),'precision':4},
             'soil_m':{'description':'volumetric soil moisture','units':'m3/m3','dims':('time','hru','soil'),'precision':4},
             'soil_w':{'description':'liquid volumetric soil moisture','units':'m3/m3','dims':('time','hru','soil'),'precision':4},
             'snow_t':{'description':'snow temperature','units':'K','dims':('time','hru','snow'),'precision':4},
             'snowh':{'description':'Snow depth','units':'m','dims':('time','hru',),'precision':4},
             'qsnow':{'description':'snowfall rate','units':'mm/s','dims':('time','hru',),'precision':8},
             'isnow':{'description':'Number of snow layers','units':'count','dims':('time','hru',),'precision':4},
             'fsno':{'description':'Snow-cover fraction on the ground','units':'fraction','dims':('time','hru',),'precision':4},
             'acsnow':{'description':'accumulated snow fall','units':'mm','dims':('time','hru',),'precision':4},
             'chv':{'description':'Exchange coefficient vegetated','units':'m/s','dims':('time','hru',),'precision':4},
             'chb':{'description':'Exchange coefficient bare','units':'m/s','dims':('time','hru',),'precision':4},
             'chleaf':{'description':'Exchange coefficient leaf','units':'m/s','dims':('time','hru',),'precision':8},
             'chuc':{'description':'Exchange coefficient bare','units':'m/s','dims':('time','hru',),'precision':4},
             'chv2':{'description':'Exchange coefficient 2-meter vegetated','units':'m/s','dims':('time','hru',),'precision':4},
             'chb2':{'description':'Exchange coefficient 2-meter bare','units':'m/s','dims':('time','hru',),'precision':4},
             'lfmass':{'description':'Leaf mass','units':'g/m2','dims':('time','hru',),'precision':4},
             'rtmass':{'description':'Mass of fine roots','units':'g/m2','dims':('time','hru',),'precision':4},
             'stmass':{'description':'Stem mass','units':'g/m2','dims':('time','hru',),'precision':8},
             'wood':{'description':'Mass of wood and woody roots','units':'g/m2','dims':('time','hru',),'precision':4},
             'grain':{'description':'Mass of grain','units':'g/m2','dims':('time','hru',),'precision':4},
             'gdd':{'description':'Growing degree days(10)','units':'','dims':('time','hru',),'precision':4},
             'stblcp':{'description':'Stable carbon in deep soil','units':'g/m2','dims':('time','hru',),'precision':4},
             'fastcp':{'description':'Short-lived carbon in shallow soil','units':'g/m2','dims':('time','hru',),'precision':4},
             'nee':{'description':'Net ecosystem exchange','units':'g/m2/s CO2','dims':('time','hru',),'precision':8},
             'gpp':{'description':'Net instantaneous assimilation','units':'g/m2/s C','dims':('time','hru',),'precision':8},
             'npp':{'description':'Net primary productivity','units':'g/m2/s C','dims':('time','hru',),'precision':8},
             'psn':{'description':'Total photosynthesis','units':'umol/m2/s','dims':('time','hru',),'precision':4},
             'apar':{'description':'Photosynthesis active energy by canopy','units':'W/m2','dims':('time','hru',),'precision':4},
             'emissi':{'description':'Land surface albedo','units':' ','dims':('time','hru',),'precision':4},
             'cosz':{'description':'Land surface albedo','units':' ','dims':('time','hru',),'precision':4},
             'qbase':{'description':'Excess runoff','units':'mm/s','dims':('time','hru',),'precision':4},
             'udrunoff':{'description':'Accumulated baseflow','units':'mm','dims':('time','hru',),'precision':4},
             'sfcrunoff':{'description':'Accumulated surface','units':'mm','dims':('time','hru',),'precision':4},
             'qsurface':{'description':'Surface runoff','units':'mm/s','dims':('time','hru',),'precision':4},
             'runoff':{'description':'Total runoff','units':'mm/s','dims':('time','hru',),'precision':4},
             'prcp':{'description':'Precipitation','units':'mm/s','dims':('time','hru',),'precision':4},
             'smc':{'description':'Soil water content','units':'m3/m3','dims':('time','hru','soil'),'precision':3},
             'swd':{'description':'Soil water deficit','units':'mm','dims':('time','hru',),'precision':4},
             'sstorage':{'description':'Surface storage','units':'mm','dims':('time','hru',),'precision':4},
             'swe':{'description':'Snow water equivalent','units':'mm','dims':('time','hru',),'precision':4},
             'qout_subsurface':{'description':'Subsurface flux','units':'m2/s','dims':('time','hru',),'precision':4},
             'qout_surface':{'description':'Surface flux','units':'m2/s','dims':('time','hru',),'precision':4},
             'zwt':{'description':'WTD','units':'m','dims':('time','hru',),'precision':2},
             'errwat':{'description':'errwat','units':'mm','dims':('time','hru',),'precision':4},
             'totsmc':{'description':'totsmc','units':'m3/m3','dims':('time','hru',),'precision':2},
             'hdiv':{'description':'hdiv','units':'mm/s','dims':('time','hru','soil'),'precision':4},

             'smc1':{'description':'Soil water content at the root zone','units':'m3/m3','dims':('time','hru',),'precision':3},
             'smc_root':{'description':'Soil water content at the root zone','units':'m3/m3','dims':('time','hru',),'precision':3},  

             # Water Management
             'demand_agric':{'description':'Irrigation water demand','units':'m','dims':('time','hru',),'precision':4},
             'deficit_agric':{'description':'Irrigation water deficit','units':'m','dims':('time','hru',),'precision':4},
             'irrig_agric':{'description':'Irrigated water volume','units':'m','dims':('time','hru',),'precision':4},
             'demand_indust':{'description':'Industrial water demand','units':'m','dims':('time','hru',)},'precision':4,
             'deficit_indust':{'description':'Industrial water deficit','units':'m','dims':('time','hru',),'precision':4},
             'alloc_indust':{'description':'Industrial water allocated','units':'m','dims':('time','hru',),'precision':4},
             'demand_domest':{'description':'Domestic demand','units':'m','dims':('time','hru',),'precision':4},
             'deficit_domest':{'description':'Domestic deficit','units':'m','dims':('time','hru',),'precision':4},
             'alloc_domest':{'description':'Domestic water allocated','units':'m','dims':('time','hru',),'precision':4},
             'demand_lstock':{'description':'Livestock demand','units':'m','dims':('time','hru',),'precision':4},
             'deficit_lstock':{'description':'Livestock deficit','units':'m','dims':('time','hru',),'precision':4},
             'alloc_lstock':{'description':'Livestock water allocated','units':'m','dims':('time','hru',),'precision':4},
             'alloc_sf':{'description':'Surface water allocated','units':'m','dims':('time','hru',),'precision':4},
             'alloc_gw':{'description':'Groundwater water allocated','units':'m','dims':('time','hru',),'precision':4},
             'Q':{'description':'Discharge','units':'m3/s','dims':('time','channel',),'precision':4},
             'Qc':{'description':'Discharge (channel)','units':'m3/s','dims':('time','channel',),'precision':4},
             'Qf':{'description':'Discharge (floodplain)','units':'m3/s','dims':('time','channel',),'precision':4},
             'reach_inundation':{'description':'Inundation height (reach level)','units':'m','dims':('time','channel','hband'),'precision':3},
             'A':{'description':'Cross section','units':'m2','dims':('time','channel',),'precision':4},
             'inundation':{'description':'Inundation height','units':'m','dims':('time','hru',),'precision':4},
             }

  #Create the dimensions
  #print('Creating the dimensions',flush=True)
  ntime = 24*3600*((self.fdate - self.idate).days)/self.dt
  self.ntime = ntime
  nhru = len(fp_in.dimensions['hru'])
  fp_out.createDimension('hru',nhru)
  fp_out.createDimension('time',ntime)
  fp_out.createDimension('soil',self.nsoil)
  fp_out.createDimension('snow',self.noahmp.nsnow)
  fp_out.createDimension('hband',self.nhband)
  if self.routing_module == 'kinematic':
   fp_out.createDimension('channel',self.routing.nchannel)

  #Create the output
  #print('Creating the data group',flush=True)
  grp = fp_out.createGroup('data')
  for var in self.metadata['output']['vars']:
   ncvar = grp.createVariable(var,'f4',metadata[var]['dims'],least_significant_digit=metadata[var]['precision'],fill_value=9999999999.0)#,zlib=True)
   ncvar.description = metadata[var]['description']
   ncvar.units = metadata[var]['units']

  #Create the routing output
  if self.routing_module == 'kinematic':
   #print('Creating the routing group',flush=True)
   grp = fp_out.createGroup('data_routing')
   for var in self.metadata['output']['routing_vars']:
    ncvar = grp.createVariable(var,'f4',metadata[var]['dims'],least_significant_digit=metadata[var]['precision'])#,zlib=True)
    ncvar.description = metadata[var]['description']
    ncvar.units = metadata[var]['units']


  #Create the metadata
  #print('Creating the metadata group',flush=True)
  grp = fp_out.createGroup('metadata')
  #Time
  grp.createVariable('time','f8',('time',))
  dates = grp.createVariable('date','f8',('time',))
  dates.units = 'hours since 1900-01-01'
  dates.calendar = 'standard'

  #HRU percentage coverage
  #print('Setting the HRU percentage coverage',flush=True)
  pcts = grp.createVariable('pct','f4',('hru',))
  pcts[:] = fp_in.groups['parameters'].variables['area_pct'][:]
  pcts.description = 'hru percentage coverage'
  pcts.units = '%'

  #HRU area
  #print('Setting the HRU areal coverage',flush=True)
  area = grp.createVariable('area','f4',('hru',))
  area[:] = fp_in.groups['parameters'].variables['area'][:]
  area.units = 'meters squared'
  #print('Defining the HRU ids',flush=True)
  hru = grp.createVariable('hru','i4',('hru',))
  hrus =[]
  for value in range(nhru):hrus.append(value)
  hru[:] = np.array(hrus)
  hru.description = 'hru ids'

  return

 def finalize(self,):

  #Create the restart directory if necessary
  os.system('mkdir -p %s' % self.metadata['restart']['dir'])

  #Save the restart file
  file_restart = '%s/%s.h5' % (self.metadata['restart']['dir'],self.fdate.strftime('%Y-%m-%d'))
  fp = h5py.File(file_restart,'w')
  #diagnostics
  fp['acc_et'] = self.acc_et
  fp['acc_etran'] = self.acc_etran
  fp['acc_edir'] = self.acc_edir
  fp['acc_ecan'] = self.acc_ecan
  fp['acc_prcp'] = self.acc_prcp
  fp['acc_q'] = self.acc_q
  fp['acc_errwat'] = self.acc_errwat
  fp['acc_erreng'] = self.acc_erreng
  fp['runtime'] = self.runtime
  #noahmp
  fp['smceq'] = self.noahmp.smceq[:]
  fp['albold'] = self.noahmp.albold[:]
  fp['albedo'] = self.noahmp.albedo[:]
  fp['sneqvo'] = self.noahmp.sneqvo[:]
  fp['tslb'] = self.noahmp.tslb[:]
  fp['sh2o'] = self.noahmp.sh2o[:]
  fp['smois'] = self.noahmp.smois[:]
  fp['tah'] = self.noahmp.tah[:]
  fp['eah'] = self.noahmp.eah[:]
  fp['fwet'] = self.noahmp.fwet[:]
  fp['canliq'] = self.noahmp.canliq[:]
  fp['canice'] = self.noahmp.canice[:]
  fp['tv'] = self.noahmp.tv[:]
  fp['tg'] = self.noahmp.tg[:]
  fp['qsfc1d'] = self.noahmp.qsfc1d[:]
  fp['qsnow'] = self.noahmp.qsnow[:]
  fp['isnow'] = self.noahmp.isnow[:]
  fp['zsnso'] = self.noahmp.zsnso[:]
  fp['snowh'] = self.noahmp.snowh[:]
  fp['snow'] = self.noahmp.snow[:]
  fp['snice'] = self.noahmp.snice[:]
  fp['snliq'] = self.noahmp.snliq[:]
  fp['zwt'] = self.noahmp.zwt[:]
  fp['wa'] = self.noahmp.wa[:]
  fp['wt'] = self.noahmp.wt[:]
  fp['wslake'] = self.noahmp.wslake[:]
  fp['lfmass'] = self.noahmp.lfmass[:]
  fp['rtmass'] = self.noahmp.rtmass[:]
  fp['stmass'] = self.noahmp.stmass[:]
  fp['wood'] = self.noahmp.wood[:]
  fp['stblcp'] = self.noahmp.stblcp[:]
  fp['fastcp'] = self.noahmp.fastcp[:]
  fp['lai'] = self.noahmp.lai[:]
  fp['sai'] = self.noahmp.sai[:]
  fp['cm'] = self.noahmp.cm[:]
  fp['ch'] = self.noahmp.ch[:]
  fp['tauss'] = self.noahmp.tauss[:]
  fp['smcwtd'] = self.noahmp.smcwtd[:]
  fp['deeprech'] = self.noahmp.deeprech[:]
  fp['rech'] = self.noahmp.rech[:]
  fp['tsno'] = self.noahmp.tsno[:]
  fp['grain'] = self.noahmp.grain[:]
  fp['gdd'] = self.noahmp.gdd[:]
  fp['fveg'] = self.noahmp.fveg[:]
  fp['fvgmax'] = self.noahmp.fvgmax[:]
  fp['acsnom'] = self.noahmp.acsnom[:]
  fp['acsnow'] = self.noahmp.acsnow[:]
  fp['qsfc'] = self.noahmp.qsfc[:]
  fp['sfcrunoff'] = self.noahmp.sfcrunoff[:]
  fp['udrunoff'] = self.noahmp.udrunoff[:]
  #routing
  if self.routing_module == 'kinematic':	
   fp['Q0'] = self.routing.Q0[:]
   fp['u0'] = self.routing.u0[:]
   fp['A0'] = self.routing.A0[:]
   fp['Q1'] = self.routing.Q1[:]
   fp['A1'] = self.routing.A1[:]
   fp['dA'] = self.routing.dA[:]
   fp['bcs'] = self.routing.bcs[:]
   fp['qss'] = self.routing.qss[:]
   #fp['hru_inundation'] = self.routing.hru_inundation[:]
   fp['hband_inundation'] = self.routing.hband_inundation[:]
   fp['reach2hband_inundation'] = self.routing.reach2hband_inundation[:]
  fp.close()
   
  #Close the LSM
  del self.noahmp
 
  #Close the routing module
  if self.routing_module == 'kinematic':
   del self.routing

  #Close the files
  self.input_fp.close()
  self.output_fp.close()

  return
