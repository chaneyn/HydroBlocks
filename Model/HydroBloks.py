import os
import numpy as np
import datetime
import time
import sys
sys.path.append('pyNoahMP')
sys.path.append('pyDTopmodel')
import scipy.sparse as sparse

def Finalize_Model(NOAH,TOPMODEL):

 #Deallocate the NOAH memory
 NOAH.finalize()

 #Delete the objects
 del NOAH
 del TOPMODEL

 return

def Initialize_Model(ncells,dt,nsoil,data,parameters,info,wbd):

 from NoahMP import model

 model.ncells = ncells
 model.nsoil = nsoil
 model.dt = dt
 model.nsnow = 3
 model.llanduse[:] = 'MODIFIED_IGBP_MODIS_NOAH'
 model.lsoil[:] = 'STAS'
 model.vegparm_file[:] = info['VEGPARM']#'data/VEGPARM.TBL'
 model.genparm_file[:] = info['GENPARM']#'data/GENPARM.TBL'
 model.soilparm_file[:] = info['SOILPARM']
 model.mptable_file[:] = info['MPTABLE']#'pyNoahMP/data/MPTABLE.TBL'
 #Read in the soil parameter file
 fp = open(info['SOILPARM'])
 iline = 0
 soils_data = {'MAXSMC':[],'DRYSMC':[],'REFSMC':[]}
 for line in fp:
  if (iline > 2) & (iline < 15):
   tmp = line.split(',')
   soils_data['MAXSMC'].append(float(tmp[4]))
   soils_data['DRYSMC'].append(float(tmp[2]))
   soils_data['REFSMC'].append(float(tmp[5]))
  iline = iline + 1
 fp.close()
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
 psoil = 10**parameters['log10soil']
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
 for hsu in data['hsu']:
  ihsu = data['hsu'].keys().index(hsu)
  model.vegtyp[ihsu] = data['hsu'][hsu]['land_cover'] #HERE
  model.soiltyp[ihsu] = data['hsu'][hsu]['soil_texture_class']#1#ihsu+1 #HERE
  model.smcmax[ihsu] = soils_data['MAXSMC'][model.soiltyp[ihsu]-1]
  model.sh2o[ihsu] = soils_data['MAXSMC'][model.soiltyp[ihsu]-1]
  model.smcref[ihsu] = soils_data['REFSMC'][model.soiltyp[ihsu]-1]
  model.smcdry[ihsu] = soils_data['DRYSMC'][model.soiltyp[ihsu]-1]
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

def Initialize_DTopmodel(ncells,dt,data,parameters):

 #sys.path.append("topmodel")
 import dynamic_topmodel

 #Initialize Dynamic Topmodel
 model = dynamic_topmodel.Dynamic_Topmodel(ncells)

 #Set model parameters
 model.dt = dt #seconds
 model.area[:] = 0.0 #meters^2
 model.dx[:] = 30.0 #meters
 model.m[:] = 0.01#1.0#10**parameters['log10m'] #HERE
 model.sdmax[:] = 1000.0#parameters['sdmax'] #HERE
 dem = []
 #Set cluster information
 for hsu in data['hsu']:
  ihsu = data['hsu'].keys().index(hsu)
  model.pct[ihsu] = data['hsu'][hsu]['area_pct']/100
  model.area[ihsu] = model.dx[ihsu]*data['hsu'][hsu]['area']
  model.T0[ihsu] = data['hsu'][hsu]['soil_parameters']['SATDK']
  model.sti[ihsu] = data['hsu'][hsu]['ti']
  model.beta[ihsu] = data['hsu'][hsu]['slope']
  model.carea[ihsu] = data['hsu'][hsu]['carea']
  model.channel[ihsu] = data['hsu'][hsu]['channel']
  model.surface_velocity[ihsu] = 1000.0/3600.0 #m/s
  dem.append(data['hsu'][hsu]['dem'])
 model.dem[:] = np.array(dem)
 model.sdmax[:] = 0.1#np.array(dem) - np.min(dem)
 model.pct = model.pct/np.sum(model.pct)
 ti_mean = np.sum(model.pct*model.sti[:])
 lnTe = parameters['lnTe']
 
 #Calculate the sti
 model.T0[:] = (model.dem - np.min(model.dem))*model.T0
 lnT0 = np.log(model.T0)
 lnTe = np.sum(model.pct*lnT0)
 model.sti = model.sti - (lnT0 - lnTe)

 #Set weight matrix
 model.outlet_hsu = int(data['outlet']['hsu'])
 hsu = model.outlet_hsu
 data['tp'][hsu,hsu] = data['tp'][model.outlet_hsu,hsu] - model.dx[hsu]**2/model.area[hsu]
 model.w = sparse.csr_matrix(data['tp'].T,dtype=np.float32)
 model.wfull[:] = data['tp'].T
 #Initialize the soil moisture deficit values
 model.si[:] = 0.0
 
 return model

def Update_Model(NOAH,TOPMODEL,ncores):
 
 #Set the partial pressure of CO2 and O2
 NOAH.co2air = 355.E-6*NOAH.psfc# ! Partial pressure of CO2 (Pa) ! From NOAH-MP-WRF
 NOAH.o2air = 0.209*NOAH.psfc# ! Partial pressure of O2 (Pa)  ! From NOAH-MP-WRF

 #Update NOAH
 NOAH.minzwt[:] = -0.25#-0.1*((TOPMODEL.dem - np.min(TOPMODEL.dem))+0.5)
 ntt = 1
 NOAH.dzwt[:] = NOAH.dzwt[:]/ntt
 dt = np.copy(NOAH.dt)
 NOAH.dt = dt/ntt
 for i in xrange(ntt):
  NOAH.run_model(ncores)
 NOAH.dt = np.copy(dt)

 #Calculate the updated soil moisture deficit
 si0 = np.copy(NOAH.si0)
 si1 = np.copy(NOAH.si1)

 #Calculate the change in deficit
 TOPMODEL.si[:] = si1
 TOPMODEL.dsi[:] = np.copy(si1 - si0)
 TOPMODEL.r[:] = -TOPMODEL.dsi[:]/TOPMODEL.dt

 #Add the surface runoff
 TOPMODEL.qsurf[:] = NOAH.runsf/1000.0

 #Update dynamic topmodel
 TOPMODEL.update(ncores)

 #Calculate the change in deficit
 TOPMODEL.sideep = TOPMODEL.sideep.astype(np.float32)
 TOPMODEL.si = TOPMODEL.si.astype(np.float32)
 dsi = np.copy(si1 - TOPMODEL.si)

 #Update the soil moisture values
 NOAH.dzwt[:] = np.copy(dsi+TOPMODEL.dt*TOPMODEL.ex-TOPMODEL.dt*TOPMODEL.r)
 TOPMODEL.ex[:] = 0.0

 return (NOAH,TOPMODEL)
