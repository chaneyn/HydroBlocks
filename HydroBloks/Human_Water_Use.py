import numpy as np
import os

class Human_Water_Use:

 def __init__(self,NOAH,ncells):
  self.ncells = ncells
  self.itime = 0
  self.well_depth = -3.0
  self.irrig_supply = np.zeros(ncells,dtype=np.float64)
  self.irrig_demand = np.zeros(ncells,dtype=np.float64)
  self.nroot_zone   = np.zeros(ncells)
    
 def Calculate_Irrigation_Deficit(self,smcref,sldpth,smc,nroot_zone_depht,ncells):
  """
  Calculates the irrigation deficit based on how much water [m] is necessary to 
  soil moisture reach field capacity
  """
  
  irrig_demand = np.zeros(ncells,dtype=np.float64)
  dsm = smcref - smc
  dsm [dsm < 0] = 0.0
 
  irrig_demand = [ np.sum(dsm[i,:nroot_zone_depht[i]]*sldpth[i,:nroot_zone_depht[i]])  for i in range(ncells) ]
  return irrig_demand


 def irrigation(self,NOAH):
   
   #Calculate irrigation demand
   self.irrig_demand[:] = self.Calculate_Irrigation_Deficit(NOAH.smcref,NOAH.sldpth,NOAH.sh2o,NOAH.root_depth,NOAH.ncells)
   self.irrig_supply[:] = 0.0
   
   # Calculate how much of the demand can actually be extracted
   self.irrig_supply[ self.irrig_demand > 0.0 and NOAH.zwt > self.well_depth ] = self.irrig_demand
   self.irrig_supply[ self.irrig_supply > np.abs(NOAH.zwt-self.well_depth) ] = np.abs(NOAH.zwt-self.well_depth)

   # Abstract water from the water table and add to precip
   NOAH.dzwt[ self.irrig_supply > 0 ] = NOAH.dzwt -self.irrig_supply
   NOAH.prcp[ self.irrig_supply > 0 ] = NOAH.prcp +self.irrig_supply*1000./NOAH.dt

   #print 'demand',self.irrig_demand
   #print 'supply',self.irrig_supply
   return

