import numpy as np

def Calculate_Field_Capacity_Deficit(smcref,smcwtd,zwt,sldpth,smc,nsoil):
 """
 Calculates the soil water deficit to have soil moisture in field capacity in meters
 """
 
 #print smcref,smcwtd,zwt,sldpth,smc,nsoil
 sldpth=sldpth[0]
 smc=smc[0]
 sif = 0.0

 # Substitute this by the root depth afterwards
 root_depht=2  # root layers
 noah_std_root_profile = [0.1,0.3,0.6,1.0]  # Depth in meters
 total_root_depht = np.sum(noah_std_root_profile[:root_depht]) # Meters
 total_soil_depht = [np.sum(sldpth[:i+1]) for i in range(nsoil)] # Meters
 nroot_depht = np.where(total_soil_depht<=total_root_depht)[0][-1]  # Soil Layers

 #Iterate though all soil layers
 for isoil in range(nroot_depht+1):
    #Calculate the empty space
    dsm=1.0*smcref-smc[isoil]
    if dsm > 0. : sif = sif + dsm*sldpth[isoil]

 #Add in the transmission zone deficit
 #if abs(zwt) > sum(sldpth) :
 #    dsm=smcref - smcwtd
 #    if dsm > 0. : sif = sif + (abs(zwt) - sum(sldpth))*dsm

 return sif


