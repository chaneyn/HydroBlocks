import numpy as np
import os
from datetime import datetime
from scipy.sparse import csc_matrix, csr_matrix
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)))#+'/pyHWU/')
import management_funcs as mgmt_funcs
#import matplotlib.pyplot as plt


class Human_Water_Use:

 def __init__(self,HB,info):

  NOAH = HB.noahmp
       
  ncells = NOAH.ncells
  self.itime = 0
  self.area = np.zeros(ncells)

  # tstep for water allocation 
  self.dta = 3600*24.0 #HB.dt #Update dta  -- solving water demand daily 
  self.ntt = int(np.copy(self.dta/HB.dt))

  self.hwu_flag = info['water_management']['hwu_flag']

  # Water Supply Flags
  self.hwu_gw_flag = info['water_management']['hwu_gw_flag']
  self.hwu_sf_flag = info['water_management']['hwu_sf_flag']
  if self.hwu_flag == True:
   if not any([self.hwu_gw_flag,self.hwu_sf_flag]): 
    exit('Error: All water supply flags are off - at least one water source is needed')

  # Water Demand Flags
  self.hwu_agric_flag  = info['water_management']['hwu_agric_flag']
  self.hwu_domest_flag = info['water_management']['hwu_domest_flag']
  self.hwu_indust_flag = info['water_management']['hwu_indust_flag']
  self.hwu_lstock_flag = info['water_management']['hwu_lstock_flag']
  demands = [self.hwu_agric_flag,self.hwu_domest_flag,self.hwu_indust_flag,self.hwu_lstock_flag]
  if self.hwu_flag == True:
   if not any(demands):
     exit('Error: At least one water demand is needed') 

  # Land Cover for Water Use -- note: it must be the same defined at Preprocessing.py
  self.wuse_index = {}
  self.water_use_land_cover = {
    'industrial': [13],
    'domestic':   [6,7,8,9,10,13,19],
    'livestock':  [6,7,8,9,10,19],
    'agriculture':[12,14],
    'surface_water':[11,17]}

  # Calc topographic index
  self.beta = HB.input_fp.groups['parameters'].variables['slope'][:] #Noemi
  m = HB.input_fp.groups['parameters'].variables['m'][:] #Depth to bedrock #Noemi
  #pct = HB.input_fp.groups['parameters'].variables['area_pct'][:]/100.0
  #af=10 #anisotropic factor
  #T0 = af*HB.input_fp.groups['parameters'].variables['SATDK'][:]*m
  #self.sti = HB.input_fp.groups['parameters'].variables['ti'][:] #Noemi
  #self.sti = self.sti -(np.log(T0)-np.sum(pct*np.log(T0)))
  self.dem = HB.input_fp.groups['parameters'].variables['dem'][:] 
 

  # GROUNDWATER supply variables
  if self.hwu_gw_flag == True:
   # Well_depht : Correct well_depht for shallow soils and add 20% buffer. Well depth = 80% of bedrock depth
   self.well_depth = np.maximum(-1.0*m*0.8,NOAH.zsoil[:,-1]*0.8)
   self.well_layer = np.argmin(abs(NOAH.zsoil-self.well_depth[:,np.newaxis]))
   self.supply_gw = np.zeros(ncells,dtype=np.float64)
   self.alloc_gw = np.zeros(ncells,dtype=np.float64)

   # Groundwater pumping where clay content < 50% OR hydraulic conductivity > 10^-6 
   self.mask_gw = np.array([ False if i<=10.0**-6 or j>=50. else True for i,j in zip(NOAH.satdk0,NOAH.clay_pct) ])
   # Maximum groundwater pumping distance (km)
   self.gw_dist_lim  = 5.0 #km
   # Maximum slope allowed for upstream pumping
   self.gw_slope_lim = -1./100. # m/m
      

  # SURFACE WATER supply variables
  if self.hwu_sf_flag == True:
   self.supply_sf = np.zeros(ncells,dtype=np.float64)
   self.alloc_sf  = np.zeros(ncells,dtype=np.float64)

   # Surface water availability where ti > 10 and slope < 0.01
   # (surface water only at the river streams -- update this later) with actual river storage
   cond1 = ( self.beta < 0.01) #( self.sti >= 10.0 ) & ( self.beta < 0.01)
   #cond1 = ( self.sti >= 10.0 )
   # surface water from HRUs with land cover as water or wetlands 
   cond2 = np.array([ True if x in self.water_use_land_cover["surface_water"] else False for x in NOAH.vegtyp ])
   self.mask_sf = np.array([ True if i==True or j==True else False for i,j in zip(cond1,cond2) ])
   # Maximum surface water abstractions distance (km)
   self.sf_dist_lim  = 10.0 # km
   # Maximum slope allowed for upstream abstraction
   self.sf_slope_lim = -1./100. # m/m 
   

  # AGRICULTURE demand variables
  if self.hwu_agric_flag == True:
   self.deficit_agric = np.zeros(ncells,dtype=np.float64)
   self.demand_agric  = np.zeros(ncells,dtype=np.float64)
   self.alloc_agric   = np.zeros(ncells,dtype=np.float64)
   self.irrigation    = np.zeros(ncells,dtype=np.float64)
   self.mask_agric    = [True if x in self.water_use_land_cover['agriculture'] else False for x in NOAH.vegtyp ]
   self.wuse_index['a'] = len(self.wuse_index)

   # Irrigation
   # 1: non-paddy irrigation, 2: paddy irrigation, 0: others
   self.irrig_land = HB.input_fp.groups['parameters'].variables['irrig_land'][:]
   self.mask_irrig = self.mask_agric & ( self.irrig_land > 0.0 )
   print('Irrig Percentage:', np.sum(self.mask_irrig)/float(np.sum(self.mask_agric)))
   
   # Crop calendar
   self.st_gscal = np.asarray(HB.input_fp.groups['parameters'].variables['start_growing_season'][:],dtype=np.int)
   self.en_gscal = np.asarray(HB.input_fp.groups['parameters'].variables['end_growing_season'][:],dtype=np.int)
   self.gscal = mgmt_funcs.calc_calendar(self,ncells)
   m = np.where(np.invert(self.mask_agric))[0]
   self.gscal[m,:] = 0.0
   # Note: need to test for the case with demand but zero irrigation


  # INDUSTRIAL demand variables
  if self.hwu_indust_flag == True:
   self.demand_indust  = np.zeros(ncells,dtype=np.float64)
   self.deficit_indust = np.zeros(ncells,dtype=np.float64)
   self.mask_indust    = [True if x in self.water_use_land_cover['industrial'] else False for x in NOAH.vegtyp ]
   self.wuse_index['i'] = len(self.wuse_index)


  # DOMESTIC demand variables
  if self.hwu_domest_flag == True:
   self.demand_domest  = np.zeros(ncells,dtype=np.float64)
   self.deficit_domest = np.zeros(ncells,dtype=np.float64)
   self.mask_domest    = [True if x in self.water_use_land_cover['domestic'] else False for x in NOAH.vegtyp ]
   self.wuse_index['d'] = len(self.wuse_index)


  # LIVESTOCK demand variables
  if self.hwu_lstock_flag == True:
   self.demand_lstock  = np.zeros(ncells,dtype=np.float64)
   self.deficit_lstock = np.zeros(ncells,dtype=np.float64)
   self.mask_lstock    = [True if x in self.water_use_land_cover['livestock'] else False for x in NOAH.vegtyp ]
   self.wuse_index['l'] = len(self.wuse_index)

  self.nwuse_index = len(self.wuse_index)



 def initialize_allocation(self,HB):
  NOAH = HB.noahmp
  
  if self.hwu_flag == True:
   ncells = NOAH.ncells

   # 1. HRU distances
   # Centroid Distances
   self.ctrd_lats = HB.input_fp.groups['parameters'].variables['centroid_lats'][:]
   self.ctrd_lons = HB.input_fp.groups['parameters'].variables['centroid_lons'][:]
   self.hrus_centroid_distances = mgmt_funcs.hrus_centroid_distance(self.ctrd_lats,self.ctrd_lons)  # km
   # Minimum Distance between a HRU centroid and it's neighboors boundary 
   self.hru_min_dist = HB.input_fp.groups['parameters'].variables['hru_min_dist'][:]
   self.hrus_boundary_distances = np.copy(self.hru_min_dist)  # km
   # Assign Distaces to Surface of Grounwater supply
   print("HRU's GW distance - mean:%f and std:%f" % (np.mean(self.hrus_centroid_distances[self.hrus_centroid_distances>0.]), np.std(self.hrus_centroid_distances[self.hrus_centroid_distances>0.])))
   print("HRU's SF distance - mean:%f and std:%f" % (np.mean(self.hrus_boundary_distances[self.hrus_boundary_distances>0.]), np.std(self.hrus_boundary_distances[self.hrus_boundary_distances>0.])))

   # HRU relative distances - Review this at some point
   #self.hrus_rel_dist = np.copy(self.hrus_distances)
   #for i in range(ncells): self.hrus_rel_dist[i,:] = self.hrus_distances[i,:]/np.sum(self.hrus_distances[i,:])

   # Calculate slope between HRUs
   self.hrus_slopes = mgmt_funcs.hrus_slope(self.dem,self.hrus_centroid_distances)  # m/m


   # 2. SURFACE WATER RATIO AND COST
   if self.hwu_sf_flag == True:
    # Surface Water Ratio 
    self.ratio_sf = np.ones((self.nwuse_index,ncells,ncells))

    # Update Ratio for Demands
    if self.hwu_agric_flag == True:
      m = np.copy(np.invert(self.mask_irrig))
      self.ratio_sf[self.wuse_index['a'],:,m] = 0.0
    if self.hwu_domest_flag == True:
      m = np.copy(np.invert(self.mask_domest))
      self.ratio_sf[self.wuse_index['d'],:,m] = 0.0
    if self.hwu_indust_flag == True:
      m = np.copy(np.invert(self.mask_indust))
      self.ratio_sf[self.wuse_index['i'],:,m] = 0.0
    if self.hwu_lstock_flag == True:
      m = np.copy(np.invert(self.mask_lstock))
      self.ratio_sf[self.wuse_index['l'],:,m] = 0.0

    # Update Ratio for Supply
    m = np.copy(np.invert(self.mask_sf)); self.ratio_sf[:,m,:] = 0.0
    m = np.copy(self.hrus_boundary_distances > self.sf_dist_lim ); self.ratio_sf[:,m] = 0.0
    m = np.copy(self.hrus_slopes < self.sf_slope_lim); self.ratio_sf[:,m] = 0.0
    # sum over the diff water sectors, then over the demand nodes
    total_ratio_sf = np.sum(np.sum(self.ratio_sf, axis=0),axis=1)
    # update mask to only simulate over the active supply nodes
    self.mask_sf = (self.mask_sf == True) & ( total_ratio_sf > 0.0)
    #print  self.ratio_sf, self.mask_sf, total_ratio_sf 
    m = np.copy(np.invert(self.mask_sf)); self.ratio_sf[:,m,:] = 0.0
 
    # Surface Water Cost
    self.cost_sf    = np.ones((self.nwuse_index,ncells,ncells))
    #self.cost_sf[:] = self.hrus_boundary_distances


   # 3. GROUNDWATER RATIO AND COST 
   if self.hwu_gw_flag == True:

    # Groundwater Ratio
    self.ratio_gw = np.ones((self.nwuse_index,ncells,ncells))*1.0

    # Update Ratio for Demands
    if self.hwu_agric_flag == True:
      m = np.copy(np.invert(self.mask_irrig))
      self.ratio_gw[self.wuse_index['a'],:,m] = 0.0
    if self.hwu_domest_flag == True:
      m = np.copy(np.invert(self.mask_domest))
      self.ratio_gw[self.wuse_index['d'],:,m] = 0.0
    if self.hwu_indust_flag == True:
      m = np.copy(np.invert(self.mask_indust))
      self.ratio_gw[self.wuse_index['i'],:,m] = 0.0
    if self.hwu_lstock_flag == True:
      m = np.copy(np.invert(self.mask_lstock))
      self.ratio_gw[self.wuse_index['l'],:,m] = 0.0

    # Update Ratio for Supply
    m = np.copy(np.invert(self.mask_gw)); self.ratio_gw[:,m,:] = 0.0
    m = np.copy(self.hrus_centroid_distances > self.gw_dist_lim ); self.ratio_gw[:,m] = 0.0
    m = np.copy(self.hrus_slopes < self.gw_slope_lim); self.ratio_gw[:,m] = 0.0
    # sum over the diff water sectors, then over the demand nodes
    total_ratio_gw = np.sum(np.sum(self.ratio_gw, axis=0),axis=1)
    # update mask to only simulate over the active supply nodes
    self.mask_gw = (self.mask_gw == True) & ( total_ratio_gw > 0.0)
    m = np.copy(np.invert(self.mask_gw)); self.ratio_gw[:,m,:] = 0.0


    # Update Ratio to relative values
    '''print self.hrus_distances
    print self.ratio_gw[0]
    for i in range(ncells):
      m = self.ratio_gw[0,i,:] > 0
      a = self.hrus_distances[i,:][m]
      b = np.sum(self.hrus_distances[i,:][m])
      c = np.divide(a, b, out=np.zeros_like(a), where=b!=0)
      print a,b, c
      self.ratio_gw[0,i,:][m] = 1.0-c
    '''

    # Groundwater Cost
    self.cost_gw     = np.ones((self.nwuse_index,ncells,ncells))*2.0
    #self.cost_gw[:]  = self.hrus_centroid_distances#*(self.sf_dist_lim/self.gw_dist_lim)


   # INITIALIZE ALLOCATION MODEL
   self.hwu_optimal_allocation_flag = True 
   if self.hwu_optimal_allocation_flag == True:

    from pywr.core import Model as wr_Model
    from pywr.core import Input as wr_Input
    from pywr.core import Output as wr_Output
    from pywr.core import Link as wr_Link

    # Define Model
    self.ntwkm = wr_Model(start="2017-01-01", end="2017-01-01", timestep=1, solver='glpk')

    # Define Network Supply and Demand Nodes    
    if self.hwu_agric_flag == True:
     a_nodes_names=[] #Agriculture
     for i in np.where(self.mask_irrig)[0]:
      dem = wr_Output(self.ntwkm, name='a%i' % i, min_flow = 0.0, max_flow = 0.0, cost=-999)
      a_nodes_names.append('a%i' % i)

    if self.hwu_indust_flag == True:
     i_nodes_names=[] #Industrial
     for i in np.where(self.mask_indust)[0]:
      dem = wr_Output(self.ntwkm, name='i%i' % i, min_flow = 0.0, max_flow = 0.0, cost=-999+5)
      i_nodes_names.append('i%i' % i)

    if self.hwu_domest_flag == True:
     d_nodes_names=[] #Domestic
     for i in np.where(self.mask_domest)[0]:
      dem = wr_Output(self.ntwkm, name='d%i' % i, min_flow = 0.0, max_flow = 0.0, cost=-999+15)
      d_nodes_names.append('d%i' % i)

    if self.hwu_lstock_flag == True:
     l_nodes_names=[] #Livestock
     for i in np.where(self.mask_lstock)[0]:
      dem = wr_Output(self.ntwkm, name='l%i' % i, min_flow = 0.0, max_flow = 0.0, cost=-999+10)
      l_nodes_names.append('l%i' % i)

    if self.hwu_sf_flag == True:
     s_nodes_names=[] #Suface Water
     for i in np.where(self.mask_sf)[0]:
      surf = wr_Input( self.ntwkm, name='s%i' % i, min_flow = 0.0, max_flow = 0.0, cost=1)
      s_nodes_names.append('s%i' % i)

    if self.hwu_gw_flag == True:
     g_nodes_names=[] #Groundwater
     for i in np.where(self.mask_gw)[0]:
      gw = wr_Input( self.ntwkm, name='g%i' % i, min_flow = 0.0, max_flow = 0.0, cost=5)
      g_nodes_names.append('g%i' % i)


    # Define Surface Nodes Connections
    if self.hwu_sf_flag == True:
     s_links_names = []
     m = np.copy(self.ratio_sf > 0.0)
     sf_valid_links = m.flatten()
     for n,i,j in zip(*np.where(m)):
      n = list(self.wuse_index.keys())[n]
      link_sf = wr_Link(self.ntwkm, name='s%i_%s%i' % (i,n,j), min_flow = 0.0, max_flow = 0.0)
      self.ntwkm.nodes['s%i'%i].connect(link_sf)
      link_sf.connect(self.ntwkm.nodes['%s%i'% (n,j)])
      s_links_names.append('s%i_%s%i' % (i,n,j))

    # Define Groundwater Nodes Connections
    if self.hwu_gw_flag == True:
     g_links_names = []
     m = np.copy(self.ratio_gw > 0.0)
     gw_valid_links = m.flatten()
     for n,i,j in zip(*np.where(m)):
      n = list(self.wuse_index.keys())[n]
      link_gw = wr_Link(self.ntwkm, name='g%i_%s%i' % (i,n,j), min_flow = 0.0, max_flow = 0.0)
      self.ntwkm.nodes['g%i'%i].connect(link_gw)
      link_gw.connect(self.ntwkm.nodes['%s%i'% (n,j)])
      g_links_names.append('g%i_%s%i' % (i,n,j))

    # Check if there is any source of water in the catchment
    self.valid_links = np.empty(ncells)
    self.valid_links[:] = False
    if self.hwu_sf_flag == True and any(sf_valid_links):
      self.valid_links = sf_valid_links
      if self.hwu_gw_flag == True and any(gw_valid_links):
        self.valid_links = np.concatenate([sf_valid_links,gw_valid_links])
    elif self.hwu_gw_flag == True or any(gw_valid_links):
        self.valid_links = gw_valid_links

    #Model Setup
    # Test for valid nodes and connections
    if any(self.valid_links):

      #self.ntwkm.check()
      self.ntwkm.setup()
      #print dir(self.ntwkm.graph.nodes)
      
      self.nodes_list = list(self.ntwkm.graph.nodes(data=True))
      #print self.nodes_list[0]
      #print dir(self.nodes_list[0][0])

      n_nodes = len(self.nodes_list)
      #print n_nodes

      #print dir(def_list[0][0])
      #print nodes_list[0][0].max_flow    
      nodes_names = [ str(self.nodes_list[i][0].name) for i in range(n_nodes) ]
      nodes_types = [ str(self.nodes_list[i][0].__class__.__name__) for i in range(n_nodes) ]

      # Order the list - List containing the ordered position of each node and link
      # Demands, Supply, Groundwater, Link_i_j
      if self.hwu_agric_flag == True:
        self.a_nodes_position = []
        for i in a_nodes_names: self.a_nodes_position.append( nodes_names.index(i) )
      if self.hwu_domest_flag == True:
        self.d_nodes_position = []
        for i in d_nodes_names: self.d_nodes_position.append( nodes_names.index(i) )
      if self.hwu_indust_flag == True:
        self.i_nodes_position = []
        for i in i_nodes_names: self.i_nodes_position.append( nodes_names.index(i) )
      if self.hwu_lstock_flag == True:
        self.l_nodes_position = []
        for i in l_nodes_names: self.l_nodes_position.append( nodes_names.index(i) )

      if self.hwu_sf_flag == True:
        self.s_nodes_position = []
        for i in s_nodes_names: self.s_nodes_position.append( nodes_names.index(i) )
        self.s_links_position = []
        for i in s_links_names: self.s_links_position.append( nodes_names.index(i) )

      if self.hwu_gw_flag == True:
        self.g_nodes_position = []
        for i in g_nodes_names: self.g_nodes_position.append( nodes_names.index(i) )
        self.g_links_position = []
        for i in g_links_names: self.g_links_position.append( nodes_names.index(i) )

      #print len(s_nodes_names), len(g_nodes_names), len(s_links_names), len(g_links_names)
      #print n_nodes
      #print d_nodes_names[:10]
      #print self.d_nodes_position[:10]
      #print [ nodes_names[i] for i in self.d_nodes_position[:10] ]

      #self.ntwkm.graph.nodes[0][0].max_flow=9876
      #print self.ntwkm.graph.nodes[0][0].max_flow
      #print "Allocation Network Done!"

    return 



 def Calc_Human_Water_Demand_Supply(self,HB,date):
   NOAH = HB.noahmp
   area = np.copy(self.area)
    
   # Only calculate water demand and update allocations at the dta time step.
   if (date.hour*3600.0)%self.dta == 0.0:
 
    if self.hwu_agric_flag  == True:
     # Calculate Agricultural Demand
     adem = np.copy(self.Agriculture_Demand(HB,date))
     self.demand_agric = adem/self.dta #m/s
     # Convert from m  to m3
     self.deficit_agric = adem*area
     #print 'def',self.deficit_agric

    if self.hwu_flag == True:
     # Convert from m/s to m3
     if self.hwu_indust_flag == True: self.deficit_indust = np.copy(self.demand_indust)*self.dta*area
     if self.hwu_domest_flag == True: self.deficit_domest = np.copy(self.demand_domest)*self.dta*area
     if self.hwu_lstock_flag == True: self.deficit_lstock = np.copy(self.demand_lstock)*self.dta*area

    if self.hwu_flag == True:
     # Calculate Supply [m]
     self.Calc_Water_Supply(HB)
     # Convert m to m3 
     if self.hwu_sf_flag == True: self.supply_sf = self.supply_sf*area
     if self.hwu_gw_flag == True: self.supply_gw = self.supply_gw*area

     # Dynamically update water costs
     self.Update_water_cost()
     
     # Allocation
     self.Optimal_Water_Allocation(NOAH)
     
    # Convert from m3 to m/s
    if self.hwu_agric_flag  == True: self.deficit_agric  = self.deficit_agric/(self.dta/area)

    
    if self.hwu_flag == True:
 
     with np.errstate(divide='ignore'): 
       # Convert from m3 to m/s  
       if self.hwu_indust_flag == True: self.deficit_indust = self.deficit_indust/(self.dta/area)
       if self.hwu_domest_flag == True: self.deficit_domest = self.deficit_domest/(self.dta/area)
       if self.hwu_lstock_flag == True: self.deficit_lstock = self.deficit_lstock/(self.dta/area)
     
       # Convert from m3 to m/tstep
       if self.hwu_agric_flag  == True: self.alloc_agric  = self.alloc_agric/(area/self.ntt)
       if self.hwu_gw_flag == True : self.alloc_gw = self.alloc_gw/(area/self.ntt)
       if self.hwu_sf_flag == True : self.alloc_sf = self.alloc_sf/area#/self.ntt #Abstract sf water just on dta

       # Convert from m3 to m (m3/m2)
       if self.hwu_sf_flag == True: self.supply_sf = self.supply_sf/area
       if self.hwu_gw_flag == True: self.supply_gw = self.supply_gw/area
       #print self.alloc_sf
   
    #print 'D', self.demand_domest, self.deficit_domest
    #print 'I',self.demand_indust, self.deficit_indust
    #print 'A',self.demand_agric , self.deficit_agric
    #print 'sf',self.alloc_sf
    #print 'gw', self.alloc_gw*self.ntt
    #print 'diff',np.sum(self.alloc_agric*area*self.ntt), np.sum(self.alloc_sf*area+(self.alloc_gw*area*self.ntt))
    #print 'diff', np.sum(self.demand_agric-self.deficit_agric)*self.dta, np.sum([self.alloc_sf/self.dta]*self.dta)
    #print "Alloc agr",np.sum(self.alloc_agric*self.ntt)
    #print 'alloc sup',np.sum(self.alloc_sf+(self.alloc_gw*self.ntt))
   return


 def Calc_Water_Supply(self,HB):
  NOAH = HB.noahmp
  ncells = NOAH.ncells
  sldpth = NOAH.sldpth
 
  if self.hwu_flag == True:

   # Groundwater Supply
   if self.hwu_gw_flag == True :
    
    dsm = NOAH.sh2o - NOAH.smcref[:,np.newaxis] #moisture content between 
    dsm[dsm < 0] = 0    
    self.supply_gw = np.array([ np.sum(dsm[i,:self.well_layer+1]*sldpth[i,:self.well_layer+1]) for i in range(ncells) ])
    self.supply_gw[np.invert(self.mask_gw)]=0.0 
    # Future: Include corrections for environmental flows
    self.supply_gw = self.supply_gw*0.8
    

   # Surface Water Supply
   if self.hwu_sf_flag == True :
    self.supply_sf = np.copy(NOAH.runsf)
    self.supply_sf = self.supply_sf*(HB.dt/1000.) # from mm/s to m
    self.supply_sf[np.invert(self.mask_sf)] = 0

    # Future Include correction for environmetal flows.
    self.supply_sf = self.supply_sf*0.8
    #print self.supply_sf, self.mask_sf

  return


 def Water_Supply_Abstraction(self, HB, date):
  NOAH = HB.noahmp
  ncells = NOAH.ncells
  sldpth = NOAH.sldpth
  nsoils = np.arange(NOAH.sh2o.shape[1])
 
  if self.hwu_flag == True:
   # Abstract from Surface
   if self.hwu_sf_flag == True:
    if (date.hour*3600.0)%self.dta == 0.0: # Only abstract surface water at the allocation time 
     m = (self.alloc_sf > 0.0)
     #print 'sf:', self.alloc_sf
     if any(t < 0 for t in NOAH.runsf):
       print('1 Negative Runoff!!!', NOAH.runsf)
     #print 'runoff', NOAH.runsf, self.alloc_sf*(1000.0/HB.dt)
     NOAH.runsf[m] = (NOAH.runsf - (self.alloc_sf*(1000.0/HB.dt)))[m] # mm/s
     #print 'runoff abs',NOAH.runsf
     if any(t < 0 for t in NOAH.runsf[m]):
       print('2 Negative Runoff!!!', NOAH.runsf[m], (self.alloc_sf*(1000./HB.dt))[m])
       print('salloc',self.alloc_sf*(1000./HB.dt))
       print('mask',self.mask_sf)
       exit()

   # Abstract from Grondwater
   if self.hwu_gw_flag == True:
    m = (self.alloc_gw > 0.0)
    #if (date.hour*3600)%self.dta == 0: print 'gw:', self.alloc_gw
    
    if HB.subsurface_module == 'dtopmodel':
      NOAH.dzwt[m] = (NOAH.dzwt-(self.alloc_gw))[m]  # m
    
    elif HB.subsurface_module == 'richards':
         
      dsm = NOAH.sh2o - NOAH.smcref[:,np.newaxis]
      dsm[dsm<0]=0
      volume_avail  = np.array([ dsm[i,:]*sldpth[i,:] for i in range(ncells) ]) #m per layer
      volume_final = np.copy(volume_avail)
      voli = np.copy(self.alloc_gw)
  
      # Check Layer UP
      if np.sum(voli)> 0:
       for l in reversed(nsoils[:self.well_layer+1]):
        if np.sum(voli) > 0:
         m1 = (voli > 0)      
         volume_final[m1,l] = volume_final[m1,l]-voli[m1]
         voli[m1] = 0.0
         m2 = (volume_final[:,l] < 0.0)
         voli[m2] = abs(volume_final[:,l])[m2]
         volume_final[m2,l]=0.0
        else:
         break
      # Check Layer DOWN
      if np.sum(voli)> 0:
       for l in nsoils[self.well_layer+1:]:
        if np.sum(voli)> 0:
         m1 = (voli > 0)
         volume_final[m1,l] = volume_final[m1,l]-voli[m1]
         voli[m1] = 0.0
         m2 = (volume_final[:,l] < 0.0)
         voli[m2] = abs(volume_final[:,l])[m2]
         volume_final[m2,l]=0.0
        else:
         break
      delta_dhiv = (volume_avail-volume_final)*1000.0/HB.dt
      m = (delta_dhiv > 0)
      HB.richards.hdiv[m] = (HB.richards.hdiv-delta_dhiv)[m]  # mm/s  
      
  return NOAH


 def Update_water_cost(self,):
  # Groundwater
  if self.hwu_gw_flag == True:  
   #vol = np.copy(self.supply_gw)[:,np.newaxis]+1.
   #ncost = np.copy(self.hrus_centroid_distances+1)#/vol
   #m = (ncost == np.inf) | (ncost == np.nan)
   #ncost[m] = 10000000000.
   #self.cost_gw[:] = ncost
   self.cost_gw[:] = np.copy(self.hrus_centroid_distances+0.1)

  if self.hwu_sf_flag == True:
   #vol = np.copy(self.supply_sf)[:,np.newaxis]+1.
   #ncost = np.copy(self.hrus_boundary_distances+1)#/vol
   #m = (ncost == np.inf) | (ncost == np.nan)
   #ncost[m] = 10000000000.
   #self.cost_sf[:] = ncost
   self.cost_sf[:] = np.copy(self.hrus_boundary_distances+0.1)

  return
     

 def Human_Water_Irrigation(self,HB,date):
  if self.hwu_flag == True:
   if self.hwu_agric_flag  == True:# and (date.hour*3600)%self.dta == 0:
        
    # Add as irrigation the amount of water that was allocated
    self.irrigation = np.copy(self.alloc_agric)*(1000.0/HB.noahmp.dt)  #from m/tstep to mm/s
    m = (self.mask_irrig == True)
    #print self.irrigation
    #print self.irrigation[m]
    #print 'irrig',np.sum(self.irrigation)*(HB.noahmp.dt/1000.0)
    HB.noahmp.prcp[m] = (HB.noahmp.prcp + self.irrigation)[m] #mm/s
    # Include irrigation efficiency later on

  return


 def Optimal_Water_Allocation(self,NOAH):  

    ncells = NOAH.ncells
   
    if any(self.valid_links.flatten()):

     p1 = datetime.now()
     #print 'P1', str(p1)

     # Update Node
     if self.hwu_agric_flag  == True:
      self.alloc_agric[:] = 0.0
      a_deficit = np.copy(self.deficit_agric)
      for p,d in zip(self.a_nodes_position,a_deficit[self.mask_irrig]):
        self.nodes_list[p][0].max_flow = d
     
     if self.hwu_indust_flag  == True:
      i_deficit = np.copy(self.deficit_indust)
      for p,d in zip(self.i_nodes_position,i_deficit[self.mask_indust]):
        self.nodes_list[p][0].max_flow = d
 
     if self.hwu_domest_flag  == True:
      d_deficit = np.copy(self.deficit_domest)
      for p,d in zip(self.d_nodes_position,d_deficit[self.mask_domest]):
        self.nodes_list[p][0].max_flow = d
  
     if self.hwu_lstock_flag  == True:
      l_deficit = np.copy(self.deficit_lstock)
      for p,d in zip(self.l_nodes_position,l_deficit[self.mask_lstock]):
        self.nodes_list[p][0].max_flow = d

     if self.hwu_sf_flag  == True:
      self.alloc_sf[:] = 0.0
      sf_supply = np.copy(self.supply_sf)
      #print "sf_supply", sf_supply
      #print "sf_supply", sf_supply[self.mask_sf]
      for p,s in zip(self.s_nodes_position,sf_supply[self.mask_sf]):
        self.nodes_list[p][0].max_flow = s
        self.nodes_list[p][0].cost = 1
        #print dir(self.nodes_list[p][0])
        #print self.nodes_list[p][0].min_flow[:]      
        if s < 1e-10: # Dont abstract from less the 1 liters over dt
           self.nodes_list[p][0].max_flow = 0.0
           #self.nodes_list[p][0].cost = 1e+20   
        

     if self.hwu_gw_flag  == True:
      self.alloc_gw[:] = 0.0
      gw_supply = np.copy(self.supply_gw)
      for p,s in zip(self.g_nodes_position,gw_supply[self.mask_gw]):
        self.nodes_list[p][0].max_flow = s
        self.nodes_list[p][0].cost = 5
        if s < 1e-10: # Dont abstract from less the 1 liters
          self.nodes_list[p][0].max_flow = 0.0
          #self.nodes_list[p][0].cost = 1e+20
      
     
     p2 = datetime.now()
     #print 'P2', str(p2-p1), str(p2-p1)

     # Update Links
     # Surface Water
     if self.hwu_sf_flag  == True:
      m = self.ratio_sf > 0.0
      cost  = np.copy(self.cost_sf[m])
      max_flow = np.copy(np.transpose(np.transpose(self.ratio_sf,(0, 2, 1))*sf_supply,(0,2,1))[m])
      #print max_flow
      for p,f,c in zip(self.s_links_position, max_flow, cost):
       self.nodes_list[p][0].max_flow = f
       self.nodes_list[p][0].cost = c  
       if f < 1e-10: # Dont abstract from less the 1 liters
          self.nodes_list[p][0].max_flow = 0.0
          #self.nodes_list[p][0].cost = 1e+20
       #self.nodes_list[p][0].min_flow = 1e-10
          
     # Groundwater
     if self.hwu_gw_flag  == True:
      m = self.ratio_gw > 0
      cost  = np.copy(self.cost_gw[m])
      max_flow =  np.copy(np.transpose(np.transpose(self.ratio_gw,(0, 2, 1))*gw_supply,(0,2,1))[m])
      for p,f,c in zip(self.g_links_position, max_flow, cost):
       self.nodes_list[p][0].max_flow = f
       self.nodes_list[p][0].cost = c
       if f < 1e-10: # Dont abstract from less the 1 liters
          self.nodes_list[p][0].max_flow = 0.0
          #self.nodes_list[p][0].cost = 1e+20

     # Optimize the network
     if self.check_setup() == True: 
    
        p3 = datetime.now()
        #print 'P3', str(p3-p1), str(p3-p2)
        #self.ntwkm.run() 
        self.ntwkm._step()
  
        p4 = datetime.now()
        #print 'P4', str(p4-p1), str(p4-p3)
        
        delta = np.zeros(ncells)     
        def2correct = 0.0
        
        if self.hwu_agric_flag == True:
         delta[:] = 0
         for p,h in zip(self.a_nodes_position,np.where(self.mask_irrig)[0]):
           delta[h] = np.copy(self.nodes_list[p][0].flow[0])
         #print "demad",a_deficit
         #print "delta",delta
         #print "defic",a_deficit-delta
         self.alloc_agric = np.copy(delta)
         self.deficit_agric = np.copy(a_deficit-self.alloc_agric)
         
         
        if self.hwu_indust_flag == True:
         delta[:] = 0     
         for p,h in zip(self.i_nodes_position,np.where(self.mask_indust)[0]):
           delta[h] = np.copy(self.nodes_list[p][0].flow[0])
         self.deficit_indust = np.copy(i_deficit - delta)
         
        if self.hwu_domest_flag == True:
         delta[:] = 0     
         for p,h in zip(self.d_nodes_position,np.where(self.mask_domest)[0]):
           delta[h] = np.copy(self.nodes_list[p][0].flow[0])
         self.deficit_domest = np.copy(d_deficit - delta)
       
        if self.hwu_lstock_flag == True:
         delta[:] = 0     
         for p,h in zip(self.l_nodes_position,np.where(self.mask_lstock)[0]):
           delta[h] = np.copy(self.nodes_list[p][0].flow[0])
         self.deficit_lstockl = np.copy(l_deficit - delta)
        
        if self.hwu_sf_flag == True:
         delta[:] = 0
         for p,h in zip(self.s_nodes_position,np.where(self.mask_sf)[0]):
           delta[h] = np.copy(self.nodes_list[p][0].flow[0])
         #print 'delta sf',delta
         self.alloc_sf = np.copy(delta)
         pos = np.where((self.alloc_sf>0.0) & (self.alloc_sf < 1e-10))[0]
         #print pos
         for p in pos:
          def2correct =+ np.copy(self.alloc_sf[p])
          self.alloc_sf[p] = 0.0
         self.supply_sf = np.copy(sf_supply - self.alloc_sf)
        

        if self.hwu_gw_flag == True:
         delta[:] = 0
         for p,h in zip(self.g_nodes_position,np.where(self.mask_gw)[0]):
           delta[h] = np.copy(self.nodes_list[p][0].flow[0])
         #print 'delta gw',delta
         self.alloc_gw = np.copy(delta)
         pos = np.where((self.alloc_gw>0.0) & (self.alloc_gw < 1e-10))[0]
         for p in pos:
          def2correct =+ np.copy(self.alloc_gw[p])
          self.alloc_gw[p] = 0.0
         self.supply_gw = np.copy(gw_supply - self.alloc_gw)

           
        if def2correct > 0 :
         if self.hwu_agric_flag == True:
          pos = np.argmax(self.alloc_agric)
          #print 'remove', pos, self.alloc_agric 
          self.alloc_agric[pos] = self.alloc_agric[pos]-def2correct
          self.deficit_agric[pos] = self.deficit_agric[pos] +def2correct
         elif self.hwu_indust_flag == True:
          pos = np.argmax(self.deficit_indust)
          self.deficit_indust[pos] = self.deficit_indust[pos] +def2correct
         elif self.hwu_domest_flag == True:
          pos = np.argmax(self.deficit_domest)
          self.deficit_domest[pos] = self.deficit_domest[pos] +def2correct
         elif self.hwu_lstock_flag == True:
          pos = np.argmax(self.deficit_lstock)
          self.deficit_lstock[pos] = self.deficit_lstock[pos] +def2correct
         else:
          exit('Warning!! Water Demands not closing!!!')


        #print np.sum(self.alloc_gw+self.alloc_sf-self.alloc_agric)

        p5 = datetime.now()
        #print 'P5', str(p5-p1), str(p5-p4)
          
     return


 def Agriculture_Demand(self,HB,date):
  NOAH = HB.noahmp
      
  demand_agric = np.copy(self.Calculate_Irrigation_Deficit(NOAH)) # m
  demand_agric = demand_agric/0.85 # Efficiency
  #print demand_agric
  # Limit Demand for the Crop Calendar
  mnt = date.month-1
  m = self.gscal[:,mnt] == 0
  demand_agric[m] = 0.0  

  #print 'cal',demand_agric 
  return demand_agric


 def Calculate_Irrigation_Deficit(self,NOAH):
  """
  Calculates the irrigation deficit based on how much water [m] is necessary to 
  soil moisture reach field capacity
  """
  ncells = np.copy(NOAH.ncells)
  nroot_zone_depht = np.copy(NOAH.root_depth)
  smc = np.copy(NOAH.sh2o)
  sldpth = np.copy(NOAH.sldpth)
  smcref = np.copy(NOAH.smcref) # Field Capacity
  wltsmc0 = np.copy(NOAH.wltsmc0) # Wilting point
  smcmax = np.copy(NOAH.smcmax) # Saturated 
  RAW = np.copy(NOAH.smcref-(NOAH.smcref-NOAH.wltsmc0))*0.7

  demand_agric = np.zeros(ncells,dtype=np.float64)
  soil_demand = np.zeros(ncells,dtype=np.float64)

  # Irrigation Demand for general crop areas
  m = self.mask_agric
  if any(m):
   dsm = smcref[:,np.newaxis] - smc
   dsm [dsm < 0] = 0.0
   soil_demand[:] = np.array([ np.sum(dsm[i,:nroot_zone_depht[i]]*sldpth[i,:nroot_zone_depht[i]]) for i in range(ncells) ])

   # Irrigation Trigger
   #mean_smc_root = np.array([ np.sum(smc[i,:nroot_zone_depht[i]]*sldpth[i,:nroot_zone_depht[i]])/np.sum(sldpth[i,:nroot_zone_depht[i]]) for i in range(ncells) ])
   #print mean_smc_root
   #print RAW
   #self.mask_trigger = (mean_smc_root > RAW)
   #print self.mask_trigger

   # Check if irrigation is triggered or not - Ozdodan et. al (2010)
   #MAi = (smc - wltsmc0[:,np.newaxis])/(smcref[:,np.newaxis]-wltsmc0[:,np.newaxis])
   #MA = np.array([ np.sum(MAi[i,:nroot_zone_depht[i]]*sldpth[i,:nroot_zone_depht[i]])/np.sum(sldpth[i,:nroot_zone_depht[i]]) for i in range(ncells) ])
   #print MA
   #self.mask_trigger = (MA > 0.5*smcref) 
   #print self.mask_trigger

   #soil_demand[self.mask_trigger] = 0.0
   demand_agric = soil_demand
   #print demand_agric 
 
  # Irrigation Demand for paddy  crop areas
  m = ( self.irrig_land == 2.0 ) & (self.mask_agric == True)
  if any(m):
   dsm = smcmax[:,np.newaxis] - smc
   dsm[dsm < 0] = 0.0
   soil_demand = np.array([ np.sum(dsm[i,:nroot_zone_depht[i]]*sldpth[i,:nroot_zone_depht[i]]) for i in range(ncells) ])
   demand_agric[m] = soil_demand[m]  
  
  

  return demand_agric

    
 def check_setup(self,):
   total_deficit = 0.0; total_sf_supply= 0.0; total_gw_supply=0.0
   if self.hwu_agric_flag == True:  total_deficit +=np.sum(self.deficit_agric[self.mask_irrig])
   if self.hwu_domest_flag == True: total_deficit +=np.sum(self.deficit_domest[self.mask_domest])
   if self.hwu_indust_flag == True: total_deficit +=np.sum(self.deficit_indust[self.mask_indust])
   if self.hwu_lstock_flag == True: total_deficit +=np.sum(self.deficit_lstock[self.mask_lstock])
   if self.hwu_sf_flag == True: total_sf_supply = np.sum(self.supply_sf[self.mask_sf])
   if self.hwu_gw_flag == True: total_gw_supply = np.sum(self.supply_gw[self.mask_gw])

   if total_deficit != 0.0 and ( total_sf_supply != 0.0 or total_gw_supply != 0.0 ):
    out = True
   else: 
    out = False
   return out
  
  

