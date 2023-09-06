import numpy as np
#import terrain_tools_fortran as ttf
from . import terrain_tools_fortran as ttf
from . import metrics 
import sklearn.cluster
import sklearn.linear_model
#import shapely
#import geopandas
import scipy.stats
import scipy.sparse
import copy
import time
import pickle
import copy
import numba
import rasterio

def surface_area(phi_a,phi_b,lambda_c,lambda_d):
    Re = 6371000.0 #m
    #Convert to radians
    phi_a = phi_a*np.pi/180.0
    phi_b = phi_b*np.pi/180.0
    lambda_c = lambda_c*np.pi/180.0
    lambda_d = lambda_d*np.pi/180.0
    #Calculate surface area
    return Re**2*np.abs(np.sin(phi_a)-np.sin(phi_b))*np.abs(lambda_c - lambda_d)

class terrain_analysis:

 def __init__(self,file,mask_file=False):
 
  #Read in the dem
  self.dem = rasterio.open(file).read(1)
  #self.demns = rasterio.open(file).read(1)
  #Determine the spatial resolution (assumed to be the same in x and y)
  fp = rasterio.open(file)
  if fp.crs == 'EPSG:4326': #This is a temporary fix
   phi_a = fp.bounds.bottom
   phi_b = fp.bounds.top
   lambda_c = fp.bounds.left
   lambda_d = fp.bounds.right
   area = surface_area(phi_a,phi_b,lambda_c,lambda_d)
   #Determine the average area per cell
   area = area/fp.width/fp.height
   self.dx = area**0.5
  else:
   self.dx = rasterio.open(file).res[0] #meters
  self.res = fp.res
  #Read in the bound of the data
  self.bounds = rasterio.open(file).bounds
  #Create mask
  if mask_file == False:
   self.mask = np.ones(self.dem.shape) 
   self.mask[self.dem == -9999] = 0
  #Basin area
  self.basin_area = self.dx**2*np.sum(self.mask)
  #Sink fill
  self.demns = ttf.remove_pits_planchon(self.dem,self.dx)
  #Slope and aspect
  res_array = np.copy(self.dem)
  res_array[:] = self.dx
  (slope,aspect) = ttf.calculate_slope_and_aspect(np.flipud(self.dem),res_array,res_array)
  self.slope = np.flipud(slope)
  self.aspect = np.flipud(aspect)
  #QC
  self.slope[self.slope > 1.0] = 1.0

  return

 def calculate_drainage_area(self,):

  (self.acc,self.fdir) = ttf.calculate_d8_acc(self.demns,self.mask,self.dx)

  return
  
 def delineate_river_network(self,):

  #Construct array of x,y
  x = np.linspace(self.bounds.bottom+self.res[0]/2,self.bounds.top-self.res[0]/2,self.acc.shape[0])
  y = np.linspace(self.bounds.left+self.res[1]/2,self.bounds.right-self.res[1]/2,self.acc.shape[1])
  (xs,ys) = np.meshgrid(x,y)
  thld = self.channel_threshold
  (channels,channels_wob,channel_topology,tmp1,crds) = ttf.calculate_channels_wocean_wprop_wcrds(self.acc,thld,thld,self.fdir,self.mask,np.flipud(xs.T),ys.T)
  #Compute and output the list of the channel positions
  lst_crds = []
  for icrd in range(crds.shape[0]):
    mcrd = crds[icrd,:,0] != -9999
    if (np.sum(mcrd) == 0):break
    crds_i = crds[icrd,mcrd,:]
    if crds_i.shape[0] > 1:
        lst_crds.append(shapely.geometry.LineString(np.fliplr(crds_i)))
    else:
        lst_crds.append(shapely.geometry.Point(np.flipud(crds_i[0,:])))
  self.topology = channel_topology
  self.channels_vector = geopandas.GeoSeries(lst_crds)
  self.channels_raster = np.copy(channels_wob)
  self.shreve_order = tmp1 #This is broken
  self.reach_length = []
  #Calculate the length of each reach
  self.channels_vector = self.channels_vector.set_crs("EPSG:4326")
  self.channels_vector_reproj = self.channels_vector.to_crs("+proj=laea +lon_0=-120.234375 +lat_0=49.5878293 +datum=WGS84 +units=m +no_defs")
  self.stream_length = 0
  for geom in self.channels_vector_reproj:
    self.stream_length += geom.length
    self.reach_length.append(geom.length)
  self.reach_length = np.array(self.reach_length)

  return

 def calculate_reach_properties(self,):

  channels = self.channels_raster
  topology = self.topology
  slope = self.slope
  dx = self.dx
  mask = self.mask
  acc = self.acc
  basins = self.basins
  shreve_order = self.shreve_order
  pscaling = {'channel_manning':1,'floodplain_manning':1,'bankfull_depth':1,'channel_width':1}

  #calculate the properties
  db_channels = calculate_channel_properties(channels,topology,slope,dx,mask,acc,
          acc,basins,shreve_order,pscaling)

  self.db_channels = db_channels

  return

 def delineate_basins(self,):

  self.basins = ttf.delineate_basins(self.channels_raster,self.mask,self.fdir)
  ubs = np.unique(self.basins)
  ubs = ubs[ubs != -9999]
  area = []
  for ub in ubs:
      area.append(self.dx**2*np.sum(self.basins == ub))
  self.basin_area = np.array(area)

  return

 def calculate_height_above_nearest_drainage(self,):

  self.hand = ttf.calculate_depth2channel(self.channels_raster,self.basins,self.fdir,self.demns)

  return 

 def discretize_hand(self,):

    ubs = np.unique(self.basins)
    ubs = ubs[ubs > 0]
    self.hbands = -9999*np.ones(self.basins.shape).astype(np.int32)
    hband = 1
    for ub in ubs:
        m = self.basins == ub
        tmp = self.hand[m]
        #create the bins
        bins = [0,]
        while np.max(bins) < np.max(tmp):
            bins.append(bins[-1]+self.dh)
        bins = np.array(bins)
        #Compress the bins
        if ((bins[-1] - np.max(tmp)) < self.dh) & (bins.size >= 3):
            bins[-2] = bins[-1]
            bins = bins[:-1]
        #assign the cells to a given hru
        for i in range(bins.size-1):
            m2 = m & (self.hand >= bins[i]) & (self.hand < bins[i+1])
            self.hbands[m2] = hband
            hband += 1

    return


def sink_fill(dem,dx):

 return ttf.remove_pits_planchon(dem,dx)

def delineate_river_network(network):

  bounds = network['bounds']
  acc = network['acc']
  thld = network['channel_threshold']
  fdir = network['fdir']
  mask = network['mask']
  dx = network['dx']
  #Construct array of x,y
  x = np.linspace(bounds.bottom+dx/2,bounds.top-dx/2,acc.shape[0])
  y = np.linspace(bounds.left+dx/2,bounds.right-dx/2,acc.shape[1])
  (xs,ys) = np.meshgrid(x,y)
  (channels,channels_wob,channel_topology,tmp1,crds) = ttf.calculate_channels_wocean_wprop_wcrds(acc,thld,thld,fdir,mask,np.flipud(xs.T),ys.T)
  #Compute and output the list of the channel positions
  lst_crds = []
  for icrd in range(crds.shape[0]):
    mcrd = crds[icrd,:,0] != -9999
    if (np.sum(mcrd) == 0):break
    crds_i = crds[icrd,mcrd,:]
    if crds_i.shape[0] > 1:
        lst_crds.append(shapely.geometry.LineString(np.fliplr(crds_i)))
    else:
        lst_crds.append(shapely.geometry.Point(np.flipud(crds_i[0,:])))
  network['channels_vector'] = geopandas.GeoSeries(lst_crds)
  network['channels_raster'] = np.copy(channels_wob)

  return network

def delineate_basins(channels,fdir,demns,mask=False):

  if mask == False:
    mask = np.copy(demns)
    mask[demns == -9999] = 0
    mask[demns != -9999] = 1
  basins = ttf.delineate_basins(channels,mask,fdir)

  return basins

def calculate_height_above_nearest_drainage(channels,basins,fdir,dem):

 hand = ttf.calculate_depth2channel(channels,basins,fdir,dem)

 return hand

def discretize_hand(hand,basins,dh):

    hand = np.ma.getdata(hand)
    ubs = np.unique(basins)
    ubs = ubs[ubs > 0]
    hbands = -9999*np.ones(basins.shape).astype(np.int32)
    hband = 1
    for ub in ubs:
        m = basins == ub
        tmp = hand[m]
        #create the bins
        bins = [0,]
        while np.max(bins) < np.max(tmp):
            bins.append(bins[-1]+dh)
        bins = np.array(bins)
        #Compress the bins
        if ((bins[-1] - np.max(tmp)) < dh) & (bins.size >= 3):
            bins[-2] = bins[-1]
            bins = bins[:-1]
        #assign the cells to a given hru
        for i in range(bins.size-1):
            m2 = m & (hand >= bins[i]) & (hand < bins[i+1])
            hbands[m2] = hband
            hband += 1

    return hbands

def calculate_distance(lat0,lat1,lon0,lon1):

 R = 6372800 #meters
 dlat = np.deg2rad(lat1-lat0)
 dlon = np.deg2rad(lon1-lon0)
 lat1 = np.deg2rad(lat1)
 lat0 = np.deg2rad(lat0)
 a = np.sin(dlat/2)**2 + np.cos(lat0)*np.cos(lat1)*np.sin(dlon/2)**2
 c = 2*np.arctan2(np.sqrt(a),np.sqrt(1-a))
 return R*c

def calculate_area(r):

 lats = np.linspace(r.miny,r.maxy,r.ny+1)#2*r.ny+1)
 lons = np.linspace(r.minx,r.maxx,r.nx+1)#2*r.nx+1)
 (lats,lons) = np.meshgrid(lats,lons)
 #r.dx = calculate_distance(lats[1::2,1::2],lats[1::2,1::2],lons[0:-2:2,0:-2:2],lons[2::2,2::2])
 r.dx = calculate_distance((lats[0:-1,0:-1]+lats[1:,1:])/2,
                           (lats[0:-1,0:-1]+lats[1:,1:])/2,
                           lons[0:-1,0:-1],lons[1:,1:]).astype(np.float32)
 #r.dy = calculate_distance(lats[0:-2:2,0:-2:2],lats[2::2,2::2],lons[1::2,1::2],lons[1::2,1::2])
 r.dy = calculate_distance(lats[0:-1,0:-1],lats[1:,1:],
                           (lons[0:-1,0:-1]+lons[1:,1:])/2,
                           (lons[0:-1,0:-1]+lons[1:,1:])/2).astype(np.float32)
 r.area = r.dx*r.dy

 return r

def frelief_inv(y,a,b):
 return (1 - (1 - y)**(1/b))**(1/a)

def frelief(x,a,b):
 return 1 - (1 - x**a)**b

def fwidth(x,a):
 return 1 + a*x

def fslope(x,a,b):
 return a + b*x

def normalize_variable(input,min,max):

 data = np.copy(input)
 m = data != -9999
 #if (np.max(data[m]) != np.min(data[m])):
 # data[m] = (data[m] - np.min(data[m]))/(np.max(data[m]) - np.min(data[m]))
 #if (max != min):
 # data[m] = (data[m] - min)/(max-min)
 #else:
 # data[m] = 0.0
 #if (np.std(data[m] != 0):
 # data[m] = (data[m] - np.mean(data[m]))/np.std(data[m])
 #else:
 # data[m] = 0
 if (np.nanstd(data[m]) != 0.0):
  data[m] = (data[m] - np.nanmean(data[m]))/np.nanstd(data[m])
 else:
  data[m] = 0.0

 return data

def cluster_data(X,nc):
 #Assemble sample list`
 minsamples = 10**5
 if X.shape[0] > minsamples:
  np.random.seed(1245)
  idx = np.random.choice(np.arange(X.shape[0]),minsamples)
  #idx = np.arange(X.shape[0])[0:minsamples]
 else:
  idx = np.arange(X.shape[0])
 
 #The number of clusters must be equal or smaller than the number of samples
 if idx.size < nc:nc = idx.size

 #Cluster the data
 if nc > 1:
  nfeatures = X.shape[1]
  ##laura: Initialize centroids
  centroids=np.empty([nc,nfeatures]) #laura
  max=np.max(X,axis=0) #laura
  min=np.min(X,axis=0) #laura
  rng=max-min #laura
  centroids[0,:]=min #laura
  centroids[-1,:]=max #laura
  delta=rng/(nc-1) #laura
  for c in range(1,centroids.shape[0]-1): #laura
   centroids[c,:]=min+c*delta #laura
   #model = sklearn.cluster.MiniBatchKMeans(n_clusters=nc,random_state=35799)
   #model = sklearn.cluster.KMeans(n_clusters=nc,random_state=35799,init='k-means++') #laura, commented out
  model = sklearn.cluster.KMeans(n_clusters=nc,init=centroids)
  #model = KMeans(nc)
  #model.compute(X)
  #exit()
  #model.fit(X[idx,:])
  p = model.fit(X[idx,:])
  #p = model.fit_predict(X)
  p = model.predict(X)
 else:
  p = np.zeros(X.shape[0])

 return p

def compute_performance_metrics(Xd,data):

 #Evaluate performance
 vals = []
 for var in Xd:
  #tmp = []
  ucs = np.unique(data)
  tmp = np.zeros(data.size)
  tmp[:] = -9999
  #Normalize data
  obs = normalize_variable(Xd[var]['d'],Xd[var]['min'],Xd[var]['max'])
  #ms = []
  for uc in ucs:
   m = data == uc
   #print var,uc,np.unique(obs[m])
   #p32 = np.percentile(obs[m],32)
   #p68 = np.percentile(obs[m],68)
   #print var,p10,p90,p90-p1
   #tmp.append(p68-p32)
   #tmp.append(2*np.std(obs[m]))
   tmp[m] = np.mean(obs[m])
   #ms.append(metrics.RMSE(obs[m],tmp[m]))
   #ms.append(np.abs(obs[m]-tmp[m])))
  #maes.append(metrics.MAE(Xd[var]['d'],tmp))
  RMSE = metrics.RMSE(obs,tmp)
  vals.append(RMSE)
  #vals.append(metrics.MAE(obs,tmp))
  #vals.append(np.max(ms))
  #print var,np.max(ms)
  #print var,np.max(tmp)
  #vals.append(np.mean(tmp))#percentile(tmp,90))

 return np.array(vals)

def compute_cluster_parameters(Xd,maxnc=1000):

 #Assemble tolerances
 tols = []
 for var in Xd:
  tols.append(Xd[var]['t'])
 tols = np.array(tols)

 #Initialize weights
 ws = np.ones(len(Xd.keys()))

 #Initialize maes
 maes0 = 2*tols

 #Initialize parameters
 nc = 1
 flag = False
 ncl = 0
 ncr = 0
 count = 0
 size = tols.size
 tcount = 0

 while flag == False:

  tcount += 1
  if tcount == 1000:return (nc,ws)
  #0.Prepare the data
  X = []
  #print Xd.keys()
  for var in Xd:
   #print Xd[var]['min'],Xd[var]['max']
   tmp = normalize_variable(np.copy(Xd[var]['d']),Xd[var]['min'],Xd[var]['max'])
   tmp = ws[Xd.keys().index(var)]*tmp#(tmp-np.min(tmp))/(np.max(tmp)-np.min(tmp))
   X.append(tmp)
  X = np.array(X).T
  #Remove columns with all zeros
  #m = np.sum(X,axis=0) != 0
  #X = X[:,m]
  #1.Cluster data
  data = cluster_data(X,nc)
  #2.Compute the performance metrics
  maes = compute_performance_metrics(Xd,data)
  #3.Determine the next step
  #print 't',tols
  #print ncl,nc,ncr,maes,ws,size,np.sum(maes > tols)
  #if (np.sum(maes <= tols) == maes.size) & (ncr - ncl == 1):
  if nc >= maxnc:
   nc = maxnc
   break
  elif (ncr - ncl == 1):
   break
  elif (np.sum(maes > tols) < size) & (nc <= 100) & (count < 20):
   rc = 0.01*(maes - tols)/tols
   ws = ws + rc
   ws[ws < 0] = 0#10**-10
   ws = ws/np.sum(ws)
   size = np.sum(ws > 0)
   count += 1
  elif (np.sum(maes <= tols) == maes.size) & (ncr == 0):
   ncr = nc
   nc = int(np.ceil((float(ncl) + float(ncr))/2))
   count = 0
  elif (np.sum(maes > tols) <= size) & (ncr == 0):
   ncl = nc
   nc = 2*nc
   count = 0
  elif (np.sum(maes <= tols) == maes.size) & (ncr != 0):
   ncr = nc
   nc = int(np.ceil((float(ncl) + float(ncr))/2))
   count = 0
  elif (np.sum(maes > tols) <= size) & (ncr != 0):
   ncl = nc
   nc = int(np.ceil((float(ncl) + float(ncr))/2))
   count = 0
  #4.Save some information for next iteration
  maes0[:] = maes[:]

 return (nc,ws)

def compute_basin_delineation_nbasins(dem,mask,res,nbasins):

 channel_threshold = 10**6
 #Calculate the d8 accumulation area and flow direction
 (area,fdir) = ttf.calculate_d8_acc(dem,res)
 area[mask == 0] = 0.0
 #Iterate until the number of basins match the desired (bisection)
 max_threshold = np.max(area) - res**2
 min_threshold = max_threshold/1000
 #Calculate number of basins for the two boundaries
 channels = ttf.calculate_channels(area,channel_threshold,max_threshold,fdir)
 min_basins = ttf.delineate_basins(channels,mask,fdir)
 min_nbasins = np.unique(min_basins)[1::].size
 #print min_nbasins
 #Min iteration
 channels = ttf.calculate_channels(area,channel_threshold,min_threshold,fdir)
 max_basins = ttf.delineate_basins(channels,mask,fdir)
 max_nbasins = np.unique(max_basins)[1::].size
 for i in xrange(10):
  #Calculate midpoint
  c = (np.log(max_threshold) + np.log(min_threshold))/2
  #Calculate the number of basins for the given threshold
  channels = ttf.calculate_channels(area,channel_threshold,np.exp(c),fdir)
  basins = ttf.delineate_basins(channels,mask,fdir)
  c_nbasins = np.unique(basins)[1::].size
  #print min_nbasins,c_nbasins,max_nbasins
  #Determine if we have found our solution
  if c_nbasins == nbasins:
   return basins
  #Create the new boundaries
  if nbasins < c_nbasins:
   min_threshold = np.exp(c)
   channels = ttf.calculate_channels(area,channel_threshold,min_threshold,fdir)
   max_basins = ttf.delineate_basins(channels,mask,fdir)
   max_nbasins = np.unique(max_basins)[1::].size
  else:
   max_threshold = np.exp(c)
   channels = ttf.calculate_channels(area,channel_threshold,max_threshold,fdir)
   min_basins = ttf.delineate_basins(channels,mask,fdir)
   min_nbasins = np.unique(min_basins)[1::].size

 #print "Did not converge. Returning the best"
 return basins

def define_hrus(basins,dem,channels):

 nbins = 10
 #Define the unique basins
 ubasins = np.unique(basins)[1::]
 tmp = np.copy(basins)
 tmp[:] = 0
 #Create the dem bins
 for basin in ubasins:
  smask = basins == basin
  #Bin the elevation data
  (hist,bins) = np.histogram(dem[smask],bins=nbins)
  #Place the data
  for ibin in xrange(nbins):
   smask = (basins == basin) & (dem >= bins[ibin]) & (dem < bins[ibin+1])
   tmp[smask] = np.mean(dem[smask])
 #import matplotlib.pyplot as plt
 #tmp = np.ma.masked_array(tmp,tmp==0)
 #plt.imshow(tmp)
 #plt.show()

 return 

def calculate_basin_properties(basins,res,latitude,longitude,fdir):

 nb = np.max(basins)
 (ah,lath,lonh,hid,nid) = ttf.calculate_basin_properties(basins,res,nb,fdir,
                               latitude,longitude)
 properties = {
               'area':ah,
               'latitude':lath,
               'longitude':lonh,
               'id':hid,
               'nid':nid
              }

 return properties

def reduce_basin_number(basins,bp,nbasins_goal):

 ids = bp['id']-1
 nids = bp['nid']-1
 area = bp['area']
 nbasins = ids.size

 while nbasins > nbasins_goal:
  #Determine the basin to add
  #To attempt similar area we will first determine which one decreases the area standard deviation
  astd = []
  #Get the smallest 10
  ibs = np.argsort(area)[0:10]
  
  for ib in ibs:
   area_cp = np.copy(area)
   #tmp = np.argmin(area[nids>=0])
   #ib = np.where(area_cp == area_cp[nids>=0][tmp])[0][0]
   area_cp[ids==nids[ib]] = area_cp[ids==nids[ib]] + area_cp[ib]
   astd.append(np.std(area_cp))
  astd = np.array(astd)
  tmp = np.argmin(astd[nids[ibs]>=0])
  ib = ibs[np.where(astd == astd[nids[ibs]>=0][tmp])[0][0]]
 
  #Add it to its next basin
  area[ids==nids[ib]] = area[ids==nids[ib]] + area[ib]
  #Wherever the original basin id was the next basin replace it
  nids[nids == ids[ib]] = nids[ib]
  #Set the basin map to new id
  basins[basins == ids[ib]+1] = nids[ib]+1
  #Remove the row of the basin
  ids = np.delete(ids,ib)
  nids = np.delete(nids,ib)
  area = np.delete(area,ib)
  #Update the number of basins
  nbasins = nbasins - 1

 #Reassign basins
 ubasins = np.unique(basins)[1::]
 for i in xrange(ubasins.size):
  basins[basins == ubasins[i]] = i+1
 
 #Reset the undefined values
 basins[basins <= 0] = -9999

 return basins

def calculate_basin_properties_updated(basins,res,cvs,vars):

 #Initialize properties dictionary
 #vars = ['latitude','longitude','dem','bid']
 properties = {}
 for var in vars:properties[var] = []
 properties['bid'] = []

 #Assemble masks
 masks = {}
 for i in range(basins.shape[0]):
  for j in range(basins.shape[1]):
   h = basins[i,j]
   if h == -9999:continue
   if h not in masks:masks[h] = []
   masks[h].append([i,j])
 for id in masks.keys():
  masks[id] = np.array(masks[id])

 #Iterate through each hillslope to calculate properties
 count = 0
 for uh in masks.keys():
  #tic = time.time()
  #imin = np.min(masks[uh][:,0])
  #imax = np.max(masks[uh][:,0])
  #jmin = np.min(masks[uh][:,1])
  #jmax = np.max(masks[uh][:,1])

  tmp = {}
  for var in vars:
   #print(masks[uh])
   tmp[var] = cvs[var][masks[uh][:,0],masks[uh][:,1]]#imin:imax+1,jmin:jmax+1]
  #tmp = {'latitude':latitude[imin:imax+1,jmin:jmax+1],
  #       'longitude':longitude[imin:imax+1,jmin:jmax+1],
  #       'dem':dem[imin:imax+1,jmin:jmax+1]}
  #tmask = masks[uh][[imin:imax+1,jmin:jmax+1]

  #Add properties to dictionary
  for var in tmp:
   '''#print(uh,var,np.unique(tmp))
   if np.sum(tmp[var] != -9999) > 0:
    properties[var].append(np.mean(tmp[var][tmp[var] != -9999]))
   else:
    properties[var].append(-9999)'''
   if var in ['shreve_order','carea','carea_log10']:
    properties[var].append(np.max(tmp[var]))
   else:
    properties[var].append(np.mean(tmp[var]))
  properties['bid'].append(uh)
  count += 1
  
 #Finalize the properties
 for var in properties:
  properties[var] = np.array(properties[var])

 return properties

def calculate_hillslope_properties_updated(hillslopes,dem,res,latitude,
    longitude,depth2channel,slope,aspect,tas,prec,cdir,uhrt,uhst,lt_uvt,ul_mask):

 #Convert aspect to cartesian coordinates
 x_aspect = np.sin(aspect)
 y_aspect = np.cos(aspect)

 #Initialize properties dictionary
 vars = ['latitude','longitude','dem','aspect','tas','prec','slope',
         'width_intercept','width_slope',
         'length','area','d2c_array','position_array','width_array',
         'relief','x_aspect','y_aspect','hid','relief_a','relief_b',
         'uhrt','uhst','lt_uvt','ul_mask']
 properties = {}
 for var in vars:properties[var] = []

 #Assemble masks
 tic = time.time()
 masks = {}
 for i in xrange(hillslopes.shape[0]):
  for j in xrange(hillslopes.shape[1]):
   h = hillslopes[i,j]
   if h == -9999:continue
   if h not in masks:masks[h] = []
   masks[h].append([i,j])
 for id in masks.keys():
  masks[id] = np.array(masks[id])

 #Iterate through each hillslope to calculate properties
 count = 0
 for uh in masks.keys():
  #tic = time.time()
  imin = np.min(masks[uh][:,0])
  imax = np.max(masks[uh][:,0])
  jmin = np.min(masks[uh][:,1])
  jmax = np.max(masks[uh][:,1])

  #Extract covariates for region 
  shs = np.copy(hillslopes[imin:imax+1,jmin:jmax+1])
  sd2c = np.copy(depth2channel[imin:imax+1,jmin:jmax+1])
  sslope = np.copy(slope[imin:imax+1,jmin:jmax+1])

  #Bin the d2c
  m = shs == uh
  sd2c[~m] = -9999
  nc = min(25,np.ceil(np.sum(m)*res**2/8100.0))
  nc = min(nc,np.unique(sd2c[m]).size)
  if nc > 1:
   tmp = np.sort(sd2c[m])
   bin_edges = tmp[np.arange(0,tmp.size,np.int(np.ceil(float(tmp.size)/(nc+1))))]
   tmp = np.digitize(sd2c[m],bin_edges)
   #X = sd2c[m]
   #X = X[:,np.newaxis]
   #model = sklearn.cluster.KMeans(n_clusters=nc,random_state=35799)
   #tmp = model.fit_predict(X)+1
  else:
   tmp = np.array([1,])
  cls = np.copy(sd2c)
  cls[m] = tmp[:]

  #Reassign the d2c and create a new hillslope
  data = {'slope':[],'d2c':[],'area':[]}
  hillslope = np.zeros(sd2c.shape).astype(np.int32)
  for cl in np.unique(tmp):
   m1 = cls == cl
   if np.sum(m1) == 0:continue
   #Calculate properties
   data['slope'].append(np.mean(sslope[m1]))
   data['d2c'].append(np.mean(sd2c[m1]))
   data['area'].append(res**2*np.sum(m1))
   #Add id
   hillslope[m1] = cl

  #Sort data
  argsort = np.argsort(data['d2c'])
  for var in data:
   data[var] = np.array(data[var])[argsort]

  #Construct position and length arryas
  s = data['slope']
  d2c = data['d2c']
  length = []
  slopes = []
  hand = []
  position = []
  ns = []
  #Ensure there are no zero slopes
  s[s == 0] = 10**-4
  #Use the d2c as the vertices
  for i in xrange(data['d2c'].size):
   #if (i == data['d2c'].size):
   # l = (d2c[i-1]-r)/s[i-1]#/2
   # slp = s[i-1]
   # hand.append(r + l*slp/2)
   # r = r + l*slp
   # slopes.append(slp)
   # pos = pos + l/2
   # position.append(pos)
   if i == 0:
    l = d2c[i]/s[i]#/2
    slp = s[i]
    hand.append(l*slp/2)
    r = l*slp
    slopes.append(slp)
    pos = l/2
    position.append(pos)
   else:
    l = (d2c[i] - r)/((s[i] + s[i-1])/2)
    slp = (s[i] + s[i-1])/2
    hand.append(r + l*slp/2)
    r = r + l*slp
    slopes.append((s[i] + s[i-1])/2)
    pos = pos + l/2
    position.append(pos)
   length.append(l)
  length = np.array(length)
  slopes = np.array(slopes)
  position = np.array(position)
  hand = np.array(hand)
  #Quality control
  if (np.min(length) == 0.0) or (np.max(hand) == 0.0):
   hand = np.array([0.5,1.5])
   length = np.array([10.0,10.0])
   slopes = np.array([0.1,0.1])
   position = np.array([5.0,15.0])
   data['area'] = np.array([900.0,900.0])
  #Place data
  data['position'] = position
  data['length'] = length
  data['slope'] = slopes
  data['d2c'] = hand
  '''length = []
  pos0 = 0
  dtmp = 0
  for i in xrange(data['d2c'].size):
   if (data['d2c'].size == 1):
    ld = d2c[i]/s[i]/2
    lu = ld
    ns.append(s[i])
   elif (i == data['d2c'].size-1):
    ld = ((d2c[i]-d2c[i-1])/((s[i]+s[i-1])/2))/2
    lu = ld
    ns.append((s[i]+s[i-1])/2)
   elif i == 0:
    ld = d2c[i]/s[i]/2
    lu = ((d2c[i+1]-d2c[i])/((s[i+1]+s[i])/2))/2
    ns.append((ld*s[i]+lu*(s[i+1]+s[i]))/(ld + lu))
   else:
    ld = ((d2c[i] - d2c[i-1])/((s[i-1]+s[i])/2))/2
    lu = ((d2c[i+1] - d2c[i])/((s[i]+s[i+1])/2))/2
    ns.append((ld*(s[i-1]+s[i]) + lu*(s[i]+s[i+1]))/(ld + lu))
   #Ensure ld/lu are not 0 (HACK)
   if ld == 0:
    ld = res#1.0 #meter
    ns.append(0.001)
   if lu == 0:
    lu = res#1.0 #meter
    ns.append(0.001)
   dtmp += (ld+lu)*s[i]
   print ld*s[i],(ld+lu)*s[i],s[i],d2c[i],dtmp

   pos = pos0 + ld
   position.append(pos)
   pos0 = pos + lu
   length.append(ld+lu)'''
  #data['position'] = np.array(position)
  #data['length'] = np.array(length)
  #data['nslope'] = np.array(ns)
  
  #Calculate width 
  data['width'] = data['area']/data['length']
  
  #Fit line to width and depth2channel (slope is derivative of the second)
  position = np.array([0.,]+list(data['position'])+[data['length'][-1]/2,])
  w = np.array([data['width'][0],]+list(data['width'])+[data['width'][-1],])
  #s = np.array([0.0,]+list(data['slope']))#+[0.0,])
  #s = np.array([0.0,]+list(data['nslope'])+[0.0,])
  #d2c = np.array([0.0,]+list(data['d2c'])+[data['d2c'][-1]+data['slope'][-1]*data['length'][-1]/2,])
  d2c = np.array([0.0,]+list(data['d2c'])+[data['d2c'][-1],])
  relief = d2c[-1]
  #print data['length']
  #print np.mean(data['nslope'])*np.sum(data['length'])
  #print 't,ap',relief,np.sum(data['slope']*data['length'])
  #Normalize position,width,d2c
  position = position/np.sum(length)
  d2c = d2c/relief
  #w[w > 20] = 20
  if d2c.size == 3:
   #Width
   fw = [0,1]
   #Slope
   #fs = [0,s[1]]
   #relief
   fr = [1.0,1.0]
  else:
   weights = np.cos(np.linspace(-np.pi/4,np.pi/4,position.size-2))
   weights = weights/np.sum(weights)
   #Width
   tmp = w/np.max(w)
   w[tmp > 100] = 100*tmp[tmp > 100] #Limit on width differences
   #popt, pcov = scipy.optimize.curve_fit(fwidth,position,w,bounds=([0.0,-1000],[10**4,1000]))
   z = np.polyfit(position[1:-1],w[1:-1],1,w=weights)
   #fw = [popt[1]/popt[0],1]
   fw = [z[0]/z[1],1]
   if fw[0] > 99:fw[0] = 99
   if fw[0] < -0.99:fw[0] = -0.99
   #Slope
   #z = np.polyfit(position[1:-1],s[1:-1],1,w=weights)
   #if z[0] < -1: z[0] = -1
   #if z[0] > 1: z[0] = 1
   #if z[1] < 0: z[1] = 0
   #if z[1] > 1: z[1] = 1
   #popt, pcov = scipy.optimize.curve_fit(fslope,position,s,bounds=([0.0,-1.0],[1.0,1.0]))
   #fs = [popt[1],popt[0]]
   #fs = [z[0],z[1]]
   #Relief
   #tic = time.time()
   if d2c[1:-1].size > 10:
    try:
     fr, pcov = scipy.optimize.curve_fit(frelief,position[1:-1],d2c[1:-1],bounds=([1.0,1.0],[5.0,5.0]))
    except:
     fr = [1.0,1.0]
   else:
    fr = [1.0,1.0]
  
  tmp = {'latitude':latitude[imin:imax+1,jmin:jmax+1],
         'longitude':longitude[imin:imax+1,jmin:jmax+1],
         'dem':dem[imin:imax+1,jmin:jmax+1],
         'aspect':aspect[imin:imax+1,jmin:jmax+1],
         'tas':tas[imin:imax+1,jmin:jmax+1],
         'prec':prec[imin:imax+1,jmin:jmax+1],
         'slope':slope[imin:imax+1,jmin:jmax+1],
         'x_aspect':x_aspect[imin:imax+1,jmin:jmax+1],
         'y_aspect':y_aspect[imin:imax+1,jmin:jmax+1],
         'uhrt':uhrt[imin:imax+1,jmin:jmax+1],
         'uhst':uhst[imin:imax+1,jmin:jmax+1],
         'lt_uvt':lt_uvt[imin:imax+1,jmin:jmax+1],
         'ul_mask':ul_mask[imin:imax+1,jmin:jmax+1]}

  #Add properties to dictionary
  for var in tmp:
   if np.sum(tmp[var] != -9999) > 0:
    properties[var].append(np.mean(tmp[var][tmp[var] != -9999]))
   else:
    properties[var].append(-9999)
  properties['width_intercept'].append(fw[1])
  #properties['slope_intercept'].append(fs[1])
  properties['width_slope'].append(fw[0])
  #properties['slope_slope'].append(fs[0])
  properties['relief_a'].append(fr[0])
  properties['relief_b'].append(fr[1])
  properties['length'].append(np.sum(data['length']))
  properties['area'].append(float(np.sum(data['area'])))
  properties['relief'].append(relief)
  properties['position_array'].append(position)
  properties['d2c_array'].append(d2c)
  properties['width_array'].append(w)
  properties['hid'].append(uh)
  length = np.sum(data['length'])
  #hslope = np.mean(tmp['slope'][m])#[tmp['slope'] != -9999])
  #fs = data['length']/length
  #print 'slope',hslope
  #print 'crelief',hslope*length
  #if count == 2:exit()
  count += 1
  
 #Finalize the properties
 for var in properties:
  #if var in ['position_array','width_array','d2c_array']:continue
  properties[var] = np.array(properties[var])

 #return properties
 #Write out output
 pickle.dump(properties,open('%s/hillslope_properties.pck' % cdir,'wb'),pickle.HIGHEST_PROTOCOL)
 return

def calculate_hillslope_properties(hillslopes,dem,basins,res,latitude,
    longitude,depth2channel,slope,aspect,cplan,cprof,channels,tas,prec):

 nh = np.max(hillslopes)+1
 (eh,ah,bh,lath,lonh,erange,hid,d2c,slope,haspect,hcplan,hcprof,hmaxd2c,hmind2c,htwidth,hbwidth,htas,hprec) = ttf.calculate_hillslope_properties(hillslopes,dem,basins,res,nh,latitude,longitude,depth2channel,slope,aspect,cplan,cprof,channels,tas,prec)
 properties = {'elevation':eh,
               'area':ah,
               'basin':bh,
               'latitude':lath,
               'longitude':lonh,
               'range':erange,
               'id':hid,
	       'd2c':d2c,
	       'slope':slope,
               #'c2n':hc2n,
               #'g2t':hg2t,
               #'maxsmc':hmaxsmc,
               'aspect':haspect,
               'cplan':hcplan,
               'cprof':hcprof,
               'mind2c':hmind2c,
               'maxd2c':hmaxd2c,
               'twidth':htwidth,
               'bwidth':hbwidth,
               'tas':htas,
               'prec':hprec,
              }

 #Remove nans
 m = np.isnan(eh) == 0
 for p in properties:
  properties[p] = properties[p][m]

 #Ensure the sloep values are not too large
 m = properties['slope'] > 0.4
 properties['slope'][m] = 0.4

 #Compute the hillslope lengths (the 30 accounts for the beginning and end)
 a = properties['maxd2c'] - properties['mind2c']
 b = a/properties['slope'] + res
 properties['length'] = b
 #properties['length'] = (a**2 + b**2)**0.5

 #Compute the ratio of width top to bottom
 properties['twidth'][properties['twidth'] == 0] = 1
 properties['bwidth'][properties['bwidth'] == 0] = 1
 properties['twidth'] = res*properties['twidth']
 properties['bwidth'] = res*properties['bwidth']
 r = properties['twidth']/properties['bwidth']
 #Restrict to 10/1
 m = r > 3
 r[m] = 3
 properties['twidth'][m] = 3*properties['bwidth'][m]
 properties['rwidth'] = r

 return properties

def create_tiles_kmeans(basins,covariates,ntiles):

 #Define the mask
 mask = basins > 0
 
 #Initialize the cluster number
 icluster = 0

 #Initialize the hru map
 hrus = np.empty(covariates[covariates.keys()[0]]['data'].shape).astype(np.int32)
 hrus[:] = -9999

 #Iterate through each hillslope making the hrus
 ub = np.unique(basins)[1::]
 for ib in ub:
  mask = basins == ib

  #Define the data and the bins
  X = []
  for var in covariates:
   X.append(covariates[var]['data'][mask])
  X = np.array(X).T

  #Normalize the data
  for i in xrange(X.shape[1]):
   X[:,i] = (X[:,i]-np.min(X[:,i]))/(np.max(X[:,i])-np.min(X[:,i]))
   
  #Subsample the array
  np.random.seed(1)
  minsamples = 10**5
  if X.shape[0] > minsamples:
   Xf = X[np.random.choice(np.arange(X.shape[0]),minsamples),:]
  else:
   Xf = X

  #Cluster the data
  init = 0.5*np.ones((ntiles,Xf.shape[1]))
  batch_size = 25*ntiles
  init_size = 3*batch_size
  clf = sklearn.cluster.MiniBatchKMeans(ntiles,random_state=1,init=init,batch_size=batch_size,init_size=init_size)
  #clf = sklearn.cluster.KMeans(ntiles,random_state=1)
  clf.fit(Xf)#
  clf_output = clf.predict(X)

  #Map the hrus
  hrus[mask] = clf_output+icluster
 
  #Update icluster
  icluster = np.max(hrus)+1

 #Clean up the hrus
 uhrus = np.unique(hrus)[1::]
 hrus_new = np.copy(hrus)
 for i in xrange(uhrus.size):
  hrus_new[hrus == uhrus[i]] = i
 hrus = hrus_new

 #Finalize hrus array
 hrus[hrus < 0] = -9999

 return hrus

def create_nd_histogram(hillslopes,covariates):

 undef = -9999.0
 #Construct the mask
 m = hillslopes != undef
 for var in covariates:
  m = m & (covariates[var]['data'] != -9999.0)
 
 #Initialize the cluster number
 icluster = -1

 #Initialize the hru map
 hrus = np.empty(covariates[covariates.keys()[0]]['data'].shape).astype(np.float32)
 hrus[:] = -9999

 #Iterate through each hillslope making the hrus
 uh = np.unique(hillslopes)
 uh = uh[uh != -9999]
 for ih in uh:
  mask = (hillslopes == ih) & m

  #Define the data and the bins
  bins,data = [],[]
  for var in covariates:
   bins.append(covariates[var]['nbins'])
   #Convert the data to percentiles if necessary
   if covariates[var]['type'] == 'p':
    tmp = np.copy(covariates[var]['data'][mask])
    argsort = np.argsort(tmp)
    tmp[argsort] = np.linspace(0,1,tmp.size)
    #Have this data replace the covariate information 
    covariates[var]['data'][mask] = tmp
   else:
    tmp = np.copy(covariates[var]['data'][mask])
   data.append(tmp)
   #data.append(covariates[var]['data'][mask])
  bins = np.array(bins)
  data = np.array(data).T

  #Create the histogram
  H,edges = np.histogramdd(data,bins=bins)
  H = H/np.sum(H) 

  #Create a dictionary of hru info
  clusters = {}
  Hfl = H.flat
  for i in xrange(H.size):
   coords = Hfl.coords
   if H[coords] > 0:
    icluster = icluster + 1
    clusters[icluster] = {'pct':H[coords]}
    clusters[icluster]['bounds'] = {}
    for var in covariates:
     key = covariates.keys().index(var)
     clusters[icluster]['bounds'][var] = [edges[key][coords[key]],edges[key][coords[key]+1]]
   Hfl.next()

  #Map the hru id to the grid
  for cid in clusters.keys():
   for id in covariates.keys():
    if covariates.keys().index(id) == 0: string = "(covariates['%s']['data'] >= clusters[%d]['bounds']['%s'][0]) & (covariates['%s']['data'] <= clusters[%d]['bounds']['%s'][1]) & mask" % (id,cid,id,id,cid,id)
    else: string = string +  " & (covariates['%s']['data'] >= clusters[%d]['bounds']['%s'][0]) & (covariates['%s']['data'] <= clusters[%d]['bounds']['%s'][1]) & mask" % (id,cid,id,id,cid,id)
   idx = eval('np.where(%s)' % string)
   hrus[idx] = cid + 1

 #Cleanup the hrus
 hrus = np.array(hrus,order='f').astype(np.int32)
 ttf.cleanup_hillslopes(hrus)
 hrus[hrus >= 0] = hrus[hrus >= 0] + 1

 return hrus

def create_hillslope_tiles(hillslopes,depth2channel,nbins,bins):

 undef = -9999.0
 #Construct the mask
 m = (hillslopes != undef) & (depth2channel != undef)

 #Define the clusters for each hillslope
 clusters = np.copy(hillslopes)
 uh = np.unique(hillslopes)
 uh = uh[uh != undef]
 for ih in uh:
  mask = (hillslopes == ih) & m
  tmp = np.copy(depth2channel[mask])
  argsort = np.argsort(tmp)
  tmp[argsort] = np.linspace(0,1,tmp.size)
  depth2channel[mask] = tmp
  (hist,bins) = np.histogram(tmp,bins=nbins[ih-1])#Change to predefined bins
  for ibin in xrange(nbins[ih-1]):
   #if ibin == 0:smask = mask & (depth2channel >= np.min(tmp)) & (depth2channel <= bins[ih-1][ibin+1])
   #elif ibin == nbins[ih-1]-1:smask = mask & (depth2channel >= bins[ih-1][ibin]) & (depth2channel <= np.max(tmp))
   #else: smask = mask & (depth2channel >= bins[ih-1][ibin]) & (depth2channel <= bins[ih-1][ibin+1])
   smask = mask & (depth2channel >= bins[ibin]) & (depth2channel <= bins[ibin+1])
   clusters[smask] = ibin+1

 #Cleanup the tiles
 clusters = np.array(clusters,order='f').astype(np.int32)
 ttf.cleanup_hillslopes(clusters)
 clusters[clusters >= 0] = clusters[clusters >= 0] + 1

 return clusters

def create_hillslope_tiles_updated(hillslopes,depth2channel,hillslopes_full,hp_in,hp):

 #Construct the lookup 
 lt = {}
 for i in xrange(hp_in['hid'].size):
  h = hp_in['hid'][i]
  lt[h] = hp_in['relief'][i]

 #Normalize the depth to channel by its relief
 lt.keys()
 nrelief = np.copy(depth2channel)
 for i in xrange(hillslopes_full.shape[0]):
  for j in xrange(hillslopes_full.shape[1]):
   h = hillslopes_full[i,j]
   if (h != -9999.0) & (nrelief[i,j] != -9999.0):
    nrelief[i,j] = nrelief[i,j]/lt[h]
 #Clump the max to 1
 nrelief[nrelief > 1.0] = 1.0

 undef = -9999.0
 #Construct the mask
 m = (hillslopes != undef) & (nrelief != undef)

 #Define the clusters for each hillslope
 clusters = np.copy(hillslopes)
 uh = np.unique(hillslopes)
 uh = uh[uh != undef]
 new_hand = np.copy(nrelief)
 for ih in uh:
  #Define the normalized hillslope relief that we want to use
  #p0 = hp['relief_p0'][ih-1]
  #p1 = hp['relief_p1'][ih-1]
  #nr = np.linspace(0,1,2*hp['nbins'][ih-1]+1)
  nr = np.linspace(0,1,2*hp['nbins'][ih-1]+1)[0::2]
  #nrs = nr[0::2]
  #nhand = []
  #for i in xrange(nr.size-1):
  # x = np.linspace(nr[i],nr[i+1],100)
  # nhand.append(np.mean(frelief(x,p0,p1)))
  nhand = hp['relief'][ih-1]*np.array(nr)#nhand)
  #Assign an id to each elevation tile
  mask = (hillslopes == ih) & m
  new_hand[mask] = hp['relief'][ih-1]*new_hand[mask]
  nbins = hp['nbins'][ih-1]
  for ibin in xrange(nbins):
   if ibin == 0:
    #smask = mask & (nrelief <= nr[0::2][ibin+1])
    smask = mask & (new_hand <= nhand[ibin+1])
   elif ibin == nbins-1:
    #smask = mask & (nrelief > nr[0::2][ibin])
    smask = mask & (new_hand > nhand[ibin])
   else:
    #smask = mask & (nrelief > nr[0::2][ibin]) & (nrelief <= nr[0::2][ibin+1])
    smask = mask & (new_hand > nhand[ibin]) & (new_hand <= nhand[ibin+1])
   clusters[smask] = ibin+1

 #Cleanup the tiles
 clusters = np.array(clusters,order='f').astype(np.int32)
 ttf.cleanup_hillslopes(clusters)
 clusters[clusters >= 0] = clusters[clusters >= 0] + 1

 return (clusters,new_hand)

def create_basin_tiles(basin_clusters,hand,basins,dh):

 new_hand = np.copy(hand)
 #Iterate per cluster
 ubcs = np.unique(basin_clusters)
 ubcs = ubcs[ubcs!=-9999]
 tiles = np.copy(hand).astype(np.int32)
 tiles[:] = -9999
 tiles_position = np.copy(hand).astype(np.int32)
 tiles_position[:] = -9999
 count = 0
 
 #Normalize basins and compute maximum relief
 #0.Assemble db of max relief
 db = {}
 for i in range(basins.shape[0]):
  for j in range(basins.shape[1]):
   if basins[i,j] == -9999:continue
   b = basins[i,j]
   if b not in db:db[b] = 0.0
   if hand[i,j] > db[b]:db[b] = hand[i,j]
 #1.Normalize each basin
 maxhand = np.copy(hand)
 for i in range(basins.shape[0]):
  for j in range(basins.shape[1]):
   if basins[i,j] == -9999:continue
   b = basins[i,j]
   hand[i,j] = hand[i,j]/db[b]
   maxhand[i,j] = db[b]
 #Check for nans
 hand[np.isnan(hand) == 1] = 0.0
 #2.Compute basin cluster average max hand
 for ubc in ubcs:
  m = basin_clusters == ubc
  val = np.mean(maxhand[m]) 
  hand[m] = val*hand[m]

 for ubc in ubcs:
  m = basin_clusters == ubc
  data = hand[m]
  #curate
  data[data == -9999] = np.max(data[data != -9999])
  hand[m & (hand == -9999)] = np.max(hand[m & (hand != -9999)])
  #compute number of bins
  nbins = int(np.ceil(np.max(data)/dh))
  #Compute the edges
  pedges = 2.0#1.0#2.5
  bin_edges = np.linspace(0.0,np.max(data)**(1.0/float(pedges)),nbins+1)**pedges
  #Add another bin edge to ensure the channel (hand = 0) is explicitly represented
  bin_edges = np.concatenate((np.array([0.0,]),bin_edges))
  #compute the binning
  #(hist,bin_edges) = np.histogram(data,bins='fd')#bins=nbins)
  #update edges
  #bin_edges[0] = 0.0
  #bin_edges[-1] = np.max(data)
  #Assign the tiles
  count2 = 0
  for i in range(bin_edges.size-1):
   if i == 0:m2 = m & (hand >= bin_edges[i]) & (hand <= bin_edges[i+1])
   else:m2 = m & (hand > bin_edges[i]) & (hand <= bin_edges[i+1])
   #print(i,np.sum(m2))
   if np.sum(m2) > 0:
    tiles[m2] = count
    tiles_position[m2] = count2
    new_hand[m2] = np.mean(hand[m2])
    count += 1
    count2 += 1

 return (tiles,new_hand,tiles_position)

def create_hrus_hydroblocks(hillslopes,htiles,covariates,nclusters,cid):#laura, add cid and flag_fd
 #Compute optimal number of clusters flag
 flag = False
 import sys
 
 #Curate the covariates
 for var in covariates: 
  val = np.mean(covariates[var]['d'][covariates[var]['d'] != -9999])
  covariates[var]['d'][covariates[var]['d'] == -9999] = val
 #import sklearn.cluster
 hrus = np.copy(hillslopes)
 hrus[:] = -9999
 #Iterate through each gru and tile and compute hrus
 uhs = np.unique(hillslopes)
 uhs = uhs[uhs != -9999]
 maxc = 1
 for uh in uhs:
  mh = hillslopes == uh
  uts = np.unique(htiles[mh])
  #Determine fractional coverage of each height band
  fct = []
  numb_pix=[] #laura, FD
  for ut in uts:
   fct.append(np.sum(htiles == ut)/np.sum(mh))
   numb_pix.append(np.sum(htiles == ut))#number of pixesl per height band laura
  fct = np.array(fct)
  numb_pix = np.array(numb_pix)#laura
  #Use fractions to determine the number of clusters
  tnc = np.ceil(uts.size*nclusters*fct).astype(np.int32)
  #Ensure the average number of clusters is consistent with the defined parameter
  while np.mean(tnc) > nclusters:
   tnc[np.argmax(tnc)] = tnc[np.argmax(tnc)] - 1
  
  #Process each tile
  for it in range(uts.size):#ut in uts:
   ut = uts[it]
   mt = mh & (htiles == ut)
   #Compute the parameters for the clustering algorithm
   ccp = {}
   for var in covariates:
    tmp = covariates[var]['d'][mt]
    ccp[var] = {'d':tmp,'t':covariates[var]['t'],
                'min':covariates[var]['min'],
                'max':covariates[var]['max']}
   if flag == True:
    (nc,ws) = compute_cluster_parameters(ccp,maxnc)
   else:
    nc = tnc[it]#nclusters
    ws = np.ones(len(ccp.keys()))
   #print 'hillslope: %d, tile: %d, nc: %d' % (uh,ut,nc)
   #Add weights to covariates
   ##laura: Add weights to covariates so land cover does not overwhem clustering
   geo=0 #laura
   soil=0 #laura
   lc=0 #laura
   other=0 #laura
   for var in covariates:
    if var in ['lats','lons']: #laura
     geo+=1 #laura
    elif var in ['clay','sand','silt','BB','DRYSMC','MAXSMC','REFSMC','SATPSI','SATDK','SATDW','WLTSMC','QTZ']: #laura
     soil+=1 #laura
    elif var in ['lc_w_now','lc_urb_nourb','lc_grass_forest','ndvi']: #laura
     lc+=1 #laura
    else: #laura
     other+=1 #laura
   w_gen=1/((np.sum(geo>0))+(np.sum(soil>0))+np.sum(lc>0)+np.sum(other>0)) #laura

   for var in covariates:
    if var in ['lats','lons']: #laura
     covariates[var]['w'] = w_gen/geo #laura
    elif var in ['clay','sand','silt','BB','DRYSMC','MAXSMC','REFSMC','SATPSI','SATDK','SATDW','WLTSMC','QTZ']: #laura
     covariates[var]['w'] = w_gen/soil #laura
    elif var in ['lc_w_now','lc_urb_nourb','lc_grass_forest','ndvi']: #laura
     covariates[var]['w'] = w_gen/lc #laura
    else: #laura
     covariates[var]['w'] = w_gen/other #laura
    #covariates[var]['w'] = ws[list(ccp).index(var)] #laura: commented out
    
   #prepare the covariate data
   X = []
   for var in covariates:
    #Normalize and apply weight
    tmp = covariates[var]['w']*normalize_variable(covariates[var]['d'][mt],covariates[var]['min'],covariates[var]['max'])
    #tmp[(np.isnan(tmp) == 1) | (np.isinf(tmp) == 1)] = 0.0
    #Convert to percentiles
    #argsort = np.argsort(tmp)
    #tmp[argsort] = np.linspace(0,1,tmp.size)
    X.append(tmp)
    #cluster the data
   X = np.array(X).T
   clusters = cluster_data(X,nc)+maxc
   
   hrus[mt] = clusters
   maxc = np.max(clusters)+1
 #Cleanup hrus
 ttf.cleanup_hillslopes(hrus)
 hrus[hrus >= 0] = hrus[hrus >= 0] + 1

 return hrus
 
 return

def create_hrus(hillslopes,htiles,covariates,nclusters,flag,maxnc,cdir):

 #Curate the covariates
 for var in covariates: 
  val = np.mean(covariates[var]['d'][covariates[var]['d'] != -9999])
  covariates[var]['d'][covariates[var]['d'] == -9999] = val
 import sklearn.cluster
 hrus = np.copy(hillslopes)
 hrus[:] = -9999
 #Iterate through each hillslope and tile and compute hrus
 uhs = np.unique(hillslopes)
 uhs = uhs[uhs != -9999]
 maxc = 1
 for uh in uhs:
  mh = hillslopes == uh
  uts = np.unique(htiles[mh])
  for ut in uts:
   mt = mh & (htiles == ut)
   #Compute the parameters for the clustering algorithm
   ccp = {}
   for var in covariates:
    tmp = covariates[var]['d'][mt]
    ccp[var] = {'d':tmp,'t':covariates[var]['t'],
                'min':covariates[var]['min'],
                'max':covariates[var]['max']}
   if flag == True:
    (nc,ws) = compute_cluster_parameters(ccp,maxnc)
   else:
    nc = nclusters
    ws = np.ones(len(ccp.keys()))
   #print 'hillslope: %d, tile: %d, nc: %d' % (uh,ut,nc)
   #Add weights to covariates
   for var in covariates:
    covariates[var]['w'] = ws[ccp.keys().index(var)]
   #prepare the covariate data
   X = []
   for var in covariates:
    #Normalize and apply weight
    tmp = covariates[var]['w']*normalize_variable(covariates[var]['d'][mt],covariates[var]['min'],covariates[var]['max'])
    #tmp[(np.isnan(tmp) == 1) | (np.isinf(tmp) == 1)] = 0.0
    #Convert to percentiles
    #argsort = np.argsort(tmp)
    #tmp[argsort] = np.linspace(0,1,tmp.size)
    X.append(tmp)
   #cluster the data
   X = np.array(X).T
   clusters = cluster_data(X,nc)+maxc
   #state = 35799
   #if (X.shape[0] >= nc):
   # model = sklearn.cluster.KMeans(n_clusters=nc,random_state=state)
   # clusters = model.fit_predict(X)+maxc
   #else:
    #clusters = np.zeros(X.shape[0])+maxc
   hrus[mt] = clusters
   maxc = np.max(clusters)+1

 #Cleanup hillslopes
 ttf.cleanup_hillslopes(hrus)
 hrus[hrus >= 0] = hrus[hrus >= 0] + 1

 #return hrus
 #Write out output
 pickle.dump(hrus,open('%s/hrus.pck' % cdir,'wb'),pickle.HIGHEST_PROTOCOL)
 
 return

def calculate_hru_properties(hillslopes,tiles,channels,res,nhillslopes,hrus,depth2channel,slope,basins,cdir):

 tmp = np.unique(hrus)
 tmp = tmp[tmp != -9999]
 nhru = tmp.size
 (wb,wt,l,hru_position,hid,tid,hru,hru_area,hru_dem,hru_slope) = ttf.calculate_hru_properties(hillslopes,tiles,channels,basins,nhru,res,nhillslopes,hrus,depth2channel,slope)
 #Curate (everyone must have info)
 hru_properties = {'width_bottom':wb,
                   'width_top':wt,
                   'hillslope_length':l,
                   'hillslope_position':hru_position,
                   'hillslope_id':hid,
                   'tile_id':tid,
                   'hru':hru,
                   'area':hru_area,
                   'slope':hru_slope,
                   'depth2channel':hru_dem,
                   'frac':hru_area/np.sum(hru_area)}

 #return hru_properties
 #Write out output
 pickle.dump(hru_properties,open('%s/hru_properties.pck' % cdir,'wb'),pickle.HIGHEST_PROTOCOL)

 return

def calculate_hru_properties_updated(hillslopes,tiles,res,hrus,depth2channel,slope,hp,cdir,nhand):

 #Get the hillslope fractions
 fs = []
 for ih in xrange(hp['hid'].size):
  hid = int(hp['hid'][ih])
  f = np.sum(hillslopes == hid)/float(hillslopes.size)
  fs.append(f)
 fs = np.array(fs)
 hp['frac'] = fs/np.sum(fs) #fix to 1

 #Assemble masks
 masks = {}
 for i in xrange(hillslopes.shape[0]):
  for j in xrange(hillslopes.shape[1]):
   hru = hrus[i,j]
   if hru == -9999.0:continue
   if hru not in masks:masks[hru] = []
   masks[hru].append([i,j])
 for hru in masks:
  masks[hru] = np.array(masks[hru])

 #Gather some general hru information
 hru_properties = {}
 vars = ['hillslope_id','tile_id','hru','area','hillslope_slope','hand_ecdf','hand_bedges']
 for var in vars:hru_properties[var] = []
 for hru in masks:
  iss = masks[hru][:,0]
  jss = masks[hru][:,1]
  hru_properties['hillslope_id'].append(int(np.mean(hillslopes[iss,jss])))
  hru_properties['tile_id'].append(int(np.mean(tiles[iss,jss])))
  hru_properties['hru'].append(int(hru))
  hru_properties['area'].append(np.float64(res**2*np.sum(iss.size)))
  hru_properties['hillslope_slope'].append(np.float64(np.mean(slope[iss,jss])))
  #Assign ecdf of hand - mean(hand) per hru
  #tmp = depth2channel[iss,jss]
  tmp = nhand[iss,jss]
  if np.sum(tmp != -9999) == 0:
   tmp[tmp == -9999] = 0.0
  else:
   tmp[tmp == -9999] = np.mean(tmp[tmp != -9999])
  #Compute the ecdf
  nbins = 10
  (hist,bin_edges) = np.histogram(tmp,bins=nbins)
  ecdf = np.cumsum(hist).astype(np.float32)
  ecdf = ecdf/ecdf[-1]
  ecdf = np.append(np.zeros(1),ecdf)
  hru_properties['hand_ecdf'].append(ecdf)
  hru_properties['hand_bedges'].append(bin_edges)
 for var in hru_properties:
  hru_properties[var] = np.array(hru_properties[var])
 hru_properties['frac'] = np.zeros(hru_properties['area'].size)
 #hru_properties['frac'] = np.zeros(hru_properties['area']/np.sum(hru_properties['area'])

 #Add fill for the other properties
 vars = ['hillslope_length','hillslope_hand','hillslope_position','hillslope_width','hillslope_frac',
         'soil_depth','depth_to_bedrock']
 for var in vars:
  hru_properties[var] = np.zeros(hru_properties['area'].size).astype(np.float64)

 #Associate the hillslope properties
 for ih in xrange(hp['hid'].size):
  hid = int(hp['hid'][ih])
  #print hid
  m = hru_properties['hillslope_id'] == hid
  #Extract the tids
  (tids,idx) = np.unique(hru_properties['tile_id'][m],return_inverse=True)
  #Compute the normalized relief
  nrelief = np.linspace(0,1,2*tids.size+1)[0::2]
  #Compute the correposponding lengths
  p0 = hp['relief_p0'][ih]
  p1 = hp['relief_p1'][ih]
  length = hp['length'][ih]*(frelief_inv(nrelief[1:],p0,p1) - frelief_inv(nrelief[0:-1],p0,p1))
  #Compute the relief for each segment
  #hand = []
  #for i in xrange(nrelief.size-1):
  # x = np.linspace(nrelief[i],nrelief[i+1],100)
  # hand.append(np.mean(frelief(x,p0,p1)))
  #hand = hp['relief'][ih]*np.array(hand)
  hand = hp['relief'][ih]*(nrelief[0:-1]+nrelief[1:])/2
  #Compute the width for each segment
  pos = frelief_inv(nrelief,p0,p1)
  p0 = hp['width_p0'][ih]
  width = (fwidth(pos[1:],p0) + fwidth(pos[0:-1],p0))/2
  #Convert all to float64
  length = length.astype(np.float64)
  width = width.astype(np.float64)
  #Correct to the actual area (This is necessary for conservation in div_it...)
  #print idx
  #print length
  #print width
  #tmp = hru_properties['area'][m][idx]
  #tmp = np.sum(length*width)*tmp/np.sum(tmp)
  #r = tmp/(length*width)
  #length = r*length
  hand = hand.astype(np.float64)
  #Compute the fractions
  frac = (width*length)/np.sum(width*length)
  #Compute the positions
  positions = np.linspace(0,1,2*tids.size+1)[1::2]
  #Place the properties
  hru_properties['hillslope_length'][m] = length[idx]
  hru_properties['hillslope_hand'][m] = hand[idx]
  hru_properties['hillslope_position'][m] = positions[idx]
  hru_properties['hillslope_width'][m] = width[idx]
  #Compute the hillslope fraction
  for it in xrange(tids.size):
   m1 = m & (hru_properties['tile_id'] == tids[it])
   f = hru_properties['area'][m1]/np.sum(hru_properties['area'][m1])
   hru_properties['hillslope_frac'][m1] = frac[it]*f
  #Set the overall fraction
  hru_properties['frac'][m] = hp['frac'][ih]*hru_properties['hillslope_frac'][m]
  #print ih,tids.size
  #Determine if hillslope is in the lowlands or uplands (per pelletier 2016)
  if hp['ul_mask'][ih] >= 1.5: #LOWLAND
   soil_thickness = 2.0
   sedimentary_thickness =  hp['lt_uvt'][ih] - soil_thickness
   if sedimentary_thickness < 0:sedimentary_thickness = 0.0
   soil_depth = soil_thickness*np.ones(tids.size)
   depth_to_bedrock = (soil_thickness + sedimentary_thickness)*np.ones(tids.size)
  elif hp['ul_mask'][ih] < 1.5: #UPLAND
   soil_thickness = np.linspace(2.0,hp['uhst'][ih],tids.size)
   regolith_thickness = np.linspace(hp['lt_uvt'][ih],hp['uhrt'][ih],tids.size)
   soil_depth = soil_thickness
   depth_to_bedrock = regolith_thickness#soil_thickness + regolith_thickness
  hru_properties['soil_depth'][m] = soil_depth[idx]
  hru_properties['depth_to_bedrock'][m] = depth_to_bedrock[idx]

 #return hru_properties
 #Write out output
 pickle.dump(hru_properties,open('%s/hru_properties.pck' % cdir,'wb'),pickle.HIGHEST_PROTOCOL)

 return
                       
#def cluster_hillslopes(hp,hillslopes,nclusters,covariates):
def cluster_hillslopes(hillslopes,covariates,hp_in,nclusters,ws):

 #Add weights to covariates
 for var in covariates:
  covariates[var]['w'] = ws[covariates.keys().index(var)]

 #import sklearn.cluster
 X = []
 for var in covariates:
  otmp = np.copy(covariates[var]['d'])
  otmp[(np.isnan(otmp) == 1) | (np.isinf(otmp) == 1)] = 0.0
  tmp = np.copy(otmp)
  #Normalize and apply weight
  tmp = covariates[var]['w']*(tmp-np.min(tmp))/(np.max(tmp)-np.min(tmp))
  #Convert to percentiles
  #argsort = np.argsort(tmp)
  #tmp[argsort] = np.linspace(0,1,tmp.size)
  #Group all the 0s together
  #tmp[otmp == 0.0] = np.mean(tmp[otmp == 0.0])
  X.append(tmp)
 X = np.array(X).T
 #state = 35799#80098
 clusters = cluster_data(X,nclusters)+1
 #model = sklearn.cluster.KMeans(n_clusters=nclusters,random_state=state)
 #clusters = model.fit_predict(X)+1
 #Clean up the hillslopes
 hillslopes = np.array(hillslopes,order='f').astype(np.int32)
 ttf.cleanup_hillslopes(hillslopes)
 #Assign the new ids to each hillslpe
 hillslopes_clusters = ttf.assign_clusters_to_hillslopes(hillslopes,clusters)
 #Determine the number of hillslopes per cluster
 uclusters = np.unique(clusters)
 nhillslopes = []
 for cluster in uclusters:
  nhillslopes.append(np.sum(clusters == cluster))
 nhillslopes = np.array(nhillslopes)

 #Compute the average value for each cluster of each property
 hp_out = {}
 hp_out['hid'] = []
 for cluster in uclusters:
  hp_out['hid'].append(cluster)
  m = clusters == cluster
  #Compute fraction dependent on area of hillslope
  frac = hp_in['area'][m]/np.sum(hp_in['area'][m])
  for var in hp_in:
   if var not in hp_out:hp_out[var] = []
   #hp_out[var].append(np.median(hp_in[var][m]))
   hp_out[var].append(np.sum(frac*hp_in[var][m]))
  #Calculate the fraction
  if 'frac' not in hp_out:hp_out['frac'] = []
  hp_out['frac'].append(np.sum(hp_in['area'][m])/np.sum(hp_in['area']))
 for var in hp_out:
  hp_out[var] = np.array(hp_out[var])
 
 return (hillslopes_clusters,nhillslopes,hp_out)

def cluster_hillslopes_updated(hillslopes,covariates,hp_in,nclusters,ws,dh,max_nbands,min_nbands):

 #Add weights to covariates
 for var in covariates:
  covariates[var]['w'] = ws[covariates.keys().index(var)]

 X = []
 for var in covariates:
  otmp = np.copy(covariates[var]['d'])
  otmp[(np.isnan(otmp) == 1) | (np.isinf(otmp) == 1)] = 0.0
  tmp = np.copy(otmp)
  #Normalize and apply weight
  tmp = covariates[var]['w']*normalize_variable(tmp,covariates[var]['min'],covariates[var]['max'])
  #tmp = covariates[var]['w']*(tmp-np.min(tmp))/(np.max(tmp)-np.min(tmp))
  X.append(tmp)
 X = np.array(X).T
 clusters = cluster_data(X,nclusters)+1
 #Clean up the hillslopes
 hillslopes = np.array(hillslopes,order='f').astype(np.int32)
 ttf.cleanup_hillslopes(hillslopes)
 #Assign the new ids to each hillslpe
 hillslopes_clusters = ttf.assign_clusters_to_hillslopes(hillslopes,clusters)
 #Determine the number of hillslopes per cluster
 uclusters = np.unique(clusters)
 #nhillslopes = []
 #for cluster in uclusters:
 # nhillslopes.append(np.sum(clusters == cluster))
 #nhillslopes = np.array(nhillslopes)

 #Compute the average value for each hillslope of each property
 hp_out = {}
 hp_out['hid'] = []
 for cluster in uclusters:
  hp_out['hid'].append(cluster)
  m = clusters == cluster
  #Compute fraction dependent on area of hillslope
  frac = hp_in['area'][m]/np.sum(hp_in['area'][m])
  for var in hp_in:
   if var in ['position_array','width_array','d2c_array','hid']:continue
   if var not in hp_out:hp_out[var] = []
   hp_out[var].append(np.sum(frac*hp_in[var][m]))
  #Calculate the fraction
  if 'frac' not in hp_out:hp_out['frac'] = []
  hp_out['frac'].append(np.sum(hp_in['area'][m])/np.sum(hp_in['area']))

 #Compute the average width and d2c function
 vars = ['relief_p0','relief_p1','width_p0','w','p','d']
 for var in vars:
  if var in ['w','p','d']:
   hp_out[var] = {}
  else:
   hp_out[var] = []
 for cluster in uclusters:
  '''ids = np.where(clusters == cluster)[0]
  d = []
  p = []
  w = []
  for id in ids:
   print id
   d = d + list(hp_in['d2c_array'][id])
   #w = w + (1 + list(hp_in['width_array'][id])
   w = w + list(1 + hp_in['position_array'][id]*hp_in['width_slope'][id])
   p = p + list(hp_in['position_array'][id])
  d = np.array(d)
  w = np.array(w)'''
  mc = np.where(clusters == cluster)[0]
  d = hp_in['d2c_array'][mc]
  #print hp_in['position_array'][mc].shape
  #print hp_in['width_slope'][mc].shape
  if len(hp_in['position_array'][mc].shape) > 1:
   w = 1 + hp_in['position_array'][mc,:]*hp_in['width_slope'][mc][:,np.newaxis]
  else:
   w = 1 + hp_in['position_array'][mc]*hp_in['width_slope'][mc]
  p = hp_in['position_array'][mc]
  p = np.concatenate(p)
  d = np.concatenate(d)
  w = np.concatenate(w)
  hp_out['p'][cluster-1] = p
  hp_out['d'][cluster-1] = d
  hp_out['w'][cluster-1] = w
  #Fit curve to d2c
  #fr, pcov = scipy.optimize.curve_fit(frelief,p,d)#,bounds=([0.0,-1000],[10**4,1000]))
  #try:
  try:
   fr, pcov = scipy.optimize.curve_fit(frelief,p,d,bounds=([1.0,1.0],[5.0,5.0]))
  except:
   fr = [1.0,1.0]
  hp_out['relief_p0'].append(fr[0])
  hp_out['relief_p1'].append(fr[1])
  #Fit line to width
  try:
   fw, pcov = scipy.optimize.curve_fit(fwidth,p,w,bounds=([-0.99,],[99,]))
  except:
   fw = [1.0,]
  hp_out['width_p0'].append(fw[0])
  #plt.plot(p,d,'bo',alpha=0.05)
  #plt.plot(p,frelief(p,fr[0],fr[1]),'ro',alpha=0.05)
  #plt.show()
  
 #Convert to arrays
 for var in hp_out:
  if var in ['p','d','w']:continue
  hp_out[var] = np.array(hp_out[var])

 #Define the number of elevation tiles per cluster
 tile_relief = dh#md['clustering']['tile_relief']
 max_ntiles = max_nbands#md['clustering']['max_ntiles']
 min_ntiles = min_nbands#md['clustering']['max_ntiles']
 nbins = np.round(hp_out['relief']/tile_relief).astype(np.int)
 nbins[nbins < min_ntiles] = min_ntiles
 nbins[nbins > max_ntiles] = max_ntiles
 hp_out['nbins'] = nbins

 #Set some constraints
 m = hp_out['length'] > 10000
 hp_out['length'][m] = 10000

 return (hillslopes_clusters,hp_out)

def cluster_basins_updated(basins,covariates,hp_in,nclusters):

 X = []
 for var in covariates:
  otmp = np.copy(covariates[var]['d'])
  otmp[(np.isnan(otmp) == 1) | (np.isinf(otmp) == 1)] = 0.0
  tmp = np.copy(otmp)
  #Normalize and apply weight
  tmp = normalize_variable(tmp,covariates[var]['min'],covariates[var]['max'])
  X.append(tmp)
 X = np.array(X).T
 clusters = cluster_data(X,nclusters)+1
 #Create the mapping
 mapping = np.zeros(np.max(hp_in['bid'])+1)
 mapping[:] = -9999
 for i in range(hp_in['bid'].size):
  #print(i,hp_in['bid'][i],clusters[i])
  mapping[hp_in['bid'][i]] = clusters[i]
 #print(mapping)
 #print(hp_in['bid'])
 #print(clusters)
 #exit()
 #Clean up the basins
 #basins = np.array(basins,order='f').astype(np.int32)
 #ttf.cleanup_hillslopes(basins)
 #Assign the new ids to each hillslpe
 basins_clusters = ttf.assign_clusters_to_hillslopes(basins,mapping)
 #Determine the number of basins per cluster
 uclusters = np.unique(clusters)

 return (basins_clusters,)#hp_out)

def curate_hru_properties(hru_properties,hp):

 #hp['length'][:] = 1000.0
 #hp['slope'][:] = 0.1
 #hp['rwidth'][:] = 1.0
 #Iterate per hillslope
 hru_properties['wspec'] = np.copy(hru_properties['slope'])
 hru_properties['wspec'][:] = 0.0
 for hid in hp['hid']:
  hid = int(hid)
  m = hru_properties['hillslope_id'] == hid
  if np.sum(m) == 0:continue
  #redo the length
  (d2c,idx) = np.unique(hru_properties['depth2channel'][m],return_inverse=True)
  #Calculate the update properties for the elevation tiles
  hlength = hp['length'][hid-1]/d2c.size*np.ones(d2c.size)
  width = np.linspace(1,hp['rwidth'][hid-1],d2c.size+1)
  w0 = (width[1:]+width[0:-1])/2
  #Compute the fractions
  f0 = hlength*w0/np.sum(hlength*w0)
  f1 = []
  tids = hru_properties['tile_id'][m]
  for tid in np.unique(tids):
   m1 = (tids == tid)
   f1.append(np.sum(hru_properties['area'][m][m1]/np.sum(hru_properties['area'][m])))
  f1 = np.array(f1)
  #Correct the width and length to ensure f1 is met
  hlength = (f1/f0)**0.5*hlength
  w1 = (f1/f0)**0.5*w0
  hpos = np.cumsum(hlength) - hlength[0]/2
  helev = hp['slope'][hid-1]*hpos
  slope = hp['slope'][hid-1]*np.ones(d2c.size)
  twidth = w1/w0*width[1:]
  bwidth = w1/w0*width[0:-1]
  #Split up the relevant information among the intra-elevation tiles
  utids = np.unique(tids)
  t1,b1 = [],[]
  for it in xrange(utids.size):
   m1 = tids == utids[it]
   f = hru_properties['area'][m][m1]/np.sum(hru_properties['area'][m][m1])
   t1 = t1 + list(twidth[it]*f)
   b1 = b1 + list(bwidth[it]*f)
  t1 = np.array(t1)
  b1 = np.array(b1)
  wspec = (t1+b1)/2

  #Place the parameters
  hru_properties['hillslope_length'][m] = hlength[idx]
  hru_properties['slope'][m] = slope[idx]
  hru_properties['depth2channel'][m] = helev[idx]
  hru_properties['hillslope_position'][m] = hpos[idx]
  hru_properties['width_top'][m] = twidth[idx]#twidth
  hru_properties['width_bottom'][m] = bwidth[idx]#bwidth
  hru_properties['wspec'][m] = wspec[:] #weighted width

 return hru_properties

def polygonize_raster(data):

 din = np.copy(data,order='F')
 dout = np.copy(data,order='F')
 dout[:] = -9999
 ttf.polygonize_raster(din,dout)

 return dout

def compute_polygon_info(polygons,clusters,res):

 #Construct x,y using resolution 
 x = np.linspace(res/2,polygons.shape[0]*res - res/2,polygons.shape[0])
 y = np.linspace(res/2,polygons.shape[1]*res - res/2,polygons.shape[1])
 (xs,ys) = np.meshgrid(y,x)
 xs = np.copy(xs,order='F')
 ys = np.copy(ys,order='F')

 #Assemble i/o arrays
 pcxy = np.zeros((np.int(np.max(polygons)) + 1,2),order='F').astype(np.float32)
 pd2o = np.zeros((4*xs.size,2),order='F').astype(np.float32)
 cd2o = np.zeros((4*xs.size,2),order='F').astype(np.float32)
 pd2o[:] = -9999
 cd2o[:] = -9999

 #Compute properties
 ttf.compute_polygon_info(polygons,clusters,xs,ys,pcxy,pd2o,cd2o)
 pd2o = pd2o[pd2o[:,0]!=-9999,:]
 cd2o = cd2o[cd2o[:,0]!=-9999,:]
 pmatrix = scipy.sparse.coo_matrix((np.ones(pd2o.shape[0]),(pd2o[:,0],pd2o[:,1])),shape=(np.int(np.max(polygons)) + 1,np.int(np.max(polygons)) + 1),dtype=np.float32)
 pmatrix = pmatrix.tocsr()

 #Compute length between centroids
 tmp = scipy.sparse.find(pmatrix)
 dist = ((pcxy[tmp[0],0] - pcxy[tmp[1],0])**2 + (pcxy[tmp[0],1] - pcxy[tmp[1],1])**2)**0.5

 #Assemble list of polygons per cluster
 ucls = np.unique(clusters)
 ucls = ucls[ucls != -9999]
 db = {}
 for cl in ucls:
  db[cl] = np.unique(polygons[clusters == cl])

 #Compute average length per cluster
 #for cl in db:
 # print(cl,np.mean(


 #Create database
 db = {'centroid':pcxy,}

 return db

def calculate_channel_properties(channels,channel_topology,slope,eares,mask,area_all,area_all_cp,basins,pscaling):

 #Compute channel properties
 (channel_slope,channel_length,channel_mannings,channel_width,channel_bankfull,channel_topology,channel_area,reach_area,floodplain_mannings) = calculate_channel_properties_workhorse(channels,channel_topology,slope,eares,mask,area_all,area_all_cp,basins,pscaling['channel_manning'],pscaling['floodplain_manning'],pscaling['bankfull_depth'],pscaling['channel_width'])

 #Assemble final database
 db_channels = {'slope':channel_slope,'length':channel_length,'manning_channel':channel_mannings,
                'width':channel_width,'bankfull':channel_bankfull,'topology':channel_topology,
                'acc':channel_area,'area':reach_area,
                'manning_floodplain':floodplain_mannings}

 return copy.deepcopy(db_channels)

@numba.jit(nopython=True,cache=True)
def calculate_channel_properties_workhorse(channels,channel_topology,slope,eares,mask,area_all,area_all_cp,basins,m_c_n,m_fp_n,m_bd,m_cw):

 #Convert to 0-index
 channel_topology[channel_topology > 0] = channel_topology[channel_topology > 0] - 1

 #Determine total number of channels
 nc = channel_topology.size
 channel_slope = np.zeros(nc)
 channel_length = np.zeros(nc)
 channel_area = np.zeros(nc)
 channel_area_cp = np.zeros(nc)
 reach_area = np.zeros(nc)
 count = np.zeros(nc)
 for i in range(channels.shape[0]):
  for j in range(channels.shape[1]):
   channel = channels[i,j]
   basin = basins[i,j]
   if ((basin == -9999) | (basin == 0)):continue
   reach_area[basin-1] += eares**2
   if ((channel == -9999) | (channel == 0)):continue
   channel_slope[channel-1] += slope[i,j]
   channel_length[channel-1] += eares
   if area_all[i,j] > channel_area[channel-1]:channel_area[channel-1] = area_all[i,j]
   if area_all_cp[i,j] > channel_area_cp[channel-1]:channel_area_cp[channel-1] = area_all_cp[i,j]
   count[channel-1] += 1

 #Compute averages
 channel_slope = channel_slope/count

 #Add others
 channel_mannings = m_c_n*0.035*np.ones(nc) #Assume natural gravel stream bed
 floodplain_mannings = m_fp_n*0.15*np.ones(nc) #Assume natural gravel stream bed
 #Compute bankfull width and height (only applicable for United States ();
 channel_width = m_cw*2.70*(channel_area_cp*10**-6)**0.352#30.0*np.ones(nc) #meter
 #channel_width[channel_width > 50] = 50 #NWC -> hack to make sure model works (03/31/23)
 #channel_width = 90*np.ones(nc) #meter
 tmp = reach_area/channel_length
 channel_width[channel_width > tmp] = tmp[channel_width > tmp]
 channel_bankfull = m_bd*0.30*(channel_area_cp*10**-6)**0.213#1.0*np.ones(nc) #meter
 #channel_bankfull = 1.0*np.ones(nc) #meter
 #Compute bankfull width and height (only applicable for Atlantic Plain of US)
 #channel_width = 2.22*(channel_area*10**-6)**0.363#30.0*np.ones(nc) #meter
 #channel_bankfull = 0.24*(channel_area*10**-6)**0.323#1.0*np.ones(nc) #meter

 return (channel_slope,channel_length,channel_mannings,
         channel_width,channel_bankfull,channel_topology,channel_area,reach_area,
         floodplain_mannings)

def calculate_inlets_oulets(channels_wob,fdir,area,mask,lats,lons,mask_all,area_all):

 (inlet_org,inlet_dst,outlet_org,outlet_dst) = calculate_inlets_oulets_workhorse(channels_wob,fdir,area,mask,lats,lons,mask_all,area_all)
 db_io = {'inlet':{'org':inlet_org,'dst':inlet_dst},
          'outlet':{'org':outlet_org,'dst':outlet_dst}}

 return db_io

@numba.jit(nopython=True,cache=True)
def calculate_inlets_oulets_workhorse(channels_wob,fdir,area,mask,lats,lons,mask_all,area_all):

 #Construct positions array
 positions = np.zeros((8,2)).astype(np.int32)
 pos = 0
 for k in range(-1,2):
  for l in range(-1,2):
   if (k == 0) & (l == 0):continue
   positions[pos,0] = k
   positions[pos,1] = l
   pos = pos + 1
 i_inlet = -1
 i_outlet = -1
 inlet_org = np.zeros((250,5))
 inlet_dst = np.zeros((250,5))
 outlet_org = np.zeros((250,5))
 outlet_dst = np.zeros((250,5))
 for i in range(channels_wob.shape[0]):
  for j in range(channels_wob.shape[1]):
   if (channels_wob[i,j] <= 0) :continue
   #Determine inlets
   for ipos in range(positions.shape[0]):
    inew = i+positions[ipos,0]
    jnew = j+positions[ipos,1]
    if ((inew < 0) | (jnew < 0) | (inew >= channels_wob.shape[0]) | (jnew >= channels_wob.shape[1])):continue
    if ((fdir[inew,jnew,0] == i+1) & (fdir[inew,jnew,1] == j+1)) & (mask[inew,jnew] == False):
     #if ((area[inew,jnew] > 10**4) & (mask_all[inew,jnew] != -9999)):
     if ((area[inew,jnew] > 10**5) & (mask_all[inew,jnew] != -9999)):
      i_inlet += 1
      #Determine origin and destination grid id, lat, and lon
      inlet_org[i_inlet,0] = mask_all[inew,jnew]
      inlet_org[i_inlet,1] = lats[inew,jnew]
      inlet_org[i_inlet,2] = lons[inew,jnew]
      inlet_org[i_inlet,3] = -9999
      inlet_org[i_inlet,4] = area_all[inew,jnew]
      inlet_dst[i_inlet,0] = mask_all[i,j]
      inlet_dst[i_inlet,1] = lats[i,j]
      inlet_dst[i_inlet,2] = lons[i,j]
      inlet_dst[i_inlet,3] = channels_wob[i,j]-1
      inlet_dst[i_inlet,4] = area_all[i,j]
   #Determine outlets
   inew = fdir[i,j,0]-1
   jnew = fdir[i,j,1]-1
   if ((area[inew,jnew] > 10**4) & (mask[inew,jnew] == False) & (mask_all[inew,jnew] != -9999)):
      #print(inew,jnew,inew.dtype,jnew.dtype,area.dtype)
      #if (area[inew,jnew] > 10**5):
      # if (mask[inew,jnew] == False):
      #  if (mask_all[inew,jnew] != -9999):
      i_outlet += 1
      #Determine origin and destination grid id, lat, and lon
      outlet_dst[i_outlet,0] = mask_all[inew,jnew]
      outlet_dst[i_outlet,1] = lats[inew,jnew]
      outlet_dst[i_outlet,2] = lons[inew,jnew]
      outlet_dst[i_outlet,3] = -9999
      outlet_dst[i_outlet,4] = area_all[inew,jnew]
      outlet_org[i_outlet,0] = mask_all[i,j]
      outlet_org[i_outlet,1] = lats[i,j]
      outlet_org[i_outlet,2] = lons[i,j]
      outlet_org[i_outlet,3] = channels_wob[i,j]-1
      outlet_org[i_outlet,4] = area_all[i,j]

 inlet_dst = inlet_dst[0:i_inlet+1,:]
 inlet_org = inlet_org[0:i_inlet+1,:]
 outlet_dst = outlet_dst[0:i_outlet+1,:]
 outlet_org = outlet_org[0:i_outlet+1,:]

 return (inlet_org,inlet_dst,outlet_org,outlet_dst)

@numba.jit(nopython=True,cache=True)
def determine_band_position(nhand,bin_edges,basin_clusters):
    #Determine the position of the height band with respect to the channel
    band_position = np.copy(nhand).astype(np.int32)
    basin_clusters = basin_clusters.astype(np.int32)
    band_position[:] = -9999
    for i in range(nhand.shape[0]):
        for j in range(nhand.shape[1]):
            if (nhand[i,j] == -9999) | (basin_clusters[i,j] == -9999):continue
            bc = basin_clusters[i,j]-1
            for k in range(bin_edges[bc,:].size-1):
                if bin_edges[bc,k+1] == -9999:break
                if (nhand[i,j] >= bin_edges[bc,k]) & (nhand[i,j] <= bin_edges[bc,k+1]):
                    band_position[i,j] = k
                    break
    return band_position

@numba.jit(nopython=True,cache=True)
def curate_band_position(band_position,basins):
    basins = basins.astype(np.int32)
    #Determine the count of each id per basin
    ubs = np.unique(basins)
    ubs = ubs[ubs != -9999]
    obdb = np.zeros((ubs.size,np.max(band_position)+1))
    for i in range(band_position.shape[0]):
        for j in range(band_position.shape[1]):
            if band_position[i,j] == -9999:continue
            basin = basins[i,j]-1
            bp = band_position[i,j]
            obdb[basin,bp] += 1
            
    #Iterate per basin count and look for problematic ones and assign a mapping to correct
    mapping = np.zeros((ubs.size,np.max(band_position)+1)).astype(np.int32)
    mapping[:] = -9999
    bdb = np.copy(obdb)
    for ib in range(bdb.shape[0]):
        flag = 0
        for l in range(bdb.shape[1]-1):
            if (bdb[ib,l] == 0) & (np.sum(bdb[ib,l:])!=0):
                flag = 1
                #Shift all down by one
                while bdb[ib,l] == 0:
                    bdb[ib,l:-1] = bdb[ib,l+1:]
                    bdb[ib,-1] = 0
        #Assemble mapping by starting from the back
        if flag == 1:
            bls = np.arange(bdb.shape[1])[::-1]
            m = bls.size-1
            for l in bls:
                if (obdb[ib,l] != 0) & (np.sum(obdb[ib,:l] == 0) > 0):#(bdb[ib,l] != obdb[ib,l]):
                    #Search for the first one in the update one that is not 0
                    if l < m:m = l
                    bls2 = np.arange(bdb.shape[1])[:m][::-1]
                    for m in bls2:
                        if bdb[ib,m] != 0:
                            if m != l:mapping[ib,l] = m
                            break
    #Correct the original map (and save correction flag)
    band_correction = np.copy(band_position)
    band_correction[:] = -9999
    band_position2 = np.copy(band_position)
    band_position2[:] = -9999
    for i in range(band_position.shape[0]):
        for j in range(band_position.shape[1]):
            if band_position[i,j] == -9999:continue
            basin = basins[i,j]-1
            bp = band_position[i,j]
            if mapping[basin,bp] != -9999:
                band_position2[i,j] = mapping[basin,bp]
                band_correction[i,j] = 1
            else:
                band_position2[i,j] =  band_position[i,j]
            
    return (band_position2,band_correction)

@numba.jit(nopython=True,cache=True)
def determine_band_hand(nhand,band_position,basin_clusters,band_correction):
    
    basin_clusters = basin_clusters.astype(np.int32)
    #1.Construct database of nhand value per band per basin cluster
    ubcs = np.unique(basin_clusters)
    ubcs = ubcs[ubcs != -9999]
    vals = 10**5*np.ones((ubcs.size,np.max(band_position)+1))
    #counts = np.zeros((ubcs.size,np.max(band_position)+1))
    for i in range(nhand.shape[0]):
        for j in range(nhand.shape[1]):
            if ((basin_clusters[i,j] == -9999) | (band_position[i,j] == -9999)):continue
            if (band_correction[i,j] == 1):continue
            bc = basin_clusters[i,j]-1
            bp = band_position[i,j]
            if nhand[i,j] < vals[bc,bp]:
             vals[bc,bp] = nhand[i,j]
    ##2.0 Place the means
    #2.0 Place the mins
    band_nhand = np.copy(nhand)
    band_nhand[:] = -9999.0
    for i in range(nhand.shape[0]):
        for j in range(nhand.shape[1]):
            if ((basin_clusters[i,j] == -9999) | (band_position[i,j] == -9999)):continue
            bc = basin_clusters[i,j]-1
            bp = band_position[i,j]
            band_nhand[i,j] = vals[bc,bp]
            
    return band_nhand

@numba.jit(nopython=True,cache=True)
def determine_band_id(band_position,basin_clusters):
    
    basin_clusters = basin_clusters.astype(np.int32)
    #Assemble a database to assign a unique id to each band position and basin cluster
    ubcs = np.unique(basin_clusters)
    ubcs = ubcs[ubcs != -9999]
    bps = np.zeros((ubcs.size,np.max(band_position)+1))
    for i in range(band_position.shape[0]):
        for j in range(band_position.shape[1]):
            if ((basin_clusters[i,j] == -9999) | (band_position[i,j] == -9999)):continue
            bc = basin_clusters[i,j]-1
            bp = band_position[i,j]
            bps[bc,bp] = 1
    #Determine the id for each site
    bids = np.zeros(bps.shape).astype(np.int32)
    count = 0
    for i in range(bps.shape[0]):
        for j in range(bps.shape[1]):
            if bps[i,j] == 1:
             bids[i,j] = count
             count += 1
    #Place the data
    band_id = np.zeros(basin_clusters.shape).astype(np.int32)
    band_id[:] = -9999
    for i in range(band_position.shape[0]):
        for j in range(band_position.shape[1]):
            if ((basin_clusters[i,j] == -9999) | (band_position[i,j] == -9999)):continue
            bc = basin_clusters[i,j]-1
            bp = band_position[i,j]
            band_id[i,j] = bids[bc,bp]
    return band_id

@numba.jit(nopython=True,cache=True)
def normalize_hand_values(hand,basins):
    #Normalize the hand values in each basin across the clusters
    nhand = np.copy(hand)
    basins = basins.astype(np.int32)
    ubs = np.unique(basins)
    ubs = ubs[ubs != -9999]
    db = np.zeros(ubs.size)
    for i in range(basins.shape[0]):
     for j in range(basins.shape[1]):
      if basins[i,j] == -9999:continue
      b = basins[i,j]-1
      if hand[i,j] > db[b]:db[b] = hand[i,j]
    #1.Normalize each basin
    maxhand = np.copy(hand)
    for i in range(basins.shape[0]):
     for j in range(basins.shape[1]):
      if basins[i,j] == -9999:continue
      b = basins[i,j]-1
      if db[b] != 0:nhand[i,j] = hand[i,j]/db[b]
      maxhand[i,j] = db[b]
    
    return (nhand,maxhand)

#@numba.jit(nopython=True,cache=True)
def cdf_match_hand_values(hand,basins,basin_clusters):

 #CDF match the hand values of all basins that belong to a given cluster
 ubcs = np.unique(basin_clusters).astype(np.int32)
 ubcs = ubcs[ubcs != -9999]
 nhand = np.copy(hand)
 npct = 250
 for ub in ubcs:
    #Determine total number of cells of non-zero hand values
    m = (basin_clusters == ub) & (hand != 0)
    ncells = np.sum(m)
    vals = np.zeros(npct)
    pcts = np.linspace(0,1,npct)
    #Compute ecdf of all non-zero hand values for each basin
    ubs = np.unique(basins[m])
    tmp = 0
    for b in ubs: 
        m1 = (basins == b) & (hand != 0)
        vals1 = hand[m1]
        #Remove the HAND derived channel depth (HACK)
        vals1 = vals1 - np.min(vals1) + 10**-10
        hand[m1] = vals1[:]
        #Construct ecdf for the basin data
        argsort = np.argsort(vals1)
        pcts1 = np.zeros(vals1.size)
        pcts1[argsort] = np.linspace(0,1,pcts1.size)
        #Interpolate to the desired pcts
        fct = np.sum(m1)/np.sum(m)
        vals += fct*np.interp(pcts,pcts1[argsort],vals1[argsort])
    #CDF match the hands of each basin
    for b in ubs:
        m1 = (basins == b) & (hand != 0)
        vals1 = hand[m1]
        #Construct ecdf for the basin data
        argsort = np.argsort(vals1)
        pcts1 = np.zeros(vals1.size)
        pcts1[argsort] = np.linspace(0,1,pcts1.size)
        #CDF match to the basin cluster cdf
        vals2 = np.copy(vals1)
        vals2[argsort] = np.interp(pcts1[argsort],pcts,vals)
        nhand[m1] = vals2[:] 
        '''if np.max(nhand[m1] > 1000):
         print(pcts1[argsort])
         print(vals1)
         print(vals2)
         print(vals)
         exit()'''
    
 return nhand

def create_basin_tiles_updated(basin_clusters,hand,basins,n_binning,cid,max_nbins=100):
    
    #print('before',max(np.unique(hand)))
    nhand = cdf_match_hand_values(hand,basins,basin_clusters)
    #print('after',max(np.unique(nhand)))
    nhand = np.copy(hand)
    
    #Compute the bin edges
    ubcs = np.unique(basin_clusters).astype(np.int32)
    ubcs = ubcs[ubcs != -9999]
    bin_edges = np.zeros((ubcs.size,100))
    bin_edges[:] = -9999
    for ubc in ubcs:
        m = basin_clusters == ubc
        #data = hand[m]
        data = nhand[m]
        #Calculate pct areal coverage of subbasin from channel cells
        a = np.sum(data == 0)/data.size
        #Calculate percentiles
        argsort = np.argsort(data)
        pcts = np.copy(data)
        pcts[argsort] = np.linspace(0,1,data.size)
        #Compute the edges
        x = np.arange(max_nbins)
        tmp1 = a*n_binning**x
        m1 = (tmp1 <= 1.0) & (tmp1 >= a)
        tmp1 = tmp1[m1]
        #Find the closest matches in the pcts
        tmp = np.copy(tmp1)
        tmp[:] = -9999
        for i in range(tmp.size):
         if i == 0:tmp[i] = 0.0
         else:
            argmin = np.argmin(np.abs(tmp1[i] - pcts))
            tmp[i] = data[argmin]
        #Min size of tmp
        if tmp.size < 2:
         tmp = np.zeros(2)
        #Enforce the last value to be the max
        tmp[-1] = np.max(data)
        #Place the data
        bin_edges[ubc-1,:tmp.size+1] = np.concatenate((np.array([0.0,]),tmp))
    #Determine band position
    band_position = determine_band_position(nhand,bin_edges,basin_clusters)
    
    #Correct band position
    (band_position2,band_correction) = curate_band_position(band_position,basins)
    
    #Compute band hand
    band_nhand = determine_band_hand(nhand,band_position2,basin_clusters,band_correction)
    
    #Compute band id
    band_id = determine_band_id(band_position2,basin_clusters)
    
    return (band_id,band_nhand,band_position)

@numba.jit(nopython=True,cache=True)
def transform_arcgis_fdir(fdir_arcgis):
    fdir = np.zeros((fdir_arcgis.shape[0],fdir_arcgis.shape[1],2)).astype(np.int32)
    fdir[:] = -9999
    for i in range(fdir_arcgis.shape[0]):
        for j in range(fdir_arcgis.shape[1]):
            if fdir_arcgis[i,j] == 1:
                fdir[i,j,0] = i
                fdir[i,j,1] = j+1
            if fdir_arcgis[i,j] == 2:
                fdir[i,j,0] = i+1
                fdir[i,j,1] = j+1
            if fdir_arcgis[i,j] == 4:
                fdir[i,j,0] = i+1
                fdir[i,j,1] = j
            if fdir_arcgis[i,j] == 8:
                fdir[i,j,0] = i+1
                fdir[i,j,1] = j-1
            if fdir_arcgis[i,j] == 16:
                fdir[i,j,0] = i
                fdir[i,j,1] = j-1
            if fdir_arcgis[i,j] == 32:
                fdir[i,j,0] = i-1
                fdir[i,j,1] = j-1
            if fdir_arcgis[i,j] == 64:
                fdir[i,j,0] = i-1
                fdir[i,j,1] = j
            if fdir_arcgis[i,j] == 128:
                fdir[i,j,0] = i-1
                fdir[i,j,1] = j+1
    #add one to all indices for fortran indexing
    for i in range(fdir_arcgis.shape[0]):
        for j in range(fdir_arcgis.shape[1]):
            if fdir[i,j,0] != -9999:
                fdir[i,j,0] = fdir[i,j,0] + 1
                fdir[i,j,1] = fdir[i,j,1] + 1

    return fdir
