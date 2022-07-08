import numpy as np
#import terrain_tools_fortran as ttf
from geospatialtools import terrain_tools_fortran as ttf
from geospatialtools import metrics 
import sklearn.cluster
import sklearn.mixture
import sklearn.linear_model
import scipy.stats
import copy
import time
import pickle

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
 if (max != min):
  data[m] = (data[m] - min)/(max-min)
 else:
  data[m] = 0.0

 data[data > 1.0] = 1.0
 data[data < 0] = 0.0

 return data

def cluster_data(X,nc):

 #Assemble sample list
 minsamples = int((10**5)/2)  # Taking very long per tile for large catchments
 if X.shape[0] > minsamples:
  np.random.seed(1245)
  idx = np.random.choice(np.arange(X.shape[0]),minsamples)
 else:
  idx = np.arange(X.shape[0])
 
 #The number of clusters must be equal or smaller than the number of samples
 if idx.size <= nc: nc = int(idx.size/10)  # Noemi
 
 #Cluster the data
 if nc > 1:
  nfeatures = X.shape[1]
  #model = sklearn.cluster.MiniBatchKMeans(n_clusters=nc,random_state=35799)
  model = sklearn.cluster.KMeans(n_clusters=nc,random_state=35799)
  p = model.fit_predict(X)
  
 else:
  p = np.zeros(X.shape[0])

 return p



def alternative_cluster_data(X,nc):

 #Assemble sample list
 minsamples = int((10**5)/2)  # Taking very long per tile for large catchments
 if X.shape[0] > minsamples:
  np.random.seed(1245)
  idx = np.random.choice(np.arange(X.shape[0]),minsamples)
 else:
  idx = np.arange(X.shape[0])

 #The number of clusters must be equal or smaller than the number of samples
 if idx.size <= nc: nc = int(idx.size/10)  # Noemi

 #Cluster the data
 if nc > 1:
  nfeatures = X.shape[1]
  #model = sklearn.mixture.GaussianMixture(n_components=nc, covariance_type='tied', random_state=35799)
  #model = sklearn.mixture.BayesianGaussianMixture(n_components=nc, covariance_type='full', random_state=35799)
  model = sklearn.cluster.SpectralClustering(n_clusters=nc, assign_labels="discretize", random_state= 35799)
  p = model.fit_predict(X)

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
   tmp = ws[Xd.keys().index(var)]*tmp
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
 (ah,lath,lonh,hid,nid) = ttf.calculate_basin_properties(basins,res,nb,fdir, latitude,longitude)
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

 ##Construct the lookup 
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

def create_basin_tiles(basin_clusters,hand,input_dem,basins,dh):

 new_hand = np.copy(hand)
 #Iterate per cluster
 ubcs = np.unique(basin_clusters)
 ubcs = ubcs[ubcs!=-9999]
 tiles = np.copy(hand).astype(np.int32)
 tiles[:] = -9999
 tiles_position = np.copy(hand).astype(np.int32)
 tiles_position[:] = -9999
 count = 0

 '''
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
 ''' 
 
 '''           
 #compute number of bins
 data_max = np.max(hand)
 pedges = 1.5
 #if data_max < 2: exit('check terrain_tools.py  sub-basins is too flat for high bands')
 if dh <= 10  and data_max > 30:   dh = 20
 if dh <= 20  and data_max > 100:  dh = 30
 if dh <= 30  and data_max > 200:  dh = 50; pedges = 2.0
 if dh <= 50  and data_max > 300:  dh = 60
 if dh <= 60  and data_max > 400:  dh = 70
 if dh <= 70  and data_max > 600:  dh = 80; pedges = 2.5
 if dh <= 80  and data_max > 800:  dh = 90
 if dh <= 90  and data_max > 1000:  dh = 100

 nbins = int(np.ceil(data_max/dh))
 bin_edges = np.linspace( 3**(1.0/float(pedges)), data_max**(1.0/float(pedges)), nbins+1)**pedges
 bin_edges = np.concatenate(([0.0,0.1],bin_edges))
 if data_max < 3:  bin_edges = np.array([0.0,0.1,data_max])
 print(['{:.2f}'.format(i) for i in bin_edges],flush=True)
 '''

 for ubc in ubcs:                
  m = basin_clusters == ubc      
  data = hand[m]

  #curate        -- can be removed since spatial imputation is now done at Preprocessing.
  #data[data == -9999] = np.max(data[data != -9999])                
  #hand[m & (hand == -9999)] = np.max(hand[m & (hand != -9999)])    

  
  dh = 50#10
  #compute number of bins
  data_max = np.max(data)
  #if data_max < 2: exit('check terrain_tools.py  sub-basins is too flat for high bands')
  #if dh <= 10  and data_max > 30:   dh = 20
  #if dh <= 30  and data_max > 100:  dh = 40
  if dh <= 40  and data_max > 200:  dh = 60
  if dh <= 50  and data_max > 300:  dh = 70
  if dh <= 60  and data_max > 400:  dh = 80
  if dh <= 70  and data_max > 500:  dh = 100
  if dh <= 100 and data_max > 750:  dh = 150
  if dh <= 100 and data_max > 1000: dh = 200
  #dh = 50  
 
  nbins = int(np.ceil(data_max/dh))
  if nbins < 1 : nbins=1 # Noemi

  #Compute the edges for the high bands
  pedges = 1.5  #2.5 # 3.0 #2.5   # Noemi
  bin_edges = np.linspace( 7**(1.0/float(pedges)), data_max**(1.0/float(pedges)), nbins+1)**pedges
  bin_edges = np.concatenate(([0.0,0.1,3],bin_edges))

  if data_max < 10:  bin_edges = np.array([0.0,0.1,3,data_max])
  if data_max < 7 :  bin_edges = np.array([0.0,0.1,data_max])
  if data_max < 0.1: bin_edges = np.array([0.0,0.1]) 
  if data_max > 30:
     it = 0
     while bin_edges[4] > 20 and it < 100:
       pedges = pedges+0.1
       bin_edges = np.linspace( 7**(1.0/float(pedges)), data_max**(1.0/float(pedges)), nbins+1)**pedges
       bin_edges = np.concatenate(([0.0,0.1,3],bin_edges))
       it = it + 1
       if pedges > 10.0 : 
         it = 0
         pedges = 1.5
         nbins = nbins+1

  #print('{:.2f}'.format(pedges),['{:.2f}'.format(i) for i in bin_edges],flush=True)
  #exit() 
 
  #compute the binning           
  #(hist,bin_edges) = np.histogram(data,bins='fd')#bins=nbins)     
  #update edges 
  #bin_edges[0] = 0.0            
  #bin_edges[-1] = np.max(data)  
  #Assign the tiles              

  # Combute edges for the dem bins
  dh_dem = 200
  dem = input_dem[m]
  delta_dem = np.max(dem)-np.min(dem)
  dem_nbins = int(np.around(delta_dem/dh_dem,0))
  if dem_nbins < 1: dem_nbins=1
  dem_edges = np.linspace(np.min(dem),np.max(dem),dem_nbins+1)
  #dem_edges = np.linspace(np.min(dem), np.max(dem), 2)

  count2 = 0
  for i in range(bin_edges.size-1):               
   if i == 0:
     m2 = m & (hand >= bin_edges[i]) & (hand <= bin_edges[i+1])      
   else:
     m2 = m & (hand >  bin_edges[i]) & (hand <= bin_edges[i+1]) 

   if np.sum(m2) > 0:

     tiles_position[m2] = count2
     if i == 0 : new_hand[m2] = 0.0
     else: new_hand[m2] = np.mean(hand[m2])
     #new_hand[m2] = np.mean(hand[m2])
     count2 += 1

     for j in range(dem_edges.size-1):
       if j == 0: m3 = m2 & (input_dem >= dem_edges[j]) & (input_dem <= dem_edges[j+1])
       else:      m3 = m2 & (input_dem > dem_edges[j]) & (input_dem <= dem_edges[j+1])
       #m3 = m2

       if np.sum(m3) > 0:
         tiles[m3] = count
         count += 1
         #new_hand[m3] = np.mean(hand[m3])

     #if np.sum(m2) > 0:            
     # tiles[m2] = count            
     # tiles_position[m2] = count2  
     # new_hand[m2] = np.mean(hand[m2])              
     # count += 1  
     # count2 += 1       

 return (tiles,new_hand,tiles_position)

def create_hrus_hydroblocks(hillslopes,htiles,channels,covariates,nclusters):

 #Compute optimal number of clusters flag
 flag = False

 #Curate the covariates
 for var in covariates:
  val = np.mean(covariates[var]['d'][covariates[var]['d'] != -9999])
  covariates[var]['d'][covariates[var]['d'] == -9999] = val
 import sklearn.cluster
 hrus = np.copy(hillslopes)
 hrus[:] = -9999
 #Iterate through each gru and tile and compute hrus
 uhs = np.unique(hillslopes)
 uhs = uhs[uhs != -9999]
  
 # define average size for hrus # Noemi
 catch_ngrids = np.sum(hillslopes >= 0)
 grid_area_km2 = (100*0.00027777777777 * 100*0.00027777777777) #1 km2
 grids_per_km2 = 1/grid_area_km2
 
 # 36 grids -> 1km
 # 1296 -> 1km2
 # 648  grids ~ 0.5 km2/hru
 
 scale_hru = 0.5*grids_per_km2  # number of grids per HRU. Set to 0.5 km2 per HRU
 if catch_ngrids/scale_hru < 10:  scale_hru = catch_ngrids/10
 if catch_ngrids/scale_hru > 400: scale_hru = catch_ngrids/400.0 
 
 # define binary channels
 is_channel = np.zeros(channels.shape)
 is_channel[channels.data > 0] = 1
 
 maxc = 1
 for uh in uhs:
  mh = hillslopes == uh  # or subbasins
  uts = np.unique(htiles[mh]) 
 
  # Normalize DEM at the subbasin level -- Noemi
  #tmp = covariates['dem']['d'][mh]
  #covariates['dem']['min'] = np.min(tmp[tmp!=-9999])
  #covariates['dem']['max'] = np.max(tmp[tmp!=-9999])  

  for ut in uts:
   mt = mh & (htiles == ut)

   is_channel_hand = np.max(is_channel[mt]) == 1.0

   #Compute the parameters for the clustering algorithm
   ccp = {}
   for var in covariates:
    tmp = covariates[var]['d'][mt] 
    ccp[var] = {'d':tmp,'t':covariates[var]['t'],
                'min':covariates[var]['min'],
                'max':covariates[var]['max']}
    if var in ['lats','lons','basins']: # Normalize lat lon at the tile level -- Noemi
      ccp[var] = {'d':tmp,'t':covariates[var]['t'],            
                 'min':np.min(tmp[tmp!=-9999]),     
                 'max':np.max(tmp[tmp!=-9999])}
    # if this hand is the hand with the channels of this subbasins, only cluster on lat, lon and dem -- Noemi
    if is_channel_hand:
      #if var in ['sand','clay','silt','aspect']:
      if var not in ['lat','lon','dem','basins']:
        ccp[var]['d'][:] = 0
        ccp[var]['min'] = 0
        ccp[var]['max'] = 0

      # Normalize DEM if in the channel level
      #if var in ['dem']:
      #  tmp = covariates['dem']['d'][mt]
      #  ccp['dem']['min'] = np.min(tmp[tmp!=-9999])
      #  ccp['dem']['max'] = np.max(tmp[tmp!=-9999])


   if flag == True:
    (nc,ws) = compute_cluster_parameters(ccp,maxnc)
   else:
    #nc = nclusters
    size_tile = np.sum((mt == 1)) 
    nc = int(np.ceil(size_tile/scale_hru)) 

    if nc < 1: nc = 1 # Noemi

    # Check if this high band is the 'bottom' one 
    #if is_channel_hand :  #if nc < 1: nc = 1  # Noemi -- add more HRUS at the river channel
    #  nc = int(nc*2)
    #  #print(uh,uts,min_hand,nc,flush=True)

    if nc > (size_tile/10): nc = int(np.ceil(size_tile/10)) # at least about 10 pixels in each cluster, otherwise small clusters disappear in hru_latlon -- maybe will be solved when dropping ea_tif inputs.
    ws = np.ones(len(ccp.keys()))
   #print 'hillslope: %d, tile: %d, nc: %d' % (uh,ut,nc)
   #Add weights to covariates
   for var in covariates:
    covariates[var]['w'] = ws[list(ccp).index(var)]
   #prepare the covariate data
   X = []
   for var in covariates:
    #Normalize and apply weight
    tmp = covariates[var]['w']*normalize_variable(covariates[var]['d'][mt],covariates[var]['min'],covariates[var]['max'])
    mask_out = ( (np.isnan(tmp) | np.isinf(tmp)) | (tmp == -9999) )  # Noemi
    #print(tmp.shape,flush=True)
    #tmp = tmp[~mask_out]
    #print(tmp.shape,flush=True)
    if sum(mask_out) > 0 : exit('%s has -9999 or invalid entries on the clustering' % var)
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

 #Cleanup hrus
 ttf.cleanup_hillslopes(hrus)
 hrus[hrus >= 0] = hrus[hrus >= 0] + 1


 return hrus
 

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

 varea = np.sum(basins >= 0) # grids #Noemi
 # 33 grids -> 1km
 # 1089 -> 1 km2 
 varea_km2 = varea/1089.
 print('basin area km2: %.2f' % varea_km2, flush=True)
 scale_hru = 5.                      # 10 km2 per subbasin         
 nclusters = int(np.round(varea_km2/scale_hru,0))       
 #print('subbasins: %.2f' % nclusters, flush=True)
 max_nclusters = 50  # Set a maximum number of sub-basins  -- Noemi
 if nclusters == 0: nclusters=1
 if nclusters > max_nclusters: nclusters=max_nclusters
 
 '''
 # use more subbasins if there is a steep topography
 if nclusters < max_nclusters:
  dhmax = covariates['dem']['max']-covariates['dem']['min']
  step_ratio = (dhmax/2.)/(varea_km2) # m/km2
  if step_ratio > 1: nclusters = int(nclusters*step_ratio)
  if nclusters > max_nclusters: nclusters=max_nclusters
 '''

 # number of existent subbasins
 nbasins = np.unique(basins)
 if -9999 in nbasins: nbasins = len(nbasins[1:] )
 if max_nclusters > nbasins: max_nclusters = nbasins

 print('subbasins: %i' % nclusters, flush=True)
 X = []
 for var in covariates:
  otmp = np.copy(covariates[var]['d'])
  otmp[(np.isnan(otmp) == 1) | (np.isinf(otmp) == 1)] = -9999.0
  tmp = np.copy(otmp) 
  #Normalize and apply weight
  tmp = normalize_variable(tmp,covariates[var]['min'],covariates[var]['max'])
  X.append(tmp)

 X = np.array(X).T
 
 # discontinuous clusters
 clusters = cluster_data(X,nclusters)+1

 # continuous clusters
 #clusters = alternative_cluster_data(X,nclusters)+1

 #Create the mapping
 mapping = np.zeros(np.max(hp_in['bid'])+1)
 mapping[:] = -9999
 for i in range(hp_in['bid'].size):
  mapping[hp_in['bid'][i]] = clusters[i]
 
 #Assign the new ids to each hillslpe
 basins_clusters = ttf.assign_clusters_to_hillslopes(basins,mapping)

 #Determine the number of basins per cluster
 #uclusters = np.unique(clusters)

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

