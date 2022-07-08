import numpy as np
#from skimage.segmentation import find_boundaries, clear_border
import random

def hrus_centroid_distance(lats,lons):
    radius = 6371 # km

    ncells = len(lats)
    distance = np.zeros((ncells,ncells))
    for i, lat1, lon1 in zip(range(ncells),lats,lons):
     for j, lat2, lon2 in zip(range(ncells),lats,lons):
      if j>i:
       dlat = np.radians(lat2-lat1)
       dlon = np.radians(lon2-lon1)
       a = np.sin(dlat/2) * np.sin(dlat/2) + np.cos(np.radians(lat1)) \
         * np.cos(np.radians(lat2)) * np.sin(dlon/2) * np.sin(dlon/2)
       c = 2 * np.math.atan2(np.sqrt(a), np.sqrt(1-a))
       distance[i,j] = radius * c
       distance[j,i] = radius * c

    return distance

def hrus_slope(elevation,dist):
    ncells = len(elevation)

    slope = np.zeros((ncells,ncells))
    distance = np.array(dist)*1000.0 # convert to m

    for i, h1 in zip(range(ncells),elevation):
     for j, h2 in zip(range(ncells),elevation):
      if j>i : 
       if distance[i,j] > 0:
        slope[i,j] = (h1-h2)/distance[i,j]
        slope[j,i] = -slope[i,j]
       else:
        slope[i,j] = 0.0
        slope[j,i] = 0.0

    return slope

def calc_calendar(self,ncells):
    grow = np.zeros((ncells,12),dtype=int)
    st_gscal = np.copy(self.st_gscal)-1
    en_gscal = np.copy(self.en_gscal)-1    

    leng = en_gscal-st_gscal+1
   
    m = leng > 0
    for i in np.where(m)[0]:
      grow[i,st_gscal[i]:en_gscal[i]+1] = 1

    m = leng < 0
    for i in np.where(m)[0]:
      grow[i,:] = 1
      grow[i,en_gscal[i]+1:st_gscal[i]] = 0
 
  
    return grow

def calculate_min_distance(hsu,nhru,cluster_ids,lats,lons,clats,clons):
  radius = 6367.0
  # Minimum distance between HRUs

  # Get lat lon from the borders
  idx = (cluster_ids == hsu)
  idx = clear_border(idx,bgval=False)
  idx = find_boundaries(idx, mode='inner')
  bd_lats = lats[idx].flatten()
  bd_lons = lons[idx].flatten()

  if len(bd_lats) < 1 :
   idx = (cluster_ids == hsu)
   idx = find_boundaries(idx, mode='inner')
   bd_lats = lats[idx].flatten()
   bd_lons = lons[idx].flatten()

  # Get unique lat,lon values and sample 50 points
  points = set(zip(bd_lats,bd_lons))
  nsamp = 1#30
  if len(points) <= nsamp: nsamp = int(len(points)/2.)
  if len(points) <= 5: nsamp = len(points)

  points = random.sample(points, nsamp)
  bd_lats = np.array(list(zip(*points))[0])
  bd_lons = np.array(list(zip(*points))[1])

  distance = np.ones(nhru)*10000000.

  #Calculate the distance of a boundary to a centroid of each hru

  for lat, lon in zip(bd_lats,bd_lons):
    dlat = np.radians(lat-clats)
    dlon = np.radians(lon-clons)
    a = np.sin(dlat/2) * np.sin(dlat/2) + np.cos(np.radians(clats)) \
      * np.cos(np.radians(lat)) * np.sin(dlon/2) * np.sin(dlon/2)
    c = np.zeros((len(a)))
    for count in range(len(a)):
      c[count] = 2 * np.math.atan2(np.sqrt(a[count]), np.sqrt(1-a[count]))
    dist = radius * c
    distance[dist < distance] = dist[dist < distance]
#  for hrs in range(nhru):
#    if hrs == hsu:
#      distance[hrs] = 0.0
#    else:
#      clat = clats[hrs]
#      clon = clons[hrs]
#
#      for lat, lon in zip(bd_lats,bd_lons):
#        dlat = np.radians(lat-clat)
#        dlon = np.radians(lon-clon)
#        a = np.sin(dlat/2) * np.sin(dlat/2) + np.cos(np.radians(clat)) \
#          * np.cos(np.radians(lat)) * np.sin(dlon/2) * np.sin(dlon/2)
#        c = 2 * np.math.atan2(np.sqrt(a), np.sqrt(1-a))
#        dist = radius * c
#        if dist < distance[hrs]: distance[hrs] = dist
#      #print hsu, hrs, dist, distance[hrs] 

  #print hsu, distance
  return distance
#


