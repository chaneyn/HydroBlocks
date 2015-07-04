import gdal
import osgeo
import os
from osgeo import osr
import numpy as np
import osgeo
import os

def extract_point_data(file,lats,lons):
 
 #Open file and get geotransformation
 ds = gdal.Open(file)
 gt = ds.GetGeoTransform()
 rb = ds.GetRasterBand(1)

 #Compute ilats and ilons
 ilons = np.round((np.array(lons) - gt[0])/gt[1]).astype(np.int)
 ilats = np.round((np.array(lats) - gt[3])/gt[5]).astype(np.int)

 #Extract data
 values = []
 for i in xrange(ilons.size):
  values.append(rb.ReadAsArray(ilons[i],ilats[i],1,1)[0])

 return np.array(values).astype(np.float)

def read_raster(file):

 #Read in the raster
 dataset = gdal.Open(file)

 #Get dimensons
 nx = dataset.RasterXSize
 ny = dataset.RasterYSize

 #Retrieve band
 band = dataset.GetRasterBand(1)

 #Convert to numpy array
 data = band.ReadAsArray(0,0,nx,ny).astype(np.float32)

 return data

def read_raster_subarea(file,metadata):

 #Read in the raster
 dataset = gdal.Open(file)

 #Get dimensons
 nx = metadata['nx']#dataset.RasterXSize
 ny = metadata['ny']#dataset.RasterYSize
 ixmin = metadata['ixmin']
 iymin = metadata['iymin']

 #Retrieve band
 band = dataset.GetRasterBand(1)

 #Convert to numpy array
 data = band.ReadAsArray(ixmin,iymin,nx,ny).astype(np.float32)

 return data

'''def write_raster(metadata,data,file):

 cols = metadata['nlon']
 rows = metadata['nlat']
 minlon = metadata['minlon']
 if minlon > 180: minlon = minlon - 360
 bands = 1
 driver = gdal.GetDriverByName('GTiff')
 #Create file
 ds = driver.Create(file,cols,rows,1,gdal.GDT_Float32)
 #Set geo information
 ds.SetGeoTransform([minlon,metadata['res'],0,metadata['maxlat'],0,-metadata['res']])
 proj = osr.SpatialReference()
 proj.SetWellKnownGeogCS("EPSG:4326")
 ds.SetProjection(proj.ExportToWkt())
 outband = ds.GetRasterBand(1)
 outband.WriteArray(np.flipud(data),0,0)
 ds = None

 return'''

def shapefile2raster(raster_in,shp_in,raster_out,workspace,field,layer):

 #Extract coordinates and projection info from the target file
 ds = gdal.Open(raster_in)
 gt = ds.GetGeoTransform()
 cols = ds.RasterXSize
 rows = ds.RasterYSize
 srs = osgeo.osr.SpatialReference()
 srs.ImportFromWkt(ds.GetProjection())
 proj4 = srs.ExportToProj4()

 #Rasterize the shapefile
 minx = gt[0]
 miny = gt[3]+rows*gt[5]
 maxx = gt[0]+cols*gt[1]
 maxy = gt[3]
 shp_out = '%s/%d' % (workspace,np.random.randint(10**5))
 #os.system('ogr2ogr -spat %.16f %.16f %.16f %.16f -overwrite -select CELLVALUE,MUKEY %s %s' % (dims['minlon'],dims['minlat'],dims['maxlon'],dims['maxlat'],shp_out,shp_in))
 os.system("ogr2ogr -spat %.16f %.16f %.16f %.16f -overwrite -t_srs '%s' %s %s" % (minx,miny,maxx,maxy,proj4,shp_out,shp_in))
 os.system('gdal_rasterize -ot Float32 -l %s -init -9999 -a %s -te %.16f %.16f %.16f %.16f -tr %.16f %.16f %s %s' % (layer,field,minx,miny,maxx,maxy,gt[1],gt[5],shp_out,raster_out))
 os.system('rm -r %s' % shp_out)

 return

def raster2raster(raster_template,raster_in,raster_out):

 #Extract coordinates and projection info from the target file
 ds = gdal.Open(raster_template)

def write_raster(file,metadata,data):

 driver = gdal.GetDriverByName('GTiff')
 ds_out = driver.Create(file,metadata['nx'],metadata['ny'],1,gdal.GDT_Float32)
 ds_out.GetRasterBand(1).WriteArray(data.astype(np.float32))
 ds_out.SetGeoTransform(metadata['gt'])
 ds_out.SetProjection(metadata['projection'])
 ds_out.GetRasterBand(1).SetNoDataValue(metadata['nodata'])
 ds_out = None

 return

def retrieve_metadata(raster):

 metadata = {}
 #Extract coordinates and projection
 ds = gdal.Open(raster)
 gt = ds.GetGeoTransform()
 cols = ds.RasterXSize
 rows = ds.RasterYSize
 srs = osgeo.osr.SpatialReference()
 srs.ImportFromWkt(ds.GetProjection())
 metadata['proj4'] = srs.ExportToProj4()
 metadata['minx'] = gt[0]
 metadata['miny'] = gt[3]+rows*gt[5]
 metadata['maxx'] = gt[0]+cols*gt[1]
 metadata['maxy'] = gt[3]
 metadata['resx'] = gt[1]
 metadata['resy'] = gt[5]
 metadata['gt'] = gt
 metadata['nx'] = cols
 metadata['ny'] = rows
 metadata['projection'] = ds.GetProjection()

 return metadata
