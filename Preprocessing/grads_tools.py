import grads
import numpy as np
import datetime
#ga = grads.GrADS(Bin='grads',Window=False,Echo=False)

def open_grads(bin):

 ga = grads.GrADS(Bin=bin,Window=False,Echo=False)
 #ga = grads.GrADS(Bin=bin,Window=False,Echo=True,Verb=2,Strict=0)#,Echo=False)

 return ga

def extract_point_data(file,lats,lons,var,type):

 #Open file
 ga("%s %s" % (type,file))

 #Extract data
 values = []
 for i in xrange(len(lons)):
  ga("set lat %f" % lats[i])
  ga("set lon %f" % lons[i])
  values.append(np.array(ga.expr(var)))

 #Close file
 ga("close 1")

 return np.array(values)

def datetime2gradstime(date):

 #Convert datetime to grads time
 str = date.strftime('%HZ%d%b%Y')

 return str

def gradstime2datetime(str):

 #Convert grads time to datetime
 date = datetime.datetime.strptime(str,'%HZ%d%b%Y')

 return date

def retrieve_metadata(ga):
 
 metadata = {}

 #Obtain the metadata
 metadata['nlat'] = ga.query('dims').ny
 metadata['nlon'] = ga.query('dims').nx
 metadata['minlat'] = ga.query('dims').lat[0]
 metadata['minlon'] = ga.query('dims').lon[0]
 metadata['maxlat'] = ga.query('dims').lat[1]
 metadata['maxlon'] = ga.query('dims').lon[1]
 vars = ga.query('file').vars
 tmp = ga.exp(vars[0])
 metadata['res'] = tmp.grid.lat[1] - tmp.grid.lat[0]
 metadata['undef'] = np.min(np.ma.getdata(tmp))

 return metadata
