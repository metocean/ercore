#!/usr/bin/env python
import os
import numpy
from ercore import ERCoreException,ERConfigException,copydoc,ncep2dt,dt2ncep
from _flib_ercore import interph,interp3d,interpz,inpoly
import datetime
import netCDF4 as nc
import glob
import re
import shelve

R2D=180./numpy.pi
D2R=1/R2D
ARAD=R2D/6367456.
DMIN=1.e-10
ALMOST_ZERO = 1.e-6


def slope_correction(p,topo,uu):
    '''Vertical velocity correction for slope of seabed'''
    frac = numpy.ones(topo.shape[0])
    ind = (topo[:,0]!=0)
    frac[ind] = numpy.minimum(abs(p[ind,2]/topo[ind,0]),1)
    return frac*(uu[:,0]*topo[:,1]+uu[:,1]*topo[:,2])


class DataException(ERCoreException):
  pass
 
class RectInterpolator(object):
  """Interpolation class for regular grids"""
  def __init__(self,x0,y0,x1,y1,dx,dy,lev=None,lat=None):
    """Constructor for interpolator
    Arguments:
      x0: Bottom left x-coordinate
      y0: Bottom left y-coordinate
      x1: Top right x-coordinate
      y1: Top right y-coordinate
      dx: Grid spacing x
      dy: Grid spacing y
      lev: Levels for 3D grid
      lat: Latitude values for geodetic grid
    """
    self.x0=x0
    self.y0=y0
    self.x1=x1
    self.y1=y1
    self.idx=1./dx
    self.idy=1./dy
    self.lev=lev
    if lat:
      self.geod=True
      self.lat=lat[:]
    
  def __call__(self,dat,p): 
    if self.lev is not None:
      return interp3d(dat,p[:,0],p[:,1],p[:,2],self.x0,self.y0,self.idx,self.idy,self.lev)
    else:
      return interph(dat,p[:,0],p[:,1],self.x0,self.y0,self.idx,self.idy)
    
  def ingrid(self,p):
    inbbox=(p[:,0]>=self.x0) & (p[:,0]<=self.x1) & (p[:,1]>=self.y0) & (p[:,1]<=self.y1)
    return inbbox
    
  def grad(self,dat):
    # compute gradient - used when called in GriddedTopo to work out vertical velocity term
    # allowing nice seabed-following particles motions near the seabed
    # 
    dhdy,dhdx=numpy.gradient(dat)
    dhdy*=self.idy
    dhdx*=self.idx
    if self.geod:
      mfy,mfx=numpy.meshgrid(numpy.tile(ARAD,dhdy.shape[1]),ARAD/numpy.cos(D2R*self.lat))
      dhdy*=mfy
      dhdx*=mfx
    return dhdy,dhdx
      
class FEInterpolator(object):
  """Interpolation class for Finite element grids"""
  def __init__(self,lon,lat,tri,lev=None,geod=True):
    """Constructor for interpolator
    Arguments:
      lon: Array of longitudes for each node
      lat: Array of latitudes for each node
      geod: Geodetic coordinates (True/False)
      lev: Levels for 3D grid
    """
    from scipy.spatial import cKDTree #quick nearest-neighbor lookup
    from scipy.spatial import ConvexHull # convex hull of grid points
    from matplotlib.path import Path # convex hull of grid points
    self.lon = lon
    self.lat = lat
    self.x0 = min(lon)
    self.x1 = max(lon)
    self.y0 = min(lat)
    self.y1 = max(lat)
    self.geod = geod
    self.lev = lev
    #tree for quick nearest-neighbor lookup
    self.tree = cKDTree(numpy.vstack((self.lon,self.lat)).T)
    # build convex hull of points for particle-in-mesh checks
    hull_obj = ConvexHull(numpy.vstack([self.lon,self.lat]).T)
    self.hull_path = Path(numpy.vstack([self.lon[hull_obj.vertices],self.lat[hull_obj.vertices]]).T)  # Matplotlib Path object
    self.hull = numpy.vstack([self.lon[hull_obj.vertices],self.lat[hull_obj.vertices]]).T
    # save triangles - these will be needed to compute the depth gradient
    self.tri = tri

  def __call__(self,dat,p):
    dist,i=self.tree.query(p[:,:2],3, n_jobs=-1) #quick nearest-neighbor lookup
    
    if i.max()>=dat.shape[-1]:
      raise DataException('Finite element interpolation out of range')
    dist[dist<DMIN]=DMIN
    fac=(1./dist)
    
    if dat.ndim==2:
      tmp=(fac*dat.take(i,1)).sum(-1)/fac.sum(-1)
      datout = interpz(tmp.astype('f'),p[:,2],self.lev)
    else:
      datout = (fac*dat.take(i)).sum(-1)/fac.sum(-1)
    # mask out of grid points - setting velocities to 0.0
    datout[numpy.where(~self.ingrid(p) )] = 0.0 


    return datout
    
  def ingrid(self,p):
    # this checks if particle positions are within the mesh extents
    # Note this is using the convex hull of points instead of the true mesh boundaries
    # *** we would need infos in the netcdf itself to have bnd infos or a more evolved trimesh package ?
    #
    # inmesh=inpoly(p[:,0],p[:,1],self.hull[:,0],self.hull[:,1]) # this doesnt work...not sure why
    # xp = [-1.0,1.0,1.0,-1.0,-1.0]
    # yp = [-1.0,-1.0,1.0,1.0,-1.0]
    # inpoly(0,0,xp,yp) yields zero ???

    in_hull =self.hull_path.contains_points(p[:,:2]) # use the Path object (scipy) defined in def init
    
    return in_hull
    
  def grad(self,dat):
    # gradient computation of a scalar field on a triangular mesh 
    # 
    # based on example here, uisng matplotlib :
    # https://matplotlib.org/examples/pylab_examples/trigradient_demo.html
    # http://www.edwilliams.org/avform.htm#Dist
    from matplotlib.tri import Triangulation,CubicTriInterpolator
    import matplotlib.cm as cm
    import math
    # need to access the scalar data e.g. depth
    # triangle connectivity
    # then apply map factors
    # then implement proper use in code to work out 
    # a vertical velocity allowing nice seabed-following motions near the seabed
    if False:
      print 'Computing bathy gradients from triangular mesh'
      # convert lon/lat to true x,y
      if False: # some test using different appraoch for the computation of map factors
        lat_rad=numpy.pi*self.lat/180 # latitude in radians
        deg_in_m_lat=111132.92-559.82 * numpy.cos(2* lat_rad)+1.175*numpy.cos(4*lat_rad)
        deg_in_m_lon=111412.84 * numpy.cos(lat_rad)- 93.5 * numpy.cos(3*lat_rad)
        #https://knowledge.safe.com/articles/725/calculating-accurate-length-in-meters-for-lat-long.html
        # see also http://www.edwilliams.org/avform.htm
        x = self.lon * deg_in_m_lon
        y = self.lat * deg_in_m_lat
        # convert grid lon/lat to cartesian coordinates for the gradient computation - see geodcalc in /materials/__init__.py
        # mfx,mfy allows converting meters to degrees
        mfy = numpy.tile(ARAD,self.lat.shape[0]) #self.mfx[:,1]=numpy.tile(ARAD,self.npmax+1)
        mfx = mfy / numpy.cos(D2R*self.lat)#  self.mfx[:,0]=self.mfx[:,1]/numpy.cos(D2R*self.pos[:,1])
        x1 = self.lon / mfx
        y1 = self.lat / mfy

      # convert grid lon/lat to cartesian coordinates for the gradient computation - see geodcalc in /materials/__init__.py
      # mfx,mfy allows converting meters to degrees
      mfy = numpy.tile(ARAD,self.lat.shape[0]) #self.mfx[:,1]=numpy.tile(ARAD,self.npmax+1)
      mfx = mfy / numpy.cos(D2R*self.lat)#  self.mfx[:,0]=self.mfx[:,1]/numpy.cos(D2R*self.pos[:,1])
      # get triangulation connectivty
      self.tri_id = self.tri - 1 # switch to index-based traingle connectivty , index=0 refers to  node number 1       
      #create the triangulation 
      triang = Triangulation(self.lon,self.lat,self.tri_id) # use triangle data connectivity from netcdf file
      # # create the interpolator
      data = CubicTriInterpolator(triang, dat)
      # # compute depth gradients at mesh nodes, but could be anywhere else
      # NOTE : here the <depth> field is already negative down; this yields positive depth gradient when the
      # seabed gets shallower, which is what we want
      # (e.g. particle goes from -5m to -3m : dz = -3 -(-5) = +2 >> dz/dx >0      
      (dhdx,dhdy) = data.gradient(self.lon,self.lat)
      # NOTE : at this stage, gradients dh/dx, dh/dy are [meter/degrees]
      # so these should be used on particle displacement in [degrees] (rather than [meters]).

      # switch to absolute gradient - consistent with approach for regular grid...
      # still not sure this is right...
      dhdx *= mfx
      dhdy *= mfy
      import pdb;pdb.set_trace()
      return dhdx,dhdy
    else:
      return 0.*dat,0.*dat ##!!!!Not correct - just for testing - but this needs to be implemented


class FieldData:
  """Base class for fields"""
  is3d=False
  res=1.0
  def __init__(self,id,vars,**options):
    """Constructor for field data
  Arguments:
    id: Unique id
    vars: List of string ids for variables
    options: optional extra keyword arguments
  """
    self.id=id
    self.vars=vars if isinstance(vars,list) else [vars]
    self.__dict__.update(options)
    
  def interp(self,p,time,age):
    """Return field values for each particle position at given time
  Arguments:
    p: List of particle positions
    time: Specified time as NCEP decimal (optional)
    age: Age of particles (optional)
    imax: Dimensions of field (i.e. 2 or 3)
  """
    pass
    
class ConstantData(FieldData):
  """Class to specify a constant field"""
  @copydoc(FieldData.__init__)
  def __init__(self,id,vars,**options):
    """  options:
    constant values for each of the variables"""
    FieldData.__init__(self,id,vars,**options)
    if self.is3d:
      if not options.has_key('levels') and not options.has_key('depth'):
        raise ConfigException('Constant 3D data must specify levels')
      if options.has_key('levels'):
        self.lev=numpy.array(map(float,options['levels']))
      elif options.has_key('depth'):
        self.lev=-1*numpy.array(map(float,options['depth']))
    self.dat=[]
    for v in self.vars:
      if not options.has_key(v):
        raise ERConfigException('Constant data for variable %s missing' % (v))
      vin=options[v] if hasattr(options[v], '__len__') else [options[v]]
      self.dat.append(numpy.array(map(float,vin)))
    if self.is3d and (len(self.lev)>1) and (self.lev[1]-self.lev[0])<0:
      self.lev=self.lev[::-1]
      for id,d in enumerate(self.dat):
        self.dat[id]=d[::-1]
      
      
  @copydoc(FieldData.interp)
  def interp(self,p,time=None,age=None,imax=3):
    datout=numpy.zeros((len(p),imax))
    for id,d in enumerate(self.dat[:imax]):
      if self.is3d:
        datout[:,id]=numpy.interp(p[:,2],self.lev,d)
      else:
        datout[:,id]=d
    return datout  
      
  def __str__(self):
    strout='Constant data %s:\n' % (self.id)
    for iv,v in enumerate(self.vars):
      strout+=' %s:%s\n' % (v,self.dat[iv])
    strout+=' levels:%s' % (self.lev if self.is3d else 'None')
    return strout

class GridData(FieldData):
  geod=True
  surfsub=False
  zinvert=False
  @copydoc(FieldData.__init__)
  def __init__(self,id,vars,**options):
    """  options:
    file: filename of gridded data file
    zcoord: direction of z-coordinate [up/down] (optional, default is 'up')
    **'up' means positive up and 'down' means positive down
    """
    #import cdms2,re,glob
    FieldData.__init__(self,id,vars,**options)
    self.keystore = []
    if 'file' not in options.keys():
      raise ConfigException('Grid data source must specify filename')
    self.filename = options.pop('file')
    self.filelist = self.load_files(self.filename)
    if len(self.filelist) == 0:
      raise ERConfigException('No files found that match template %s' % (filetmpl))

    self.filelist.sort()
    self.buf0={}
    self.vars=vars if isinstance(vars,list) else [vars]
    self.time=[]
    ncfile = nc.Dataset(self.filelist[0])
    cfile=ncfile.variables
    
    if cfile.has_key('time'): # time-varying data
      self.time=[]
      self.buf1={}
      self.bufstore={}
    else: # No time variable : tidal constituent file
      self.time=None
    

    for v in self.vars:
      if not cfile[v]:raise DataException('Variable %s not found in grid file %s' % (v,self.filename))
      if self.time is not None and cfile[v].shape[0] == cfile['time'].shape[0]:
        self.buf0[v]=None
        self.buf1[v]=None    
      else:
        # load variable data
        self.buf0[v]=cfile[v][:]

    # HORIZONTAL COORDINATES
    if self.geod:
      for vlon,vlat in [('lon','lat'),('longitude','latitude')]:
        if vlon in cfile.keys() and vlat in cfile.keys():
          lon, lat = cfile[vlon], cfile[vlat]      
          break
      if vlon not in cfile.keys() or vlat not in cfile.keys(): # if the allocation above didnt happen during the loop
        raise DataException('Dataset needs (lon, lat) or (longitude,latitude) variables if geod = True')
    else:
      for vlon,vlat in [('x','y'),('easting','northing'),('eastings','northings')]:
        if vlon in cfile.keys() and vlat in cfile.keys():
          x, y = cfile[vlon], cfile[vlat]   # 
          # would be better to differentiate name so there are no confusion ?
          # then need to split for case geod = true or false for interpolator definition below 
          break
      if vlon not in cfile.keys() or vlat not in cfile.keys(): # if the allocation above didnt happen during the loop
        raise DataException('Dataset needs (x, y) or (easting,northing) or (eastings,northings) variables if geod = False')

    # Check is it is a finite elements grid - look for triangulation variables
    # self.nv = True if 'nv' in cfile.keys() or 'elem' in cfile.keys() else None
    for var in ['nv', 'elem', 'ele']:
      if var in cfile.keys():
        self.nv = True
        self.tri = cfile[var][:,:]
        break
      else:
        self.nv = False

    for var in ['zlevels', 'lev', 'levels', 'level','nz']:
      #if var in cfile.keys(): >> This will not be correct where both 2d and 3d data are read from a single netcdf file (e.g. reading 3d current and 2D topo, or 2D tidal currents)
      if var in cfile[vars[0]].dimensions:
        try : 
          self.lev = cfile[var][:] # assumes that the "level" dimension and variable have the same name - not always true
          self.is3d = True
          break
        except : #"level" dimension and variable have different names - maybe we should simply modify the SELFE post-processing script ?
          #reloop through common names
          for var in ['zlevels', 'lev', 'levels', 'level','nz']:
            try :
              self.lev = cfile[var][:]
              self.is3d = True
              break
            except: 
              pass
      else:
        self.lev = None
        self.is3d = False

    if options.pop('zcoord','up')=='down':
      # maybe force input of zcoord to reduce confusion ?
      print 'Initialization - Inverting vertical levels of %s -  %s'  % (self.id,self.file)
      if self.is3d: self.lev=-self.lev
      self.zinvert=True

    if self.nv:
      # Is finite element grid
      lon=lon[:]
      lat=lat[:]
      if cfile.has_key('mask'):
        self.mask=numpy.where(cfile['mask'][:] != False )[0]
        lon=lon.take(self.mask)
        lat=lat.take(self.mask)
      else:
        self.mask=None
      self.interpolator=FEInterpolator(lon,lat,self.tri,self.lev,self.geod)
    elif lat and lon:
      self.mask=None # Mask not implemented for standard grids yet
      self.interpolator=RectInterpolator(lon[0],lat[0],lon[-1],lat[-1],
                                        (lon[1]-lon[0]) if len(lon)>1 else 0,
                                        (lat[1]-lat[0]) if len(lat)>1 else 0,
                                        self.lev,
                                        lat if self.geod else None)
    else:
       raise DataException('Gridded file %s structure not understood' % (self.filename)) 

    self.res=min(self.interpolator.x1-self.interpolator.x0,self.interpolator.y1-self.interpolator.y0)

    self.__dict__.update(options)

    ncfile.close()

    if self.time is None:
      bfile = nc.Dataset(self.filelist[0])
      self.files=[bfile]
    else: # Loop through files to get the full time vector
      self.timeindex=[0]
      self.files=[]
      self.flen=[]
      for filepath in self.filelist:
        bfile = nc.Dataset(filepath)
        self.files.append(bfile) #Open all the files
        # convert time in netcdf, to time convention used in ERcore i.e. NCEP/CF time, with fractions of days
        # A typical netcdf following the CF convention is expected to have time saved as :
        # 
        # double time(time) ;
        #   time:long_name = "time" ;
        #   time:units = "days since 1990-1-1 0:0:0" ;
        # 
        # Note the reference date "since XXXX" could be any date.

        # get start time string
        start_time_str = re.search('(?<=\s)\d.+$', bfile.variables['time'].units).group()   # get the file start time from the units attribute 
        # convert reference date to datetime instance
        try:
          start_time = datetime.datetime.strptime(start_time_str,'%Y-%m-%d %H:%M:%S')
        except:
          try:
            start_time = datetime.datetime.strptime(start_time_str,'%Y-%m-%d %H:%M')
          except:
            # import pdb;pdb.set_trace()
            start_time = datetime.datetime.strptime(start_time_str,'%Y-%m-%d')  # internal MSL, UDS-formatted file use 1-1-1
              

          # add more template here if required 
                
        # identify time units,:  seconds,hours,days, and convert to NCEP/CF fraction of days as used in ERcore
        
        if 'seconds' in bfile.variables['time'].units:
          deltas = [datetime.timedelta(seconds=float(t)) for t in bfile.variables['time'][:]] # deltas is incremental number of sedconds since file start               
        elif 'hours' in bfile.variables['time'].units:
          deltas = [datetime.timedelta(hours=float(t)) for t in bfile.variables['time'][:]] # deltas is incremental number of sedconds since file start    
        elif 'days' in bfile.variables['time'].units:
          deltas = [datetime.timedelta(days=float(t)) for t in bfile.variables['time'][:]] # deltas is incremental number of sedconds since file start    
          # there may be a bug when netcdf files are already ion fraction of days - to double check
        time0 = [ dt2ncep(start_time+delta) for delta in deltas ] # convert time to fraction of days - CF-compliant
          # if time is already as fraction of days - since 1-1-1     as in uds    
          # time0=bfile.variables['time'][:]

        if (len(self.time)>0) and (time0[0]<self.time[-1]):raise DataException('For templated time files times must be increasing - time in file %s less than preceeding file' % (bfile.filepath())) 
        self.time.extend(time0) #Add times in file to time list
        self.flen.append(len(time0))
        self.timeindex.append(len(self.time)) #Add start time index of next file
        #self.load_keystore(filepath) # probably related to cryptage so not needed ?
      self.flen.append(self.flen[-1])
      self.reset()
  

  def load_files(self, cfile):
    files = []
    if isinstance(cfile, list):
      for filename in cfile:
        filetmpl = re.sub('%Y|%m|%d|%H|%M','*', filename)
        files.extend(glob.glob(filetmpl))
    elif isinstance(cfile, (str,unicode)):
      filetmpl = re.sub('%Y|%m|%d|%H|%M','*', cfile)
      files.extend(glob.glob(filetmpl))
    files.sort()
    return files   
    
  def reset(self):
    """Reset time counter in file"""
    self.buftime=-1
    self.tind=0
    self.t0=self.time[0]
    self.t1=self.time[0]
    self.fileind0=0
    self.fileind1=0
    self.file0=self.files[0]
    self.file1=self.files[0]
    self.buftime=0
    

  def get(self,time): #Time is current model time
    """Get data slab for specified time
    Arguments:
      time: Time as NCEP decimal days
    """
    start = datetime.datetime.now()
    if time is None or self.time is None:
      return [self.buf0[v] for v in self.vars]
    if time==self.buftime:
      return [self.bufstore[v] for v in self.vars]
    readfile=(time<=self.time[0])
    while (time>self.t1) and (self.tind+1)<len(self.time):
      self.tind+=1
      self.t0=self.t1
      self.t1=self.time[self.tind]
      self.fileind0=self.fileind1
      if (self.tind>=self.timeindex[self.fileind0+1]):
        self.fileind1+=1
      readfile=True
    # for backwards in time
    while (time < self.t0) and (self.tind>1):
      self.tind -= 1  # index of t1
      self.t1=self.t0
      self.t0=self.time[self.tind-1]
      self.fileind1=self.fileind0
      if (self.tind <= self.timeindex[self.fileind0]):
        self.fileind0 -= 1
      readfile=True
    #print 'time,t0,t1,tind,fileind0,fileind1,readfile = %10.3f, %10.3f, %10.3f, %3i, %3i, %3i, %s' % (time, self.t0, self.t1, self.tind, self.fileind0, self.fileind1, readfile)
    if readfile:
      for v in self.vars:
        ind0=max(0,self.tind-self.timeindex[self.fileind0]-1)
        ind1=min(self.tind-self.timeindex[self.fileind1],self.flen[self.fileind1]-1)
        #print '%s %d %d %d %d %s %s' % (v,self.fileind0,self.fileind1,ind0,ind1, self.t0, self.t1)
        ncfile1 = self.files[self.fileind0]
        ncfile2 = self.files[self.fileind1]
        var1 = ncfile1.variables[v]
        var2 = ncfile2.variables[v]
        arr1 = var1[ind0][:]
        arr2 = var2[ind1][:]
        self.buf0[v]= arr1 
        self.buf1[v]= arr2
        if numpy.any(self.mask): # mask specified in netcdf 
          self.buf0[v]=self.files[self.fileind0][v][ind0].filled()
          self.buf1[v]=self.files[self.fileind1][v][ind1].filled()
        else:  # if no specific mask input - set "bad data" to 0.0
          id0 = numpy.where(numpy.abs(self.buf0[v][:])>1e10) or numpy.where(numpy.abs(self.buf0[v][:]) == 999) #or self.buf0[v][:] == ncfile2.variables[v].missing_value
          self.buf0[v][id0] = 0.0
          id1 = numpy.where(numpy.abs(self.buf1[v][:])>1e10) or numpy.where(numpy.abs(self.buf1[v][:]) == 999) 
          self.buf1[v][id1] = 0.0


    if (self.tind==0) and (time<self.time[0]):
      print 'Warning: model time (%s) before start time (%s) of data %s' % (ncep2dt(time), ncep2dt(self.time[0]), self.id)
    elif (self.tind==len(self.time)-1) and (time>self.time[-1]):
      print 'Warning: model time (%s) after end time (%s) of data %s' % (ncep2dt(time), ncep2dt(self.time[-1]), self.id)

    if self.t0==self.t1:
      tfac=0.
    else:
      tfac=min(max(0,(time-self.t0)/(self.t1-self.t0)),1)
    out=[]
    for v in self.vars:
      dat=(1.-tfac)*self.buf0[v]+tfac*self.buf1[v]
      if self.is3d and self.surfsub:
        if self.nv:
          dat=dat[:,:,:]-dat[0,:,:]
        else:
          dat=dat[:,:,:]-dat[0,:,:]
      out.append(dat)
      self.bufstore[v]=dat
    self.buftime=time
    return out
  
  @copydoc(FieldData.interp)
  def interp(self,p,time=None,imax=1):
    dat=self.get(time)
    datout=numpy.zeros((len(p),imax))
    for id,d in enumerate(dat[:imax]):
      datout[:,id]=self.interpolator(d,p)
    # Note on masking bad data : for now, if no specific mask are avaialble in the netcdf
    # bad data is set to 0.0 in the self.get function  
    # there may be a more elegant way to do this ? ideally a dynamic mask including 
    # depth and elevation (when available?)
    return datout
    
  def __str__(self):
    timestr=self.time[0].strftime()+' to '+self.time[-1].strftime() if self.time else 'Static'
    return """Grid data %s:
  File:%s
  Variables:%s
  Time:%s
  Levels:%s""" % (self.id,self.file,','.join(self.vars),timestr,self.lev if self.is3d else 'None')


class ConstantTide(ConstantData):
  def __init__(self,id,vars,**options):
    """  options:
    constant values for each of the variables"""
    FieldData.__init__(self,id,vars,**options)
    if self.is3d and not options.has_key('levels'):
      raise ConfigException('Constant 3D data must specify levels')
      self.lev=numpy.array(map(float,options['levels']))
      if options.get('levels','depth')=='depth':self.lev=-self.lev
    self.dat=[0.,0.,0.]
    if not options.has_key('cons'):
      raise ConfigException('Constituent ids missing')
    self.cons=options['cons']
    for iv,v in enumerate(self.vars):
      vamp=v+'_amp'
      vpha=v+'_pha'
      if not options.has_key(vamp):
        raise ConfigException('Amplitude data for variable %s missing' % (v))
      if not options.has_key(vpha):
        raise ConfigException('Phase data for variable %s missing' % (v))
      amp=options[vamp] if isinstance(options[vamp],list) else [options[vamp]]
      pha=options[vpha] if isinstance(options[vpha],list) else [options[vpha]]
      if (len(amp)<>len(self.cons)) or (len(pha)<>len(self.cons)):
        raise ConfigException('Number of amplitudes and phases must match number of constituents')
      self.tidestr[iv]=TideStr(amp,pha,icons=self.cons,lat=options.get('lat',0.))
  
  def interp(self,p,time,imax=2):
    self.dat=[self.tidestr[tide].ts([time]) for tide in self.tidestr]
    ConstantData.interp(self,time,imax=imax)
    
  
class GriddedTide(GridData): # Tidal consituent grid
  def __init__(self,id,vars,**options):
    from ercore.lib.tide import TideStr
    import datetime
    
    GridData.__init__(self,id,[vars[0]+'_amp'],**options)
    self.vars=vars
    fcons=[cons.tostring().rstrip().rstrip('\x00').upper() for cons in self.files[0].variables['cons'][:]] # ercore_nc
    
    consindex=[]
    self.cons=[]
    consreq=options.get('cons',fcons)
    for icons,cons in enumerate(fcons):
      if cons in consreq:
         consindex.append(icons)
         self.cons.append(cons)
    self.ncons=len(consindex)
    self.amp={}
    self.phac={}
    self.phas={}
    #import pdb;pdb.set_trace() 
    #lat0 = self.files[0].lat0 # this fails if lat0 is not an attribute of netcdf file
    lat0 = numpy.mean(self.files[0].variables['lat'][:]) # do mean of latitude instead, is that what it is meant to be ? 
    # if lat==0.0  the nodal phase modulation u, and nodal amplitude correction f returned by tide.py will be 0
    # what about t0 ?
    self.tidestr=TideStr(numpy.zeros((self.ncons,1)),numpy.zeros((self.ncons,1)),self.cons,options.get('t0',datetime.datetime.now()),lat0)
    for iv,v in enumerate(self.vars):
      vamp=v+'_amp'
      vpha=v+'_pha'
      if not self.files[0].variables[vamp]:
        raise ConfigException('Amplitude data for variable %s missing' % (v))
      if not self.files[0].variables[vpha]:
        raise ConfigException('Phase data for variable %s missing' % (v))
      ncfile = self.files[0]
      varamp = ncfile.variables[vamp]
      varpha = ncfile.variables[vpha]
      arramp = varamp[consindex,:]
      arrpha = varpha[consindex,:]
      #keyamp = self.get_key(0, vamp)
      #keypha = self.get_key(0, vpha)
      #arramp = decrypt_var(varamp, keyamp, arramp)
      #arrpha = decrypt_var(varpha, keypha, arrpha)
      if numpy.any(self.mask):
        arramp = arramp.take(self.mask, -1)
        arrpha = arrpha.take(self.mask, -1)
      self.amp[v]=arramp
      self.phac[v]=numpy.cos(arrpha)
      self.phas[v]=numpy.sin(arrpha) 
      if (self.amp[v].shape[0]<>self.ncons) or (self.phac[v].shape[0]<>self.ncons):
        raise ConfigException('First dimension of amplitudes and phases must match number of constituents')
  
  def interp(self,p,time,imax):
      self.tidestr.amp=numpy.zeros((len(self.cons),len(p)))
      phac=numpy.zeros((len(self.cons),len(p)))
      phas=numpy.zeros((len(self.cons),len(p)))
      datout=numpy.zeros((len(p),imax))
      for i,v in enumerate(self.vars):
        if i>imax-1:break
        v=self.vars[i]
        for ic,a in enumerate(self.amp[v]):
          self.tidestr.amp[ic]=self.interpolator(a,p)
          phac[ic]=self.interpolator(self.phac[v][ic],p)
          phas[ic]=self.interpolator(self.phas[v][ic],p)
        self.tidestr.pha=numpy.arctan2(phas,phac)
        datout[:,i]=(self.tidestr.ts(ncep2dt(time)))
      return datout

class ConstantTopo(ConstantData):
  pass

class ConstantMover(ConstantData):
  topo=None
  z0=0.0
  def interp(self,p,time=None,imax=2):
    uu=ConstantData.interp(self,p,time,imax)
    if self.topo:
      topo=self.topo.interp(p,None,3)
      if (not self.is3d) and (self.z0>0):#Log profile for 2D case
        uu[:,:2]=uu[:,:2]*(numpy.log(dz/self.z0))/(numpy.log(topo/self.z0)-1)
      if (imax==3):  #Vertical velocity correction for slope
        w1=slope_correction(p,topo,uu)
        # uu[:,2]+=w1
        # the below is probably not right...need to include dhdx/ dhdy see infos here
        # C:\metocean\Google Drive\R&D\ercore_dev\bathys_gradients
        uu[:,2]=0.0
    return uu

class ConstantReactor(ConstantData):
  pass

class ConstantDiffuser(ConstantData):
  pass

#Variable diffusion based on Okubo 4/3 law
class VariableDiffuser(FieldData):
  alpha=0.005
  mindiffH=0.001
  maxdiffH=10
  #maxdistH=5000
  geod=True
  @copydoc(FieldData.__init__)
  def __init__(self,id,vars,**options):
    FieldData.__init__(self,id,vars,**options)
    if not hasattr(self,'P0'):
      raise ConfigException('Variable diffusion must have release origin specified')
    self.a1=86400.**1.3333*0.0001*self.alpha
    #self.mfa=numpy.array([[numpy.cos(D2R*self.P0[1])/ARAD,1./ARAD]]) if self.geod else 1.
  
  @copydoc(FieldData.interp)
  def interp(self,p,time,age,imax=2):
    diff=numpy.zeros((len(p),imax))
    #l2=numpy.minimum(((self.mfa*(p[:,:2]-self.P0[:2]))**2).sum(1),self.maxdistH**2)[:,None]
    diff=numpy.tile(self.a1*age[:,None]**1.3333,(1,imax))
    diff[:,:2]=numpy.minimum(numpy.maximum(diff[:,:2],self.mindiffH),self.maxdiffH)
    if imax==3:
      diff[:,2]=self.diffz
    return diff
    
 
class GriddedTopo(GridData):
  refloat=0.
  @copydoc(GridData.__init__)
  def __init__(self,id,vars,**options):
    GridData.__init__(self,id,vars,**options)
    v=self.vars[0]
    if self.zinvert:self.buf0[v]=-self.buf0[v]
    # we should probably be consistent with GriddedMover and use zcoords rather than zinvert ?
    dhdy,dhdx=self.interpolator.grad(self.buf0[v])
    self.buf0['dhdx']=dhdx
    self.vars.append('dhdx')
    self.buf0['dhdy']=dhdy
    self.vars.append('dhdy')

  def intersect(self,pos,post,state,t1,t2):
    """Test for intersection of particles with topo
    Arguments:
      pos: previous particle positions
      post: new particle positions
      #   state: particle state
    Returns:
      New particle positions after intersection (state modified in place)
    """
    dep1=self.interp(pos,imax=1)[:,0]  # particle depths at t1
    dep2=self.interp(post,imax=1)[:,0] # particle depths at t2
    # -1* to ensure that outputs False if deps are actually the same.
    # identify active particles (state=1) that will reach the seabed at t2 (post)
    
    ind=((post[:,2]-dep2 < -1*ALMOST_ZERO) & (state>0)) # index of particles that would go below seabed
    pout=post[:,:]
    if ind.sum(): # some particles would go below seabed from t1 to t2
    # find where they would intersect
      #---not necessary anymore----------------------------------
      # first check if there were any particles that were ALREADY below seabed at t1
      #  
      # id_check1 = (pos[:,2]-dep1 < -1*ALMOST_ZERO)
      # pout[id_check1,:]=pos[id_check1,:] # if they were, they couldnt have moved - set back to pos
      # pout[id_check1,2]=pos[id_check1,2] # reset particle depth as dep1
      #-------------------------------------
      
      denom=(dep2[ind]-dep1[ind]+pos[ind,2]-post[ind,2]) #  use Thales theorem to work out intersection point
      f=(pos[ind,2]-dep1[ind])/denom 
      f[denom == 0.0] = 0.0 # stay at same position
      pout[ind,:]=pos[ind,:]+f[:,None]*(post[ind,:]-pos[ind,:]) 
      state[ind]+=1 # state set to 2 to identify it is on seabed
      print '%s particle(s) intersected seabed' %(ind.sum())
      # import pdb;pdb.set_trace()
      pout[ind,2] = numpy.maximum(dep2[ind],pout[ind,2]) # ensure that particle depth is not below seabed
      # the line pos[ind,:]+f[:,None]*(post[ind,:]-pos[ind,:]) may not work if
      # post[ind,2] == pos[ind,2] as for the case of passive tracer, with no diffz
      # 
      # This issue should be fixed when we have terms accounting for bathy gradients
    return pout

def intersect_free_surface(self,pos,post,state,t1,t2):
  """Test for intersection of particles with free surface
  Arguments:
    pos: previous particle positions
    post: new particle positions
  Returns:
    New particle positions after intersection (state of particles that intersected the surface
    is set back to 1 
  """
  elev1=self.interp(pos,t1,imax=1)[:,0]
  elev2=self.interp(post,t2,imax=1)[:,0]
  # to ensure that outputs False if elevs are actually the same.  
  ind=((post[:,2]-elev2 > ALMOST_ZERO) & (state>0))
  pout=post[:,:]
  if ind.sum():
    # using https://github.com/metocean/ercore/blob/ercore_nc/ercore/fields.py
    # import pdb;pdb.set_trace()  
    pout[ind,2] = elev2[ind] # simply set new elevation to elev2
    state[ind]=1 # state set to 1 i.e. keep active
    
    # note we could use something of the type below - looking for the actual intersection, see here
    #  but this is probably overkill given the small slopes expected for the sea surface
    # https://github.com/metocean/ercore/blob/ercore_nc_tests/ercore/fields.py
    #   denom=(elev2[ind]-elev1[ind]+pos[ind,2]-post[ind,2])
    #   f=(pos[ind,2]-elev1[ind])/denom
    #   pout[ind,:]=pos[ind,:]+f[:,None]*(post[ind,:]-pos[ind,:])     
 

  return pout


class ConstantElevation(ConstantData):
  @copydoc(ConstantData.__init__)
  def intersect(self,pos,post,state,t1,t2):
    return intersect_free_surface(self,pos,post,state,t1,t2)

class GriddedElevation(GridData):
  @copydoc(GridData.__init__)
  def intersect(self,pos,post,state,t1,t2):
    return intersect_free_surface(self,pos,post,state,t1,t2)


class TidalElevation(GriddedTide):
  @copydoc(GridData.__init__)
  def intersect(self,pos,post,state,t1,t2):
    return intersect_free_surface(self,pos,post,state,t1,t2)


class GriddedMover(GridData):
  z0=0.001 # setting default z0 to 0.001 m - must be >0 otherwise the log profile is not used  
  topo=None
  @copydoc(FieldData.interp)
  def interp(self,p,time=None,imax=2):
    
    uu=GridData.interp(self,p,time,imax)
    #import pdb;pdb.set_trace()
    #print '%s' % (uu[-1,:])
    #print '%s' % (time)
    #print '%s' % (p[-1,:])
    #Correct for vertical motion   
    if self.topo:
      topo=self.topo.interp(p,None,3)
      dz=numpy.maximum(p[:,2]-topo[:,0],self.z0) # particle depth above seabed (positive)
      if (not self.is3d) and (self.z0>0.0):#Log profile for 2D depth averaged currents
        # use of standard log profile to extrapolate depth-averaged current to particls depths p[:,2]
        uu[:,:2]=uu[:,:2]*(numpy.log(dz[:,numpy.newaxis]/self.z0))/(numpy.log(abs(topo[:,0:1])/self.z0)-1)
      if (imax==3):  #Vertical velocity correction for bathymetric gradients
        # the below is probably not right...need to include dhdx/ dhdy see infos here
        # C:\metocean\Google Drive\R&D\ercore_dev\bathys_gradients
        w1=slope_correction(p,topo,uu)
        # uu[:,2]+=w1
        uu[:,2]=0.0
      # we should probably apply a log profile for the 3D case for the region from the last wet bin to bottom ? 
    return uu

class TidalMover(GriddedTide):
  z0=0.001 # default z0
  topo=None
  @copydoc(FieldData.interp)
  def interp(self,p,time=None,imax=2):
    uu=GriddedTide.interp(self,p,time,imax)
    #Correct for vertical motion
    if self.topo:
      topo=self.topo.interp(p,None,3)
      dz=numpy.maximum(p[:,2]-topo[:,0],self.z0)
      if (not self.is3d) and (self.z0>0.0):#Log profile for 2D case
        uu[:,:2]=uu[:,:2]*(numpy.log(dz[:,numpy.newaxis]/self.z0))/(numpy.log(abs(topo[:,0:1])/self.z0)-1)
      if (imax==3): #Vertical velocity correction for slope
        w1=slope_correction(p,topo,uu)
        uu[:,2]+=w1
    return uu

class GriddedReactor(GridData):
  pass

class GriddedDiffuser(GridData):
  pass

class GriddedSticker(GridData):
  pass

class GriddedDataGroup(FieldData):
  def __init__(self,id,vars,members,**options):
    FieldData.__init__(self,id,vars,**options)
    self.members=members if isinstance(members,list) else [members]
  
  def interp(self,p,time=None,imax=1):
    if len(self.members)==1:#Trivial case of only one member
      return self.members[0].interp(p,time,imax)
    datout=numpy.zeros((len(p),imax))
    pmask=self.members[0].interpolator.ingrid(p)
    if pmask.any():
      datout[pmask,:]=self.members[0].interp(p[pmask,:],time,imax)
    for member in self.members[1:]:
      nmask=member.interpolator.ingrid(p)
      nmask=nmask&~pmask
      if not nmask.any():continue
      datout[nmask,:]=member.interp(p[nmask,:],time,imax)
      pmask=pmask|nmask
    return datout
  
  def intersect(self,pos,post,state,t1,t2):
    if len(self.members)==1:#Trivial case of only one member
      return self.members[0].intersect(pos,post,state,t1,t2)
    pout=post[:,:]
    pmask=self.members[0].interpolator.ingrid(post)
    if pmask.any():
      pout[pmask,:]=self.members[0].intersect(pos[pmask,:],post[pmask,:],state[pmask],t1,t2)
    for member in self.members[1:]:
      nmask=member.interpolator.ingrid(post)
      nmask=nmask&~pmask
      if not nmask.any():continue
      pout[nmask,:]=member.intersect(pos[nmask,:],post[nmask,:],state[nmask],t1,t2)
      pmask=pmask|nmask
    return pout
    
  
    
   
if __name__=="__main__":
  import sys
  print ConstantData.__init__.__doc__
  import yaml
  t=yaml.load("""
ConstantTopo:
  id:test
  vars:topo
  topo:10""")
  t=ConstantTopo('test','topo',topo=10)
  print t
  
  

