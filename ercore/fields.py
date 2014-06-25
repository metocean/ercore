#!/usr/bin/env python
import numpy
from ercore import ERCoreException,ERConfigException,copydoc,ncep2dt
from scipy.io.netcdf import netcdf_file
from _flib_ercore import interph,interp3d,interpz,inpoly

R2D=180./numpy.pi
D2R=1/R2D
ARAD=R2D/6367456.
DMIN=1.e-10

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
  def __init__(self,lon,lat,lev=None,geod=True,mask=None):
    """Constructor for interpolator
    Arguments:
      lon: Array of longitudes for each node
      lat: Array of latitudes for each node
      geod: Geodetic coordinates (True/False)
      lev: Levels for 3D grid
    """
    from scipy.spatial import cKDTree
    self.lon=lon[~mask] if mask else lon
    self.lat=lat[~mask] if mask else lat
    self.mask=~mask if mask else None
    self.x0=min(lon)
    self.x1=max(lon)
    self.y0=min(lat)
    self.y1=max(lat)
    self.geod=geod
    self.lev=lev
    self.tree=cKDTree(numpy.vstack((self.lon,self.lat)).T)
    
  def __call__(self,dat,p):
    if self.mask:
      dat=dat[self.mask]
    dist,i=self.tree.query(p[:,:2],3)
    dist[dist<DMIN]=DMIN
    if dat.ndim==2:
      tmp=dist.sum(-1)*((1./dist)*dat.take(i,1)).sum(-1)
      return interpz(tmp.filled().astype('f'),p[:,2],self.lev)
    else:
      return dist.sum(-1)*((1./dist)*dat.take(i)).sum(-1)
    
  def ingrid(self,p):
    inbbox=(p[:,0]>=self.x0) & (p[:,0]<=self.x1) & (p[:,1]>=self.y0) & (p[:,1]<=self.y1)
    #inmesh=inpoly(p[inbbox,0],p[inbbox,1],self.bndx,self.bndy)
    #inbbox[inbbox]=inmesh
    return inbbox
    
  def grad(self,dat):
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
      if not options.has_key('levels'):raise ConfigException('Constant 3D data must specify levels')
      self.lev=numpy.array(map(float,options['levels']))
      if options.get('levels','depth')=='depth':self.lev=-self.lev
    self.dat=[]
    for v in self.vars:
      if not options.has_key(v):
        raise ERConfigException('Constant data for variable %s missing' % (v))
      vin=options[v] if isinstance(options[v],list) else [options[v]]
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
    zcoord: direction of z-coordinate [up/down] (optional)
"""
    import cdms2,re,glob
    FieldData.__init__(self,id,vars,**options)
    if 'file' not in options.keys():
      raise ConfigException('Grid data source must specify filename')
    self.filename=options.pop('file')
    filetmpl=re.sub('%Y|%m|%d|%H','*',self.filename)
    self.filelist=glob.glob(filetmpl)
    if len(self.filelist)==0:
      raise ERConfigException('No files found that match template %s' % (filetmpl))
    self.filelist.sort()
    self.buf0={}
    self.vars=vars if isinstance(vars,list) else [vars]
    self.time=[]
    cfile=cdms2.open(self.filelist[0])
    if cfile['time']:
      self.time=[]
      self.buf1={}
      self.bufstore={}
    else:
      self.time=None
    for v in self.vars:
      if not cfile[v]:raise DataException('Variable %s not found in grid file %s' % (v,self.filename))
      if cfile[v].getTime():
        self.buf0[v]=None
        self.buf1[v]=None    
      else:
        self.buf0[v]=cfile[v][:].asma()
    lon=cfile[v].getLongitude()
    lat=cfile[v].getLatitude()
    nv=cfile['nv']
    self.lev=None
    self.is3d=((lon and (len(cfile[v].shape)==4)) or (nv and (len(cfile[v].shape)==3)))
    if self.is3d:
      self.lev=cfile[v].getAxis(1)[:]
      if self.lev is None:
        raise DataException('3D data file must specify levels')
    if options.pop('zcoord','up')=='down':
      if self.is3d:self.lev=-self.lev
      self.zinvert=True
    if nv:
      lon=cfile['lon'][:] if cfile['lon'] else cfile['longitude'][:]
      lat=cfile['lat'][:] if cfile['lat'] else cfile['latitude'][:]
      self.interpolator=FEInterpolator(lon.asma(),lat.asma(),self.lev,self.geod,cfile['mask'])
    elif lon:
      self.interpolator=RectInterpolator(lon[0],lat[0],lon[-1],lat[-1],(lon[1]-lon[0]) if len(lon)>1 else 0,(lat[1]-lat[0]) if len(lat)>1 else 0,self.lev,lat if self.geod else None)
    else:
      raise DataException('Gridded file %s structure not understood' % (self.filename))
    self.res=min(self.interpolator.x1-self.interpolator.x0,self.interpolator.y1-self.interpolator.y0)
    self.__dict__.update(options)
    cfile.close()
    if self.time is None:
      self.files=[cdms2.open(self.filelist[0])]
    else:
      self.timeindex=[0]
      self.files=[]
      self.flen=[]
      for cfile in self.filelist:
        self.files.append(cdms2.open(cfile)) #Open all the files
        time0=[t.torelative('days since 1-1-1').value for t in self.files[-1]['time'].asRelativeTime()]
        if (len(self.time)>0) and (time0[0]<self.time[-1]):raise DataException('For templated time files times must be increasing - time in file %s less than preceeding file' % (cfile))
        self.time.extend(time0) #Add times in file to time list
        self.flen.append(len(time0))
        self.timeindex.append(len(self.time)) #Add start time index of next file
      self.flen.append(self.flen[-1])
      self.reset()
      
    
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
      time: Time as NCEP decimal
    """
    if time is None or self.time is None:
      return [self.buf0[v] for v in self.vars]
    if time==self.buftime:
      return [self.bufstore[v] for v in self.vars]
    readfile=False
    while (time>self.t1) and (self.tind+1)<len(self.time):
      self.tind+=1
      self.t0=self.t1
      self.t1=self.time[self.tind]
      self.fileind0=self.fileind1
      if (self.tind>=self.timeindex[self.fileind0+1]):
        self.fileind1+=1
      readfile=True
    if readfile:
      for v in self.vars:
        ind0=self.tind-self.timeindex[self.fileind0]
        ind1=min(self.tind-self.timeindex[self.fileind1]+1,self.flen[self.fileind1]-1)
        #print '%s %d %d %d %d' % (v,self.fileind0,self.fileind1,ind0,ind1)
        self.buf0[v]=self.files[self.fileind0][v][ind0]
        self.buf1[v]=self.files[self.fileind1][v][ind1]
    if (self.tind==0) and (time<self.time[0]):
      print 'Warning: model time before start time of data %s' % self.id
    elif (self.tind==len(self.time)-1) and (time>self.time[-1]):
      print 'Warning: model time after end time of data %s' % self.id
    tfac=min(max(0,(time-self.t0)/(self.t1-self.t0)),1)
    out=[]
    for v in self.vars:
      dat=(1.-tfac)*self.buf0[v]+tfac*self.buf1[v]
      if self.is3d and self.surfsub:
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
    
  
class GriddedTide(GridData):
  def __init__(self,id,vars,**options):
    from ercore.lib.tide import TideStr
    import datetime
    GridData.__init__(self,id,[vars[0]+'_amp'],**options)
    self.vars=vars
    fcons=[cons.tostring().rstrip().rstrip('\x00').upper() for cons in self.files[0]['cons'][:]]
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
    self.tidestr=TideStr(numpy.zeros((self.ncons,1)),numpy.zeros((self.ncons,1)),self.cons,options.get('t0',datetime.datetime.now()),options.get('lat',0.0))
    for iv,v in enumerate(self.vars):
      vamp=v+'_amp'
      vpha=v+'_pha'
      if not self.files[0][vamp]:
        raise ConfigException('Amplitude data for variable %s missing' % (v))
      if not self.files[0][vpha]:
        raise ConfigException('Phase data for variable %s missing' % (v))
      self.amp[v]=self.files[0][vamp][:].asma()[consindex]
      self.phac[v]=numpy.cos(self.files[0][vpha][:].asma()[consindex])
      self.phas[v]=numpy.sin(self.files[0][vpha][:].asma()[consindex])
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
  z0=0
  def interp(self,p,time=None,imax=2):
    uu=ConstantData.interp(self,p,time,imax)
    if self.topo:
      topo=self.topo.interp(p,None,3)
      if (not self.is3d) and (self.z0>0):#Log profile for 2D case
        uu[:,:2]=uu[:,:2]*(numpy.log(dz/self.z0))/(numpy.log(topo/self.z0)-1)
      if (imax==3):  #Vertical velocity correction for slope
        w1=p[:,2]/topo[:,0]*(uu[:,0]*topo[:,1]+uu[:,1]*topo[:,2])
        uu[:,2]+=w1
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
    dhdy,dhdx=self.interpolator.grad(self.buf0[v])
    self.buf0['dhdx']=dhdx
    self.vars.append('dhdx')
    self.buf0['dhdy']=dhdy
    self.vars.append('dhdy')
  
  def intersect(self,pos,post,state):
    """Test for intersection of particles with topo
    Arguments:
      pos: previous particle positions
      post: new particle positions
      state: particle state
    Returns:
      New particle positions after intersection (state modified in place)
    """
    dep1=self.interp(pos,imax=1)[:,0]
    dep2=self.interp(post,imax=1)[:,0]
    ind=((post[:,2]<dep2) & (state>0))
    pout=post[:,:]
    if ind.sum():
      denom=(dep2[ind]-dep1[ind]+pos[ind,2]-post[ind,2])
      f=(pos[ind,2]-dep1[ind])/denom
      pout[ind,:]=pos[ind,:]+f[:,None]*(post[ind,:]-pos[ind,:])
      state[ind]+=1
    pout[:,2]=numpy.minimum(pout[:,2],0)
    return pout

class GriddedMover(GridData):
  z0=0.
  topo=None
  @copydoc(FieldData.interp)
  def interp(self,p,time=None,imax=2):
    uu=GridData.interp(self,p,time,imax)
    #Correct for vertical motion
    if self.topo:
      topo=self.topo.interp(p,None,3)
      if (not self.is3d) and (self.z0>0):#Log profile for 2D case
        uu[:,:2]=uu[:,:2]*(numpy.log(dz/self.z0))/(numpy.log(topo/self.z0)-1)
      if (imax==3):  #Vertical velocity correction for slope
        w1=(p[:,2]/topo[:,0])*(uu[:,0]*topo[:,1]+uu[:,1]*topo[:,2])
        uu[:,2]+=w1
    return uu

class TidalMover(GriddedTide):
  z0=0.
  topo=None
  @copydoc(FieldData.interp)
  def interp(self,p,time=None,imax=2):
    uu=GriddedTide.interp(self,p,time,imax)
    #Correct for vertical motion
    if self.topo:
      topo=self.topo.interp(p,None,3)
      dz=p[:,2]-topo[:,0]
      if (not self.is3d) and (self.z0>0):#Log profile for 2D case
        uu[:,:2]=uu[:,:2]*(numpy.log(dz[:,numpy.newaxis]/self.z0))/(numpy.log(abs(topo[:,0:1])/self.z0)-1)
      if (imax==3): #Vertical velocity correction for slope
        uu[:,2]+=(1.-dz/topo[:,0])*(uu[:,0]*topo[:,1]+uu[:,1]*topo[:,2])
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
  
  def intersect(self,pos,post,state):
    if len(self.members)==1:#Trivial case of only one member
      return self.members[0].intersect(pos,post,state)
    pout=post[:,:]
    pmask=self.members[0].interpolator.ingrid(post)
    if pmask.any():
      pout[pmask,:]=self.members[0].intersect(pos[pmask,:],post[pmask,:],state[pmask])
    for member in self.members[1:]:
      nmask=member.interpolator.ingrid(post)
      nmask=nmask&~pmask
      if not nmask.any():continue
      pout[nmask,:]=member.intersect(pos[nmask,:],post[nmask,:],state[nmask])
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
  
  

