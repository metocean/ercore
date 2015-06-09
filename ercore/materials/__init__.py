#!/usr/bin/env python
import numpy
import copy
from ercore import dt2ncep,parsetime,ObjectList,ERRuntimeException
from ercore.lib import pres,temppot,dens
from ercore._flib_ercore import slipvel

R2D=180./numpy.pi
D2R=1/R2D
ARAD=R2D/6367456.
R1=1/8.31
SPI=0.5235987755982 #pi/6
PI2=2*numpy.pi

def eqnstate(P,T,Mg,Z=1.0):
  return R1*P*Mg/Z/(T+273)

#Base class for all materials - all materials must inherit from this class 
class _Material:
  """Initialization:
    <MaterialClass>(id,nbuff,movers,reactors,stickers,diffusers,tstart,tend,outfile,**properties)
    Arguments:
      id: Unique id
      nbuff: Total number of particles in buffer
      movers: List of mover id strings
      reactors: List of reactor id strings
      stickers: List of sticker id strings
      diffusers: List of diffuser id strings
      tstart: Starting time for release
      tend: Ending time for release
      tstep: Timestep of release
      outfile: Filename of output file
      P0: Initial position of release
      spawn: Number of spawned particles (per day)
      reln: Number of particles per release
      R0: Total release of material
      Q0: Flux of material (per day) 
      **properties: Optional keyword arguments specifying additional properties
      
    Properties:
    """
  is3d=True
  geod=False
  default_props={}
  status_codes={0:'Not released',1:'Released and active',-1:'Stuck to shoreline or bottom',-2:'Dead'}
  def __init__(self,id,nbuff,movers=[],reactors=[],stickers=[],diffusers=[],tstart=None,tend=None,tstep=0.,outfile=None,P0=[0,0,0],spawn=1,reln=0,R0=1.,Q0=1.,unstick=0.,**prop):
    self.id=id
    self.np=0
    self.ninc=1 #Counter for unique numbering
    self.tstep=prop.get('tstep',0.)
    self.npmax=nbuff
    self.tstart=parsetime(tstart) if tstart else 0.
    self.tend=parsetime(tend) if tend else 1.e10
    self.state=numpy.zeros((nbuff+1),'i') #State of particle
    self.age=numpy.zeros((nbuff+1))
    self.mass=numpy.ones((nbuff+1)) #Mass is used for particle weighting - could also represent volume or counts
    self.nid=numpy.zeros((nbuff+1)) #Array for unique numbering
    self.props=copy.copy(self.default_props)
    self.props['P0']=P0
    self.props['spawn']=spawn
    self.props.update(prop)
    self.pos=numpy.array(P0)*numpy.ones((nbuff+1,3))
    self.post=numpy.array(P0)*numpy.ones((nbuff+1,3))
    self.reln=reln #Particles per release
    self.R=R0 #Total release of material
    self.Q=Q0 # Flux of material per day
    self.unstick=unstick # Can become unstuck - number is halflife
    self.movers=movers
    self.reactors=reactors
    self.stickers=stickers
    self.diffusers=diffusers
    self.mfx=numpy.ones((nbuff+1,3))
    self.arrays=[]
    self.children={}
    self.outfile=outfile if outfile else 'ercore.'+self.id+'.out'
    self.relsumt=0.
    if self.tstart!=self.tend:self._npt=1.*self.reln/abs(self.tend-self.tstart)
    
  def yamlstr(self):
    str="""
  class: %s
  #Unique id for this material
  id: %s
  #[List of mover id strings]
  movers:
  #[List of reactor id strings]
  reactors:
  #[List of sticker id strings]
  stickers:
  #[List of diffuser id strings]
  diffusers:
  #Start time of release as timestring
  tstart:
  #End time of release as timestring
  tend:
  #Timestep of material dynamics
  tstep:0.
  #Output filename
  outfile:
  #Initial position of release
  P0:[0,0,0]
  #Spawning rate
  spawn:1
  #Release number of particles (for complete release)
  reln:0
  #Total release mass
  R0:1.
  #Flux of material per day (optional)
  Q0:0       
  """ % (self.__class__.__name__,self.__class__.__name__+'1')

  def initialize(self,t1,t2):
    """Initialize material at time t1  """
    self.tcum=self.tstep
    if not self.tstart:self.tstart=t1
    if not self.tend:self.tend=self.tstart
    if self.Q<=0.:self.Q=self.R/(abs(self.tend-self.tstart) if self.tend!=self.tstart else 1.) #See release below - this makes self.Q=self.R at release time
  
  def geodcalc(self,init=False):
    """Calculate the map factors for geodetic coordinates"""
    if init:
      self.mfx[:,1]=numpy.tile(ARAD,self.npmax+1)
      self.mfx[:,0]=self.mfx[:,1]/numpy.cos(D2R*self.pos[:,1])
      self.geod=True
    self.mfx[:self.np,0]=self.mfx[:self.np,1]/numpy.cos(D2R*self.pos[:self.np,1])
     
  def __str__(self):
    """String representation of material"""
    if self.np>0:
      minn=numpy.minimum(self.pos[:self.np],0)
      maxx=numpy.maximum(self.pos[:self.np],0)
    else:
      minn=maxx=numpy.nan
    props='\n'.join([p+' '+self.props[p] for p in self.props])
    return """%s of %d particles (buffer %d)
  Properties:
  %s
  Spread:
  %f to %f x-direction
  %f to %f y-direction
  %f to %f z-direction""" % (self.type,self.np,self.pos.shape[0],tuple(props),minn[0],maxx[0],minn[1],maxx[1],minn[2],maxx[2])
  
    
  def __index__(self,ind):
    """Return particle at index"""
    if ind>self.np:return None
    return Particle(self.pos[ind,:],self.state[ind],self.age[ind],self.mass[ind],self.props)
    
  def _reset(self,i0):
    """Reset and shuffle arrays after particles removed"""
    if len(self.arrays)==0:
      for a in dir(self):
        if isinstance(getattr(self,a),numpy.ndarray):
          self.arrays.append(getattr(self,a))
    if isinstance(i0,int):
      for a in self.arrays:
        a[:-i0]=a[i0:]
        a[-i0:]=a[-1]
      self.np-=i0
      return True
    else:
      nind=i0.sum()
    if nind==0:return False
    self.np-=nind
    for a in self.arrays:
      a[:-nind]=a[~i0]
      a[self.np:self.np+nind]=a[-1]
    return True
    
  def fheader(self):
    """Return file header for output"""
    return 'Time\tid\tx\ty\tz\tstate\tage\tmass\n'
   
  def sfprint(self,t):
    """Return string dump of all particles at specified timestep"""
    str=''
    for i in range(0,self.np):
      str+="%f\t%d\t%f\t%f\t%f\t%d\t%f\t%f\n" % ((t,self.nid[i])+tuple(self.pos[i])+(self.state[i],self.age[i],self.mass[i]))
    return str
    
  def release(self,t1,t2,**k):
    """Release all particles between time t1 and t2"""
    if t2>t1:
      if (self.tstart>t2) or (self.tend<t1):return
      dt1=t2-self.tstart
      dt2=self.tend-t1
    else:
      if (t2>self.tstart) or (t1<self.tend):return
      dt1=self.tstart-t2
      dt2=t1-self.tend
    dt=min(abs(t2-t1),dt1,dt2,abs(self.tend-self.tstart))
    if k.has_key('nprel'): #Prescribed release size 
      np=k['nprel']
    elif dt==0: #Start and end time the same
      np=k.get('nprel',self.reln)
      dt=1.
    else: #Incremental release
      self.relsumt+=abs(dt)
      np=k.get('nprel',abs(int(self.relsumt*self._npt)))
      if np>0:self.relsumt-=1.0*np/self._npt
    if np==0:return
    np=int(np)
    nmax=self.npmax-self.np   
    if np>nmax:
      print 'Warning: particles exhausted for '+self.id
      np=nmax
      if np==0:return 0
    np1=self.np+np
    self.mass[self.np:np1]=self.Q*abs(dt)/np
    self.mass0=self.mass[self.np]
    self.state[self.np:np1]=1 #Released
    self.age[self.np:np1]=0.
    self.nid[self.np:np1]=range(self.ninc,self.ninc+np)
    self.ninc+=np
    for key in k:
      if hasattr(self,key):
        array=getattr(self,key)
        if isinstance(array,numpy.ndarray):
          array[self.np:np1]=k[key][:np]
        else:
          array=k[key]
    self.np=np1
    return np
    
  def advect(self,t1,t2,order=4):
    """Do advection for time t1 to t2"""
    pass
  
  def diffuse(self,t1,t2):
    """Do diffusion for t1 to t2"""
    pass
  
  def react(self,t1,t2):
    """Do reaction for t1 to t2"""
    pass
  
  def spawn(self,t1,t2):
    """Do spawning for t1 to t2"""
    return None
  
  def stick(self,t1,t2):
    """Do sticking for t1 to t2"""
    if self.np<1:return
    np=self.np
    posi=numpy.where(self.state[:np,None]<0,self.pos[:np,:],self.post[:np,:])
    for sticker in self.stickers:
      self.post[:self.np,:]=sticker.intersect(self.pos[:self.np,:],posi,self.state[:self.np])
      if self.unstick<=0.:
        self.state[self.state>1]=-1
  
  def die(self,t1,t2):
    """Kill particles between times t1 and t2 and remove from simulation"""
    if self.np==0:return
    dead=(self.age>self.props.get('maxage',1.e20)) | (self.mass<=0.0001*self.mass0)
    self.state[dead]=-2
    
  
class PassiveTracer(_Material):
  __doc__=_Material.__doc__
  #Movers:currents0,currents1...
  #Diffusers:diffh,diffv...
  def advect(self,t1,t2,order=4):
    np=self.np
    if np==0:return
    imax=3 if self.is3d else 2
    if len(self.movers)==0:
      self.post[:np,:imax]=self.pos[:np,:imax]
      return
    dt=86400.*(t2-t1)
    try:
#4th order RungeKutta advection
      kx1=self.movers[0].interp(self.pos[:np,:],t1,imax=imax)
      for mover in self.movers[1:]:
        kx1+=mover.interp(self.pos[:np,:],t1,imax=imax)
      kx1*=dt*self.mfx[:np,:imax]
      if order==1:
        self.post[:np,:imax]=self.pos[:np,:imax]+kx1
        return
      pxt=self.pos[:np,:imax]+0.5*kx1
      t12=t1+0.5*(t2-t1)
      kx2=0.5*(self.movers[0].interp(pxt,t1,imax=imax)+self.movers[0].interp(pxt,t12,imax=imax))
      for mover in self.movers[1:]:
        kx2+=0.5*(mover.interp(pxt,t1,imax=imax)+mover.interp(pxt,t12,imax=imax))
      kx2*=dt*self.mfx[:np,:imax]
      if order==2:
        self.post=self.pos[:np,:imax]+kx2
        return
      if order==3:
        pxt=self.pos[:np,:imax]+2*kx2-kx1
      elif order==4:
        pxt=self.pos[:np,:imax]+0.5*kx2
      kx3=0.5*(self.movers[0].interp(pxt,t1,imax=imax)+self.movers[0].interp(pxt,t12,imax=imax))
      for mover in self.movers[1:]:
        kx3+=0.5*(mover.interp(pxt,t1,imax=imax)+mover.interp(pxt,t12,imax=imax))
      kx3*=dt*self.mfx[:np,:imax]
      if order==3:
        self.post[:np,:imax]=self.pos[:np,:imax]+(1/6.)*(kx1+4*kx2+kx3)
      elif order==4:
        pxt=self.pos[:np,:imax]+kx3
        kx4=self.movers[0].interp(pxt,t2,imax=imax)
        for mover in self.movers[1:]:
          kx4+=mover.interp(pxt,t2,imax=imax)
        kx4*=dt*self.mfx[:np,:imax]
        self.post[:np,:imax]=self.pos[:np,:imax]+(1/6.)*(kx1+2*(kx2+kx3)+kx4)
    except:
      import traceback
      raise ERRuntimeException(traceback.format_exc())
      
  def diffuse(self,t1,t2):
    if (len(self.diffusers)==0) or self.np==0:return
    dt=86400.*abs(t2-t1)
    np=self.np
    imax=3 if self.is3d else 2
    for diffuser in self.diffusers:
      # Diffusion references:
      # Garcia-Martinez and Tovar, 1999 - Computer Modeling of Oil Spill Trajectories With a High Accuracy Method
      # Lonin, S.A., 1999. Lagrangian model for oil spill diffusion at sea. Spill Science and Technology Bulletin, 5(5): 331-336 
      diff=(6*dt*diffuser.interp(self.pos[:np],t1,self.age[:np],imax=imax))**0.5
      self.post[:np,:imax]+=numpy.random.uniform(-diff,diff,size=(np,imax))*self.mfx[:np,:imax] #self.mfx=map factors i.e. meters to lat/lon
      # correction for vertical diffusion resulting in above sea-surface Zlevel
      self.post[:np,2]=numpy.minimum(self.post[:np,2],0.0)
    
    
class BuoyantTracer(PassiveTracer):
  __doc__=PassiveTracer.__doc__+"""
    w0: Rise velocity (m/s) [-ve for sinking]
  """
  default_props={'w0':0.0}
  
  def initialize(self,t1,t2):
    PassiveTracer.initialize(self,t1,t2)
    self.w0=numpy.tile(self.props.get('w0',0.),self.npmax+1)
      
  def advect(self,t1,t2,order=4):
    if self.np==0:return
    PassiveTracer.advect(self,t1,t2,order)
    self.post[:self.np,2]+=86400.*(t2-t1)*self.w0[:self.np]
    
    
class Drifter(PassiveTracer):
  #Reactors:wind
  default_props={'dw_min':0.01,'dw_max':0.04,'cw_min':0.0,'cw_max':0.0}
  __doc__=PassiveTracer.__doc__+"""
    dw_min: Minimum downwind windage factor
    dw_max: Maximum downwind windage factor
    cw_min: Minimum downwind windage factor
    cw_max: Maximum downwind windage factor"""
  is3d=False
  def advect(self,t1,t2,order=4):
    if (self.np==0): return
    PassiveTracer.advect(self,t1,t2)
    if len(self.reactors)==0:return
    np=self.np
    dt=86400.*(t2-t1)
#Windage with cross and down-wind ranges
    wvec=0.5*(self.reactors[0].interp(self.pos[:np],t1,imax=2)+self.reactors[0].interp(self.post[:np],t2,imax=2))
    wsp=((wvec**2).sum(1))**0.5
    wdir=numpy.arctan2(wvec[:,1],wvec[:,0])
    wcos=numpy.cos(wdir)
    wsin=numpy.sin(wdir)
    dwr=numpy.random.uniform(self.props['dw_min'],self.props['dw_max'],np)
    if self.props['cw_max']>0.0:
      cwr=numpy.random.uniform(self.props['cw_min'],self.props['cw_max'],np)
      cwr=numpy.where(numpy.random.randint(0,2,np),cwr,-cwr)
    else:
      cwr=0.
    self.post[:np,0]+=dt*self.mfx[:np,0]*wsp*(dwr*wcos-cwr*wsin)
    self.post[:np,1]+=dt*self.mfx[:np,1]*wsp*(dwr*wsin+cwr*wcos)
    
  
class BDTracer(BuoyantTracer):
   #Reactors:temp,salt - note that reactor ids must start with salt or temp
   #Default: air in water
  default_props={'IFT':0.0728,'temp':20,'visco':0,'Mg':0.02897,'nmol':0,'db':0.005}
  __doc__=BuoyantTracer.__doc__+"""
    IFT: Interfacial tension
    temp: Initial temperature
    visco: Ambient viscosity
    Mg: Molar mass
    nmol: Number of moles
    db: Bubble diameter (m)
  """
  
  def initialize(self,t1,t2):
    BuoyantTracer.initialize(self,t1,t2)
    if self.props['nmol']==0:
      temp=self.reactors['temp'].interp(self.pos[:1,:],t1)[:,0]
      salt=self.reactors['salt'].interp(self.pos[:1,:],t1)[:,0]
      Pw=-pres(self.pos[0,2],self.pos[0,1])
      T=temppot(salt, temp, 0, Pw)
      gdens=self.eqnstate(101300+9.81*Pw*dens(salt,temp,Pw),T)
      self.props['nmol']=gdens*SPI*self.props['db']**3/self.props['Mg']
      self.dens=numpy.tile(gdens,self.npmax+1)
    self.db=(self.props['nmol']*self.props['Mg']/SPI/self.dens)**0.33333
  
  def eqnstate(self,P,T):
    return eqnstate(P,T,self.props['Mg'])
      
  def react(self,t1,t2):
    np=self.np
    if np==0:return
    temp=self.reactors['temp'].interp(self.pos[:np,:],t1)[:,0]
    salt=self.reactors['salt'].interp(self.pos[:np,:],t1)[:,0]
    Pw=-pres(self.pos[:np,2],self.pos[:np,1] if self.geod else 0.)
    temp=temppot(salt, temp, 0, Pw) #Convert to absolute temperature
    adens=dens(salt, temp, Pw)
    self.dens[:np]=self.eqnstate(101300+9.81*Pw*adens,temp)
    self.db[:np]=(self.props['nmol']*self.props['Mg']/SPI/self.dens[:np])**0.33333
    self.w0[:np]=slipvel.bubble_slip(self.db[:np],adens,adens-self.dens[:np],temp,self.props['visco'],self.props['IFT'])
  


