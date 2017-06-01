#!/usr/bin/env python
import numpy
import copy
from ercore import dt2ncep,parsetime,ObjectList,ERRuntimeException
from ercore.lib import pres,temppot,dens
from ercore._flib_ercore import slipvel
from shapely.geometry import Point,Polygon

R2D=180./numpy.pi
D2R=1/R2D
ARAD=R2D/6367456.
R1=1/8.31
SPI=0.5235987755982 #pi/6
PI2=2*numpy.pi

def get_random_point_in_polygon(nbuff,poly):
     (minx, miny, maxx, maxy) = poly.bounds
     matrix=numpy.zeros(shape=(nbuff,2))
     for np in range(0,nbuff):
        while matrix[np,0]==0:
                Xp=(maxx-minx)*numpy.random.random(1)+minx
                Yp=(maxy-miny)*numpy.random.random(1)+miny
                p = Point(Xp,Yp)
                if poly.contains(p):
                        matrix[np,0]=Xp
                        matrix[np,1]=Yp
     return matrix

def get_random_point_in_circle(nbuff,P0,radius):
  """
  define 'nbuff' random points within a circle of 
  center 'P0' [1x3] ([lon,lat,z] as input in Material class) 
  'radius' in meters
  """      
  #1 deg lat = 110 km, 1 deg lon= 111.32km*cos(lat) 
  # Note this will not make a perfect circle, but approximation is likely good enough.
  matrix=numpy.zeros(shape=(nbuff,2))
  lat_rad=numpy.pi*P0[1]/180 # latitude in radians
  deg_in_m_lat=111132.92-559.82 * numpy.cos(2* lat_rad)+1.175*numpy.cos(4*lat_rad)
  deg_in_m_lon=111412.84 * numpy.cos(lat_rad)- 93.5 * numpy.cos(3*lat_rad)
  #https://knowledge.safe.com/articles/725/calculating-accurate-length-in-meters-for-lat-long.html
  radius_lon=radius/deg_in_m_lon
  radius_lat=radius/deg_in_m_lat
  rand_angle=numpy.random.uniform(0.0,1.0,nbuff)
  rand_radius=numpy.random.uniform(0.0,1.0,nbuff)    
  matrix[:,0]=P0[0]+radius_lon*rand_radius*numpy.ones((1,nbuff))*numpy.cos(rand_angle*2*numpy.pi)
  matrix[:,1]=P0[1]+radius_lat*rand_radius*numpy.ones((1,nbuff))*numpy.sin(rand_angle*2*numpy.pi)
  return matrix

def eqnstate(P,T,Mg,Z=1.0):
  return R1*P*Mg/Z/(T+273)

#Base class for all materials - all materials must inherit from this class 
class _Material(object):
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
      P0: Initial position of release - Note the convention of particle vertical level Z is negative downards, where sea surface=0m i.e. -10 s 10 m below sea surface
      spawn: Number of spawned particles (per day)
      reln: Number of particles per release
      R0: Total release of material
      Q0: Flux of material (per day)
      unstick - flag 0/1 to specfiy if the material can unstick from its stickers, can be an array of length(stickers)
      ischild : flag 0/1 to identify if material is a child of another material (e.g. "spawned")  
      **properties: Optional keyword arguments specifying additional properties
      
    Properties:
    """
  is3d=True
  geod=False
  default_props={}
  status_codes={0:'Not released',1:'Released and active',-1:'Stuck to shoreline or bottom',-2:'Dead'}
  def __init__(self,id,nbuff,movers=[],reactors=[],stickers=[],diffusers=[],tstart=None,tend=None,tstep=0.,tstep_release=0.,outfile=None,P0=[0,0,0],spawn=1,reln=0,R0=1.,Q0=1.,unstick=0.,**prop):
    super(_Material, self).__init__()
    self.id=id
    self.np=0
    self.ninc=1 #Counter for unique numbering
    self.tstep=prop.get('tstep',0.)
    self.tstep_release=tstep_release
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
    if not hasattr(self.props,'ischild'):self.props['ischild']=0
    self.props.update(prop)
    # import pdb;pdb.set_trace()
    # Release Options
    #
    # Initialization of particle position vector 
    if numpy.size(P0[2])==1:
    # P0:single position x,y,z
    # pre-allocate with a single release position
      self.pos=numpy.array(P0)*numpy.ones((nbuff+1,3))
      self.post=numpy.array(P0)*numpy.ones((nbuff+1,3))
      self.dep = numpy.ones((nbuff+1))*-999.
      self.elev = numpy.zeros((nbuff+1)) 
    elif numpy.size(P0[2])==2:
    # P0:single x,y position but range of vertical z level
    # random position allocated within that range
      # Initialize variables - with lowest 
      P01=[P0[0],P0[1],min(P0[2])]
      self.pos=numpy.array(P01)*numpy.ones((nbuff+1,3))
      self.post=numpy.array(P01)*numpy.ones((nbuff+1,3))
      self.dep = numpy.ones((nbuff+1))*-999.
      self.elev = numpy.zeros((nbuff+1)) 
      #thickness of z layer
      dz=abs(P0[2][0]-P0[2][1])
      # depth are supposed to be negative
      zz=numpy.random.random(nbuff+1)
      self.pos[:,2]=min(P0[2])+dz*zz
      self.post[:,2]=min(P0[2])+dz*zz
      
    if "circular_radius" in self.props:
      #release in a circle rather than at a single X,Y point location
      point_in_circle = get_random_point_in_circle(nbuff+1,self.props['P0'],self.props['circular_radius'])
      self.pos[:,0]=point_in_circle[:,0]
      self.pos[:,1]=point_in_circle[:,1]
      self.post[:,0]=self.pos[:,0]
      self.post[:,1]=self.pos[:,1]
      self.dep = numpy.ones((nbuff+1))*-999.
      self.elev = numpy.zeros((nbuff+1))  
      # Particles release depths may need to be updated according to water depths at new locations within the circle
      for sticker in stickers:
        if 'GriddedTopo' in sticker.__class__.__name__:
           topo=sticker.interp(self.pos,imax=1)# get depths new particles locations within the release circle         
           self.pos[:,2] = numpy.maximum.reduce([self.pos[:,2],topo[:,0] +0.1])
           self.post[:,2] = self.pos[:,2]
           print 'Updating intial particles depths within release circle based on GriddedTopo'
           if (self.pos[:,2]<topo[:,0]).any() : import pdb;pdb.set_trace() 
               
    if "polygon" in self.props:
      # release in a polygon shape
      self.polygon=Polygon(self.props['polygon'])
      point_in_poly = get_random_point_in_polygon(nbuff+1,self.polygon)
      self.pos[:,0]=point_in_poly[:,0]
      self.pos[:,1]=point_in_poly[:,1]
      self.post[:,0]=self.pos[:,0]
      self.post[:,1]=self.pos[:,1]
      # Particles release depths may need to be updated according to water depths at new locations within the polygons
      for sticker in stickers:
        if 'GriddedTopo' in sticker.__class__.__name__:
           topo=sticker.interp(self.pos,imax=1)# get depths new particles locations within the release circle         
           self.pos[:,2] = numpy.maximum.reduce([self.pos[:,2],topo[:,0] +0.1])
           self.post[:,2] = self.pos[:,2]
           print 'Updating intial particles depths within release polygon based on GriddedTopo' 
           if (self.pos[:,2]<topo[:,0]).any() : import pdb;pdb.set_trace() 
 
    
    if "variable_reln" in self.props:
      # time-varying number of particles to release - prescribed in a file
      # load file  - two columns [time(CF-compliant days since 1-1-1) nb_part_to_release]
      self.variable_reln=numpy.loadtxt(self.props['variable_reln'])


    if "variable_poly" in self.props:
      # time-varying polygons to use for particle release - prescribed in a file
      # load file - the number of columns depends on the max poly length  [time(CF-compliant days since 1-1-1) x1,y1,x2,y2 etc...]
      self.variable_poly=numpy.loadtxt(self.props['variable_poly'])


    #Could add same code for range of X,Y ?
    # e.g. if numpy.size(P0[0])==2 & numpy.size(P0[1])==2 
    # then release along a line [X1,Y1] - [X2,Y2]
    #
    # End of Release Options
    self.reln=reln #Particles per release
    self.R=R0 #Total release of material
    self.Q=Q0 # Flux of material per day
    # switch to allow the particle to unstick from sticker (0:cannot unstick/1:can unstick)
    
    if not isinstance(unstick,list):unstick=[unstick] # if only a single element input, not as a list
    if len(unstick)==1:
      self.unstick=[unstick[0] for _ in xrange(len(stickers))] # if unstick is a single element, replicate to fit number of stickers
    else:
      self.unstick=unstick

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
    """Reset and shuffle arrays after particles removed
    i0 is a boolean array which is True for active particles / False otherwise
    """
    if len(self.arrays)==0:
      for a in dir(self):
        if isinstance(getattr(self,a),numpy.ndarray):
          self.arrays.append(getattr(self,a))
    # i0 is an integer array
    if isinstance(i0,int):
      for a in self.arrays:
        a[:-i0]=a[i0:]
        a[-i0:]=a[-1]
      self.np-=i0
      return True
    else:
    # i0 is a boolean array
      nind=i0.sum()
    if nind==0:return False
    
    # update the number of active particles. nind = nb of part with state<0, which need to be removed
    self.np-=nind
    self.np=max(self.np,0) #make sure np does not become <0
    # fill locations of removed particles with new ones
    # all the random position generation is now done at release time
    for a in self.arrays:
      if len(a)!=self.npmax+1: continue # in case of variable reln or poly
      # shuffle arrays removing the dead particles 
      a[:-nind]=a[~i0]
      a[self.np:self.np+nind]=a[-1]  # backfill the locations of removed particles with new ones - here last array item


    # ####OLD BLOCK##########################################################################################
    # # define indices of initial particles position/depth to use to backfill the shuffled array
    # # TO MOVE TO THE RELEASE FUNCTION -
    # # generation of random position should happen at release time, rather than reset time >> more intuitive approach 
    # if numpy.size(self.props['P0'][2])==2 or "circular_radius" in self.props or "polygon" in self.props or hasattr(self, 'polygon'): #then release depth is random within a range
    #   # generate nind array indices, picked randomly within the range [self.np+nind:end]
    #   fill_id=numpy.random.randint(self.np+nind, len(self.pos[:,0]), nind) 
    # else:
    #   # generate nind array indices, picked within the range [self.np+nind:end]
    #   #fill_id=numpy.arange(len(self.pos[:,0])-nind+1,len(self.pos[:,0]),1,'int')
    #   fill_id=-1 # using -1 replicate the last particle position/depth of the arrays - i.e. a[-1]
    # for a in self.arrays:
    #   if len(a)!=self.npmax+1: continue # in case of variable reln or poly
    #   # shuffle arrays removing the dead particles 
    #   a[:-nind]=a[~i0]
    #   a[self.np:self.np+nind]=a[fill_id]  # fill the locations of removed particles with new ones - indices fill_id defined above depending on release type (fixed or within vertical range)    
    #   #a[self.np:self.np+nind]=a[-1]  # in ercore_nc branch
    #   #a[self.np:self.np+nind]=a[-nind-1:-1] # in master branch
    #   # these dont work properly in case of random release within a poly and/or vertical range
    # # print self.pos[self.np:self.np+nind,2]
    #  ##############################################################################################
    return True
    
  def fheader(self):
    """Return file header for output"""
    #return 'Time\tid\tx\ty\tz\tstate\tage\tmass\n'
    return 'Time\tid\tx\ty\tz\tstate\tage\tmass\tzbottom\telev\n'
   
  def sfprint(self,t):
    """Return string dump of all particles at specified timestep"""
    str=''
    for i in range(0,self.np):
      #str+="%f\t%d\t%f\t%f\t%f\t%d\t%f\t%f\n" % ((t,self.nid[i])+tuple(self.pos[i])+(self.state[i],self.age[i],self.mass[i]))
      str+="%f\t%d\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\n" % ((t,self.nid[i])+tuple(self.pos[i])+(self.state[i],self.age[i],self.mass[i],self.dep[i],self.elev[i]))
    return str
    
  def release(self,t1,t2,**k):
    
    """Release all particles between time t1 and t2
    **k can be used to pass some array information from a parent material
    e.g. in the case of a BuoyantPlume becoming a BuoyantTracer"""  

    if t2>t1:
      if (self.tstart>t2) or (self.tend<t1):return
      dt1=t2-self.tstart
      dt2=self.tend-t1
    else:
      if (t2>self.tstart) or (t1<self.tend):return
      dt1=self.tstart-t2
      dt2=t1-self.tend
    dt=min(abs(t2-t1),dt1,dt2,abs(self.tend-self.tstart))

    #release types
    #
    if k.has_key('nprel'): #Prescribed release size 
      np=k['nprel']
    elif dt==0: #Start and end time the same - release all at once
      np=k.get('nprel',self.reln)
      dt=1.
    elif hasattr(self,'variable_reln'): # number of released particle is defined from a file
      id_reln=numpy.where(numpy.abs(self.variable_reln[:,0]-t2)<=1e-6) # find correct time step
      if not id_reln[0]: 
        print 'No release number provided for t= %s - using last reln' % (t2)
        id_reln=numpy.where(numpy.abs(self.variable_reln[:,0]==self.variable_reln[-1,0]))
      np=self.variable_reln[id_reln,1]                                 # find number of particles to release at t2
      np=float(np)
    else: #Incremental release number of particle released at each time step is reln/(total_duration/timestep)
      self.relsumt+=abs(dt)
      np=k.get('nprel',abs(int(self.relsumt*self._npt)))
      if np>0:self.relsumt-=1.0*np/self._npt
    # staged release 
    if self.tstep_release>0.0:
      dt1=t2-self.tstart #time since start of model start      
      if abs(((dt1*24)/self.tstep_release)-round(((dt1*24)/self.tstep_release)))<1e-3:
        #checks if the runtime to date is a true multiple of the release time 
        #if yes release particles , if no then no release
        nb_rel=int(round(self.tstep_release/(dt*24))) 
        #nb release=nb of release that would have occurred over tstep_release if continuous, so that reln is still relevant
        #so np should be
        np=nb_rel*np
      else: #then no release
        #import pdb;pdb.set_trace()
        np=0
      #print 'Releasing %s' % (np)
      if np==0:return 0
    #

    if np==0:return
    np=int(np)
    nmax=self.npmax-self.np

    #nmax is max number of particles that can be released before reaching self.npmax
    nmax=self.npmax-self.np   
    if np>nmax:
      print 'Warning: particles exhausted for '+self.id
      np=nmax
      if np==0:return 0

    np1=self.np+np
    self.mass[self.np:np1]=self.Q*abs(dt)/np
    self.mass0=self.mass[self.np]
    self.state[self.np:np1]=1 #set particle state to 1 = released
    self.age[self.np:np1]=0.
    self.nid[self.np:np1]=range(self.ninc,self.ninc+np)
    self.ninc+=np

    # if **k is provided - replace content of self arrays with that of **k
    # this happens when the material is a child from another parent material (BuoyantPlume to BuoyantTracer) 
    for key in k:
      if hasattr(self,key):
        array=getattr(self,key) # fill property 'key' with array in **k
        if isinstance(array,numpy.ndarray):
          array[self.np:np1]=k[key][:np] # update only the position of particle "freshly" released i.e. [np:np1]
        else:
          array=k[key]

    # generate random positions - depending on release types
    # 
    # random release within a polygon
    # need to get another random position because if particles die, it recycles positions
    if hasattr(self, 'polygon'):
      # release in a polygon shape
      point_in_poly = get_random_point_in_polygon(np,self.polygon)
      self.pos[self.np:np1,0] = point_in_poly[:,0]
      self.pos[self.np:np1,1] = point_in_poly[:,1]
      self.post[self.np:np1,0]=self.pos[self.np:np1,0]
      self.post[self.np:np1,1]=self.pos[self.np:np1,1]
      # not correcting for topo here because sticker function will do that afterwards
    
    # random release within time-varying polygon
    if hasattr(self, 'variable_poly'):
      # release in a time-varying polygon shape is defined from a file 
      #e.g. variable_poly='variable_poly.txt'
      id_poly=numpy.where(numpy.abs(self.variable_poly[:,0]-t2)<=1e-6) # find correct time step
      # format the polygon coordinates into [[x1,y1],[x2,y2], ...] for shapeluy Polygon
      if not id_poly[0]: 
        print 'No polygon provided for t= %s - using last available polygon' % (t2)
        id_poly=numpy.where(numpy.abs(self.variable_poly[:,0]==self.variable_poly[-1,0]))

      poly_tmp=self.variable_poly[id_poly,1:].squeeze()
      poly=[[ poly_tmp[0] , poly_tmp[1] ]]
      for ii in range(2,poly_tmp.shape[0],2):
        poly+=[ [poly_tmp[ii],poly_tmp[ii+1]]  ]
      self.polygon=Polygon(poly)
      point_in_poly = get_random_point_in_polygon(np,self.polygon)
      self.pos[self.np:np1,0] = point_in_poly[:,0]
      self.pos[self.np:np1,1] = point_in_poly[:,1]
      self.post[self.np:np1,0]=self.pos[self.np:np1,0]
      self.post[self.np:np1,1]=self.pos[self.np:np1,1]
      # not correcting for topo here because sticker function will do that afterwards
    
    # random release within a circle
    if "circular_radius" in self.props:
      #release in a circle rather than at a single X,Y point location
      point_in_circle = get_random_point_in_circle(np,self.props['P0'],self.props['circular_radius'])
      self.pos[self.np:np1,0]=point_in_circle[:,0]
      self.pos[self.np:np1,1]=point_in_circle[:,1]
      self.post[self.np:np1,0]=self.pos[self.np:np1,0]
      self.post[self.np:np1,1]=self.pos[self.np:np1,1]     

    #random release with depth band      
    if numpy.size(self.props['P0'][2])==2:
      #thickness of z layer
      dz=abs(self.props['P0'][2][0]-self.props['P0'][2][1])
      # depth are supposed to be negative
      zz=numpy.random.random(np) 
      # generate random levels within the specified band
      self.pos[self.np:np1,2]=min(self.props['P0'][2])+dz*zz
      self.post[self.np:np1,2]=min(self.props['P0'][2])+dz*zz
    # update total number of particles
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
    # check for intersection with stickers
    for cnt,sticker in enumerate(self.stickers): # shoreline or boundary sticker
      if ('Shoreline' in  sticker.__class__.__name__) or ('Boundary' in  sticker.__class__.__name__):
        posi[:self.np,:]=sticker.intersect(self.pos[:self.np,:],posi,self.state[:self.np])
      else: # 2D sticker - GriddedTopo, GriddedElevation
        posi[:self.np,:]=sticker.intersect(self.pos[:self.np,:],posi,self.state[:self.np],t1,t2)
      # posi is the matrix of intersection positions
      # particles that intersected the sticker will be flagged with self.state==2

      # additional checks for GriddedTopo and Elevation cases
      if 'GriddedTopo' in sticker.__class__.__name__:
        self.dep[:self.np]=sticker.interp(posi[:self.np,:],imax=1)[:,0]
      if 'Elevation' in sticker.__class__.__name__:
        self.elev[:self.np]=sticker.interp(posi[:self.np,:],t2,imax=1)[:,0]

      # check is material should unstick from sticker
      if self.unstick[cnt]<=0.: # default, particle cannot unstick  > unstick= 0.0
        self.state[self.state>1]=-1 # this way particles will be removed from computation
      elif (self.state>1).any(): # particle can unstick  , unstick= 1.0
        # posi is the position of intersection with shoreline
        posi[numpy.where(self.state>1),:]=self.pos[numpy.where(self.state>1),:] # set posi back to the position particles were before sticking
        self.state[numpy.where(self.state>1)]=1 # set unstuck particles back to active state=1

      # update particle position 
      self.post[:self.np,:]=posi[:self.np,:]

    #if self.unstick<=0.:
    # self.state[self.state>1]=-1
    #self.post[:self.np,:]=posi[:self.np,:] 
  
  def die(self,t1,t2):
    """Kill particles between times t1 and t2 and remove from simulation"""
    if self.np==0:return
    #dead=(self.age>self.props.get('maxage',1.e20)) | (self.mass<=0.0001*self.mass0)
    #self.state[dead]=-2
    dead=(self.age[:self.np]>self.props.get('maxage',1.e20)) | (self.mass[:self.np]<=0.0001*self.mass0)
    self.state[:self.np][dead]=-2
    
  
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
      #import pdb;pdb.set_trace()
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
    w0: Rise velocity (m/s) [ negative w0  sinking, positive w0 > buoyant]
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
    ul#Windage with cross and down-wind ranges
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
  


