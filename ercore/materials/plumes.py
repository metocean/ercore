#!/usr/bin/env python
from ercore.lib import pres,temppot,dens
from ercore.materials import _Material
from ercore._flib_ercore import plume
from ercore import ERCoreException
import numpy,sys

plume.ES=2.0
plume.EF=0.0
PI=numpy.pi
PI2=2*PI

#Core plume class
#This contains an internal lagrangian plume model (a la CorMix)
class Plume(_Material):
  default_props={
    'B0':1.0,
    'V0':[0,0,1],
    'Vb':0.12,
    'D0':1027.0,
    'C0':1.0,
    'E':2,
  }
  __doc__=_Material.__doc__+"""
    B0: Initial jet radius <float>
    V0: Initial jet velocity <list> of u,v,w components
    Vb: Bubble/dropplet terminal velocity (m/s) <float>
    D0: Jet intitial density (kg/m^3) <float>
    C0: Jet initial concentration <float>
    E: Entrainment constant <float>
  """
  
  def fheader(self):
    return "Time\tx\ty\tz\tu\tv\tw\tb\th\tmass\n"
  
  def sfprint(self,t,deadonly=False):
    str=''
    for i in range(0,self.np):
      str+="%f\t%.10f\t%.10f\t%.10f\t%f\t%f\t%f\t%f\t%f\t%f\n" % ((t,)+tuple(self.post[i])+tuple(self.u[i])+(self.b[i],self.h[i],self.mass[i]))
      # using %.10f to correctly output small-scale advection of plumes geographical coordinates
    return str
  
  def __str__(self):
    str='Plume with %d elements:\n' % (self.np)
    for i in range(self.np):
      str+='%g\t%g\t%g\t%g\t%g\t%g\t%g\n' % (tuple(self.post[i,:])+(self.vmod[i],self.b[i],self.h[i],self.conc[i]))
    return str
  
  def _get_ambient(self,t):
    np=self.np-1
    # get ambient current velocity
    vel=self.movers[0].interp(self.post[np:np+1,:],t,3)[0]
    for mover in self.movers[1:]:
      vel+=mover.interp(self.post[np:np+1,:],t,3)[0]
    self.ambients=vel,0,0,self.props['D0']
    
  def initialize(self,t1,t2):
    _Material.initialize(self,t1,t2)
    self.V0=numpy.sqrt((numpy.array(self.props['V0'])**2).sum()) # jet velocity magnitude
    self.dt0=10.*self.props['B0']/self.V0 # criteria on time step to use
    self.h=numpy.zeros(self.state.shape) # 
    self.b=self.props['B0']+self.h
    self.u=self.props['V0']*numpy.ones((len(self.state),1))
    self.vmod=numpy.sqrt((self.u**2).sum(1))
    self.ddens=0.+self.h
    self.conc=self.props['C0']+self.h
    self.M0=self.V0*numpy.array(self.props['B0'])**2
    self.Vb=self.props.get('Vb',0.12) #Terminal velocity
      
  def _get_entrainment(self,type='jirka'):
    #Entrainment calculation - Yapa and Zheng, 1997
    if False:
      if type=='jirka':
        return plume.entrainmentJ(self.u,self.ambients[0][:self.np,:],self.b,self.h,self.ddens,self.np)
      else:
        return plume.entrainmentF(self.u,self.ambients[0][:self.np,:],self.b,self.h,self.ddens,self.np)
    else:
      np=self.np-1
      U=self.ambients[0][:]   # ambient current velocity vector
      u1=self.u[np,:]-U       # relative jet current velocity vector
      modv=(u1**2).sum()**0.5 # magnitude of relative jet current velocity
      Ua=(U**2).sum()**0.5    # ambient current velocity vector
      # angle between flow and jet - we take the flow axis as x-axis i.e. relative angle
      anguu=numpy.arctan2(self.u[np,1],self.u[np,0])-numpy.arctan2(U[1],U[0]) 
      # angles :
      # theta is the angle between the x-axis and jet in the x-z plane (i.e. horizontal plane) [deg]
      # sigma is the angle between the x-axis and jet in the x-z plane (i.e vertical plane) [deg]
      sin_theta=self.u[np,2]/self.vmod[np] 
      cos_theta=numpy.sqrt(self.u[np,0]**2+self.u[np,1]**2)/self.vmod[np]
      sin_sigma=numpy.sin(anguu)
      cos_sigma=numpy.cos(anguu)
      cos_thetasig=numpy.abs(cos_sigma*cos_theta)
  
      pib2=1.414*PI*self.b[np] #2PIb/sqrt(2) 
      g1=9.81*abs(self.ddens[np])
      
      Qp=min(0.83,0.6*sin_theta*g1*self.b[np]/(modv*modv)) # shear-related entrainment
      Qw=abs(0.055*Ua*cos_thetasig/(modv+Ua))              # shear-related entrainment
      Qt=0.5*(1-cos_thetasig**2)**0.5                      # drag-related entrainment
      return pib2*self.h[np]*(modv*(0.055+Qp+Qw)+Ua*Qt)    # dV : volume of water entrained in plume
      
    
  def _mix(self,dt):
    Q=self._get_entrainment('jirka') # dV : volume of water entrained in plume
    dml=dt*self.ambients[3]*Q #dM : mass of water entrained in plume (self.ambients[3] is ambient water density)
    ml1=self.mass[self.np-1]+dml #New mass at t+1
    return dml,ml1
  
  def _add_momentum(self,vnew,dt): #Hook for additional momentum terms
    return vnew
    
  def release(self,t1,t2):
    self.np=0
    self.children={}
    if t2<self.tstart or t1>self.tend:return 0
    nt=numpy.ceil(86400*(t2-t1)/self.dt0)
    self.dt=86400.*(t2-t1)/nt
    self.nt=int(nt)
    if self.nt>self.npmax:
      print 'Warning: particles exhausted for '+self.id
      self.nt=self.npmax
    h0=self.V0*self.dt
    self.h[0]=h0
    self.mass[0]=PI*h0*self.b[0]**2*self.dens[0]
    self.state[:nt]=1
    self.np=nt
    self.post[:]=self.pos[0]
    self.age[:]=0.
    return nt
    
  #This is currently set up so that the plume should be able to evolve and terminate inside a single master timestep
  def advect(self,t1,t2,order=4):
    for np in range(1,self.nt+1):
      self.np=np # np is t+1, np-1 is t
      if not ((np-1) % 1000):self._get_ambient(t1)
      dml,ml1=self._mix(self.dt) # change of mass of plume element, and new mass at t+1
      vstar=(self.mass[np-1]*self.u[np-1]+dml*self.ambients[0])/ml1 # jet velocity vector at t+1
      if numpy.isnan(vstar).any():
        raise ERCoreException("Velocity is NAN")
      vnew=self._add_momentum(vstar,ml1,self.dt) # add buoyancy term in the vertical velocity component of jet element
      self.u[np,:]=vnew                          # jet velocity vector at t+1
      self.vmod[np]=numpy.sqrt((vnew**2).sum())  # jet velocity magnitude at t+1
      self.h[np]=self.h[np-1]*self.vmod[np]/self.vmod[np-1]  #length variation of cylindrical plume element - h   
      self.mass[np]=ml1                                      # new mass of cylindrical plume element at t+1 
      self.b[np]=numpy.sqrt(ml1/(PI*self.dens[np]*self.h[np]))  #radius b of cylindrical plume element at t+1 (using mass conservation)
      self.post[np,:]=self.post[np-1,:]+self.dt*vnew*self.mfx[0] # update position of center of plume element at t+1
      # at this stage particles positions is the center of the plume elements as it mixes with ambient water and current 
      if self.post[np,2]>0:self.post[np,2]=0 #Plume transition to free droplets at sea surface
      self.age[np]=self.age[np-1]+self.dt/86400.

      # print '[np,x,y,z,ujet,vjet,wjet,conc,dens,temp,salt] = [ %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s]' % (np,self.post[np,0],self.post[np,1],self.post[np,2],vnew[0],vnew[1],vnew[2],self.conc[np],self.dens[np],self.temp[np],self.salt[np])
      # print '[np,mass] = [ %s,%s]' % (np,self.mass[np])
      # import pdb;pdb.set_trace()
 
      if self.terminate():
        print 'Plume submodel %s terminated after %d seconds' % (self.id,86400.*self.age[np])
        return
    print "Warning: plume submodel %s not terminated - increase master time step" % (self.id)

  def terminate(self):
    # check if we can terminate the plume submodel
    U=self.ambients[0]
    modv=((self.u[self.np,:]-U)**2).sum()**0.5
    # terminate if the jet velocity becomes smaller than a given terminal velocity Vb (set to 0.12 by default),
    # or the plume element reaches the surface or the seabed
    return (modv<=self.Vb) or (self.post[self.np,2]>-0.01)
      
  def randcyl(self,ind,np,surface=False):
    u=self.u[ind,:]
    n=u.shape[0]
    if surface:
      r0=numpy.tile(self.b[ind][:,None],(1,np))
      ang0=numpy.pi*(numpy.random.random((n,np))+0.5)+numpy.tile(numpy.arctan2(u[:,1:2],u[:,0:1]),(1,np))
    else:
      r0=numpy.sqrt(numpy.random.random((n,np))*self.b[ind][:,None]**2)
      ang0=PI2*numpy.random.random((n,np))
    x0=r0*numpy.cos(ang0)
    y0=r0*numpy.sin(ang0)
    z=self.h[ind][:,None]*(numpy.random.random((n,np))-0.5)
    #x-rotation
    A=(u[:,1]**2+u[:,2]**2)**0.5
    sin=-u[:,1]/A
    cos=u[:,2]/A
    y=y0*cos[:,None]-z*sin[:,None]
    z0=y0*sin[:,None]+z*cos[:,None]
    #y-rotation
    A=(u[:,0]**2+u[:,2]**2)**0.5
    sin=u[:,0]/A
    cos=u[:,2]/A
    x=x0*cos[:,None]+z0*sin[:,None]
    z=z0*cos[:,None]-x0*sin[:,None]
    return numpy.tile(self.pos[ind,:],(np,1))+numpy.vstack((x.T.flatten()*self.mfx[0,0],y.T.flatten()*self.mfx[0,0],z.T.flatten())).T
      
#Generic buoyant plume class
class BuoyantPlume(Plume):
  __doc__=Plume.__doc__+"""
    T0: Initial temperature (C) <float>
    S0: Initial salinity (PSU) <float>
    D0: Initial density (kg/m^3) <float>
    Cpl: Specific heat capacity <float>
  """
  def initialize(self,t,t2):
    Plume.initialize(self,t,t2)
    self.temp=self.h+self.props['T0']+273. #Initial temperature
    self.salt=self.h+self.props.get('S0',0.) #Initial salinity
    self.conc=self.h+self.props.get('C0',1.) #Initial concentration
    P=-pres(self.pos[:,2],self.pos[:,1] if self.geod else 0.)
    D=dens(self.salt, self.temp, P)
    self.dens=self.conc*self._densfunc(0)+(1.-self.conc)*D
    self.ddens=self.h+0.
    self.Cpr=self.props.get('Cpa',3.93)/self.props['Cpl']
  
  def _get_ambient(self,t):
    Plume._get_ambient(self,t)
    np=self.np-1 
    temp=self.reactors[0].interp(self.post[np:np+1,:],t)[0]  # outputs array with imax length
    salt=self.reactors[1].interp(self.post[np:np+1,:],t)[0]
    P=pres(self.post[np,2],self.post[np:np+1,1] if self.geod else 0.)
    temp=temppot(salt, temp, 0, P) #Convert to absolute temperature
    den=dens(salt, temp, P)
    # Note from Rosa: added a fix for when imax>0. not tested in all conditions.
    # print 'salt = ', salt
    temp = temp[0]
    salt = salt[0]
    den  = den[0]
    self.ambients=self.ambients[0],temp,salt,den
    #print str(self.post[np,2])+str(self.ambients[0])
  
  def _densfunc(self,np):
    return self.props['D0'] # if np==0 else numpy.tile(self.props['D0'],np)
    
  def _mix(self,dt):
    np=self.np
    V,T,S,D=self.ambients
    Q=self._get_entrainment('jirka') #dV : volume of water entrained in plume element
    dml=dt*D*Q #,numpy.newaxis] #dM : mass of water entrained in plume element (D is ambient water density)
    ml1=self.mass[np-1]+dml #New mass of plume element
    self.conc[np]=(self.mass[np-1]*self.conc[np-1])/ml1                 #Concentration at t+1
    self.salt[:np]=(S*dml+self.mass[:np]*self.salt[:np])/ml1            #Salinity at t+1
    self.temp[np]=(self.Cpr*T*dml+self.mass[np-1]*self.temp[np-1])/ml1  #Temperature at t+1
    self.dens[np]=self.conc[np]*self._densfunc(np)+(1.-self.conc[np])*D #Density at t+1
    self.ddens[np]=(self.ambients[3]-self.dens[np])/self.dens[np]       #relative density at t+1
    return dml,ml1
  
  def _add_momentum(self,vstar,ml1,dt):
    vstar[2]+=dt*9.81*self.ddens[self.np] #Buoyancy term
    return vstar
      
  # Buoyant plume class
  # based on Lee, J.H.W. and Cheung, V. (1990) Generalized Lagrangian model for buoyant jets in current. 
  # Journal of Environmental Engineering, ASCE,116(6), pp. 1085-1105.
  # can probably merge with Buoyant Plume at some stage but keeping separate to not break previous tests
class BuoyantPlume_JETLAG(Plume):
  default_props={
    'B0':1.0,
    'V0':[0,0,1],
    'Vb':0.05,
    'D0':1.0,
    'C0':1.0,
    'E':2,
  }
  __doc__=Plume.__doc__+"""
    From Plume subclass
    B0: Initial jet radius / Port diameter <float>
    V0: Initial jet velocity <list> of u,v,w components
    D0: Jet intitial density (kg/m^3) <float>
    C0: Jet initial concentration <float>
    **not used
    **Vb: Bubble/dropplet terminal velocity (m/s) <float>
    **E: Entrainment constant <float>
    - specific to BuoyantPlume_JETLAG
    T0: Initial jet temperature (C) <float>
    S0: Initial jet salinity (PSU) <float>
    **Cpl: Specific heat capacity <float>
  """
  def initialize(self,t1,t2):
    _Material.initialize(self,t1,t2)
    self.V0=numpy.sqrt((numpy.array(self.props['V0'])**2).sum()) # jet velocity magnitude
    self.dt0=1.0*self.props['B0']/self.V0                        # criteria on time step to use as recommend in Lee Cheung 8d.
    self.b=self.props['B0'] *numpy.ones((len(self.state),1))     # initial radius of cylindrical plume element 8b
    self.h=self.props['B0'] *numpy.ones((len(self.state),1))     # initial height of cylindrical plume element
    self.u=self.props['V0'] *numpy.ones((len(self.state),1))     # jet velocity vector
    self.vmod=numpy.sqrt((self.u**2).sum(1)) # jet velocity magnitude
    self.conc=self.props['C0']*numpy.ones((len(self.state),1))
    self.dens=self.props['D0']*numpy.ones((len(self.state),1))      # initial jet density
    self.ddens=0.*numpy.ones((len(self.state),1)) # density difference
    # P=-pres(self.pos[:,2],self.pos[:,1] if self.geod else 0.)
    # D=dens(self.salt, self.temp, P)
    # self.dens=self.conc*self._densfunc(0)+(1.-self.conc)*D  
    self.temp=self.props['T0']  * numpy.ones((len(self.state),1))#+273. #Initial temperature deg
    self.salt=self.props.get('S0',0.) * numpy.ones((len(self.state),1)) #Initial salinity ppt
    # self.Cpr=self.props.get('Cpa',3.93)/self.props['Cpl']
    # initial dimensions of plume element > done in release function
    h0=self.V0*self.dt0
    self.h[0]=h0 #
    self.mass[0]=PI*h0*self.b[0]**2*self.dens[0]
    self.mass0=self.mass[0]
    #jet angle quantities
    self.sin_phi=0.*numpy.ones((len(self.state),1)) 
    self.cos_phi=0.*numpy.ones((len(self.state),1)) 
    self.sin_theta=0.*numpy.ones((len(self.state),1)) 
    self.cos_theta=0.*numpy.ones((len(self.state),1)) 
    self.M0=self.V0*numpy.array(self.props['B0'])**2
    self.Vb=self.props.get('Vb',0.12) #Terminal velocity
    
  def release(self,t1,t2):
    # define the timestep of plume submodel and release initial plume elements
    self.np=0
    self.children={}
    if t2<self.tstart or t1>self.tend:return 0
    nt=numpy.ceil(86400*(t2-t1)/self.dt0) # number of time steps for plume model during 1 master dt 
    # import pdb;pdb.set_trace()
    self.dt=86400.*(t2-t1)/nt 
    self.nt=int(nt)
    if self.nt>self.npmax:  
      # the number of plume timesteps is larger than the nbuff input for Material
      # the plume model will often end before reaching the self.nt-th timestep 
      print 'Warning: particles exhausted for '+self.id
      self.nt=self.npmax # reduce the number of plume model timesteps to npmax
    h0=self.V0*self.dt
    self.h[0]=h0 
    self.mass[0]=PI*h0*self.b[0]**2*self.dens[0]
    self.state[:nt]=1
    self.np=nt
    self.post[:]=self.pos[0]
    self.age[:]=0.
    return nt

  def _get_ambient(self,t):
    
    Plume._get_ambient(self,t)  #this yields  [vel,0,0,self.props['D0']]
    np=self.np-1
    # not super generic - this assumes that first reactor is temp, second is salt....
    temp=self.reactors[0].interp(self.post[np:np+1,:],t)[0]  # outputs array with imax length
    salt=self.reactors[1].interp(self.post[np:np+1,:],t)[0]
    # compute pressure, potential temperature, and density based on ambient conditions
    # using UNESCO 1983 algorithms
    P=pres(self.post[np,2],self.post[np:np+1,1] if self.geod else 0.) 
    temp=temppot(salt, temp, 0, P) #Convert to absolute temperature
    den=dens(salt, temp, P)
    # Note from Rosa: added a fix for when imax>0. not tested in all conditions.
    #print 'salt = ', salt
    temp = temp[0]
    salt = salt[0]
    den  = den[0]
    self.ambients=self.ambients[0],temp,salt,den
    #[U-vectori T S D]
    #print str(self.post[np,2])+str(self.ambients[0])

  def _get_entrainment(self,type='jetlag'):
    #Entrainment calculation - Yapa and Zheng, 1997
    if type=='jetlag':
      
      np=self.np-1 # time t
      U=self.ambients[0][:]   # ambient current velocity vector
      # u1=self.u[np,:]-U       # relative jet current velocity vector
      # modv=(u1**2).sum()**0.5 # magnitude of relative jet current velocity
      Ua=(U**2).sum()**0.5    # ambient current velocity magnitude
      # angle between flow and jet - we take the flow axis as x-axis i.e. relative angle
      anguu=numpy.arctan2(self.u[np,1],self.u[np,0])-numpy.arctan2(U[1],U[0]) 
      # jet angles : using sames names as in Lee Cheung for consistency :
      # theta is the angle between the x-axis and jet in the x-y plane (i.e. horizontal plane) [deg]
      # phi is the angle between the x-axis and jet in the x-z plane (i.e vertical plane) [deg]
      self.sin_phi[np]=self.u[np,2]/self.vmod[np] 
      self.cos_phi[np]=numpy.sqrt(self.u[np,0]**2+self.u[np,1]**2)/self.vmod[np]
      self.sin_theta[np]=numpy.sin(anguu) # here theta is taken as the angle between flow axis and jet axis
      self.cos_theta[np]=numpy.cos(anguu) # i.e. this means we use the flow axis as reference axis
      # shear-entrainment e.q 10-11
      vjet_mod=self.vmod[np]-Ua*self.cos_phi[np]*self.cos_theta[np]  # relative jet velocity 
      #Note Ua*cos_thetasig is the ambient current component projected onto the jet axis
      alpha=1.0 #proportionality constant - 1.0 used by spearman in TASS for example
      rho_a=self.ambients[3] #ambient seawater density
      d_rho=rho_a-self.dens[np]    
      froude=numpy.abs( alpha*vjet_mod/( (9.81*(d_rho/rho_a)*self.b[np])**0.5) ) # Densimetric Froude Number
      #entrainment coefficient
      E=numpy.sqrt(2.)* ( 0.057+(0.554*self.sin_phi[np]/(froude**2)) ) / ( 1 + 5*Ua*self.cos_phi[np]*self.cos_theta[np] / vjet_mod ) 
      # >> this will be done in the mix function :
      # dM_shear=rho_a.*2*pi.*self.b*self.h*E*vjet_mod*dt;
      dV_shear=2.*numpy.pi*self.b[np]*self.h[np]*E*vjet_mod # volume of water entrained in plume element due to shear entrainment
      # forced entrainment eq. 20
      mf1=2.*numpy.sqrt(self.sin_phi[np]**2 + self.sin_theta[np]**2 - (self.sin_phi[np]*self.sin_theta[np])**2 )  #entrainment due to the projected plume area normal to the cross flow
      if np>1:
        db=self.b[np]-self.b[np-1]
        if db<0: import pdb;pdb.set_trace() # 
        ds=self.vmod[np]*self.dt0# jet arc length over dt 
        mf2=numpy.pi* (db/ds) * self.cos_phi[np] * self.cos_theta[np] #correction due to the growth of plume radius
        mf3= (numpy.pi/2) * self.b[np] * ( self.cos_phi[np]*self.cos_theta[np] - self.cos_phi[np-1]*self.cos_theta[np-1] ) / ds #correction due to the curvature of the trajectory.
      else:
        mf2=0
        mf3=0
      #dM_forced(k)=rho_a.*Ua.*hk(k).*bk(k).*(mf1(k) + mf2(k) +mf3(k)) .* dt(k);
      dV_forced=Ua * self.h[np] * self.b[np] * (mf1 + mf2 +mf3) # volume of water entrained in plume element due to forced entrainment

      return numpy.maximum(dV_shear,dV_forced)    # eq 21

  def _densfunc(self,np):
    return self.props['D0'] # if np==0 else numpy.tile(self.props['D0'],np)
  
  def advect(self,t1,t2,order=4):
    Plume.advect(self,t1,t2,order=4)
    # For now this advect particles following the center of the plume element ([x,y,z])
    # we could consider releasing within each of the plume element rathen than at its center
    # or we could release only at the end of the plume model i.e. when we go into far field
  
  def spawn(self,t1,t2):
    # the computed nearfield plume dynamics are used to seed the model with particles for far-field dispersion

    # this is handled using the "spawn" function and the creation of a "child" of the BuoyantPlume_JETLAG Material
    # the child name must be specified in BuoyantPlume_JETLAG function call spawn_class='xxxx'
  
    self.state[:]=-2 # inactivate plume particles
    # for now release particles at last point of plume
    # child material release position
    #for now release at the center of the last plume element
    # could be along the plume path ? or within last plume element ? etc...
    final_pos=numpy.tile(self.post[self.np-1,:], (int(self.npmax+1),1) ) # self.post[self.np-1,:] is position of the plume at end of nearfield dynamics
    # set the mass as the last concentration computed in the nearfield plume subplume
    # can be used on to infer dilution after nearfield dynamics
    final_conc=numpy.tile(self.conc[self.np-1], (int(self.npmax+1),1) ) 
    # update positions and masses of the plume's child class
    # this should not overwrite positions of particle that are already suspended - ok : see release function   
    self.children[self.props['spawn_class']]={'pos':final_pos,'post':final_pos,'mass':final_conc }
    

  def terminate(self):
    # check if we can terminate the plume submodel
    U=self.ambients[0]
    modv=((self.u[self.np,:]-U)**2).sum()**0.5
    # terminate if the jet velocity becomes smaller than a given terminal velocity Vb (set to 0.12 by default),
    # or the plume element reaches the surface or the seabed
    water_depth=self._get_depth()
    # print water_depth
    # import pdb; pdb.set_trace()
    # print 'pos %.10f %.10f %.1f' % (self.post[self.np,0],self.post[self.np,1],self.post[self.np,2])
    return (modv<=self.Vb) or (self.post[self.np,2]>-0.01) or (self.post[self.np,2]<water_depth) 
    
    # >> add case for surface release i.e. check on seabed impact -     water_depth= look if there is a GriddedTopo field ??
    # >> look how the release is actually handled - released in last plume element ? or within each element through nearfield dilution? 

  def _mix(self,dt):
    np=self.np
    V,T,S,D=self.ambients
    Q=self._get_entrainment('jetlag') #dV : volume of water entrained in plume element
    dml=dt*D*Q #,numpy.newaxis] #dM : convert to mass of water entrained in plume element over dt (D is ambient water density)
    ml1=self.mass[np-1]+dml #New mass of plume element
    self.conc[np]=(self.mass[np-1]*self.conc[np-1] )/ml1                 #Concentration at t+1 eq. 2d (this assumes ambient conc Ca=0)
    self.salt[np]=(self.mass[np-1]*self.salt[np-1] + S*dml)/ml1           #Salinity at t+1 eq. 2a
    self.temp[np]=(self.mass[np-1]*self.temp[np-1] + T*dml)/ml1          #Temperature at t+1 eq. 2b
    self.dens[np]=self.conc[np]*self._densfunc(np)+(1.-self.conc[np])*D  #Density at t+1 , alternatively could use UNESCO function dens=f(Tj,Sj) 
    self.ddens[np]=(D-self.dens[np])/D                                   #relative density at t+1 (rho_a-rho_jet)/rho_a
    return dml,ml1
  
  def _add_momentum(self,vstar,ml1,dt):
    vstar[2]+=dt*9.81*self.ddens[self.np] #Buoyancy term
    return vstar

  def _get_depth(self):
    # get depth at plume element position
    np=self.np
    for mover in self.movers[0:]:
      if mover.topo:
        depth= mover.topo.interp(numpy.array([self.post[np,:]]) ,None,imax=3)[:,0]
        return depth