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
    'D0':1.0,
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
      str+="%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % ((t,)+tuple(self.pos[i])+tuple(self.u[i])+(self.b[i],self.h[i],self.mass[i]))
    return str
  
  def __str__(self):
    str='Plume with %d elements:\n' % (self.np)
    for i in range(self.np):
      str+='%g\t%g\t%g\t%g\t%g\t%g\t%g\n' % (tuple(self.pos[i,:])+(self.vmod[i],self.b[i],self.h[i],self.conc[i]))
    return str
  
  def _get_ambient(self,t):
    np=self.np-1
    vel=self.movers[0].interp(self.post[np:np+1,:],t,3)[0]
    for mover in self.movers[1:]:
      vel+=mover.interp(self.post[np:np+1,:],t,3)[0]
    self.ambients=vel,0,0,self.props['D0']
    
  def initialize(self,t1,t2):
    _Material.initialize(self,t1,t2)
    self.V0=numpy.sqrt((numpy.array(self.props['V0'])**2).sum())
    self.dt0=10.*self.props['B0']/self.V0
    self.h=numpy.zeros(self.state.shape)
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
      U=self.ambients[0][:]
      u1=self.u[np,:]-U
      modv=(u1**2).sum()**0.5
      Ua=(U**2).sum()**0.5  
      anguu=numpy.arctan2(self.u[np,1],self.u[np,0])-numpy.arctan2(U[1],U[0])
    
      sin_theta=self.u[np,2]/self.vmod[np]
      cos_theta=numpy.sqrt(self.u[np,0]**2+self.u[np,1]**2)/self.vmod[np]
      sin_sigma=numpy.sin(anguu)
      cos_sigma=numpy.cos(anguu)
      cos_thetasig=numpy.abs(cos_sigma*cos_theta)
  
      pib2=1.414*PI*self.b[np] #2PIb/sqrt(2) 
      g1=9.81*abs(self.ddens[np])
      
      Qp=min(0.83,0.6*sin_theta*g1*self.b[np]/(modv*modv))
      Qw=abs(0.055*Ua*cos_thetasig/(modv+Ua))
      Qt=0.5*(1-cos_thetasig**2)**0.5
      return pib2*self.h[np]*(modv*(0.055+Qp+Qw)+Ua*Qt)
      
    
  def _mix(self,dt):
    Q=self._get_entrainment('jirka')
    dml=dt*self.ambients[3]*Q #Change of liquid mass
    ml1=self.mass[self.np-1]+dml #New mass
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
    h0=self.V0*self.dt
    if self.nt>self.npmax:
      print 'Warning: particles exhausted for '+self.id
      self.nt=self.npmax
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
      self.np=np
      if not ((np-1) % 1000):self._get_ambient(t1)
      dml,ml1=self._mix(self.dt)
      vstar=(self.mass[np-1]*self.u[np-1]+dml*self.ambients[0])/ml1
      if numpy.isnan(vstar).any():
        raise ERCoreException("Velocity is NAN")
      vnew=self._add_momentum(vstar,ml1,self.dt)
      self.vmod[np]=numpy.sqrt((vnew**2).sum())
      #length variation of cylindrical plume element - h
      self.h[np]=self.h[np-1]*self.vmod[np]/numpy.sqrt((self.u[np-1,:]**2).sum())     
      self.u[np,:]=vnew
      self.mass[np]=ml1
      # work out radius b variation of cylindrical plume element based on mass conservation law 
      self.b[np]=numpy.sqrt(ml1/(PI*self.dens[np]*self.h[np])) 

      self.post[np,:]=self.post[np-1,:]+self.dt*vnew*self.mfx[0]
      if self.post[np,2]>0:self.post[np,2]=0
      #Plume transition to free droplets
      self.age[np]=self.age[np-1]+self.dt/86400.

      print 'height of cylindrical plume element %s' % (self.h[np])
      print 'radius of cylindrical plume element %s' % (self.b[np])
      print 'height of cylindrical plume element %s' % (self.post[np,2])

      
      import pdb;pdb.set_trace()
 
      if self.terminate():
        print 'Plume submodel %s terminated after %d seconds' % (self.id,86400.*self.age[np])
        return
    print "Warning: plume submodel %s not terminated" % (self.id)
      
  def terminate(self):
    U=self.ambients[0]
    modv=((self.u[self.np,:]-U)**2).sum()**0.5
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
    print 'salt = ', salt
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
    Q=self._get_entrainment('jirka')
    dml=dt*D*Q #,numpy.newaxis] #Change of liquid mass
    ml1=self.mass[np-1]+dml #New mass
    self.conc[np]=(self.mass[np-1]*self.conc[np-1])/ml1 #Concentration
    self.salt[:np]=(S*dml+self.mass[:np]*self.salt[:np])/ml1 #Salinity
    self.temp[np]=(self.Cpr*T*dml+self.mass[np-1]*self.temp[np-1])/ml1 #Temperature
    self.dens[np]=self.conc[np]*self._densfunc(np)+(1.-self.conc[np])*D
    self.ddens[np]=(self.ambients[3]-self.dens[np])/self.dens[np]
    return dml,ml1
  
  def _add_momentum(self,vstar,ml1,dt):
    vstar[2]+=dt*9.81*self.ddens[self.np] #Buoyancy term
    return vstar
      
  
      
          
