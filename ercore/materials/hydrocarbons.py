#!/usr/bin/env python
from plumes import BuoyantPlume
from ercore.materials import *
import numpy

#Hydrocarbon class to hold all chemical processes
class HydroCarbon(object):
  pass
  
#Subsurface hydrocarbon plume
class HCPlume(BuoyantPlume,HydroCarbon):
  __doc__=BuoyantPlume.__doc__+"""
    db: Bubble size <float>
    wb: Vertical velocity of bubbles (fixed) <float>
    GOR: Gas/Oil ratio at outlet <float>
    Mg: Molar mass of gas  <float>
  """
  
  def initialize(self,t1,t2):
    self.hasgas=((self.props.get('db',0)>0 or self.props.get('wb',0)>0) and self.props.get('GOR',0)>0)
    J=numpy.sqrt((self.props['V0']**2).sum())*numpy.pi*self.props['B0']**2
    if self.hasgas:
      self.np=1
      self._get_ambient(t1)
      self.gdens=numpy.tile(eqnstate(101000-9.81*self.pos[0,2]*self.ambients[3],self.props['T0'],self.props['Mg']),self.npmax+1)
      fg=self.props.get('GOR')/(1.+self.props.get('GOR')) #Initial fraction of gas by volume
      self.Jg=fg*J*self.gdens[0] #Mass flux of gas
      self.db=numpy.tile(self.props.get('db',0.005),self.npmax+1) #Bubble diameter
      self.nmol=self.gdens[0]*SPI*self.db[0]**3/self.props['Mg']
      self.mq=numpy.zeros(self.npmax+1) #Mass flux of escaping bubbles/CV/s
      self.lcf=numpy.zeros(self.npmax+1) #Fractional offset of bubble plume WRT radius
      self.acf=numpy.ones(self.npmax+1) #Fractional area overlap
      self.wb=numpy.tile(self.props.get('wb',slipvel.bubble_slip(self.db[0],self.ambients[3],self.ambients[3]-self.gdens[0],self.props['visc'],self.ambients[1],self.props['IFT'])),self.npmax+1) #Fixed vertical velocity of gas bubbles
    else:
      fg=0
    self.Jl=(1-fg)*J*self.props['D0'] #Mass flux of liquid
    BuoyantPlume.initialize(self,t1,t2) #Mass variable is total mass of all states, dens is composite density.
    
      
  def _densfunc(self,np):
    if self.hasgas:
      mg=self.Jg*self.acf[np-1]
      return (mg+self.Jl)/(mg/self.gdens[np-1]+self.Jl/self.props['D0']) #Composite density
    else:
      return self.props['D0']
    
    
  def _add_momentum(self,vstar,ml1,dt):
    np=self.np
    vstar[2]+=dt*9.81*self.ddens[np]
    if self.hasgas:
      self.gdens[np]=eqnstate(101400-9.81*self.pos[np-1,2]*self.ambients[3],self.temp[np],self.props['Mg'])
      if not self.props.has_key('wb'):
        ddens=self.ambients[3]-self.gdens[np]
        self.db[np]=(self.nmol*self.props['Mg']/SPI/self.gdens[np])**0.33333
        wbnew=slipvel.bubble_slip(self.db[np],self.ambients[3],ddens,self.ambients[1],self.props['visc'],self.props['IFT'])
        dwb=(wbnew-self.wb[np-1])/dt
        self.wb[np]=wbnew
      else:
        dwb=0.
      vstar[2]-=self.nmol*self.props['Mg']*dwb/ml1 #Correction term for change of gas slip velocity
      self.mq[np]=0
      if self.props.get('gsep',True):
        if (self.acf[np-1]>0):
          #print str(self.acf[np-1])+' '+str(np)
          Ua=numpy.sqrt((self.ambients[0]**2).sum())
          Ub=self.ambients[0]
          Uup=Ub[2]+self.wb[np]
          angb=numpy.arctan2(Uup,numpy.sqrt(Ub[0]**2+Ub[1]**2))
          angp=numpy.arctan2(self.u[np-1,2],numpy.sqrt(self.u[np-1,0]**2+self.u[np-1,1]**2))
          if (self.lcf[np-1]>0) | ((self.post[np-1,2]-self.props['P0'][2]>self.M0/Ua) & (angb>angp)):
            lcfnew=min(2.,self.lcf[np-1]+dt*self.wb[np]*numpy.sin(angp)/self.b[np-1])
            Anew=max(0,(2*numpy.arccos(0.5*lcfnew)-0.5*lcfnew*numpy.sqrt(4-lcfnew**2))/numpy.pi)
            self.mq[np]=self.Jg*(self.acf[np-1]-Anew)
            self.lcf[np]=lcfnew
            self.acf[np]=Anew
          vstar-=(dt*self.mq[np]*self.u[np-1,:]/ml1)
        else:
          self.acf[np]=0.
    return vstar
  
  #Become free floating oil droplets and gas bubbles
  def spawn(self,t1,t2):
    out={}
    nsplit=10
    npar=nsplit*(max(int((t2-t1)*self.props['spawn'])/nsplit,1))
    np=self.np
    #Initialize oil droplets
    pos0=self.randcyl(slice(np,np+1),npar)
    drops=HCDroplets('droplets',npar,is3d=True,movers=self.movers,reactors=self.reactors,diffusers=self.diffusers,P0=self.pos[np,:],D0=self.props['D0'])
    drops.initialize(t1,t2)
    if self.geod:drops.geodcalc(True)
    drops.pos[:npar,:]=pos0
    #Initialize gas
    if self.hasgas:
      nn=0
      gas=HCGas('gas',npar,is3d=True,movers=self.movers,reactors=self.reactors,diffusers=self.diffusers,P0=self.pos[np,:],Mg=self.props['Mg'])
      gas.initialize(t1,t2)
      if self.geod:gas.geodcalc(True)
      if (self.acf[np]>0.):
        nn=max(int(self.acf[np]*npar),nsplit)
        pos0=self.randcyl(slice(np,np+1),nn)
      self.state[np]=-1
      if self.props.get('gsep',True):
        sind=numpy.where(self.mq[:np]>0)[0]
        nsind=len(sind)
        if nsind>0:
          np0=max(int((npar-nn)/nsind),1)
          pos0=numpy.vstack((pos0,self.randcyl(sind,np0,surface=True)))
      numpy.random.shuffle(pos0)
      nn=min(len(pos0),npar)
      gas.pos[:nn,:]=pos0[:nn,:]
      ngas=nsplit*(nn/nsplit)
      gas.mass[:nn]=self.Jg*86400.*(t2-t1)/ngas
    tt2=t1+(t2-t1)/nsplit
    #Do near field advection/diffusion
    for it in range(0,nsplit):
      drops.release(t1,tt2,nprel=npar/nsplit,age=self.age[np:np+1])
      drops.react(t1,tt2)
      drops.advect(t1,tt2)
      drops.diffuse(t1,tt2)
      drops.pos[:drops.np]=drops.post[:drops.np]
      if self.hasgas:
        gas.release(t1,tt2,nprel=ngas/nsplit,age=self.age[np:np+1])
        gas.react(t1,tt2)
        gas.advect(t1,tt2)
        gas.diffuse(t1,tt2)
        gas.die(t1,tt2)
        gas.pos[:gas.np]=gas.post[:gas.np]
    self.children['HCDroplets']={'pos':drops.pos,'mass':numpy.repeat(self.Jl*86400.*(t2-t1)/npar,npar),'nprel':npar}
    if self.hasgas:self.children['HCGas']={'pos':gas.pos,'mass':numpy.repeat(self.Jg*86400.*(t2-t1)/ngas,gas.np),'nprel':gas.np}
    return True
  
#Subsurface hydrocarbon gas 
class HCGas(BDTracer):
  __doc__=BDTracer.__doc__
  status_code=BDTracer.status_codes.update({-3:'Gas at surface'})
  def stick(self,t1,t2):
    if self.np==0:return
    ind=(self.post[:,2]>=0.)
    self.state[ind]=-3 #Gas reaches surface
    self.post[ind,2]=0.
  
#Subsurface hydrocarbon droplets
class HCDroplets(BDTracer):
  __doc__=BDTracer.__doc__+"""
    IFT: Interfacial tension of droplets <float>
    D0: Droplet density (kg/m^3) <float>
    db: Bubble diameter (m)
    Cpl: Specific heat capacity <float>
  """
  status_codes=BDTracer.status_codes.update({-3:'Droplets at surface'})
  default_props={'IFT':0.04,'temp':20,'D0':870,'Mg':400,'nmol':0,'db':0.005,'visco':0}
  def spawn(self,t1,t2):
    if self.np==0:return {}
    ind=(self.post[:,2]>=0.) & (self.state>0)
    nind=ind.sum()
    if nind:
      self.post[ind,2]=0.
      self.children['HCSlick']={'pos':self.post[ind,:],'mass':self.mass[ind],'nprel':nind}
      self.state[ind]=-3 #Droplets reach surface
      return True
    else:
      self.children={}
      return False
      
  def eqnstate(self,P,T):
    return self.props['D0']
  
#Surface hydrocarbon slick
class HCSlick(Drifter):
  __doc__=Drifter.__doc__
  status_code=Drifter.status_codes.update({-2:'Fully weathered'})
  #reactors:wind,sst,hs
  def react(self,t1,t2):
    pass #Weathering processes here

