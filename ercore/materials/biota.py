#!/usr/bin/env python
import numpy
from ercore.materials import BuoyantTracer
from ercore.lib import Sun

#Class for planktonic biota
class Plankton(BuoyantTracer):
  default_props={'spwn':1,
                 'spawnclass':None,
                 'spawnage':0,
                 'spawnratio':1,
                 'mortality':0,
                 'vposday':0,
                 'vposnight':0,
                 'vspeed':0,
                 'tempmin':None,
                 'tempmax':None,
                 'temptol':1,
                 'temptaxis':0,
                 'saltmin':None,
                 'saltmax':None,
                 'salttol':1,
                 'salttaxis':0}
  __doc__=BuoyantTracer.__doc__+"""
    spwn: Number of offspring
    spawnclass: Class of offspring
    mortality: Mortality rate
    vposday: Day time vertical position (m)
    vposnight: Night time vertical position (m)
    vspeed: Vertical migration rate (m/s)
    tempmin: Minimum temperature tolerated (C)
    tempmax: Maximum temperature tolerated (C)
    temptaxis: Temperature taxis rate (m/s per C/m)
    saltmin: Minimum salinity tolerated (PSU)
    saltmax: Maximum salinity tolerated (PSU)
    salttaxis: Salinity taxis rate (m/s per PSU/m)
  """
  
  def __init__(self,id,nbuff,**k):
    BuoyantTracer.__init__(self,id,nbuff,**k)
    self.props['vposday']=-abs(self.props['vposday'])
    self.props['vposnight']=-abs(self.props['vposnight'])
    if self.props['mortality']>0:self.props['maxage']=min(self.props['maxage'],7./self.props['mortality']) #this will remove particles where mass<1/1000th of original value
    
  
  def react(self,t1,t2):
    np=self.np
    if np==0:return
    #Mortality - note killing after maxage already handled in material base class
    if self.props['mortality']>0:self.mass[:np]-=self.props['mortality']*(t2-t1)*self.mass[:np]
    
    if (self.props['tempmin'] is not None) or (self.props['tempmax'] is not None):
      temp=self.reactors[0].interp(self.pos[:np,:],t1)[:,0]
      if (self.props['tempmax'] is not None):
        m=(temp-self.props['tempmax']+self.props['temptol'])/self.props['temptol']
        self.mass[:np]-=numpy.maximum(numpy.minimum(m,1),0)*self.mass[:np]
      if (self.props['tempmin'] is not None):
        m=(self.props['tempmin']+self.props['temptol']-temp)/self.props['temptol']
        self.mass[:np]-=numpy.maximum(numpy.minimum(m,1),0)*self.mass[:np]
    if (self.props['saltmin'] is not None) or (self.props['saltmax'] is not None):
      salt=self.reactors[1].interp(self.pos[:np,:],t1)[:,0]
      if (self.props['saltmax'] is not None):
        m=(salt-self.props['saltmax']+self.props['salttol'])/self.props['salttol']
        self.mass[:np]-=numpy.maximum(numpy.minimum(m,1),0)*self.mass[:np]
      if (self.props['saltmin'] is not None):
        m=(self.props['saltmin']+self.props['salttol']-salt)/self.props['salttol']
        self.mass[:np]-=numpy.maximum(numpy.minimum(m,1),0)*self.mass[:np]  
    
    #Vertical migration
    if self.props['vspeed']>0:
      sun=Sun(self.pos[:np,0],self.pos[:np,1])
      stimes=sun.getTimes(t1)
      #Go down if day time
      if (t1>=stimes['Dawn']) and ((stime['Dawn']>stime['SunsetStart']) or (t1<=stime['SunsetStart'])):
        self.pos[:np,2]-=(t2-t1)*self.props['vspeed']
        self.pos[:np,2][self.pos[:np,2]<self.props['vposday']]=self.props['vposday']
      #Go up if night time
      if (t1>=stimes['Sunset Start']) and ((stime['SunsetStart']>stime['Dawn']) or (t1<=stime['Dawn'])):
        self.pos[:np,2]+=(t2-t1)*self.props['vspeed']
        self.pos[:np,2][self.pos[:np,2]>self.props['vposnight']]=self.props['vposnight']
    #Thermotaxis
    #Halotaxis
    
  
  def spawn(self,t1,t2):#Become something else or spawn
    if self.props['spawnclass']:
      ind=(self.age>=self.props['spawnage'])
      nind=ind.sum()
      if nind:
        self.children[self.props['spawnclass']]={'pos':self.post[ind,:],'mass':self.props['spawnratio']*self.mass[ind],'nprel':nind*self.props['spawnratio']}
        if self.props['maxage']<=self.props['spawnage']:
          self.state[ind]=-2 #Has become next life stage
        return True
      else:
        self.children={}
        return False
    
  
