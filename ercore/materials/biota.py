#!/usr/bin/env python
import numpy
from ercore.materials import BuoyantTracer
from ercore.lib import Sun

#Class for planktonic biota
class Plankton(BuoyantTracer):
  default_props={'spwn':1,
                 'spawnclass':None,
                 'spawnage':1,
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
    spawnage: Age at which spwaning occurs (days)
    spawnclass: Class of offspring
    spawnratio: Ratio of offspring numbers/mass to parent
    mortality: Mortality rate (per day)
    vposday: Day time vertical position (m)
    vposnight: Night time vertical position (m)
    vspeed: Vertical migration rate (m/s)
    tempmin: Minimum temperature tolerated (C)
    tempmax: Maximum temperature tolerated (C)
    temptaxis: Temperature taxis rate (m/s per C/m)
    saltmin: Minimum salinity tolerated (PSU)
    saltmax: Maximum salinity tolerated (PSU)
    salttaxis: Salinity taxis rate (m/s per PSU/m)
    
    If tempmin/tempmax is specified a reactor with id starting with 'temp' must be provided
    If saltmin/saltmax is specified a reactor with id starting with 'salt' must be provided
  """
  status_codes=BuoyantTracer.status_codes.update({-3:'Transition to next life stage'})
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
    
    if 'temp' in self.reactors:
      if (self.props['tempmin'] is not None) or (self.props['tempmax'] is not None):
        temp=self.reactors['temp'].interp(self.pos[:np,:],t1)[:,0]
        if (self.props['tempmax'] is not None):
          m=(temp-self.props['tempmax']+self.props['temptol'])/self.props['temptol']
          if (m>0).any():print 'Temperature maximum reached for some %s' % (self.id)
          self.mass[:np]-=numpy.maximum(numpy.minimum(m,1),0)*self.mass[:np]
        if (self.props['tempmin'] is not None):
          m=(self.props['tempmin']+self.props['temptol']-temp)/self.props['temptol']
          self.mass[:np]-=numpy.maximum(numpy.minimum(m,1),0)*self.mass[:np]
          if (m>0).any():print 'Temperature minimum reached for some %s' % (self.id)
    if 'salt' in self.reactors:
      if (self.props['saltmin'] is not None) or (self.props['saltmax'] is not None):
        salt=self.reactors['salt'].interp(self.pos[:np,:],t1)[:,0]
        if (self.props['saltmax'] is not None):
          m=(salt-self.props['saltmax']+self.props['salttol'])/self.props['salttol']
          if (m>0).any():print 'Salinity maximum reached for some %s' % (self.id)
          self.mass[:np]-=numpy.maximum(numpy.minimum(m,1),0)*self.mass[:np]
        if (self.props['saltmin'] is not None):
          m=(self.props['saltmin']+self.props['salttol']-salt)/self.props['salttol']
          if (m>0).any():print 'Salinity minimum reached for some %s' % (self.id)
          self.mass[:np]-=numpy.maximum(numpy.minimum(m,1),0)*self.mass[:np]  
    
    #Vertical migration
    if self.props['vspeed']>0:
      sun=Sun(self.pos[:np,0],self.pos[:np,1])
      stimes=sun.getTimes(t1)
      #Go down if day time
      day=(t1>=stimes['Dawn']) & ((stimes['Dawn']>stimes['SunsetStart']) | (t1<=stimes['SunsetStart']))
      self.pos[:np,2]=numpy.where(day,numpy.maximum(self.pos[:np,2]-(t2-t1)*self.props['vspeed'],self.props['vposday']),self.pos[:np,2])
      #Go up if night time
      night=(t1>=stimes['SunsetStart']) & ((stimes['SunsetStart']>stimes['Dawn']) | (t1<=stimes['Dawn']))
      self.pos[:np,2]=numpy.where(night,numpy.minimum(self.pos[:np,2]+(t2-t1)*self.props['vspeed'],self.props['vposnight']),self.pos[:np,2])
    #Thermotaxis
    #Halotaxis
    
  
  def spawn(self,t1,t2):#Become something else or spawn
    if self.props['spawnclass']:
      ind=(self.state) & (self.age>=self.props['spawnage'])
      nind=ind.sum()
      if nind:
        if self.props['maxage']<=self.props['spawnage']:#Has become next life stage
          self.children[self.props['spawnclass']]={'pos':self.post[ind,:],'mass':self.props['spawnratio']*self.mass[ind],'nprel':nind}
          self.state[ind]=-2 
        else:#Has spawned offspring
          self.children[self.props['spawnclass']]={'pos':self.post[ind,:],'mass':self.props['spawnratio']*self.mass[ind],'nprel':nind}
        return True
      else:
        self.children={}
        return False
    
  
