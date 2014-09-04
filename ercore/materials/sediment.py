#!/usr/bin/env python
from ercore.material import *
from ercore.lib import dens

#Class for sediment
#Suspended sediment is just a (negatively) buoyant tracer
#Resuspension and bedload could easily be incorporated - see my Master thesis ("A particle tracking model for sediment transport in the Nearshore Zone", Johnson 2000)
class Sediment(BuoyantTracer):
  __doc__=PassiveTracer.__doc__+"""
    d0:   Sediment grain size (mm) - If specified will over-ride w0
    rho_s:Sediment density (kg/m^3) 
  """
  default_props={'d0':None,'rho_s':2650}
  def initialize(self,t1,t2):
    BuoyantTracer.initialize(self,t1,t2)
    if self.props.get('d0'):
      w0=0.0005*(props.get('rho_s')-dens(35,20))*props.get('d0')**2 #This is calculate for sea water with 35 PPT and 20C
      self.w0=numpy.tile(w0,self.npmax+1)
      
  def react(self,t1,t2):
    pass
    #Could more sophisticated calculations of fall velocity based on temp/salinity
    #Also resuspension could go in here dependant on waves and currents
    
  def diffuse(self,t1,t2):
    BuoyantTracer.diffuse(self,t1,t2)
    #Heavy particle diffusion mods could go in here
  
  
