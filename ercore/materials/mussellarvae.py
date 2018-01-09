#!/usr/bin/env python
import numpy
from ercore.materials import PassiveTracer

# The MusselLarvae class inherits from the Passive tracer base class
# sticking beahaviour with respect to the "stickers" will change depending on the 
# age of the larvae
# assuming a critical age of 20 days, and a max age of 30 days :
# age = 0-20 days - larvae is re-entrained if they contact the shoreline or seabed 
# age = 20-30 days - larvae remain in place if they contact the shoreline
# age > 30 days larvae remain in position


class MusselLarvae(PassiveTracer):
  __doc__=PassiveTracer.__doc__+"""
    w0: Rise velocity (m/s) [ negative w0  sinking, positive w0 > buoyant]
    critical_age : age at which behaviour with respect to "stickers" changes from non-sticky to non-sticky, in days
    maxage : age at which the larvae permanently sticks, in days
  """
  default_props={'w0':0.0,'critical_age' : 20.0}
  
  def initialize(self,t1,t2):
    PassiveTracer.initialize(self,t1,t2)
    self.w0=numpy.tile(self.props.get('w0',0.),self.npmax+1)
      
  def advect(self,t1,t2,order=4):
    if self.np==0:return
    PassiveTracer.advect(self,t1,t2,order)
    self.post[:self.np,2]+=86400.*(t2-t1)*self.w0[:self.np]

  def stick(self,t1,t2):
    """Do sticking for t1 to t2 - accounting for critical_age - edited from the Material class in __init__.py """
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
        # print self.dep[:self.np]
        # import pdb;pdb.set_trace()
      if 'GriddedElevation' in sticker.__class__.__name__:
        self.elev[:self.np]=sticker.interp(posi[:self.np,:],t2,imax=1)[:,0] 

      
      # check is material should unstick from sticker based on its age
      # 
      # if particles are younger than the critical age - they will always be recirculated,
      # if they are older than the critical age, they will stick , and be removed from the model.
      # import pdb;pdb.set_trace()      
      id_to_unstick = numpy.where(self.age[:self.np]<=self.props.get('critical_age'))
      self.state[id_to_unstick] = 1 # set state of particles we want to unstick, back to an active state=1

      # update particle position 
      self.post[:self.np,:]=posi[:self.np,:]

      # note : posi was already back to the position particles were before sticking so no need to re-change it

      # print 'after'
      # print numpy.min(self.pos[:self.np,2])
      # print numpy.min(self.post[:self.np,2])  
      # if numpy.min(self.pos[:self.np,2])<-1.8: import pdb;pdb.set_trace()
    
    # additional checks on wetting/drying if applicable
    if 'GriddedTopo' in self.stickers and 'GriddedElevation' in self.stickers:
      updated_pos_depth = check_wet_dry(self,self.pos[:self.np,:],self.post[:self.np,:],self.state,t1,t2)
      self.post[:self.np,:] = updated_pos_depth

  def react(self,t1,t2):
    """ Here we can add reactions with temp/salt etc.. if applicable"""
    pass
  
  def die(self,t1,t2):
    """Kill particles between times t1 and t2 and remove from simulation"""
    if self.np==0:return
    #dead=(self.age>self.props.get('maxage',1.e20)) | (self.mass<=0.0001*self.mass0)
    #self.state[dead]=-2
    dead = ( self.age[:self.np] > self.props.get('maxage',1.e20) )
    self.state[:self.np][dead]=-2

  # def diffuse(self,t1,t2):
  #   """Do diffusion for t1 to t2"""
  #   pass

  # def react(self,t1,t2):
  #   """Do reaction for t1 to t2"""
  #   pass

  # def spawn(self,t1,t2):
  #   """Do spawning for t1 to t2"""
  #   pass

  # def release(self,t1,t2,**k):
  #   """Release all particles between time t1 and t2
  #   **k can be used to pass some array information from a parent material
  #   e.g. in the case of a BuoyantPlume becoming a BuoyantTracer""" 
  #   pass 