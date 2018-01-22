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

      # update particle positions to intersection points 
      self.post[:self.np,:]=posi[:self.np,:]


    # After the sticker loop was done, check if particles should unstick from "sticker" based on its age
    # 
    # if particles are younger than the critical age - they will always be recirculated,
    # if they are older than the critical age, they will stick , and be removed from the model.
    # 
    if (self.state>1).any(): # some particles intersected land - need to check if we recirculate or remove
      # import pdb;pdb.set_trace()

      # self.pos[self.state>1] # last position before intersection
      # self.post[self.state>1] # position on land 
      # posi[self.state[:self.np]>1] # intersection point // note len(posi)=self.np, while len(self.pos) = nbuff 
      
      # identify particles that were set to be removed (i.e. state==2), but that should be recirculated i.e. age<critical age
      id_to_unstick = numpy.where( (self.state > 1) & (self.age<=self.props.get('critical_age'))  ) 
      
      # print 'ID to unstick %s ' % (id_to_unstick)
      # print 'recirculating %s particles' % (len(id_to_unstick[0]))
      # print 'ages : %s ' % (self.age[id_to_unstick[0]])
      # print 'MAX age : %s ' % (max(self.age))

      # set state of those particles we want to unstick, back to an active state=1
      self.state[id_to_unstick] = 1 
      # move those particles back to the last position before sticking (rather than intertsection points in Material class)
      self.post[id_to_unstick,:] = self.pos[id_to_unstick,:]
      # print self.post[id_to_unstick,:]
      # print self.pos[id_to_unstick,:]
      id_to_keep_stuck = numpy.where( (self.state > 1) & (self.age>self.props.get('critical_age'))  )
      self.state[id_to_keep_stuck] = -1 # this would be done anyways in the main loop , but we explicitely do it here for clarity 
      if self.state[id_to_keep_stuck].any():
        #import pdb;pdb.set_trace()
        print 'Removing %s particles with age>critical_age' % (len(id_to_keep_stuck[0]))

     
    # For the particles whose state was set to 2 , but are not younger than critical_age : no change, they will be removed

    
    # # additional checks on wetting/drying if applicable
    # if 'GriddedTopo' in self.stickers and 'GriddedElevation' in self.stickers:
    #   updated_pos_depth = check_wet_dry(self,self.pos[:self.np,:],self.post[:self.np,:],self.state,t1,t2)
    #   self.post[:self.np,:] = updated_pos_depth

  def react(self,t1,t2):
    """ Here we can add reactions with temp/salt etc.. if applicable"""
    pass
  
  def die(self,t1,t2):
    """Kill particles between times t1 and t2 and remove from simulation"""
    if self.np==0:return
    # print max( self.age[:self.np])
    # if ( self.age[:self.np] > self.props.get('maxage',1.e20) ).any() :
      # import pdb;pdb.set_trace()

    dead = ( self.age[:self.np] > self.props.get('maxage',1.e20) )
    self.state[:self.np][dead] = -2

    # Not sure this is working ....check 


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