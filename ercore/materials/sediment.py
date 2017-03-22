#!/usr/bin/env python
from ercore.materials import *
from ercore.lib import dens

#Class for sediment
#Suspended sediment is just a (negatively) buoyant tracer
#Resuspension and bedload could easily be incorporated - see my Master thesis ("A particle tracking model for sediment transport in the Nearshore Zone", Johnson 2000)
#Could add more sophisticated calculations of fall velocity based on temp/salinity


#from  https://github.com/csherwood-usgs/crspy/blob/master/crspy.py
def qkhfs( w, h ):
    """
    Quick iterative calculation of kh in gravity-wave dispersion relationship
    kh = qkhfs(w, h )
    
    Input
        w - angular wave frequency = 2*pi/T where T = wave period [1/s]
        h - water depth [m]
    Returns
        kh - wavenumber * depth [ ]
    Orbital velocities from kh are accurate to 3e-12 !
    RL Soulsby (2006) "Simplified calculation of wave orbital velocities"
    HR Wallingford Report TR 155, February 2006
    Eqns. 12a - 14
    """
    g = 9.81
    x = w**2.0 *h/g
    y = numpy.sqrt(x) * (x<1.) + x *(x>=1.)
    # is this faster than a loop?
    t = numpy.tanh( y )
    y = y-( (y*t -x)/(t+y*(1.0-t**2.0)))
    t = numpy.tanh( y )
    y = y-( (y*t -x)/(t+y*(1.0-t**2.0)))
    t = numpy.tanh( y )
    y = y-( (y*t -x)/(t+y*(1.0-t**2.0)))
    kh = y
    return kh

class Sediment(BuoyantTracer):
  __doc__=PassiveTracer.__doc__+"""
    d0:   Sediment grain size (mm) - If specified will over-ride w0
    rho_s:Sediment density (kg/m^3)
    tau_crit_eros: critical bed shear stress for erosion in N/m2 (default=0.2)
    tau_crit_deposition: critical bed shear stress for deposition in N/m2 (default=1000 - deposition always possible)
  """
  default_props={'d0':None,'rho_s':2650,'tau_crit_eros':0.2,'tau_crit_depos':1000.0}
  def initialize(self,t1,t2):
    #import pdb;pdb.set_trace()
    BuoyantTracer.initialize(self,t1,t2)
    if self.props.get('d0'):
      w0=0.0005*(props.get('rho_s')-dens(35,20))*props.get('d0')**2 #This is calculate for sea water with 35 PPT and 20C
      self.w0=numpy.tile(w0,self.npmax+1)
    #new array to facilitate bookkeeping of successive particle deposition/erosion  
    self.on_seabed=numpy.tile(0,self.npmax+1)
    #self.on_seabed is true when a particle is deposited on the seabed, false otherwise
    self.seabedlayer_thick=0.2 # thickness of seabed layer used for resuspension
  
  def react(self,t1,t2):
    pass
   
  def stick(self,t1,t2):
    """Do sticking for t1 to t2"""
    # in the main RUN call - "react" happens before "stick" so the check on re-suspension needs to happen here
    if self.np<1:return
    np=self.np
    posi=numpy.where(self.state[:np,None]<0,self.pos[:np,:],self.post[:np,:])
    #import pdb;pdb.set_trace()
    for sticker in self.stickers: 
      # check for intersection with stickers
      if 'Shoreline' in  sticker.__class__.__name__: # shoreline sticker 
        # Returns matrix of intersection positions, how does it affect state ? 
        posi[:self.np,:]=sticker.intersect(self.pos[:self.np,:],posi,self.state[:self.np])
        if self.unstick<=0.: # by default unstick is 0.0
          self.state[self.state>1]=-1 # this way particles will be removed from computation
        else:
          self.state[self.state>1]=1 # this way particles will be set back to active
          # posi is the position of intersection with shoreline
          posi[self.state>1,:]=self.pos[self.state>1,:] # set posi back to the position particles were before sticking
        self.post[:self.np,:]=posi[:self.np,:]
      else: # 2D sticker (e.g. GriddedTopo)
        posi[:self.np,:]=sticker.intersect(self.pos[:self.np,:],posi,self.state[:self.np],t1,t2)

      # additional checks for GriddedTopo and Elevation cases
      if 'GriddedTopo' in sticker.__class__.__name__:
        self.dep[:self.np]=sticker.interp(posi[:self.np,:],imax=1)[:,0] # get depths at particles that touched seabed        
        #DEPOSITION
        if (self.state[:np]==2).sum():
          #import pdb;pdb.set_trace()
          state_tmp=self.state[:np].copy() # the .copy() is to avoid some weird behaviour, there may be a better way to do this
          on_seabed_tmp=self.on_seabed[:np].copy()
          ind=numpy.where(state_tmp==2) # id of particles that touched the seabed 
          state_tmp[ind]=1 # set considered particle's state back to 1 so that they remain in the computational pool
          # At this stage,particle that touched the seabed were already re-suspended to a small height above seabed  (in posi)
          # Get bed shear stresses at these locations i.e. self.post, and t2
          tau_cur,tau_cw,tau_max,topo = self.bedshearstress_cw(self.post[:np,:],t2) #at all active particles
          tau = numpy.maximum.reduce([tau_cur,tau_cw]) # maybe we should use taumax in presence of waves?    
          # deposition on the seabed if state==2 and tau<=tau_crit_depos and tau<=tau_crit_eros and 
          id_depos=numpy.where(numpy.logical_and.reduce(( tau<=self.props.get('tau_crit_depos') , tau<=self.props.get('tau_crit_eros') , self.state[:np]==2 )))
          id_no_depos=numpy.where( numpy.logical_and( numpy.logical_or( tau>self.props.get('tau_crit_depos') , tau>=self.props.get('tau_crit_eros') ), self.state[:np]==2 ) )
          # flag deposited particles in on_seabed array - this will allows later identification for stopping advection/settling
          on_seabed_tmp[id_depos]=1 
          #if deposition occurs set particle depths to seabed depth
          posi[id_depos,2]=topo[id_depos,0] 
          # if not deposition occurs, resuspend within a bottom layer of thickness self.seabedlayer_thick (0.2 m for now) (normally distributed)
          posi[id_no_depos,2]=topo[id_no_depos,0]+self.seabedlayer_thick*numpy.random.uniform(0.0,1.0,numpy.size(topo[id_no_depos,0]))
          #update top-copy arrays
          self.state[:np]=state_tmp 
          self.on_seabed[:np]=on_seabed_tmp 
          self.post[:self.np,:]=posi[:self.np,:]
          if (self.post[:self.np,2]>0).any():
            import pdb;pdb.set_trace()

          # print 'DEPOSITION check start'
          # print 'nb part that touched bottom and deposited=%s' % (len(id_depos[0]))
          # print 'nb part that touched bottom and were resuspended =%s' % (len(id_no_depos[0]))
          # print 'id_depos = %s' % (id_depos)
          # print 'tau = %s' % (tau)
          # print 'tau at deposition location = %s' % (tau[id_depos])
          # print 'tau at resuspension location = %s' % (tau[id_no_depos])               
          # print 'tau_depos_crit < %s , tau_eros_crit < %s' % (self.props.get('tau_crit_depos'),self.props.get('tau_crit_eros'))
          # print 'id_no_depos = %s' % (id_no_depos)
          # print 'on_seabed_tmp after check = %s' % on_seabed_tmp
          # # print 'state_tmp after check= %s' % state_tmp
          # print 'depth POS= %s' % self.pos[:np,2]
          # print 'depth POST= %s' % self.post[:np,2]
          # print 'depth POSI= %s' % posi[:np,2]
          # print 'depth after check= %s' % posi[:,2]
          # print 'TOPO= %s' % topo[:,0]
          # print 'DEPOSITION check end'
          #import pdb;pdb.set_trace()

          # if (self.post[:self.np,2]<topo[:,0]).any():
          #   import pdb;pdb.set_trace()

        # RE-SUSPENSION
        if (self.on_seabed[:np]==1).sum():
          state_tmp=self.state[:np].copy()
          on_seabed_tmp=self.on_seabed[:np].copy()
          ind=numpy.where(on_seabed_tmp==1) # id of particles previously deposited on seabed
          # Get bed shear stresses at these locations i.e. self.post, and t2
          tau_cur,tau_cw,taumax,topo = self.bedshearstress_cw(self.post[:np,:],t2) #at all active particles
          tau = numpy.maximum.reduce([tau_cur,tau_cw]) # maybe we should use taumax in presence of waves? 

          # erosion/resuspension if tau>tau_crit_eros and on_seabed=1
          id_eros=numpy.where(numpy.logical_and( tau>=self.props.get('tau_crit_eros') , self.on_seabed[:np]==1 ))
          id_no_eros=numpy.where(numpy.logical_and( tau<self.props.get('tau_crit_eros') , self.on_seabed[:np]==1 ))
          on_seabed_tmp[id_eros]=0 # particles are not deposited on the seabed anymore
          # if  resuspension occurs, resuspend within a bottom layer of thickness self.seabedlayer_thick (0.2 m for now) (normally distributed)
          posi[id_eros,2]=topo[id_eros,0]+self.seabedlayer_thick*numpy.random.uniform(0.0,1.0,numpy.size(topo[id_eros,0])) # random number between 0-1
          # ensure previously deposited particles stay on seabed
          posi[id_no_eros,2]=topo[id_no_eros,0]

          self.on_seabed[:np]=on_seabed_tmp
          self.state[:np]=state_tmp # update the top state array
          self.post[:self.np,:]=posi[:self.np,:]
          if (self.post[:self.np,2]>0).any():
            import pdb;pdb.set_trace()

          # print 'RESUSPENSION check start'
          # print 'nb part that were deposited and were NOT resuspended=%s' % (len(id_no_eros[0]))
          # print 'nb part that  were deposited and were resuspended =%s' % ( len(id_eros[0]) )
          # print 'tau = %s' % (tau)
          # print 'tau where resuspension = %s' % (tau[id_eros])
          # print 'tau where NO resuspension = %s' % (tau[id_no_eros])          
          # print 'tau_eros_crit > %s' % (self.props.get('tau_crit_eros'))
          # print 'id_no_eros = %s' % (id_no_eros)
          # print 'id_eros = %s' % (id_eros)
          # print 'on_seabed_tmp after check = %s' % on_seabed_tmp
          # print 'state_tmp after check= %s' % state_tmp
          # print 'RESUSPENSION check end'
          #import pdb;pdb.set_trace()

          # if (self.post[:self.np,2]<topo[:,0]).any():
          #   import pdb;pdb.set_trace()

        # if (self.state[:np]==2).any():
        #   import pdb;pdb.set_trace() #there should NOT be any state==2

      # if (self.post[:self.np,2]<self.dep[:self.np]).any():
      #    import pdb;pdb.set_trace()
        
      if 'Elevation' in sticker.__class__.__name__:
        self.elev[:self.np]=sticker.interp(posi[:self.np,:],t2,imax=1)[:,0]
        self.post[:self.np,:]=posi[:self.np,:]

    #if self.unstick<=0.:
    # print self.state[self.state>1]
    #  self.state[self.state>1]=-1c
    #self.post[:self.np,:]=posi[:self.np,:] 

  def advect(self,t1,t2,order=4):
    if self.np==0:return

    PassiveTracer.advect(self,t1,t2,order) # general advection of all particles
    if hasattr(self, 'on_seabed'):
      if self.on_seabed.sum(): #there are some particles that were previously deposited
        # "cancel" advection of deposited particles (set self.post to self.pos)
        self.post[numpy.nonzero(self.on_seabed),:]=self.pos[numpy.nonzero(self.on_seabed),:]
        # add vertical settling to suspended particles only, deposited particles stay at the same depth
        self.post[:self.np,2]=numpy.where(self.on_seabed[:self.np]==1,self.pos[:self.np,2],self.post[:self.np,2]+86400.*(t2-t1)*self.w0[:self.np])
      else:
        # if no particles deposited - same as in BuoyantTracer.advect
        self.post[:self.np,2]+=86400.*(t2-t1)*self.w0[:self.np]

  
  def diffuse(self,t1,t2):
    BuoyantTracer.diffuse(self,t1,t2)
    #Heavy particle diffusion mods could go in here

  def bedshearstress_cw(self,p,time=None,imax=2):
    """Computation of bed shear stress due to current and waves
    current-related stress is computed following a drag-coefficient approach
    wave-related stress is computed following Van Rijn approach
    combined wave-current mean and max stresses are computed followin Soulsby(1995) approach
    Arguments:
      self : Material object, expected to include fields movers, reactors (if input)
      p: particle positions array (Nx3)
      time: 
    Returns:
      tau_cur : current-related bed shear stress tau_cur
      tau_cw : combined mean current-wave bed shear stress 
      tau_max: combined max current-wave bed shear stress tau_max
      topo : water depth at particle positions (of first mover) 
    """
    rhow=1027 # default volumic mass for seawater
    tau_cur=numpy.tile(0.0,numpy.size(p,0)) # allocate
    tau_cw=numpy.tile(0.0,numpy.size(p,0)) # allocate
    tau_max=numpy.tile(0.0,numpy.size(p,0)) # allocate
    # current-related bed shear stress (sum of all movers)
    for mover in self.movers[0:]:
      if mover.topo: # topo needed to define bedshear stress
        topo=mover.topo.interp(p,None,3)     
        if (not mover.is3d) and (mover.z0>0): # mover is a 2D-depth averaged current
          # temporarily set mover.z0 to 0.0 so that mover.interp yields the un-corrected depth-averaged current (direct interpolation, no log profile )
          # not super elegant, probably a better way to do this - anyway to access GriddedTide or GriddedData from here ?
          z0_tmp=copy.copy(mover.z0)
          mover.z0=0.0
          u2dhim=mover.interp(p,time,imax)   #u2dhim=mover.interp(self,p,time,imax)        
          mover.z0=z0_tmp
          u2dhim_mag=(u2dhim[:,0]**2+u2dhim[:,1]**2)**0.5
          # Drag coefficient for 2D case using water depth and z0 (see COHERENS manual eq.7.2, or Delft3d)
          Cdrag=( 0.4 /(numpy.log(abs(topo[:,0] /mover.z0))-1) )**2
          #Now compute the bed shear stress [N/m2] 
          tau_cur+=rhow*Cdrag*u2dhim_mag**2             
        elif (mover.is3d) and (mover.z0>0):   # mover is a 3D current field
          #import pdb;pdb.set_trace()
          # Assume the first grid point above the bed is assumed to be the top of the logarithmic boundary layer
          # the log profile extends from than last wet bin level, to the bottom
          # see COHERENS manual eq 7.1/7.2
  
          # find closest "wet" vertical levels at each particle locations
          bin_lev=numpy.zeros(len(p[:,0]))
          for lev in mover.lev:
            bin_lev[topo[:,0]<=lev]=lev
          #vertical height from last wet vertical bin to seabed
          zb=bin_lev-topo[:,0]
          #current computed at last wet vertical bin
          uub=mover.interp(numpy.vstack((p[:,0],p[:,1],bin_lev)).T,time,imax)
          uub_mag=(uub[:,0]**2+uub[:,1]**2)**0.5
          # Drag coefficient for 3D case using zb and z0 (see COHERENS manual eq.7.2, or Delft3d)
          Cdrag=( 0.4 /(numpy.log(abs(zb /mover.z0))-1) )**2 
         #Now compute the bed shear stress [N/m2]
          tau_cur+=rhow*Cdrag*uub_mag**2

    # wave-related bed shear stress
    if len(self.reactors)>0: # check if wave forcing is included 
      # for now assume that if reactors exist, they will be correctly input
      # computation of wave-related and combined bed shear stresses based on code from 
      #https://svn.oss.deltares.nl/repos/openearthtools/trunk/matlab/general/phys_fun/bedshearstresses.m 
      #https://svn.oss.deltares.nl/repos/openearthtools/trunk/matlab/general/phys_fun/sandandmudtransport.m
      hs=self.reactors['hs'].interp(p[:1,:],time)[:,0] #wave height
      tp=self.reactors['tp'].interp(p[:1,:],time)[:,0] #wave period
      # wave-related roughness
      # vanRijn 2007 suggests same equations than for current-related roughness where 20* d50 <ksw<150*d50
      # here we are using nikuradse for consistency with the use of z0 in the mover class for now
      ksw=30*self.movers[0].z0  
      topo= self.movers[0].topo.interp(p,None,3)
      w=2*numpy.pi/tp
      kh = qkhfs( w, topo[:,0] ) # dispersion relationship
      Adelta = hs/(2.*numpy.sinh(k*topo[:,0])) # peak wave orbital excursion
      Udelta = (numpy.pi*hs)/(tp*numpy.sinh(kh))  # peak wave orbital velocity
      fw = numpy.exp(-5.977+5.213*(Adelta/ksw)**-0.194)  # wave-related friction coefficient (van Rijn)
      fw = numpy.min(fw,0.3)
      tau_wave = 0.25 * rhow * fw * (Udelta)**2 # wave-related bed shear stress
      #cycle mean bed shear stress according to Soulsby,1995
      tau_cw=tau_cur*[1+1.2*(tau_wave/(tau_cur+tau_wave))**3.2]
      # max bed shear stress during wave cycle
      tau_max=[tau_wave**2+tau_cur**2]**0.5
    else:
      tau_max=tau_cur
      tau_cw=tau_cur
    # if (tau_cur==0).any():
    #   import pdb;pdb.set_trace()
    return tau_cur,tau_cw,tau_max,topo