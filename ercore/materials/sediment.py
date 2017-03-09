#!/usr/bin/env python
from ercore.materials import *
from ercore.lib import dens

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



#Class for sediment
#Suspended sediment is just a (negatively) buoyant tracer
#Resuspension and bedload could easily be incorporated - see my Master thesis ("A particle tracking model for sediment transport in the Nearshore Zone", Johnson 2000)
#Could add more sophisticated calculations of fall velocity based on temp/salinity


class Sediment(BuoyantTracer):
  __doc__=PassiveTracer.__doc__+"""
    d0:   Sediment grain size (mm) - If specified will over-ride w0
    rho_s:Sediment density (kg/m^3)
    tau_crit_eros: critical bed shear stress for erosion in N/m2 (default=0.2)
    tau_crit_deposition: critical bed shear stress for deposition in N/m2 (default=1000 - deposition always possible)
  """
  default_props={'d0':None,'rho_s':2650,'tau_crit_eros':0.2,'tau_crit_depos':1000}
  def initialize(self,t1,t2):
    #import pdb;pdb.set_trace()
    BuoyantTracer.initialize(self,t1,t2)
    if self.props.get('d0'):
      w0=0.0005*(props.get('rho_s')-dens(35,20))*props.get('d0')**2 #This is calculate for sea water with 35 PPT and 20C
      self.w0=numpy.tile(w0,self.npmax+1)
    #new array to facilitate bookkeeping of successive deposition/erosion  
    self.bot_layer=numpy.tile(0,self.npmax+1)
    #self.bot_layer is true when a particle is deposited on the seabed, false otherwise
  
  def react(self,t1,t2):
    pass

  def bedshearstress(self,p,time=None,imax=2):
    """Computation of bed shear stress
    Arguments:
      self :  
      p: particle positions
      uses interp function from GridData
    Returns:
      Bed shear stress (current only or current+wave if wave reactors are input)
    """
    rhow=1027
    #Current-related bed shear stresses 
    tau_cur,topo=self.movers[0].bedshearstress(p,time,imax=2) # 
    for mover in self.movers[1:]: # sum contributions from all movers if applicable
      tau_cur+=mover.bedshearstress(p,t2,imax=2)[0]
    # add wave contribution of applicable
    # check if wave forcing is included    
    if len(self.reactors)>0:
      # for now assume that if reactors exist, they will be correctly input
      # based on code from 
      #https://svn.oss.deltares.nl/repos/openearthtools/trunk/matlab/general/phys_fun/bedshearstresses.m 
      #https://svn.oss.deltares.nl/repos/openearthtools/trunk/matlab/general/phys_fun/sandandmudtransport.m

      hs=self.reactors['hs'].interp(self.pos[:1,:],t1)[:,0] #wave height
      tp=self.reactors['tp'].interp(self.pos[:1,:],t1)[:,0] #wave period
      # wave related roughness
      # vanRijn 2007 suggests same equations than for current-related roughness where 20* d50 <ksw<150*d50
      # here we are using nikuradse for consistency with the use of z0 for now
      ksw=30*self.movers[0].z0  
      topo= self.movers[0].topo.interp(self.pos,None,3)
      w=2*numpy.pi/tp
      kh = qkhfs( w, topo[:,0] )
      Adelta = hs/(2.*numpy.sinh(k*topo[:,0])) # peak wave orbital excursion
      Udelta = (numpy.pi*hs)/(tp*numpy.sinh(kh))  # peak wave orbital velocity
      fw = numpy.exp(-5.977+5.213*(Adelta/ksw)**-0.194)  # wave-related friction coefficient (van Rijn)
      fw = numpy.min(fw,0.3)
      tau_wave = 0.25 * rhow * fw * (Udelta)**2 # wave-related bed shear stress
      #cycle mean bed shear stress according to Soulsby,1995
      tau_cw=tau_cur*[1+1.2*(tau_wave/(tau_cur+tau_wave))**3.2]
      # max bed shear stress during wave cycle
      taumax=[tau_wave**2+tau_cur**2]**0.5
    else:
      taumax=tau_cur
      tau_cw=tau_cur
    return tau_cur,tau_cw,taumax,topo
    # compared bed shear stresses to critical stress and apply erosion/deposition accordingly
    # 
    #consider using a "bed layer" where particles can be become deposited (but not removed from the model), or re-eroded
    #similar approach used in delft3d part
    #https://svn.oss.deltares.nl/repos/delft3d/trunk/src/engines_gpl/part/packages/kernel_f/src/part10.f90
    
    # identify particles that touched the bottom
    #import pdb;pdb.set_trace()
    #ind=(self.post[:,2]<=topo[:,0])
    #status_codes={0:'Not released',1:'Released and active',-1:'Stuck to shoreline or bottom',-2:'Dead'}
    # in the main RUN call - "react" happens before "sticking"     
    #self.props.get('tau_crit_eros')

  def stick(self,t1,t2):
    """Do sticking for t1 to t2"""
    if self.np<1:return
    np=self.np
    posi=numpy.where(self.state[:np,None]<0,self.pos[:np,:],self.post[:np,:])
    #import pdb;pdb.set_trace()
    for sticker in self.stickers: 
      # check for interestion with stickers
      if 'Shoreline' in  sticker.__class__.__name__: # shoreline sticker 
        # Returns matrix of intersection positions, how does it affect state ? 
        posi[:self.np,:]=sticker.intersect(self.pos[:self.np,:],posi,self.state[:self.np])
      else: # 2D sticker 
        posi[:self.np,:]=sticker.intersect(self.pos[:self.np,:],posi,self.state[:self.np],t1,t2)

      if 'GriddedTopo' in sticker.__class__.__name__:
        self.dep[:self.np]=sticker.interp(posi[:self.np,:],imax=1)[:,0] # get depths at particles that touched seabed        
        
        #DEPOSITION
        ind=(self.state[:np]==2)
        state_tmp=self.state[:np]
        bot_layer_tmp=self.bot_layer[:np]
        if ind.sum():
          
          # At this stage,particle that touched the seabed were already re-suspended to a small height above seabed  (in posi)
          # Get bed shear stress at these locations i.e. self.post, and t2
          tau_cur,tau_cw,taumax,topo = self.bedshearstress(self.post[:np,:],t2) #at all active particles
          tau = numpy.maximum.reduce([tau_cur,tau_cw]) # maybe we should use taumax in presence of waves?    
          # deposition on the seabed if tau<=tau_crit_depos 
          id_depos=(tau<=self.props.get('tau_crit_depos')) & ind
          posi[ind & id_depos,2]=topo[ind & id_depos,0]+0.001 # set particle depths to seabed depth+1 mm -the other ones are left suspended in water column
          state_tmp[ind & id_depos]=0 # de-activate these particles these should stay in place, no advection/settling 
          bot_layer_tmp[ind & id_depos]=1 # flag deposited particle in bot_layer arrray
          state_tmp[ind & ~id_depos]=1    # the other particles are left re-suspended in water column (already done in sticker.intersect), and state is set back from 2 to 1 - 
          #update top-copy array
          self.state[:np]=state_tmp # update the top state array
          self.bot_layer[:np]=bot_layer_tmp # flag de-activated particles following deposition on seabed
          self.post[:self.np,:]=posi[:self.np,:] # update top post array

          print tau
          print id_depos
          print 'topo = %s' % (topo[ind & id_depos,0])
          print 'posi = %s' % (posi[ind & id_depos,2])
          print 'bot_layer_tmp = %s' % bot_layer_tmp
          print 'state_tmp = %s' % state_tmp

        # RE-SUSPENSION
        ind=(self.bot_layer[:np]==1)
        state_tmp=self.state[:np]
        bot_layer_tmp=self.bot_layer[:np]
        if ind.sum():
          #import pdb;pdb.set_trace() 
          # Get bed shear stress at these locations i.e. self.post, and t2
          tau_cur,tau_cw,taumax,topo = self.bedshearstress(self.post[:np,:],t2) #at all active particles
          tau = numpy.maximum.reduce([tau_cur,tau_cw]) # maybe we should use taumax in presence of waves? 
          id_eros=(tau>=self.props.get('tau_crit_eros')) & ind
          state_tmp[ind & id_eros]=1 # re-activate particles
          print 'bot_layer_tmp = %s' % bot_layer_tmp
          bot_layer_tmp[ind & id_eros]=0 # not in the bottom layer anymore
          # re-suspend in the water column - within a 1m thick layer above the seabed (normally distributed)
          posi[ind & id_eros,2]=posi[ind & id_eros,2]+numpy.random.uniform(0.0,1.0,numpy.size(posi[ind & id_eros,2])) # random number between 0-1
          
          self.bot_layer[:np]=bot_layer_tmp
          self.state[:np]=state_tmp # update the top state array
          self.post[:self.np,:]=posi[:self.np,:]

          print tau
          print self.props.get('tau_crit_eros')
          print id_eros
          print 'topo = %s' % (topo)
          print 'posi = %s' %  posi[:self.np,:]
          print 'bot_layer_tmp = %s' % bot_layer_tmp
          print 'state_tmp = %s' % state_tmp

        if (self.state[:np]==2).any():
          import pdb;pdb.set_trace() #there should NOT be any state==2

      if 'Elevation' in sticker.__class__.__name__:
        self.elev[:self.np]=sticker.interp(posi[:self.np,:],t2,imax=1)[:,0]

    #if self.unstick<=0.:
    #  self.state[self.state>1]=-1
    #self.post[:self.np,:]=posi[:self.np,:] 
  
  def diffuse(self,t1,t2):
    BuoyantTracer.diffuse(self,t1,t2)
    #Heavy particle diffusion mods could go in here