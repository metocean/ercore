#!/usr/bin/env python
import datetime,sys,os
sys.path.append(os.path.join(os.path.dirname(__file__),'..','..','src'))
from pylab import *
from ercore._flib_ercore import interp3d,interph,slipvel
from ercore import ERcore
from ercore.materials import PassiveTracer,BuoyantTracer
from ercore.fields import GriddedMover,ConstantDiffuser,VariableDiffuser,GriddedTopo
from ercore.shoreline import Shoreline

dep=GriddedTopo('depth',['dep'],file='uds_gsb_test.nc')
currents=GriddedMover('cur',['uo','vo'],is3d=True,file='uds_gsb_test.nc',topo=dep)
tide=GriddedMover('tide',['ut','vt'],is3d=False,file='uds_gsb_test.nc')
#diff=ConstantDiffuser('diff',['diffx','diffy','diffz'],diffx=1,diffy=1,diffz=0.001)
diff=VariableDiffuser('diff',['diffx','diffy','diffz'],diffz=0.001,P0=[170.5,-46,0])
shoreline=Shoreline('shoreline','shoreline2.bnd')

cloud=PassiveTracer('cloud',10000,movers=[currents,tide],diffusers=[diff],stickers=[shoreline,dep],reln=1000,P0=[170.5,-46,0],dw_min=0.1,dw_max=0.2,tstart=datetime.datetime(2009,1,1),tend=datetime.datetime(2009,1,1,12))

particles=BuoyantTracer('particles',10000,movers=[currents,tide],diffusers=[diff],stickers=[shoreline,dep],reln=1000,P0=[170.5,-46,-40],w0=-0.1)
particles2=BuoyantTracer('particles2',10000,movers=[currents,tide],diffusers=[diff],stickers=[shoreline,dep],reln=1000,P0=[170.5,-46,-40],w0=-0.002)

ercore=ERcore(tout=3600.)
ercore.materials=[cloud,particles,particles2]
ercore.run(datetime.datetime(2009,1,1),datetime.datetime(2009,1,1,12),900)

figure()
plot(cloud.props['P0'][0],cloud.props['P0'][1],'r+')
plot(cloud.pos[:,0],cloud.pos[:,1],'.')
for i,j in enumerate(shoreline.polyi):
  n=shoreline.polyn[i]
  plot(shoreline.slx[j:j+n-1],shoreline.sly[j:j+n-1],'k')
show()
figure()
plot(cloud.pos[:,0],cloud.pos[:,2],'.')
plot(particles.pos[:,0],particles.pos[:,2],'.')
plot(particles2.pos[:,0],particles2.pos[:,2],'.')
plot(cloud.pos[:,0],dep.interp(cloud.pos))
show()
  

  

  
