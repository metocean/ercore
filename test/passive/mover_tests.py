#!/usr/bin/env python
import datetime,sys,os
sys.path.insert(0,os.path.join(os.path.dirname(__file__),'..','..'))
from pylab import *
from ercore._flib_ercore import interp3d,interph,slipvel
from ercore import ERcore
from ercore.materials import PassiveTracer,BuoyantTracer
from ercore.fields import GriddedMover,ConstantDiffuser,VariableDiffuser,GriddedTopo,ConstantMover,GridDataGroup
from ercore.shoreline import Shoreline

P0=[170.5,-46,-10]

dep=GriddedTopo('depth',['dep'],file='uds_gsb_test.nc',zinvert=True)
currents=GriddedMover('cur',['uo','vo'],is3d=True,file='uds_gsb_test.nc',topo=dep)
tide=GriddedMover('tide',['ut','vt'],is3d=False,file='uds_gsb_test.nc')
constantsurf=ConstantMover('constantsurf',['u','v'],u=1.0,v=0,topo=dep)
constant3D=ConstantMover('constant3D',['u','v'],is3d=True,levels=[0,-10,-30,-40],u=[1.0,0.5,0.0,-0.5],v=[0,0,0,0],topo=dep)
constant3Dsub=ConstantMover('constant3Dsub',['u','v'],is3d=True,levels=[0,-10,-30,-40],u=[0.0,-0.5,-1.0,-1.5],v=[0,0,0,0],surfsub=True,topo=dep)
diff=VariableDiffuser('diff',['diffx','diffy','diffz'],diffz=0.001,P0=P0)
shoreline=Shoreline('shoreline','shoreline2.bnd')
movergroup=GridDataGroup('group',['uo','vo'],[currents])

particles=BuoyantTracer('particlesA',10000,movers=[movergroup],stickers=[shoreline,dep],reln=10,P0=P0,w0=0.0)
particles3=BuoyantTracer('particlesC',10000,movers=[constant3D],stickers=[shoreline,dep],reln=10,P0=P0,w0=0.0)
particles2=BuoyantTracer('particlesB',10000,movers=[constantsurf,constant3Dsub],stickers=[shoreline,dep],reln=10,P0=P0,w0=0.0)

ercore=ERcore(tout=900.,rkorder=4)
ercore.materials=[particles]
ercore.run(datetime.datetime(2009,1,1),datetime.datetime(2009,1,3),900)
ercore.rkorder=3
ercore.materials=[particles3]
ercore.run(datetime.datetime(2009,1,1),datetime.datetime(2009,1,3),900)

#figure()
#for i,j in enumerate(shoreline.polyi):
#  n=shoreline.polyn[i]
#  plot(shoreline.slx[j:j+n-1],shoreline.sly[j:j+n-1],'k')
#show()
figure()
plot(particles.pos[:,0],particles.pos[:,2],'.')
plot(particles3.pos[:,0],particles3.pos[:,2],'.')
#plot(particles2.pos[:,0],dep.interp(particles2.pos),'.')
show()
  

  

  
