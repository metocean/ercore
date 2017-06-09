#!/usr/bin/env python
import numpy,datetime,sys,os,time
sys.path.insert(1,os.path.join(os.path.dirname(__file__),'..','..'))
#import ercore._flib_ercore as flib
from pylab import *
import matplotlib
from _flib_ercore import interp3d,interph,slipvel
from ercore import ERcore
from ercore.fields import *
from ercore.shoreline import Shoreline
from ercore.materials.hydrocarbons import HCPlume,HCGas,HCDroplets,HCSlick

test=int(sys.argv[1])
import pdb;pdb.set_trace()

if test==1:
#Test 1 - Gas plume tank scale Chen and Yapa (2004)
  current=ConstantMover('current',['u','v','w'],u=-0.1,v=0,w=0)
  tide=ConstantMover('tide',['u','v','w'],u=0.,v=0,w=0)
  wind=ConstantMover('wind',['ugrd10m','vgrd10m'],ugrd10m=0.,vgrd10m=0)
  surfmover=ConstantMover('vel',['u','v'],u=-0.1,v=0)
  temp=ConstantReactor('temp',['temp'],temp=20)
  salt=ConstantReactor('salt',['salt'],salt=0)
  diff=ConstantDiffuser('diff',['diffx','diffy','diffz'],diffx=0.01,diffy=0.01,diffz=0.001)
  surfdiff=ConstantDiffuser('surfdiff',['diffx','diffy'],diffx=0.01,diffy=0.01)
  d=0.0007 #nozzle diameter
  GOR=1.0
  F=500*(1+GOR) #ml/min
  A=0.25*numpy.pi*d*d
  V0=0.000001*F/60/A
  Z0=-0.6
  import pdb;pdb.set_trace()

  plume1=HCPlume('plume1',100000,is3d=True,movers=[current],
                                          reactors=[temp,salt],
                                          V0=numpy.array([0,0,V0]),
                                          P0=numpy.array([0,0,Z0]),
                                          T0=80,D0=870,wb=0.27,
                                          B0=d/2,GOR=GOR,Mg=28.97e-3,
                                          Cpl=1.8,spwn=1000,visc=0.033,
                                          IFT=0.02)#,wb=0.27)#,tstep=3600.)

  # plume2=HCPlume('plume2',10000,is3d=True,movers=[current],reactors=[temp,salt],V0=numpy.array([0,0,V0]),P0=numpy.array([0,0,Z0]),T0=80,D0=870,wb=0.27,B0=d/2,GOR=GOR,Mg=28.97e-3,Cpl=1.8,gsep=False)

#Test 2 - Idealised gas plume
elif test==2:
  current=ConstantMover('current',['u','v','w'],u=0.16,v=0,w=0)
  tide=ConstantMover('tide',['u','v','w'],u=0.,v=0,w=0)
  surfmover=ConstantMover('vel',['u','v'],u=0.16,v=0)
  wind=ConstantMover('wind',['ugrd10m','vgrd10m'],ugrd10m=-10,vgrd10m=3)
  temp=ConstantReactor('temp',['temp'],temp=15)
  salt=ConstantReactor('salt',['salt'],salt=35)
  diff=ConstantDiffuser('diff',['diffx','diffy','diffz'],diffx=0.2,diffy=0.2,diffz=0.001)
  surfdiff=ConstantDiffuser('surfdiff',['diffx','diffy'],diffx=1,diffy=1)
  d=0.0889 #nozzle diameter
  F=0.0184 #Volume flux liquid
  A=0.25*numpy.pi*d*d
  V0=F*(1+5.345)/A
  Z0=-400
  GOR=5.345
  import pdb;pdb.set_trace()
  plume1=HCPlume('plume1',20000,is3d=True,movers=[current],reactors=[temp,salt],diffusers=[diff],V0=numpy.array([0,0,V0]),P0=numpy.array([0,0,Z0]),T0=80,D0=893,db=0.006,B0=d,GOR=GOR,Mg=28.97e-3,Cpl=1.8,spwn=10,visc=0.033,IFT=0.02,wb=0.27,tstep=3600.)
  plume2=HCPlume('plume2',20000,is3d=True,movers=[current],reactors=[temp,salt],diffusers=[diff],V0=numpy.array([0,0,V0]),P0=numpy.array([0,0,Z0]),T0=80,D0=893,db=0.006,B0=d,GOR=GOR,Mg=28.97e-3,Cpl=1.8,spwn=10,visc=0.033,IFT=0.02,wb=0.27,gsep=False)
elif test==3:
  current=GriddedMover('ocean',['uo','vo'],file='../passive/uds_gsb_test.nc')
  tide=GriddedMover('tide',['ut','vt'],file='../passive/uds_gsb_test.nc')
  surfmover=GriddedMover('surf',['ut','vt'],file='../passive/uds_gsb_test.nc')
  wind=ConstantReactor('wind',['ugrd10m','vgrd10m'],ugrd10m=-15,vgrd10m=0)
  temp=GriddedReactor('temp',['temp'],file='../passive/uds_gsb_test.nc')
  salt=GriddedReactor('salt',['salt'],file='../passive/uds_gsb_test.nc')
  diff=ConstantDiffuser('diff',['diffx','diffy','diffz'],diffx=0.2,diffy=0.2,diffz=0.001)
  surfdiff=ConstantDiffuser('diff',['diffx','diffy'],diffx=1,diffy=1)
  shoreline=Shoreline('shore',file='../passive/shoreline2.bnd')
  d=0.0889 #nozzle diameter
  F=0.0184 #Volume flux liquid
  A=0.25*numpy.pi*d*d
  V0=F*(1+5.345)/A
  Z0=-389
  GOR=5.345
  plume1=HCPlume('plume1',20000,is3d=True,movers=[current,tide],reactors=[temp,salt],diffusers=[diff],V0=numpy.array([0,0,V0]),P0=numpy.array([170.94,-46.08,Z0]),T0=80,D0=893,db=0.006,B0=d,GOR=GOR,Mg=28.97e-3,Cpl=1.8,spwn=10,visc=0.033,IFT=0.02,wb=0.27,tstep=3600)
  plume2=HCPlume('plume2',20000,is3d=True,movers=[current,tide],reactors=[temp,salt],diffusers=[diff],V0=numpy.array([0,0,V0]),P0=numpy.array([170.94,-46.08,Z0]),T0=80,D0=893,db=0.006,B0=d,GOR=GOR,Mg=28.97e-3,Cpl=1.8,spwn=10,visc=0.033,IFT=0.02,wb=0.27,gsep=False)

gas=HCGas('gas',100000,is3d=True,movers=[current,tide],reactors=[temp,salt],diffusers=[diff],Mg=28.97e-3,P0=numpy.array([0,0,Z0]),reln=0)
drops=HCDroplets('droplets',100000,is3d=True,movers=[current,tide],reactors=[temp,salt],diffusers=[diff],P0=numpy.array([0.,0.,Z0]),reln=0)
slick=HCSlick('slick',100000,movers=[surfmover],reactors=[wind],diffusers=[surfdiff]) #,stickers=[shoreline])

ercore=ERcore(geod=True,outpath='./hydrocarbonplume_outputs')
ercore.materials=[plume1,gas,drops,slick]

def run(obj):
  t0=time.time()
  obj.run(733410,733410+48*3600.0/86400.,900.0)
  print time.time()-t0
#import profile
#profile.run('run(ercore)')


run(ercore) 

figure()
subplot(121)
for plume in ercore.materials:
  plot(plume.pos[:plume.np,0],plume.pos[:plume.np,2],'.')
  if plume.__class__.name=='HCPlume':
    plume.b[plume.np:]=0.
    ang=numpy.arctan2(plume.u[:,2],plume.u[:,0])
    #x-rotation
    A=(plume.u[:,0]**2+plume.u[:,2]**2)**0.5
    sin=-plume.u[:,0]/A
    cos=plume.u[:,2]/A
    plot(plume.pos[:,0]+plume.b*cos*plume.mfx[:,0],plume.pos[:,2]+plume.b*sin,'.')
    plot(plume.pos[:,0]-plume.b*cos*plume.mfx[:,0],plume.pos[:,2]-plume.b*sin,'.')
if not ercore.geod:axis('equal')
subplot(122)
for plume in ercore.materials:
  plot(plume.pos[:plume.np,0],plume.pos[:plume.np,1],'.')
  if plume.__class__.name=='HCPlume':
    ang=numpy.arctan2(plume.u[:,1],plume.u[:,0])
    #x-rotation
    A=(plume.u[:,0]**2+plume.u[:,1]**2)**0.5
    sin=-plume.u[:,0]/A
    cos=plume.u[:,1]/A
    plot(plume.pos[:,0]+plume.b*cos*plume.mfx[:,0],plume.pos[:,1]+plume.b*sin*plume.mfx[:,1],'.')
    plot(plume.pos[:,0]-plume.b*cos*plume.mfx[:,0],plume.pos[:,1]-plume.b*sin*plume.mfx[:,1],'.')
if test==3:
  for i,j in enumerate(shoreline.polyi):
    n=shoreline.polyn[i]
    plot(shoreline.slx[j:j+n-1],shoreline.sly[j:j+n-1],'k')
axis('equal')
show()


  
