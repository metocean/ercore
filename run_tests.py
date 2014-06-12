#!/usr/bin/env python
import numpy,datetime
#import ercore._flib_ercore as flib
from pylab import *
from _flib_ercore import interp3d,interph,slipvel
from ercore import ERcore

if True:
  zz,yy,xx=numpy.mgrid[0:-10:-1,100:201:2,0:101:10]
  d=(xx*0.5+yy)*(1+0.1*zz)
  x=numpy.array([0,0.5,5,10,23,55])
  y=numpy.array([0,100,105,123,155,178])
  z=numpy.array([0,-1,-4,-5.5,-8.9,-20])
  vout=interp3d(d,x,y,z,0,100,10,2,zz[:,0,0])
  print vout
  print (numpy.maximum(x,0)*0.5+numpy.maximum(y,100))*(1+0.1*numpy.maximum(z,-9))

if False:#Test - bubble/droplet rise velocity  ZHeng and Yapa (2000)
  d=10**numpy.arange(-1,2,0.1)
  ones=numpy.ones(d.shape)
  visc=[1.002e-3,1.002e-3,0.18,0.957e-3]
  ift=[0.0728,0.0465,0.0638,0.0725]
  bdens=[1.205,1.587,1.82,1.82]
  dens=[998.2,998.2,1234,997.7]
  temp=[20,20,25,22]
  
  for i in range(0,4):
    U=flib.slipvel.bubble_slip(0.001*d,dens[i]*ones,(dens[i]-bdens[i])*ones,visc[i]*ones,temp[i]*ones,ift)
    print U
    plot(d,U)
    #xscale('log')
    #yscale('log')
    show()
    

#Test 2 - Oil slick advection
if False:
  from ercore.materials.hydrocarbons import HCSlick
  from ercore.fields import GriddedMover,GriddedReactor,ConstantDiffuser
  from ercore.shoreline import ShorelineData
  currents=GriddedMover('cur',['water_u','water_v'],is3d=False,file='../test/slick/uds_out.gnome')
  winds=GriddedReactor('wind',['air_u','air_v'],is3d=False,file='../test/slick/uds_out.gnome')
  diff=ConstantDiffuser('diff',['diffx','diffy','diffz'],diffx=1,diffy=1,diffz=1)
  shoreline=ShorelineData('shoreline','../test/slick/shoreline.bnd')
  slick=HCSlick('slick',10000,movers=[currents],reactors=[winds],diffusers=[diff],stickers=[shoreline],reln=1000,P0=[173.9,-39.1,0],dw_min=0.1,dw_max=0.2)
  
  ercore=ERcore()
  ercore.materials=[slick]
  ercore.run(datetime.datetime(2013,3,26),datetime.datetime(2013,3,27),900)
  
  plot(slick.props['P0'][0],slick.props['P0'][1],'r+')
  plot(slick.pos[:,0],slick.pos[:,1],'.')
  for i,j in enumerate(shoreline.polyi):
    n=shoreline.polyn[i]
    plot(shoreline.slx[j:j+n-1],shoreline.sly[j:j+n-1],'k')
  show()
  #Checks out OK against GNOME
  
#Test 3 - 3D drifter advection
if True:
  from ercore.materials import PassiveTracer
  from ercore.fields import GriddedMover,ConstantDiffuser,GriddedTopo
  from ercore.shoreline import ShorelineData
  dep=GriddedTopo('depth',['dep'],file='../test/passive/uds_gsb_test.nc')
  currents=GriddedMover('cur',['uo','vo'],is3d=False,file='../test/passive/uds_gsb_test.nc',topo=dep)
  tide=GriddedMover('tide',['ut','vt'],is3d=False,file='../test/passive/uds_gsb_test.nc')
  #diff=ConstantDiffuser('diff',['diffx','diffy','diffz'],diffx=1,diffy=1,diffz=0.001)
  diff=VariableDiffuser('diff',['diffx','diffy','diffz'],diffz=0.001)
  shoreline=ShorelineData('shoreline','../test/passive/shoreline2.bnd')
  
  cloud=PassiveTracer('cloud',10000,movers=[currents,tide],diffusers=[diff],stickers=[shoreline,dep],reln=1000,P0=[170.5,-46,0],dw_min=0.1,dw_max=0.2)
  
  ercore=ERcore()
  ercore.materials=[cloud]
  ercore.run(datetime.datetime(2009,1,1),datetime.datetime(2009,1,2),900)
  
  figure()
  plot(cloud.props['P0'][0],cloud.props['P0'][1],'r+')
  plot(cloud.pos[:,0],cloud.pos[:,1],'.')
  for i,j in enumerate(shoreline.polyi):
    n=shoreline.polyn[i]
    plot(shoreline.slx[j:j+n-1],shoreline.sly[j:j+n-1],'k')
  show()
  figure()
  plot(cloud.pos[:,0],cloud.pos[:,2],'.')
  show()
  
if True:
  from ercore.fields import ConstantMover,ConstantReactor,ConstantDiffuser
  from ercore.materials.hydrocarbons import HCPlume,HCGas,HCDroplets,HCSlick
  if False:
  #Test 3 - Gas plume tank scale Chen and Yapa (2004)
    mover=ConstantMover('vel',['u','v','w'],u=-0.1,v=0,w=0)
    surfmover=ConstantMover('vel',['u','v'],u=-0.1,v=0)
    temp=ConstantReactor('temp',['temp'],temp=20)
    salt=ConstantReactor('salt',['salt'],salt=0)
    diff=ConstantDiffuser('diff',['diffx','diffy','diffz'],diffx=0.01,diffy=0.01,diffz=0.001)
    surfdiff=ConstantDiffuser('diff',['diffx','diffy'],diffx=0.01,diffy=0.01)
    d=0.0007 #nozzle diameter
    GOR=1.0
    F=500*(1+GOR) #ml/min
    A=0.25*numpy.pi*d*d
    V0=0.000001*F/60/A
    Z0=-0.6
    plume1=HCPlume('plume',100000,is3d=True,movers=[mover],reactors=[temp,salt],V0=numpy.array([0,0,V0]),P0=numpy.array([0,0,Z0]),T0=80,D0=870,wb=0.27,B0=d/2,GOR=GOR,Mg=28.97e-3,Cpl=1.8,spwn=1000)
    plume2=HCPlume('plume2',10000,is3d=True,movers=[mover],reactors=[temp,salt],V0=numpy.array([0,0,V0]),P0=numpy.array([0,0,Z0]),T0=80,D0=870,wb=0.27,B0=d/2,GOR=GOR,Mg=28.97e-3,Cpl=1.8,gsep=False)
  #Test 4 - Gas plume GOM 
  else:
    mover=ConstantMover('vel',['u','v','w'],u=0.16,v=0,w=0)
    surfmover=ConstantMover('vel',['u','v'],u=0.16,v=0)
    temp=ConstantReactor('temp',['temp'],temp=15)
    salt=ConstantReactor('salt',['salt'],salt=35)
    diff=ConstantDiffuser('diff',['diffx','diffy','diffz'],diffx=0.1,diffy=0.1,diffz=0.001)
    surfdiff=ConstantDiffuser('diff',['diffx','diffy'],diffx=0.1,diffy=0.1)
    d=0.0889 #nozzle diameter
    F=0.0184 #Volume flux liquid
    A=0.25*numpy.pi*d*d
    V0=F*(1+5.345)/A
    Z0=-400
    plume1=HCPlume('plume',20000,is3d=True,movers=[mover],reactors=[temp,salt],V0=numpy.array([0,0,V0]),P0=numpy.array([0,0,Z0]),T0=80,D0=893,db=0.006,B0=d,GOR=5.345,Mg=28.97e-3,Cpl=1.8,spwn=10,visc=0.033,IFT=0.02)
  
  gas=HCGas('gas',100000,is3d=True,movers=[mover],reactors=[temp,salt],diffusers=[diff],Mg=28.97e-3,P0=numpy.array([0,0,Z0]),reln=0)
  drops=HCDroplets('droplets',100000,is3d=True,movers=[mover],reactors=[temp,salt],diffusers=[diff],P0=numpy.array([0.,0.,Z0]),reln=0)
  slick=HCSlick('slick',10000,movers=[surfmover],diffusers=[surfdiff])
  
  ercore=ERcore(geod=False)
  ercore.materials=[plume1,gas,drops,slick]
  def run(obj):
    obj.run(0,800.0/86400.,800.0)
  #import profile
  #profile.run('run(ercore)')
  run(ercore) 
  
  for plume in ercore.materials:
    plot(plume.pos[:plume.np,0],plume.pos[:plume.np,2],'.')
    if plume.__class__.name=='HCPlume':
      ang=numpy.arctan2(plume.u[:,2],plume.u[:,0])
      #x-rotation
      A=(plume.u[:,0]**2+plume.u[:,2]**2)**0.5
      sin=-plume.u[:,0]/A
      cos=plume.u[:,2]/A
      plot(plume.pos[:,0]+plume.b*cos,plume.pos[:,2]+plume.b*sin,'.')
      plot(plume.pos[:,0]-plume.b*cos,plume.pos[:,2]-plume.b*sin,'.')
  axis('equal')
  show()
  

  
