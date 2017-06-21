#!/usr/bin/env python
import numpy,datetime,sys,os,time
sys.path.insert(1,os.path.join(os.path.dirname(__file__),'..','..'))
#import ercore._flib_ercore as flib
from pylab import *
import matplotlib
from _flib_ercore import interp3d,interph,slipvel
from ercore import ERcore
from ercore.fields import *
from ercore.shoreline import Shoreline,Boundary
from ercore.materials.plumes import BuoyantPlume,BuoyantPlume_JETLAG,BuoyantPlume_DensityCurrent
from ercore.materials import BuoyantTracer

test=int(sys.argv[1])


# PLUME PARAMETERS
d=0.2 #nozzle diameter in meter
flux=0.1 ##Volume flux liquid m3/s
area=numpy.pi*(d/2)**2
V0=flux/area #Initial jet velocity <list> of u,v,w components
Z0=-10.0
# angle between the x-axis and jet in the x-y plane (i.e. horizontal plane) [deg]
theta_jet = 0 
# angle between the x-axis and jet in the x-z plane (i.e vertical plane) [deg]
phi_jet = 90
#
theta_jet = numpy.deg2rad(theta_jet);
phi_jet   = numpy.deg2rad(phi_jet)

uj0=V0*numpy.cos(phi_jet)*numpy.cos(theta_jet)
vj0=V0*numpy.cos(phi_jet)*numpy.sin(theta_jet)
wj0=V0*numpy.sin(phi_jet)

# PLUME INFOS 
B0=d/2 #: Initial jet radius <float>
V0=numpy.array([uj0,vj0,wj0])#: Initial jet velocity <list> of u,v,w components
Vb=0.0#: Bubble/dropplet terminal velocity (m/s) <float>
C0=1.0#: Jet initial concentration <float>
E=2.0#: Entrainment constant <float>
# BUOYANT PLUME INFOS
T0=80#: Initial temperature (C) of jet <float>
S0=35#: Initial salinity (PSU) of jet <float>
D0=800# Initial density (kg/m^3) of jet <float>
Cpl=1.8#: Specific heat capacity <float>


# *** Firt reactor should be temp, second should be salt
temp=ConstantReactor('temp',['temp'],temp=12)
salt=ConstantReactor('salt',['salt'],salt=35)
spawnclass='material_farfield'


if test==1:
  # test the BuoyantPlume_JETLAG turning into a Buoyant Tracer class in constant flows
  topo=ConstantTopo('depth',['topo'],topo=-10.)
  current=ConstantMover('vel',['u','v','w'],u=+0.2,v=0,w=0,topo=topo)
  diff=ConstantDiffuser('diff',['diffx','diffy','diffz'],diffx=0.2,diffy=0.2,diffz=0.001)
  # SPAWN CLASS INFO
  spawnclass='material_farfield1'
  reln=48*3600/900 #nb part released over the entire run period > 1 per release for now

  # Defintion of plume Material
  buoyantplume_jetlag=BuoyantPlume_JETLAG(id='plume1',
                                          nbuff=1e5,
                                          is3d=True,
                                          movers=[current],
                                          reactors=[temp,salt],
                                          V0=V0,
                                          P0=numpy.array([0,0,Z0]),
                                          T0=80,
                                          S0=10,
                                          B0=B0,
                                          Vb=Vb,
                                          D0=D0,
                                          C0=1.0,
                                          E=2,
                                          Cpl=1.8,
                                          spawn_class=spawnclass,
                                          spawn_type='surface')
  # Definition of Material that will be "created" after the nearfiled dynamics of the plume
  # ***no need to define a P0
  material_farfield=BuoyantTracer(id='material_farfield1',
                                  nbuff=10000,
                                  movers=[current],
                                  reactors=[temp,salt],
                                  reln=reln,
                                  w0=-0.001,
                                  unstick=0.0,
                                  maxage=3.0,
                                  tstart=733410.,
                                  tend=733410.+2.,
                                  ischild=1)
  # Run
  ercore=ERcore(geod=True,tout=900,outpath='./outputs')
  ercore.materials=[buoyantplume_jetlag,material_farfield]
  # ercore.run(tstart,tend,pdt)
  ercore.run(733410.,733410.+48.*3600./86400.,900.0)
if test==2: 
  # outfall type simulation
  # positively buoyant plume - smaller salinity and larger temperature than ambient 
  P0=[170.5,-46,-10]
  dep = GriddedTopo('depth',['dep'],file='../passive/uds_gsb_test.nc',zinvert=True)
  Z0=dep.interp(numpy.array([[P0[0],P0[1],0]]))
  shoreline  = Shoreline('shoreline','../passive/shoreline2.bnd')
  current = GriddedMover('cur',['uo','vo'],is3d=True,file='.././passive/uds_gsb_test.nc',topo=dep)
  # tide     = GriddedMover('tide',['ut','vt'],is3d=False,file='.././passive/uds_gsb_test.nc') #
  diff=ConstantDiffuser('diff',['diffx','diffy','diffz'],diffx=0.1,diffy=0.1,diffz=0.001)
  tstart = datetime.datetime(2009,1,1)
  tend   = datetime.datetime(2009,1,3)
  tstep  = 900.
  tout   = 900.

  # SPAWN CLASS INFO
  reln=48*3600/tstep #nb part released over the entire run period > 1 per release for now
  spawnclass='material_farfield2'
   # Defintion of plume Material
  buoyantplume_jetlag=BuoyantPlume_JETLAG(id='plume2',
                                          nbuff=1e6,
                                          is3d=True,
                                          movers=[current],
                                          reactors=[temp,salt],
                                          V0=V0,
                                          P0=numpy.array([P0[0],P0[1],Z0+.1]),
                                          T0=80,
                                          S0=10,
                                          B0=B0,
                                          Vb=Vb,
                                          D0=D0,
                                          C0=1.0,
                                          E=2,
                                          Cpl=1.8,
                                          spawn_class=spawnclass)
  # Defintion of Material that will be "created" after the nearfiled dynamics of the plume
  # ***no need to define a P0
  material_farfield=BuoyantTracer(id='material_farfield2',
                                  nbuff=10000,
                                  movers=[current],
                                  reactors=[temp,salt],
                                  stickers=[dep,shoreline],
                                  reln=reln,
                                  w0=-0.001,
                                  unstick=0.0,
                                  maxage=3.0,
                                  tstart=tstart,
                                  tend=tend,
                                  ischild=1)

  # Run
  ercore=ERcore(geod=True,tout=900,outpath='./outputs')
  ercore.materials=[buoyantplume_jetlag,material_farfield]
  ercore.run(tstart,tend,tstep)

if test==3:
  # Using Lyttelton Harbour flow model as test case
  # bottom release - outfall like - at C1
  # SITE------------------------------------------------------------------------------------------------------
  site_num=1
  #site=sys.argv[1]

  spawnclass='material_farfield3'
  site_name=['C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C2_2']
  x=[172.7211,172.7388,172.7555,172.7690,172.7829,172.7972,172.8104,172.8227,172.8361,172.8478,172.8599,172.7383]
  y=[-43.6131,-43.6133,-43.6090,-43.6073,-43.6056,-43.6039,-43.5989,-43.5939,-43.5885,-43.5838,-43.5790,-43.6111]
  
  tstart=datetime.datetime(2014,1,18,8,0,0)#tstart simulation
  tend=tstart+datetime.timedelta(days=2.)
  tstep  = 900.
  tout   = 900.


  with_resusp=0
  # BOUNDARY--------------------------------------------------------------------------------------------------
  bnd=Boundary('bnd',[172.6, 173.92, -43.7, -43.5])
  shoreline=Shoreline('shoreline','../resuspension/southnz_shoreline.bnd')
  # FALL VELOCITY AND TIMESTEPS-------------------------------------------------------------------------------
  w0=-1.0e-3
  # SOURCE DATA FILE
  DataFileLyt='/home/simon/0201_ercore_lyttelton/NEW_RUNS_DISPO/channel_runs/lyt_Original_cons.nc'
  #DEPTH-----------------------------------------------------------------
  depth_lyt_hbr=GriddedTopo('depth_lyt',['dep'],file=DataFileLyt,zinvert=True)
  Z0=depth_lyt_hbr.interp(numpy.array([[x[site_num-1],y[site_num-1],0]]))
  #DIFFUSION-----------------------------------------------------------
  diff=ConstantDiffuser('diff',['diffx','diffy','diffz'],diffx=0.1,diffy=0.1,diffz=0.0000)
  # TIDE DATA----------------------------------------------------------
  current=TidalMover('lyt_tide_hbr',['u','v'],file=DataFileLyt,topo=depth_lyt_hbr,zinvert=True,z0=0.0001)
  #-----------------------------------------------------------------------
  #RELEASE-----------------------------------------------------------
  reln=1*48*3600/tstep #nb part released over the entire run period
  ## MATERIALS
  # Defintion of plume Material
  buoyantplume_jetlag=BuoyantPlume_JETLAG(id='plume3',
                                          nbuff=1e5,
                                          is3d=True,
                                          movers=[current],
                                          reactors=[temp,salt],
                                          V0=V0,
                                          P0=numpy.array([x[site_num-1],y[site_num-1],Z0+.1]),
                                          T0=80,
                                          S0=10,
                                          B0=B0,
                                          Vb=Vb,
                                          D0=D0,
                                          C0=1.0,
                                          E=2,
                                          Cpl=1.8,
                                          spawn_class=spawnclass)
  # Defintion of Material that will be "created" after the nearfiled dynamics of the plume
  # ***no need to define a P0
  material_farfield=BuoyantTracer(id='material_farfield3',
                                  nbuff=10000,
                                  movers=[current],
                                  reactors=[temp,salt],
                                  stickers=[depth_lyt_hbr,shoreline],
                                  reln=reln,
                                  w0=w0,
                                  unstick=0.0,
                                  maxage=3.0,
                                  tstart=tstart,
                                  tend=tend,
                                  ischild=1)


  ercore=ERcore(geod=True,tout=900,outpath='./outputs')
  ercore.materials=[buoyantplume_jetlag,material_farfield]
  ercore.run(tstart,tend,tstep)



#
if test==4:
  # disposal type simulation - look for TASS addition ?
  # Using Lyttelton Harbour flow model as test case
  # near surface release of a dense sediment mixture - heavier than water  
  # SITE------------------------------------------------------------------------------------------------------
  site_num=1
    #site=sys.argv[1]
  # PLUME PARAMETERS
  d=1.0 #nozzle diameter in meter
  flux=1.0 ##Volume flux liquid m3/s
  area=numpy.pi*(d/2)**2
  V0=flux/area #Initial jet velocity <list> of u,v,w components
  Z0=-1.0 # release 1m below surface - i.e. hull of dredger
  # angle between the x-axis and jet in the x-y plane (i.e. horizontal plane) [deg]
  theta_jet = 0 
  # angle between the x-axis and jet in the x-z plane (i.e vertical plane) [deg]
  phi_jet = -90 # going down - do a sketch
  #
  theta_jet = numpy.deg2rad(theta_jet);
  phi_jet   = numpy.deg2rad(phi_jet)

  uj0=V0*numpy.cos(phi_jet)*numpy.cos(theta_jet)
  vj0=V0*numpy.cos(phi_jet)*numpy.sin(theta_jet)
  wj0=V0*numpy.sin(phi_jet)

  # PLUME INFOS 
  B0=d/2 #: Initial jet radius <float>
  V0=numpy.array([uj0,vj0,wj0])#: Initial jet velocity <list> of u,v,w components
  Vb=0.0#: Bubble/dropplet terminal velocity (m/s) <float>
  C0=1.0#: Jet initial concentration <float>
  E=2.0#: Entrainment constant <float>
  # BUOYANT PLUME INFOS
  T0=12#: Initial temperature (C) of jet <float>
  S0=35#: Initial salinity (PSU) of jet <float>
  D0=2000# Initial density (kg/m^3) of jet <float>
  Cpl=1.8#: Specific heat capacity <float>


  spawnclass='material_farfield4'
  site_name=['C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C2_2']
  x=[172.7211,172.7388,172.7555,172.7690,172.7829,172.7972,172.8104,172.8227,172.8361,172.8478,172.8599,172.7383]
  y=[-43.6131,-43.6133,-43.6090,-43.6073,-43.6056,-43.6039,-43.5989,-43.5939,-43.5885,-43.5838,-43.5790,-43.6111]
  
  tstart=datetime.datetime(2014,1,18,8,0,0)#tstart simulation
  tend=tstart+datetime.timedelta(days=2.)
  tstep  = 900.
  tout   = 900.

  with_resusp=0
  # BOUNDARY--------------------------------------------------------------------------------------------------
  bnd=Boundary('bnd',[172.6, 173.92, -43.7, -43.5])
  shoreline=Shoreline('shoreline','../resuspension/southnz_shoreline.bnd')
  # FALL VELOCITY AND TIMESTEPS-------------------------------------------------------------------------------
  w0=-1.0e-3
  # SOURCE DATA FILE
  DataFileLyt='/home/simon/0201_ercore_lyttelton/NEW_RUNS_DISPO/channel_runs/lyt_Original_cons.nc'
  #DEPTH-----------------------------------------------------------------
  depth_lyt_hbr=GriddedTopo('depth_lyt',['dep'],file=DataFileLyt,zinvert=True)
  #DIFFUSION-----------------------------------------------------------
  diff=ConstantDiffuser('diff',['diffx','diffy','diffz'],diffx=0.1,diffy=0.1,diffz=0.0000)
  # TIDE DATA----------------------------------------------------------
  current=TidalMover('lyt_tide_hbr',['u','v'],file=DataFileLyt,topo=depth_lyt_hbr,zinvert=True,z0=0.0001)
  # can add a constant mover to account for dredger speed
  # current=ConstantMover('vel',['u','v','w'],u=+0.2,v=0,w=0,topo=topo)

  #-----------------------------------------------------------------------
  #RELEASE-----------------------------------------------------------
  reln=1*48*3600/tstep #nb part released over the entire run period
  ## MATERIALS
  # Defintion of plume Material
  buoyantplume_jetlag=BuoyantPlume_JETLAG(id='plume4',
                                          nbuff=1e5,
                                          is3d=True,
                                          movers=[current],
                                          reactors=[temp,salt],
                                          V0=V0,
                                          P0=numpy.array([x[site_num-1],y[site_num-1],Z0]),
                                          T0=T0,
                                          S0=S0,
                                          B0=B0,
                                          Vb=Vb,
                                          D0=D0,
                                          C0=C0,
                                          E=2,
                                          Cpl=1.8,
                                          spawn_class=spawnclass)
  # Defintion of Material that will be "created" after the nearfiled dynamics of the plume
  # ***no need to define a P0
  material_farfield=BuoyantTracer(id='material_farfield4',
                                  nbuff=10000,
                                  movers=[current],
                                  reactors=[temp,salt],
                                  stickers=[depth_lyt_hbr,shoreline],
                                  reln=reln,
                                  w0=w0,
                                  unstick=0.0,
                                  maxage=3.0,
                                  tstart=tstart,
                                  tend=tend,
                                  ischild=1)


  ercore=ERcore(geod=True,tout=900,outpath='./outputs')
  ercore.materials=[buoyantplume_jetlag,material_farfield]
  ercore.run(tstart,tend,tstep)


if test==5:
    # SAME AS test==1 but with YAML input
  import pdb;pdb.set_trace()
  tstart=733410.
  tend= 733410.+48.*3600./86400.
  dt  = 900.0
  tout   = 900.
  ercore=ERcore(geod=True,tout=tout,outpath='.',save_summary=True)
  ercore.readYAML('plume_test_1.yaml',globals())
  ercore.run(tstart,tend,dt)


if test==6:
  # disposal type simulation - using BuoyantPlume_DensityCurrent which allows modelling
  # the descent of the dense sediment mixture following disposal or overflow, and the creation
  # of a density current following its collapase on the seabed
  # 
  # Using Lyttelton Harbour flow model as test case
  # near surface release of a dense sediment mixture - heavier than water  
  # SITE------------------------------------------------------------------------------------------------------
  site_num=1
    #site=sys.argv[1]
  # PLUME PARAMETERS
  d=10.0 #nozzle diameter in meter
  flux=100.0 ##Volume flux liquid m3/s
  area=numpy.pi*(d/2)**2
  V0=flux/area #Initial jet velocity <list> of u,v,w components
  Z0=-1.0 # release 1m below surface - i.e. hull of dredger
  # angle between the x-axis and jet in the x-y plane (i.e. horizontal plane) [deg]
  theta_jet = 0 
  # angle between the x-axis and jet in the x-z plane (i.e vertical plane) [deg]
  phi_jet = -90 # going down - do a sketch
  #
  theta_jet = numpy.deg2rad(theta_jet);
  phi_jet   = numpy.deg2rad(phi_jet)

  uj0=V0*numpy.cos(phi_jet)*numpy.cos(theta_jet)
  vj0=V0*numpy.cos(phi_jet)*numpy.sin(theta_jet)
  wj0=V0*numpy.sin(phi_jet)

  # PLUME INFOS 
  B0=d/2 #: Initial jet radius <float>
  V0=numpy.array([uj0,vj0,wj0])#: Initial jet velocity <list> of u,v,w components
  Vb=0.0#: Bubble/dropplet terminal velocity (m/s) <float>
  C0=1.0#: Jet initial concentration <float>
  E=2.0#: Entrainment constant <float>
  # BUOYANT PLUME INFOS
  T0=12#: Initial temperature (C) of jet <float>
  S0=35#: Initial salinity (PSU) of jet <float>
  D0=2000# Initial density (kg/m^3) of jet <float>
  Cpl=1.8#: Specific heat capacity <float>


  site_name=['C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C2_2']
  x=[172.7211,172.7388,172.7555,172.7690,172.7829,172.7972,172.8104,172.8227,172.8361,172.8478,172.8599,172.7383]
  y=[-43.6131,-43.6133,-43.6090,-43.6073,-43.6056,-43.6039,-43.5989,-43.5939,-43.5885,-43.5838,-43.5790,-43.6111]
  
  tstart=datetime.datetime(2014,1,18,8,0,0)#tstart simulation
  tend=tstart+datetime.timedelta(days=2.)
  tstep  = 900.
  tout   = 900.

  with_resusp=0
  # BOUNDARY--------------------------------------------------------------------------------------------------
  bnd=Boundary('bnd',[172.6, 173.92, -43.7, -43.5])
  shoreline=Shoreline('shoreline','../resuspension/southnz_shoreline.bnd')
  # FALL VELOCITY AND TIMESTEPS-------------------------------------------------------------------------------
  w0=-0.001
  # SOURCE DATA FILE
  DataFileLyt='/home/simon/0201_ercore_lyttelton/NEW_RUNS_DISPO/channel_runs/lyt_Original_cons.nc'
  #DEPTH-----------------------------------------------------------------
  depth_lyt_hbr=GriddedTopo('depth_lyt',['dep'],file=DataFileLyt,zinvert=True)
  #DIFFUSION-----------------------------------------------------------
  diff=ConstantDiffuser('diff',['diffx','diffy','diffz'],diffx=0.1,diffy=0.1,diffz=0.0000)
  # TIDE DATA----------------------------------------------------------
  current=TidalMover('lyt_tide_hbr',['u','v'],file=DataFileLyt,topo=depth_lyt_hbr,zinvert=True,z0=0.0001)
  # can add a constant mover to account for dredger speed

  # topo=ConstantTopo('depth',['topo'],topo=-15.)
  # current=ConstantMover('vel',['u','v','w'],u=0,v=0,w=0,topo=topo)

  #-----------------------------------------------------------------------
  #RELEASE-----------------------------------------------------------
  reln=1*48*3600/tstep #nb part released over the entire run period
  ## MATERIALS
  # Defintion of plume Material
  buoyantplume_jetlag=BuoyantPlume_DensityCurrent(id='plume6_densitycurrent',
                                          nbuff=1e5,
                                          is3d=True,
                                          movers=[current],
                                          reactors=[temp,salt],
                                          V0=V0,
                                          P0=numpy.array([x[site_num-1],y[site_num-1],Z0]),
                                          T0=T0,
                                          S0=S0,
                                          B0=B0,
                                          Vb=Vb,
                                          D0=2000.,
                                          C0=C0,
                                          E=2,
                                          Cpl=1.8,
                                          w0=0.001,
                                          rho_sed_dry=1900,
                                          tau_crit=0.3,
                                          z0=0.001,
                                          spawn_class='material_farfield6_densitycurrent',
                                          spawn_type='surface',
                                          formulation='tass')
  # Defintion of Material that will be "created" after the nearfiled dynamics of the plume
  # ***no need to define a P0
  material_farfield=BuoyantTracer(id='material_farfield6_densitycurrent',
                                  nbuff=1e5,
                                  movers=[current],
                                  reactors=[temp,salt],
                                  stickers=[depth_lyt_hbr,shoreline],
                                  reln=1000,
                                  w0=w0,
                                  unstick=0.0,
                                  maxage=3.0,
                                  tstart=tstart,
                                  tend=tstart,
                                  ischild=1)


  ercore=ERcore(geod=True,tout=900,outpath='./outputs')
  ercore.materials=[buoyantplume_jetlag,material_farfield]
  # ercore.run(tstart,tstart+datetime.timedelta(hours=6.),tstep)
  ercore.run(tstart,tstart+datetime.timedelta(seconds=901.),tstep)


if test==7:
  # test different spawning type using the test 1 as reference
  # test the BuoyantPlume_JETLAG turning into a Buoyant Tracer class in constant flows
  # import pdb;pdb.set_trace()
  topo=ConstantTopo('depth',['topo'],topo=-10.)
  current=ConstantMover('vel',['u','v','w'],u=+0.2,v=0,w=0,topo=topo)
  # diff=ConstantDiffuser('diff',['diffx','diffy','diffz'],diffx=0.2,diffy=0.2,diffz=0.001)
  # SPAWN CLASS INFO
  spawnclass='material_farfield1'
  reln=48*3600/900 #nb part released over the entire run period > 1 per release for now

  # Defintion of plume Material
  buoyantplume_jetlag_center=BuoyantPlume_JETLAG(id='plume1_spawn_center',
                                          nbuff=1e5,
                                          is3d=True,
                                          movers=[current],
                                          reactors=[temp,salt],
                                          V0=V0,
                                          P0=numpy.array([0,0,Z0]),
                                          T0=80,
                                          S0=10,
                                          B0=B0,
                                          Vb=Vb,
                                          D0=D0,
                                          C0=1.0,
                                          E=2,
                                          Cpl=1.8,
                                          spawn_class='material_farfield1_center',
                                          spawn_type='center')

  buoyantplume_jetlag_surface=BuoyantPlume_JETLAG(id='plume1_spawn_surface',
                                          nbuff=1e5,
                                          is3d=True,
                                          movers=[current],
                                          reactors=[temp,salt],
                                          V0=V0,
                                          P0=numpy.array([0,0,Z0]),
                                          T0=80,
                                          S0=10,
                                          B0=B0,
                                          Vb=Vb,
                                          D0=D0,
                                          C0=1.0,
                                          E=2,
                                          Cpl=1.8,
                                          spawn_class='material_farfield1_surface',
                                          spawn_type='surface')

  buoyantplume_jetlag_cylinder=BuoyantPlume_JETLAG(id='plume1_spawn_cylinder',
                                          nbuff=1e5,
                                          is3d=True,
                                          movers=[current],
                                          reactors=[temp,salt],
                                          V0=V0,
                                          P0=numpy.array([0,0,Z0]),
                                          T0=80,
                                          S0=10,
                                          B0=B0,
                                          Vb=Vb,
                                          D0=D0,
                                          C0=1.0,
                                          E=2,
                                          Cpl=1.8,
                                          spawn_class='material_farfield1_cylinder',
                                          spawn_type='cylinder')
  # only change id and spawn details



  # Definition of Material that will be "created" after the nearfiled dynamics of the plume
  # ***no need to define a P0
  material_farfield_center=BuoyantTracer(id='material_farfield1_center',
                                  nbuff=10000,
                                  movers=[current],
                                  reactors=[temp,salt],
                                  reln=1000,
                                  w0=-0.001,
                                  unstick=0.0,
                                  maxage=3.0,
                                  tstart=733410.,
                                  tend=733410.,
                                  ischild=1)
  material_farfield_surface=BuoyantTracer(id='material_farfield1_surface',
                                  nbuff=10000,
                                  movers=[current],
                                  reactors=[temp,salt],
                                  reln=1000,
                                  w0=-0.001,
                                  unstick=0.0,
                                  maxage=3.0,
                                  tstart=733410.,
                                  tend=733410.,
                                  ischild=1)
  material_farfield_cylinder=BuoyantTracer(id='material_farfield1_cylinder',
                                  nbuff=10000,
                                  movers=[current],
                                  reactors=[temp,salt],
                                  reln=1000,
                                  w0=-0.001,
                                  unstick=0.0,
                                  maxage=3.0,
                                  tstart=733410.,
                                  tend=733410.,
                                  ischild=1)


  # Run

  import pdb;pdb.set_trace()

  ercore=ERcore(geod=True,tout=900,outpath='./outputs')
  ercore.materials=[buoyantplume_jetlag_center,material_farfield_center]
  # ercore.run(tstart,tend,pdt)
  ercore.run(733410.-900/86400,733410.+2.*3600./86400.,900.0)  
 # Run
  ercore=ERcore(geod=True,tout=900,outpath='./outputs')
  ercore.materials=[buoyantplume_jetlag_surface,material_farfield_surface]
  # ercore.run(tstart,tend,pdt)
  ercore.run(733410.-900/86400,733410.+2.*3600./86400.,900.0)  

 # # Run
  ercore=ERcore(geod=True,tout=900,outpath='./outputs')
  ercore.materials=[buoyantplume_jetlag_cylinder,material_farfield_cylinder]
  # ercore.run(tstart,tend,pdt)
  ercore.run(733410.-900/86400,733410.+2.*3600./86400.,900.0)  
















# def run(obj):
#   t0=time.time()
#   obj.run(733410,733410+48*3600.0/86400.,900.0)
#   print time.time()-t0
# run(ercore) 

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


  
