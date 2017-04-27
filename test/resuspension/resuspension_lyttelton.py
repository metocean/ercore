#!/usr/bin/env python
import numpy,datetime,sys,os,time
#sys.path.append(os.path.join('..'))
from ercore import ERcore,dt2ncep
from ercore.fields import *
from ercore.materials import BuoyantTracer
from ercore.shoreline import Boundary
from ercore.shoreline import Shoreline
from ercore.materials.sediment import Sediment

# SITE------------------------------------------------------------------------------------------------------
site_name=['C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C2_2']
x=[172.7211,172.7388,172.7555,172.7690,172.7829,172.7972,172.8104,172.8227,172.8361,172.8478,172.8599,172.7383]
y=[-43.6131,-43.6133,-43.6090,-43.6073,-43.6056,-43.6039,-43.5989,-43.5939,-43.5885,-43.5838,-43.5790,-43.6111]
#site=sys.argv[1]
site=int(sys.argv[1])
with_resusp=int(sys.argv[2])
# BOUNDARY--------------------------------------------------------------------------------------------------
bnd=Boundary('bnd',[172.6, 173.92, -43.7, -43.5])
shoreline=Shoreline('shoreline','southnz_shoreline.bnd')
# FALL VELOCITY AND TIMESTEPS-------------------------------------------------------------------------------
pw0=[1.0e-3,1.3659e-3,2.7559e-3,8.493e-3]
grainsize=['class1','class2','class3','class4']
pw0=[1.0e-3]
grainsize=['class1']
# SOURCE DATA FILE
DataFileLyt='/home/simon/0201_ercore_lyttelton/NEW_RUNS_DISPO/channel_runs/lyt_Original_cons.nc'
#DEPTH-----------------------------------------------------------------
depth_lyt_hbr=GriddedTopo('depth_lyt',['dep'],file=DataFileLyt,zinvert=True)
Z0=depth_lyt_hbr.interp(numpy.array([[x[site-1],y[site-1],0]]))

#DIFFUSION-----------------------------------------------------------
diff=ConstantDiffuser('diff',['diffx','diffy','diffz'],diffx=0.0,diffy=0.0,diffz=0.0000)
# based on ELDER using depth and mean speed at site Klg~0.2-0.3 Klat ~0.02-0.0.3 > mean ~0.1
#Note : to have effective Z diffusion 10 times smaller than diff XY the coeff  has to be 10e2 smaller 
# TIDE DATA----------------------------------------------------------
#lyt_flow_hbr=GriddedTide('lyt_tide_hbr',['u','v'],file=DataFileLyt,topo=depth_lyt_hbr,zinvert=True,z0=0.001)
lyt_flow_hbr=TidalMover('lyt_tide_hbr',['u','v'],file=DataFileLyt,topo=depth_lyt_hbr,zinvert=True,z0=0.0001)
#-----------------------------------------------------------------------
ercore=ERcore(geod=True,tout=360.)
#RELEASE-----------------------------------------------------------
reln=[1200] #nb part release over the entire run period > 1 month was 300,000
pdt=[360,360,180,60]


Z=[[0,Z0],[Z0+2,Z0]]
zlabel=['random_watercolumn','random_bottom_cylinder_r100']
#neap=19-Jan-2014 09:57:19
#spring 03-Jan-2014 15:48:00
# RUN LOOP----------------------------------------------------------
for ip,w0 in enumerate(pw0):
    for iz,zlev in enumerate(Z):
        # tstart=datetime.datetime(2014,1,1,0,0,0)#tstart simulation
        # tend=tstart+datetime.timedelta(days=30) # allow 2 days so that some particle are released already
        tstart=datetime.datetime(2014,1,18,8,0,0)#tstart simulation
        tstartp=tstart+datetime.timedelta(hours=0)# tstart release- 2 days each
        tend=tstart+datetime.timedelta(days=2) # allow 2 days so that some particle are released already
        tend=tstart+datetime.timedelta(hours=14)
        tendp=tend#tend release
        if with_resusp:
          if iz==1: 
            sediment=Sediment('lyt_test_site%s_%s_%s' %         (site,zlabel[iz],grainsize[ip]),5000,movers=[lyt_flow_hbr],diffusers=[diff],stickers=[depth_lyt_hbr,shoreline],reln=reln[0],P0=[x[site-1],y[site-1],Z[iz]],w0=-pw0[ip],tstart=tstartp,tend=tendp,tau_crit_eros=0.2,tau_crit_depos=1000.0,unstick=0.0,maxage=3.0)
          else:
            sediment=Sediment('lyt_test_site%s_%s_%s' %         (site,zlabel[iz],grainsize[ip]),5000,movers=[lyt_flow_hbr],diffusers=[diff],stickers=[depth_lyt_hbr,shoreline],reln=reln[0],P0=[x[site-1],y[site-1],Z[iz]],w0=-pw0[ip],tstart=tstartp,tend=tendp,tau_crit_eros=0.2,tau_crit_depos=1000.0,unstick=0.0,maxage=3.0,circular_radius=100.00)
        else:
          if iz==1: 
            sediment=BuoyantTracer('lyt_test_site%s_%s_%s' %         (site,zlabel[iz],grainsize[ip]),5000,movers=[lyt_flow_hbr],diffusers=[diff],stickers=[depth_lyt_hbr,shoreline],reln=reln[0],P0=[x[site-1],y[site-1],Z[iz]],w0=-pw0[ip],tstart=tstartp,tend=tendp,unstick=0.0,maxage=3.0,circular_radius=100.00)
          else:
            sediment=BuoyantTracer('lyt_test_site%s_%s_%s' %         (site,zlabel[iz],grainsize[ip]),5000,movers=[lyt_flow_hbr],diffusers=[diff],stickers=[depth_lyt_hbr,shoreline],reln=reln[0],P0=[x[site-1],y[site-1],Z[iz]],w0=-pw0[ip],tstart=tstartp,tend=tendp,unstick=0.0,maxage=3.0)
          
        ercore.materials=[sediment]

        ercore.run(tstart,tend,pdt[ip])