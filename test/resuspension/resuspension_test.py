#!/usr/bin/env python
import datetime,sys,os
#sys.path.insert(0,os.path.join(os.path.dirname(__file__),'..','..'))
#sys.path.insert(1,'/home/rosa/git/ercore/')
from pylab import *
from ercore._flib_ercore import interp3d,interph,slipvel
from ercore import ERcore,ncep2dt,dt2ncep
from ercore.materials import PassiveTracer,BuoyantTracer
from ercore.materials.sediment import Sediment
from ercore.fields import GriddedMover,ConstantDiffuser,VariableDiffuser,GriddedTopo,ConstantMover,GriddedDataGroup
from ercore.shoreline import Shoreline

P0=[170.5,-46,-10]
#import pdb;pdb.set_trace()
dep=GriddedTopo('depth',['dep'],file='../passive/uds_gsb_test.nc',zinvert=True) # here we set zinvert=true because depth are positive downward while ercore conventon is negative downward
currents=GriddedMover('cur',['uo','vo'],is3d=True,file='../passive/uds_gsb_test.nc',topo=dep,z0=0.001,zcoord='down') 
# in roms/poms etc.. levels are given positive down i.e. 5 is below surface, this needs to be inverted in ercore to fit convention
# tide=GriddedMover('tide',['ut','vt'],is3d=False,file='../passive/uds_gsb_test.nc') #
#shoreline=Shoreline('shoreline','../passive/shoreline2.bnd')
shoreline=Shoreline('shoreline','otago_shoreline2.bnd')
diff=ConstantDiffuser('diff',['diffx','diffy','diffz'],diffx=0.1,diffy=0.1,diffz=1.0e-6)

#import pdb; pdb.set_trace()

tstart = datetime.datetime(2009,1,1)
tend   = datetime.datetime(2009,1,3)
tstep  = 900.
tout   = 900.
tau_crit_eros = 0.1  # for now need to be explicitely specified - could be included in the function based on d50
tau_crit_depos = 1000.0  # for now need to be explicitely specified - could be included in the function based on d50
n_part=2

# CASE 1 : discrete release - low critical shear stress for erosion so that resuspension always occurs
particles0=Sediment('particles0',10000,movers=[currents],stickers=[dep,shoreline],diffusers=[],reln=n_part,P0=P0,w0=-10e-3,tstart=tstart,tend=tstart,tau_crit_eros=.01,unstick=0.0)
# CASE 2 : discrete release - high critical shear stress for erosion so that resuspension never occurs
particles1=Sediment('particles1',10000,movers=[currents],stickers=[dep,shoreline],diffusers=[],reln=n_part,P0=P0,w0=-10e-3,tstart=tstart,tend=tstart,tau_crit_eros=10.0,unstick=0.0)
# CASE 3 : continuous release - high critical shear stress for erosion so that resuspension never occurs
tstart = datetime.datetime(2009,1,5)
tend   = datetime.datetime(2009,1,15)
P0=[170.5,-46,-1.0]
particles2=Sediment('particles2',10000,movers=[currents],stickers=[dep,shoreline],diffusers=[diff],reln=480*2,P0=P0,w0=-1e-3,tstart=tstart,tend=tend,tau_crit_eros=10.0,unstick=0.0)
# CASE 4 : continuous release - low critical shear stress for erosion so that resuspension always occurs
particles3=Sediment('particles3',10000,movers=[currents],stickers=[dep,shoreline],diffusers=[diff],reln=480*2,P0=P0,w0=-1e-3,tstart=tstart,tend=tend,tau_crit_eros=0.15,unstick=0.0)

ercore=ERcore(tout=tout,rkorder=4)

###### start tests #####
currents.reset()
ercore.materials=[particles2]
ercore.run(tstart,tend,tstep)

currents.reset()
ercore.materials=[particles3]
ercore.run(tstart,tend,tstep)
import pdb;pdb.set_trace()

currents.reset()
ercore.materials=[particles1]
ercore.run(tstart,tend,tstep)

currents.reset()
ercore.materials=[particles2]
ercore.run(tstart,tend,tstep)

