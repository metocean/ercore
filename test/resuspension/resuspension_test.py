#!/usr/bin/env python
import datetime,sys,os
#sys.path.insert(0,os.path.join(os.path.dirname(__file__),'..','..'))
#sys.path.insert(1,'/home/rosa/git/ercore/')
from pylab import *
from ercore._flib_ercore import interp3d,interph,slipvel
from ercore import ERcore,ncep2dt,dt2ncep
from ercore.materials import PassiveTracer,BuoyantTracer
from ercore.fields import GriddedMover,ConstantDiffuser,VariableDiffuser,GriddedTopo,ConstantMover,GriddedDataGroup
from ercore.shoreline import Shoreline

P0=[170.5,-46,-10]
#import pdb;pdb.set_trace()
dep=GriddedTopo('depth',['dep'],file='../passive/uds_gsb_test.nc',zinvert=True)
currents=GriddedMover('cur',['uo','vo'],is3d=True,file='../passive/uds_gsb_test.nc',topo=dep)
tide=GriddedMover('tide',['ut','vt'],is3d=False,file='../passive/uds_gsb_test.nc') #
shoreline=Shoreline('shoreline','../passive/shoreline2.bnd')


tstart = datetime.datetime(2009,1,1)
tend   = datetime.datetime(2009,1,3)
tstep  = 900.
tout   = 3600.

particles=BuoyantTracer('particles',10000,movers=[currents,tide],stickers=[shoreline,dep],reln=10,P0=P0,w0=-1e-3,tstart=tstart,tend=tstart)

ercore=ERcore(tout=tout,rkorder=4)

###### start tests #####

ercore.materials=[particles]
ercore.run(tstart,tend,tstep)
