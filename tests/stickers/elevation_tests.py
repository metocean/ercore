#!/usr/bin/env python
import sys; sys.path += ["/source/ercore"]


import numpy,datetime
from ercore import ERcore
from ercore.fields import ConstantMover,ConstantElevation
from ercore.materials import BuoyantTracer,PassiveTracer
from ercore import dt2ncep


############################
## Rising moving particle ##
############################

t1 = datetime.datetime(2000,1,1)
t2 = datetime.datetime(2000,1,1,1)

h0 = 1.1

current = ConstantMover('cur',['u','v'],u=1./3600., v=0)
elev = ConstantElevation('elev', ['elev'], elev=h0)

# rising particle
p1 = BuoyantTracer('p1', nbuff=1000,geod=False,
                    movers=[current], stickers=[elev],unstick=1,
                    tstart=t1,tend=t1, tstep=0., 
                    w0=0.5/3600.,
                    reln=2,P0=[0,0,h0-2],outfile='elev1.out')

# released at surface
p2 = PassiveTracer('p2', nbuff=1000,geod=False,
                    movers=[current], stickers=[elev],unstick=1,
                    tstart=t1,tend=t1, tstep=0., 
                    reln=2,P0=[0,0,h0],outfile='elev2.out')


# sinking particle, released above surface
p3 = BuoyantTracer('p3', nbuff=1000,geod=False,
                    movers=[current], stickers=[elev],unstick=1,
                    tstart=t1,tend=t1, tstep=0., 
                    w0=-0.5/3600.,
                    reln=2,P0=[0,0,h0+2],outfile='elev3.out')

ercore=ERcore(geod=False)
ercore.materials=[p1,p2,p3]
ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)



