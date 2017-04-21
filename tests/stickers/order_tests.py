#!/usr/bin/env python
import sys; sys.path += ["/source/ercore"]

import numpy,datetime
from ercore import ERcore
from ercore.fields import ConstantMover,GriddedTopo
from ercore.materials import PassiveTracer
from ercore.shoreline import Shoreline

t1 = datetime.datetime(2000,1,1)
t2 = datetime.datetime(2000,1,1,1)


from utils import make_dep    
x = numpy.arange(4)
y = numpy.arange(4)
xx,yy = numpy.meshgrid(x,y)
dep = numpy.ones(xx.shape)
dep[:,:2] = 20.
dep[:,2:] = 10.
print dep
filedep = 'dep1.nc'
make_dep(x,y,dep,filedep)

##########################################
## meet the shoreline before the bottom ##
##########################################

fileshore = 'shore1.bnd'
f = open(fileshore, 'w+')
f.write('"Map Bounds","1",4\n')
f.write('%i,%i\n' % (0,0))
f.write('%i,%i\n' % (0,3))
f.write('%i,%i\n' % (3,3))
f.write('%i,%i\n' % (3,0))
f.write('"1","1",5\n')
f.write('%.1f,%i\n' % (0.5,0))
f.write('%.1f,%i\n' % (0.5,3))
f.write('%i,%i\n' % (3,3))
f.write('%i,%i\n' % (3,0))
f.write('%.1f,%i\n' % (0.5,0))
f.close()


dep=GriddedTopo('depth',['dep'], file=filedep, zinvert=True)
shore = Shoreline(id='shore', file=fileshore)
current = ConstantMover('cur',['uo','vo'],uo=0.1/3600.,vo=0.0)


# start above bottom
p1 = PassiveTracer('p1', nbuff=1000,geod=False,
                    movers=[current], stickers=[dep,shore],unstick=1,
                    tstart=t1,tend=t1, tstep=0., 
                    reln=2,P0=[0,2,-16],outfile='order_dep_shore1_p1.out')

# start above bottom
p2 = PassiveTracer('p2', nbuff=1000,geod=False,
                    movers=[current], stickers=[shore,dep],unstick=1,
                    tstart=t1,tend=t1, tstep=0., 
                    reln=2,P0=[0,2,-16],outfile='order_shore1_dep_p1.out')

ercore=ERcore(geod=False)
ercore.materials=[p1,p2]
ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)
