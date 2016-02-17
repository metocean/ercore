#!/usr/bin/env python
import numpy,datetime
from ercore import ERcore
from ercore.fields import ConstantMover,GriddedTopo,GriddedMover
from ercore.shoreline import Shoreline
from ercore.materials import PassiveTracer
from ercore import dt2ncep

test = 0
######################################
## Just bounding box + depth of 2 m ##
######################################

if test < 1:

    fileshore = 'shore_x2.bnd'
    f = open(fileshore, 'w+')
    f.write('"Map Bounds","1",4\n')
    f.write('%i,%i\n' % (0,0))
    f.write('%i,%i\n' % (0,3))
    f.write('%i,%i\n' % (3,3))
    f.write('%i,%i\n' % (3,0))
    f.write('"1","1",5\n')
    f.write('%.1f,%i\n' % (1.1,0))
    f.write('%.1f,%i\n' % (1.1,3))
    f.write('%i,%i\n' % (3,3))
    f.write('%i,%i\n' % (3,0))
    f.write('%.1f,%i\n' % (1.1,0))
    f.close()

    from utils import make_dep    
    x = numpy.arange(4)
    y = numpy.arange(4)
    xx,yy = numpy.meshgrid(x,y)
    dep = numpy.ones(xx.shape)
    dep[:,:2] = 0.5
    dep[:,2:] = -0.5
    print dep
    filedep = 'shore_dep_map.nc'
    make_dep(x,y,dep,filedep)

    shore = Shoreline(id='shore_map', file=fileshore) #, refloat=1)
    dep=GriddedTopo('depth',['dep'], file=filedep, zinvert=True)

    t1 = datetime.datetime(2000,1,1)
    t2 = datetime.datetime(2000,1,1,1)
    current = ConstantMover('cur',['uo','vo'],uo=1./3600.,vo=0.0)
    bombs = PassiveTracer('bombs', nbuff=1000,geod=False,
                                movers=[current], stickers=[shore,dep],unstick=1,
                                tstart=t1,tend=t1, tstep=0.,
                                reln=1,P0=[0,0,-0.25],outfile='shore_dep_map.out')
    ercore=ERcore(geod=False)
    ercore.materials=[bombs]
    ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)

