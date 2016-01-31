#!/usr/bin/env python
import numpy,datetime
import cdms2
from ercore import ERcore
from ercore.fields import ConstantMover,GriddedTopo,GriddedMover
from ercore.shoreline import Shoreline
from ercore.materials import PassiveTracer
from ercore import dt2ncep

#######################
## Just bounding box ##
#######################

# fileshore = 'shore_map.bnd'
# f = open(fileshore, 'w+')
# f.write('"Map Bounds","1",4\n')
# f.write('%i,%i\n' % (0,0))
# f.write('%i,%i\n' % (0,3))
# f.write('%i,%i\n' % (3,3))
# f.write('%i,%i\n' % (3,0))
# f.close()
#
# shore = Shoreline(id='shore_map', file=fileshore)
# t1 = datetime.datetime(2000,1,1)
# t2 = datetime.datetime(2000,1,1,1)
# current = ConstantMover('cur',['uo','vo'],uo=1./3600.,vo=0.0)
# bombs = PassiveTracer('bombs', nbuff=1000,geod=False,
#                             movers=[current], stickers=[shore],unstick=1,
#                             tstart=t1,tend=t1, tstep=0.,
#                             reln=2,P0=[0,0,0],outfile='shore_map.out')
# ercore=ERcore(geod=False)
# ercore.materials=[bombs]
# ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)
#
# current = ConstantMover('cur',['uo','vo'],uo=0,vo=-1./3600.)
# bombs = PassiveTracer('bombs', nbuff=1000,geod=False,
#                             movers=[current], stickers=[shore],unstick=1,
#                             tstart=t1,tend=t1, tstep=0.,
#                             reln=2,P0=[0,3,0],outfile='shore_map2.out')
# ercore=ERcore(geod=False)
# ercore.materials=[bombs]
# ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)

############################
## 3x3 wih shoreline on 2 ##
############################

# fileshore = 'shore_x2.bnd'
# f = open(fileshore, 'w+')
# f.write('"Map Bounds","1",4\n')
# f.write('%i,%i\n' % (0,0))
# f.write('%i,%i\n' % (0,3))
# f.write('%i,%i\n' % (3,3))
# f.write('%i,%i\n' % (3,0))
# f.write('"1","1",5\n')
# f.write('%.1f,%i\n' % (1.9,0))
# f.write('%.1f,%i\n' % (1.9,3))
# f.write('%i,%i\n' % (3,3))
# f.write('%i,%i\n' % (3,0))
# f.write('%.1f,%i\n' % (1.9,0))
# f.close()
#
# shore = Shoreline(id='shore', file=fileshore)
# t1 = datetime.datetime(2000,1,1)
# t2 = datetime.datetime(2000,1,1,1)
# current = ConstantMover('cur',['uo','vo'],uo=1./3600.,vo=0.0)
# bombs = PassiveTracer('bombs', nbuff=1000,geod=False,
#                             movers=[current], stickers=[shore],unstick=1,
#                             tstart=t1,tend=t1, tstep=0.,
#                             reln=2,P0=[0,1,3],outfile='shore_x2.out')
# ercore=ERcore(geod=False)
# ercore.materials=[bombs]
# ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)
#

###########
## Bride ##
###########

#   0   1   2   3   4   5   6
# 6 -------------------------
# 5         |   |
# 4         -----
# 3
# 2         -----
# 1         |   |
# 0 --------------------------

fileshore = 'bridge.bnd'
f = open(fileshore, 'w+')
f.write('"Map Bounds","1",4\n')
f.write('%i,%i\n' % (0,0))
f.write('%i,%i\n' % (0,6))
f.write('%i,%i\n' % (6,6))
f.write('%i,%i\n' % (6,0))
f.write('"0","1",5\n')  # bottom
f.write('%.1f,%.1f\n' % (2,0))
f.write('%.1f,%.1f\n' % (3,0))
f.write('%.1f,%.1f\n' % (3,2))
f.write('%.1f,%.1f\n' % (2,2))
f.write('%.1f,%.1f\n' % (2,0))
f.write('"1","1",5\n')
f.write('%.1f,%.1f\n' % (2,4))
f.write('%.1f,%.1f\n' % (3,4))
f.write('%.1f,%.1f\n' % (3,6))
f.write('%.1f,%.1f\n' % (2,6))
f.write('%.1f,%.1f\n' % (2,4))
f.close()

# sticks to bridge
# but passes bridge if velocity to high (uo=2/3600)
# doesn't pass if too too big (uo=1/3600, dt=7200)

shore = Shoreline(id='shore', file=fileshore)
t1 = datetime.datetime(2000,1,1)
t2 = datetime.datetime(2000,1,1,1)
current = ConstantMover('cur',['uo','vo'],uo=2./3600.,vo=0.0)
bottom = PassiveTracer('top', nbuff=1000,geod=False,
                            movers=[current], stickers=[shore],unstick=1,
                            tstart=t1,tend=t1, tstep=0.,
                            reln=2,P0=[0,1,3],outfile='bridge_bottom.out')

mid = PassiveTracer('bombs', nbuff=1000,geod=False,
                            movers=[current], stickers=[shore],unstick=1,
                            tstart=t1,tend=t1, tstep=0.,
                            reln=2,P0=[0,3,3],outfile='bridge_middle.out')

top = PassiveTracer('bottom', nbuff=1000,geod=False,
                            movers=[current], stickers=[shore],unstick=1,
                            tstart=t1,tend=t1, tstep=0.,
                            reln=2,P0=[0,5,3],outfile='bridge_top.out')

ercore=ERcore(geod=False)
ercore.materials=[top,mid,bottom]
ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)
