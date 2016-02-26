#!/usr/bin/env python
import sys; sys.path += ["/source/ercore"]

import numpy,datetime
# import cdms2
from ercore import ERcore
from ercore.fields import ConstantMover
from ercore.shoreline import Shoreline
from ercore.materials import PassiveTracer
# from ercore import dt2ncep

t1 = datetime.datetime(2000,1,1)
t2 = datetime.datetime(2000,1,1,1)

#######################
## Just bounding box ##
#######################

fileshore = 'bbox.bnd'
f = open(fileshore, 'w+')
f.write('"Map Bounds","1",4\n')
f.write('%i,%i\n' % (0,0))
f.write('%i,%i\n' % (0,3))
f.write('%i,%i\n' % (3,3))
f.write('%i,%i\n' % (3,0))
f.close()

shore = Shoreline(id='shore', file=fileshore)

# ## released WEST  ##
# current = ConstantMover('cur',['uo','vo'],uo=0.5/3600.,vo=0.0)

# # unstick
# p1 = PassiveTracer('p1', nbuff=1000,geod=False,
#                     movers=[current], stickers=[shore],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0,0,0],outfile='bbox_W1.out')


# p2 = PassiveTracer('p2', nbuff=1000,geod=False,
#                     movers=[current], stickers=[shore],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0.1,0,0],outfile='bbox_W2.out')

# # stick
# p1s = PassiveTracer('p1s', nbuff=1000,geod=False,
#                     movers=[current], stickers=[shore],unstick=0,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0,0,0],outfile='bbox_W1_stick.out')

# p2s = PassiveTracer('p2s', nbuff=1000,geod=False,
#                     movers=[current], stickers=[shore],unstick=0,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0.1,0,0],outfile='bbox_W2_stick.out')


# ercore=ERcore(geod=False)
# ercore.materials=[p1,p2,p1s,p2s]
# ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)


# ## released EAST side ##

# current = ConstantMover('cur',['uo','vo'],uo=-0.5/3600.,vo=0.0)

# # unstick
# p1 = PassiveTracer('p1', nbuff=1000,geod=False,
#                     movers=[current], stickers=[shore],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[3,0,0],outfile='bbox_E1.out')


# p2 = PassiveTracer('p2', nbuff=1000,geod=False,
#                     movers=[current], stickers=[shore],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[2.9,0,0],outfile='bbox_E2.out')


# ercore=ERcore(geod=False)
# ercore.materials=[p1,p2] #,p1s,p2s]
# ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)

# ## released SOUTH side ##

# current = ConstantMover('cur',['uo','vo'],uo=0, vo=0.5/3600.)

# # unstick
# p1 = PassiveTracer('p1', nbuff=1000,geod=False,
#                     movers=[current], stickers=[shore],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0,0,0],outfile='bbox_S1.out')


# p2 = PassiveTracer('p2', nbuff=1000,geod=False,
#                     movers=[current], stickers=[shore],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0,0.1,0],outfile='bbox_S2.out')


# ercore=ERcore(geod=False)
# ercore.materials=[p1,p2] #,p1s,p2s]
# ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)

# ## released NORTH side ##

current = ConstantMover('cur',['uo','vo'],uo=0, vo=-0.5/3600.)

# unstick
p1 = PassiveTracer('p1', nbuff=1000,geod=False,
                    movers=[current], stickers=[shore],unstick=1,
                    tstart=t1,tend=t1, tstep=0., 
                    reln=2,P0=[0,3,0],outfile='bbox_N1.out')


p2 = PassiveTracer('p2', nbuff=1000,geod=False,
                    movers=[current], stickers=[shore],unstick=1,
                    tstart=t1,tend=t1, tstep=0., 
                    reln=2,P0=[0,2.9,0],outfile='bbox_N2.out')

# Start OUTSIDE of BBOX
p3 = PassiveTracer('p3', nbuff=1000,geod=False,
                    movers=[current], stickers=[shore],unstick=1,
                    tstart=t1,tend=t1, tstep=0., 
                    reln=2,P0=[0,3.2,0],outfile='bbox_N3.out')


ercore=ERcore(geod=False)
ercore.materials=[p1,p2,p3]
ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)



############################
## 3x3 wih shoreline on 2 ##
############################

fileshore = 'shore_x2.bnd'
f = open(fileshore, 'w+')
f.write('"Map Bounds","1",4\n')
f.write('%i,%i\n' % (0,0))
f.write('%i,%i\n' % (0,3))
f.write('%i,%i\n' % (3,3))
f.write('%i,%i\n' % (3,0))
f.write('"1","1",5\n')
f.write('%.1f,%i\n' % (1.9,0))
f.write('%.1f,%i\n' % (1.9,3))
f.write('%i,%i\n' % (3,3))
f.write('%i,%i\n' % (3,0))
f.write('%.1f,%i\n' % (1.9,0))
f.close()

shore = Shoreline(id='shore', file=fileshore)
current = ConstantMover('cur',['uo','vo'],uo=0.5/3600.,vo=0.0)

# unstick
p1 = PassiveTracer('p1', nbuff=1000,geod=False,
                    movers=[current], stickers=[shore],unstick=1,
                    tstart=t1,tend=t1, tstep=0., 
                    reln=2,P0=[0,0,0],outfile='shore_x2.out')

# stick
p1s = PassiveTracer('p2', nbuff=1000,geod=False,
                    movers=[current], stickers=[shore],unstick=0,
                    tstart=t1,tend=t1, tstep=0., 
                    reln=2,P0=[0,0,0],outfile='shore_x2_stick.out')


ercore=ERcore(geod=False)
ercore.materials=[p1,p1s]
ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)

# ###########
# ## Bride ##
# ###########


# #   0   1   2   3   4   5   6
# # 6 -------------------------
# # 5         |   |
# # 4         -----
# # 3
# # 2         -----
# # 1         |   |
# # 0 --------------------------

# if test < 3:


#     fileshore = 'bridge.bnd'
#     f = open(fileshore, 'w+')
#     f.write('"Map Bounds","1",4\n')
#     f.write('%i,%i\n' % (0,0))
#     f.write('%i,%i\n' % (0,6))
#     f.write('%i,%i\n' % (6,6))
#     f.write('%i,%i\n' % (6,0))
#     f.write('"0","1",5\n')  # bottom
#     f.write('%.1f,%.1f\n' % (2,0))
#     f.write('%.1f,%.1f\n' % (3,0))
#     f.write('%.1f,%.1f\n' % (3,2))
#     f.write('%.1f,%.1f\n' % (2,2))
#     f.write('%.1f,%.1f\n' % (2,0))
#     f.write('"1","1",5\n')
#     f.write('%.1f,%.1f\n' % (2,4))
#     f.write('%.1f,%.1f\n' % (3,4))
#     f.write('%.1f,%.1f\n' % (3,6))
#     f.write('%.1f,%.1f\n' % (2,6))
#     f.write('%.1f,%.1f\n' % (2,4))
#     f.close()

#     # sticks to bridge
#     # but passes bridge if velocity to high (uo=2/3600)
#     # doesn't pass if too too big (uo=1/3600, dt=7200)

#     shore = Shoreline(id='shore', file=fileshore)
#     t1 = datetime.datetime(2000,1,1)
#     t2 = datetime.datetime(2000,1,1,1)
#     current = ConstantMover('cur',['uo','vo'],uo=1./3600.,vo=0.0)
#     bottom = PassiveTracer('top', nbuff=1000,geod=False,
#                                 movers=[current], stickers=[shore],unstick=0,
#                                 tstart=t1,tend=t1, tstep=0.,
#                                 reln=2,P0=[0,1,3],outfile='bridge_bottom.out')

#     mid = PassiveTracer('bombs', nbuff=1000,geod=False,
#                                 movers=[current], stickers=[shore],unstick=1,
#                                 tstart=t1,tend=t1, tstep=0.,
#                                 reln=2,P0=[0,3,3],outfile='bridge_middle.out')

#     top = PassiveTracer('bottom', nbuff=1000,geod=False,
#                                 movers=[current], stickers=[shore],unstick=1,
#                                 tstart=t1,tend=t1, tstep=0.,
#                                 reln=2,P0=[0,5,3],outfile='bridge_top.out')

#     ercore=ERcore(geod=False)
#     ercore.materials=[top,mid,bottom]
#     ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)
