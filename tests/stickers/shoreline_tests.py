#!/usr/bin/env python
import sys; sys.path += ["/source/ercore"]

import numpy,datetime
# import cdms2
from ercore import ERcore
from ercore.fields import ConstantMover, GriddedMover
from ercore.shoreline import Shoreline
from ercore.materials import PassiveTracer
# from ercore import dt2ncep

t1 = datetime.datetime(2000,1,1)
t2 = datetime.datetime(2000,1,1,1)

#######################
## Just bounding box ##
#######################

# fileshore = 'bbox.bnd'
# f = open(fileshore, 'w+')
# f.write('"Map Bounds","1",4\n')
# f.write('%i,%i\n' % (0,0))
# f.write('%i,%i\n' % (0,3))
# f.write('%i,%i\n' % (3,3))
# f.write('%i,%i\n' % (3,0))
# f.close()

# shore = Shoreline(id='shore', file=fileshore)

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

# current = ConstantMover('cur',['uo','vo'],uo=0, vo=-0.5/3600.)

# # unstick
# p1 = PassiveTracer('p1', nbuff=1000,geod=False,
#                     movers=[current], stickers=[shore],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0,3,0],outfile='bbox_N1.out')


# p2 = PassiveTracer('p2', nbuff=1000,geod=False,
#                     movers=[current], stickers=[shore],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0,2.9,0],outfile='bbox_N2.out')

# # Start OUTSIDE of BBOX
# p3 = PassiveTracer('p3', nbuff=1000,geod=False,
#                     movers=[current], stickers=[shore],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0,3.2,0],outfile='bbox_N3.out')


# ercore=ERcore(geod=False)
# ercore.materials=[p1,p2,p3]
# ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)



# ############################
# ## 3x3 wih shoreline on 2 ##
# ############################

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

# shore = Shoreline(id='shore', file=fileshore)
# current = ConstantMover('cur',['uo','vo'],uo=0.5/3600.,vo=0.0)

# # unstick
# p1 = PassiveTracer('p1', nbuff=1000,geod=False,
#                     movers=[current], stickers=[shore],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0,0,0],outfile='shore_x2.out')

# # stick
# p1s = PassiveTracer('p1s', nbuff=1000,geod=False,
#                     movers=[current], stickers=[shore],unstick=0,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0,0,0],outfile='shore_x2_stick.out')


# ercore=ERcore(geod=False)
# ercore.materials=[p1,p1s]
# ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)

# #########################
# ## See if particle     ##
# ## keeps moving        ##
# #########################
# # NOT MOVING AFTER HITING SHORE - PROBLEM?

# current = ConstantMover('cur',['uo','vo'],uo=0.1/3600.,vo=0.1/3600.)

# # unstick
# p1 = PassiveTracer('p1', nbuff=1000,geod=False,
#                     movers=[current], stickers=[shore],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[1.5,1,0],outfile='shore_x3_p2.out')

# # stick
# p1s = PassiveTracer('p1s', nbuff=1000,geod=False,
#                     movers=[current], stickers=[shore],unstick=0,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[1.5,1,0],outfile='shore_x3_p2_stick.out')


# ercore=ERcore(geod=False)
# ercore.materials=[p1,p1s]
# ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)


#################################
## tidal                       ##
## to see if it keeps moving   ##
#################################


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



from utils import make_dep,make_curr2d
x = numpy.arange(4)
y = numpy.arange(4)
# xx,yy = numpy.meshgrid(x,y)
# dep = numpy.ones(xx.shape)
# dep[:,:2] = 20.
# dep[:,2:] = 10.
# print dep
# filedep = 'dep1.nc'
# make_dep(x,y,dep,filedep)

tstart = datetime.datetime(2000,1,1)
tmid   = datetime.datetime(2000,1,1,6)
tend   = datetime.datetime(2000,1,1,12)
times = [0,6*3600,13*3600]
units = tstart.strftime('seconds since %Y-%m-%d %H:%M:%S')
u = numpy.ones((len(times),len(x),len(y)))
v = numpy.ones((len(times),len(x),len(y)))
u[0,:,:] = 1./3600.
u[1,:,:] = 0.
u[2,:,:] = -1./3600.
v[:,:,:] = 0.01/3600.
filecur='cur.nc'
make_curr2d (x,y,u,v,times,units,filecur)

#dep=GriddedTopo('depth',['dep'], file=filedep, zinvert=True)
shore = Shoreline(id='shore', file=fileshore)
current=GriddedMover ('cur',['u','v'],file=filecur)


p1 = PassiveTracer('p1', nbuff=1000,geod=False,
                    movers=[current], stickers=[shore],unstick=1,
                    tstart=t1,tend=t1, tstep=0., 
                    reln=2,P0=[0,0,-16],outfile='shore_x2_tide.out')

p1s = PassiveTracer('p1s', nbuff=1000,geod=False,
                    movers=[current], stickers=[shore],unstick=0,
                    tstart=t1,tend=t1, tstep=0., 
                    reln=2,P0=[0,0,-16],outfile='shore_x2_tide_stick.out')


ercore=ERcore(geod=False)
ercore.materials=[p1,p1s]
ercore.run(t=tstart,tend=tend,dt=3600)




# current=GriddedMover ('cur',['uo','vo'],file='../passive/uds_gsb_test.nc')

############
## Bridge ##
############


#   0   1   2   3   4   5   6
# 6 -------------------------
# 5         |   |
# 4         -----
# 3
# 2         -----
# 1         |   |
# 0 --------------------------

# fileshore = 'bridge.bnd'
# f = open(fileshore, 'w+')
# f.write('"Map Bounds","1",4\n')
# f.write('%i,%i\n' % (0,0))
# f.write('%i,%i\n' % (0,6))
# f.write('%i,%i\n' % (6,6))
# f.write('%i,%i\n' % (6,0))
# f.write('"0","1",5\n')  # south leg
# f.write('%.1f,%.1f\n' % (2,0))
# f.write('%.1f,%.1f\n' % (3,0))
# f.write('%.1f,%.1f\n' % (3,2))
# f.write('%.1f,%.1f\n' % (2,2))
# f.write('%.1f,%.1f\n' % (2,0))
# f.write('"1","1",5\n') # north leg
# f.write('%.1f,%.1f\n' % (2,4))
# f.write('%.1f,%.1f\n' % (3,4))
# f.write('%.1f,%.1f\n' % (3,6))
# f.write('%.1f,%.1f\n' % (2,6))
# f.write('%.1f,%.1f\n' % (2,4))
# f.close()

# # sticks to bridge
# # passes bridge if velocity to high (uo=2/3600)

# # e.g. passes bridge if L < uo*dt
# # where L is the width of the bridge leg

# shore = Shoreline(id='shore', file=fileshore)

# ## brigde slow ##
# current = ConstantMover('cur',['uo','vo'],uo=0.5/3600.,vo=0.0)

# # north leg
# p1 = PassiveTracer('p1', nbuff=1000,geod=False,
#                     movers=[current], stickers=[shore],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0,5,0],outfile='bridge_N.out')


# p2 = PassiveTracer('p2', nbuff=1000,geod=False,
#                     movers=[current], stickers=[shore],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0,3,0],outfile='bridge_mid.out')

# p3 = PassiveTracer('p3', nbuff=1000,geod=False,
#                     movers=[current], stickers=[shore],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0,1,0],outfile='bridge_S.out')


# ercore=ERcore(geod=False)
# ercore.materials=[p1,p2,p3]
# ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)

# ## brigde fast ##
# current = ConstantMover('cur',['uo','vo'],uo=2./3600.,vo=0.0)

# # north leg
# p1 = PassiveTracer('p1', nbuff=1000,geod=False,
#                     movers=[current], stickers=[shore],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0,5,0],outfile='bridge2_N.out')


# p2 = PassiveTracer('p2', nbuff=1000,geod=False,
#                     movers=[current], stickers=[shore],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0,3,0],outfile='bridge2_mid.out')

# p3 = PassiveTracer('p3', nbuff=1000,geod=False,
#                     movers=[current], stickers=[shore],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0,1,0],outfile='bridge2_S.out')


# ercore=ERcore(geod=False)
# ercore.materials=[p1,p2,p3]
# ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)
