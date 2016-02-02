#!/usr/bin/env python
import numpy,datetime
import cdms2
from ercore import ERcore
from ercore.fields import ConstantMover,GriddedTopo,GriddedMover, ConstantElevation
from ercore.shoreline import Shoreline
from ercore.materials import PassiveTracer, BuoyantTracer
from ercore import dt2ncep

t1 = datetime.datetime(2000,1,1)
t2 = datetime.datetime(2000,1,1,1)

############################
## Rising moving particle ##
############################

current = ConstantMover('cur',['u','v'],u=1./3600., v=0)
elev = ConstantElevation('elev', ['elev'], elev=1.)

p1 = BuoyantTracer('p1', nbuff=1000,geod=False,
                    movers=[current], stickers=[elev],
                    tstart=t1,tend=t1, tstep=0., w0=0.5/3600.,
                    reln=2,P0=[0,0,0],outfile='p1.out')

ercore=ERcore(geod=False)
ercore.materials=[p1]
ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)

#
#
#
# def make_dep (deplevels=[2,1],filedep='dep1.nc'):
#     x = numpy.arange(4)
#     y = numpy.arange(3)
#     xx,yy = numpy.meshgrid(x,y)
#     dep = numpy.ones(xx.shape)
#     print dep.shape
#     dep[:,:2] = deplevels[0]
#     dep[:,2:] = deplevels[1]
#     print 'x = ', x
#     print 'y = ', y
#     print 'dep = ', dep
#     dep = dep.astype('float32')
#     cdms2.setNetcdfDeflateFlag(0)
#     cdms2.setNetcdfShuffleFlag(0)
#     nc = cdms2.open(filedep,'w+')
#     xax = nc.createAxis('lon', x)
#     yax = nc.createAxis('lat', y)
#     dvar = nc.createVariable('dep','f',[yax,xax])
#     dvar[:,:] = dep[:,:]
#     nc.close()
#

#
# ## Tank with [0,-1]
# filedep='dep.nc'
# unstick = 1
# deplevels=[1,0]
# h = 2  # z
# P0 = [0,1,0.5]
#
# make_dep(deplevels=[0,-1], filedep=filedep)
# dep=GriddedTopo('depth',['dep'], file=filedep, zinvert=True)
# u = numpy.ones(len(deplevels))/3600.
# v = numpy.zeros(len(deplevels))
# current = ConstantMover('cur',['u','v'],u=u, v=v, is3d=True, depth=deplevels, topo=dep)
# # All wet
#
# elev = ConstantElevation('elev', ['elev'], elev=h)
#
# # will meet bottom
# p1 = PassiveTracer('p1', nbuff=1000,geod=False,
#                             movers=[current], stickers=[elev,dep],
#                             unstick=unstick,
#                             tstart=t1,tend=t1, tstep=0., w0=1.,
#                             reln=2,P0=[0,1,0.5],outfile='p1.out')
# #
# # # will always be underwater
# # p2 = PassiveTracer('p2', nbuff=1000,geod=False,
# #                             movers=[current], stickers=[elev,dep],
# #                             unstick=unstick,
# #                             tstart=t1,tend=t1, tstep=0.,
# #                             reln=2,P0=[0,1,1.5],outfile='p2.out')
# #
# # # staring above free surface
# # p3 = PassiveTracer('p3', nbuff=1000,geod=False,
# #                             movers=[current], stickers=[elev,dep],
# #                             unstick=unstick,
# #                             tstart=t1,tend=t1, tstep=0.,
# #                             reln=2,P0=[0,1,2.5],outfile='p3.out')
# #
# # # starting below surface but moving up
# # p4 = PassiveTracer('p4', nbuff=1000,geod=False,
# #                             movers=[current], stickers=[elev,dep],
# #                             unstick=unstick,
# #                             tstart=t1,tend=t1, tstep=0., w0=1.,
# #                             reln=2,P0=[0,1,0],outfile='p4.out')
# ercore=ERcore(geod=False)
# ercore.materials=[p1,p2,p3,p4]
# ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)
