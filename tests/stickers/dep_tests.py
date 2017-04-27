#!/usr/bin/env python
import sys; sys.path += ["/source/ercore"]

import numpy,datetime
from ercore import ERcore
from ercore.fields import ConstantMover,GriddedTopo,GriddedMover
from ercore.materials import PassiveTracer,BuoyantTracer

t1 = datetime.datetime(2000,1,1)
t2 = datetime.datetime(2000,1,1,1)

# ###########################
# ## below to above bottom ##
# ###########################

# from utils import make_dep    
# x = numpy.arange(4)
# y = numpy.arange(4)
# xx,yy = numpy.meshgrid(x,y)
# dep = numpy.ones(xx.shape)
# dep[:,:2] = -0.5
# dep[:,2:] = 0.5
# print dep
# filedep = 'dep1.nc'
# make_dep(x,y,dep,filedep)

# dep=GriddedTopo('depth',['dep'], file=filedep)#, zinvert=True)
# current = ConstantMover('cur',['uo','vo'],uo=0.5/3600.,vo=0)

# # rising particle, start below bottom

# p0 = BuoyantTracer('p0', nbuff=1000,geod=False,
#                     movers=[current], stickers=[dep],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     w0=0.5/3600.,
#                     reln=2,P0=[0,0,-22],outfile='dep_p0.out')


# ercore=ERcore(geod=False)
# ercore.materials=[p0]
# ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)


###########################
## positive to negative z ##
###########################

from utils import make_dep    
x = numpy.arange(4)
y = numpy.arange(4)
xx,yy = numpy.meshgrid(x,y)
dep = numpy.ones(xx.shape)
dep[:,:2] = 0.5
dep[:,2:] = -0.5
print dep
filedep = 'dep2.nc'
make_dep(x,y,dep,filedep)

dep=GriddedTopo('depth',['dep'], file=filedep, zinvert=True)
current = ConstantMover('cur',['uo','vo'],uo=-0.5/3600.,vo=0)

p0 = PassiveTracer('p0', nbuff=1000,geod=False,
                    movers=[current], stickers=[dep],unstick=1,
                    tstart=t1,tend=t1, tstep=0., 
                    reln=2,P0=[2.5,0,1],outfile='dep2_p0.out')


ercore=ERcore(geod=False)
ercore.materials=[p0]
ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)
ddd

# #####################
# ## meet the bottom ##
# #####################

# from utils import make_dep    
# x = numpy.arange(4)
# y = numpy.arange(4)
# xx,yy = numpy.meshgrid(x,y)
# dep = numpy.ones(xx.shape)
# dep[:,:2] = 20.
# dep[:,2:] = 10.
# print dep
# filedep = 'dep1.nc'
# make_dep(x,y,dep,filedep)

# dep=GriddedTopo('depth',['dep'], file=filedep, zinvert=True)
# current = ConstantMover('cur',['uo','vo'],uo=0.5/3600.,vo=0.0)

# # start below bottom
# p1 = PassiveTracer('p1', nbuff=1000,geod=False,
#                     movers=[current], stickers=[dep],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0,2,-30],outfile='dep1_p1.out')

# p1s = PassiveTracer('p1s', nbuff=1000,geod=False,
#                     movers=[current], stickers=[dep],unstick=0,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0,2,-30],outfile='dep1_p1stick.out')

# # start at bottom
# p2 = PassiveTracer('p2', nbuff=1000,geod=False,
#                     movers=[current], stickers=[dep],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0,2,-20],outfile='dep1_p2.out')

# p2s = PassiveTracer('p2s', nbuff=1000,geod=False,
#                     movers=[current], stickers=[dep],unstick=0,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0,2,-20],outfile='dep1_p2stick.out')


# # start above bottom
# p3 = PassiveTracer('p3', nbuff=1000,geod=False,
#                     movers=[current], stickers=[dep],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0,2,-16],outfile='dep1_p3.out')

# p3s = PassiveTracer('p3s', nbuff=1000,geod=False,
#                     movers=[current], stickers=[dep],unstick=0,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[0,2,-16],outfile='dep1_p3stick.out')


# ercore=ERcore(geod=False)
# ercore.materials=[p1,p2,p3, p1s,p2s,p3s]
# ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)

# ###########################
# ## start at bottom slope ##
# ###########################

# dep=GriddedTopo('depth',['dep'], file=filedep, zinvert=True)
# current = ConstantMover('cur',['uo','vo'],uo=-0.5/3600.,vo=0.0)


# p4 = PassiveTracer('p4', nbuff=1000,geod=False,
#                     movers=[current], stickers=[dep],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[1.5,2,-15],outfile='dep1_p4.out')


# p4s = PassiveTracer('p4s', nbuff=1000,geod=False,
#                     movers=[current], stickers=[dep],unstick=0,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[1.5,2,-15],outfile='dep1_p4s.out')



# ercore=ERcore(geod=False)
# ercore.materials=[p4,p4s]
# ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)

#################################
## meet the bottom in diagonal ##
## to see if it keeps moving   ##
#################################

# from utils import make_dep    
# x = numpy.arange(4)
# y = numpy.arange(4)
# xx,yy = numpy.meshgrid(x,y)
# dep = numpy.ones(xx.shape)
# dep[:,:2] = 20.
# dep[:,2:] = 10.
# print dep
# filedep = 'dep1.nc'
# make_dep(x,y,dep,filedep)
 
# dep=GriddedTopo('depth',['dep'], file=filedep, zinvert=True)
# current = ConstantMover('cur',['uo','vo'],uo=0.1/3600.,vo=0.1/3600.)


# p1 = PassiveTracer('p1', nbuff=1000,geod=False,
#                     movers=[current], stickers=[dep],unstick=1,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[1,2,-16],outfile='dep1_diag_p1.out')

# p1s = PassiveTracer('p1s', nbuff=1000,geod=False,
#                     movers=[current], stickers=[dep],unstick=0,
#                     tstart=t1,tend=t1, tstep=0., 
#                     reln=2,P0=[1,2,-16],outfile='dep1_diag_p1stick.out')


# ercore=ERcore(geod=False)
# ercore.materials=[p1,p1s]
# ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)


#################################
## tidal                       ##
## to see if it keeps moving   ##
#################################

from utils import make_dep,make_curr2d
x = numpy.arange(4)
y = numpy.arange(4)
xx,yy = numpy.meshgrid(x,y)
dep = numpy.ones(xx.shape)
dep[:,:2] = 20.
dep[:,2:] = 10.
print dep
filedep = 'dep1.nc'
make_dep(x,y,dep,filedep)

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

dep=GriddedTopo('depth',['dep'], file=filedep, zinvert=True)
current=GriddedMover ('cur',['u','v'],file=filecur)


p1 = PassiveTracer('p1', nbuff=1000,geod=False,
                    movers=[current], stickers=[dep],unstick=1,
                    tstart=t1,tend=t1, tstep=0., 
                    reln=2,P0=[0,0,-16],outfile='dep1_diag2_p1.out')

p1s = PassiveTracer('p1s', nbuff=1000,geod=False,
                    movers=[current], stickers=[dep],unstick=0,
                    tstart=t1,tend=t1, tstep=0., 
                    reln=2,P0=[0,0,-16],outfile='dep1_diag2_p1stick.out')


ercore=ERcore(geod=False)
ercore.materials=[p1,p1s]
ercore.run(t=tstart,tend=tend,dt=3600)




