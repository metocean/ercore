#!/usr/bin/env python
import numpy,datetime
#sys.path.insert(1,os.path.join(os.path.dirname(__file__),'..','..'))
#import ercore._flib_ercore as flib
from pylab import *
from ercore import ERcore
from ercore.fields import ConstantMover,GriddedTopo
from ercore.materials import PassiveTracer

x = numpy.arange(5)
y = numpy.arange(4)
xx,yy = numpy.meshgrid(x,y)
z = -1*(4.-xx)
print z
z = z.astype('float32')

import cdms2
cdms2.setNetcdfDeflateFlag(0)
cdms2.setNetcdfShuffleFlag(0)
nc = cdms2.open('dep.nc','w+')
xax = nc.createAxis('x', x)
yax = nc.createAxis('y', y)
dep = nc.createVariable('dep','f',[yax,xax])
dep[:,:] = z[:,:]
nc.close()


dep=GriddedTopo('depth',['dep'], file='dep.nc',zinvert=True,geod=True)
current = ConstantMover('cur',['uo','vo'],uo=1./3600.,vo=0.0)
t1 = datetime.datetime(2000,1,1)
t2 = datetime.datetime(2000,1,1,1)

particles = PassiveTracer('bombs', nbuff=1000,
                            movers=[current], stickers=[dep], geod=False,
                            tstart=t1,tend=t1, tstep=0.,
                            reln=1,
                            P0=[0,1,-2])

# import pdb; pdb.set_trace()
ercore=ERcore(geod=False)
ercore.materials=[particles]
ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)

#
# colors=['b+','r+','g+','m+']
# figure()
# for rk in range(4,0,-1):
#   part=copy.deepcopy(particles)
#   ercore.materials=[part]
#   if not os.path.isdir('rk'+str(rk)):os.mkdir('rk'+str(rk))
#   ercore.rkorder=rk
#   ercore.outpath='rk'+str(rk)
#   ercore.run(datetime.datetime(2009,1,1),datetime.datetime(2009,1,2),900)
#   plot(part.pos[:,0],part.pos[:,2],colors[rk-1], marker='o')
#
# show()
#
#
#
#
#
#
#
# particles=BuoyantTracer('particles',10000,geod=False,movers=[current],diffusers=[],
# stickers=[],
# reln=1000,P0=[0,0,0],w0=-0.1,tstart=datetime.datetime(2009,1,1),tend=datetime.datetime(2009,1,1))
#
#
# self,id,nbuff,movers=[],reactors=[],stickers=[],diffusers=[],
# tstart=None,tend=None,tstep=0.,outfile=None,P0=[0,1,0],spawn=1,reln=0,R0=1.,Q0=1.,unstick=0.,**prop):
