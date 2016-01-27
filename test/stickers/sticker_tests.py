#!/usr/bin/env python
import numpy,datetime
import cdms2
from ercore import ERcore
from ercore.fields import ConstantMover,GriddedTopo,GriddedMover
from ercore.materials import PassiveTracer
from ercore import dt2ncep

t1 = datetime.datetime(2000,1,1)
t2 = datetime.datetime(2000,1,1,1)

### NO BOTTOM ###
# current = ConstantMover('cur',['uo','vo'],uo=1./3600.,vo=0.0)
# bombs = PassiveTracer('bombs', nbuff=1000,
#                             movers=[current], #stickers=[dep],
#                             geod=False,
#                             tstart=t1,tend=t1, tstep=0.,
#                             reln=1,
#                             P0=[0,1,0],
#                             outfile='surface_no_bottom.out')
#
#
# ercore=ERcore(geod=False)
# ercore.materials=[bombs]
 # ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)

### BOTTOM ###
def make_dep (deplevels=[2,1],filedep='dep1.nc'):
    x = numpy.arange(4)
    y = numpy.arange(3)
    xx,yy = numpy.meshgrid(x,y)
    dep = numpy.ones(xx.shape)
    print dep.shape
    dep[:,:2] = deplevels[0]
    dep[:,2:] = deplevels[1]
    print 'x = ', x
    print 'y = ', y
    print 'dep = ', dep
    dep = dep.astype('float32')
    cdms2.setNetcdfDeflateFlag(0)
    cdms2.setNetcdfShuffleFlag(0)
    nc = cdms2.open(filedep,'w+')
    xax = nc.createAxis('lon', x)
    yax = nc.createAxis('lat', y)
    dvar = nc.createVariable('dep','f',[yax,xax])
    dvar[:,:] = dep[:,:]
    nc.close()

def make_currents (t1,t2,u,fileout):
    x = numpy.arange(4)
    y = numpy.arange(3)
    xx,yy = numpy.meshgrid(x,y)
    dt = t2 - t1
    dtsec = dt.days*24.*3600. + dt.seconds
    uu = numpy.zeros((3,len(y),len(x)))
    vv = numpy.zeros((3,len(y),len(x)))
    ww = numpy.zeros((3,len(y),len(x)))
    uu[0,...] = u
    uu[2,...] = -1*u
    print 'x = ', x
    print 'y = ', y
    print 'uu = ', uu
    cdms2.setNetcdfDeflateFlag(0)
    cdms2.setNetcdfShuffleFlag(0)
    nc = cdms2.open(fileout,'w+')
    xax = nc.createAxis('lon', x)
    yax = nc.createAxis('lat', y)
    tax = nc.createAxis('time', numpy.array(map(numpy.float32,[0,dtsec/2.,dtsec])))
    tax.units = t1.strftime('seconds since %Y-%m-%d %H:%M:%S')
    uvar = nc.createVariable('u','f',[tax,yax,xax])
    vvar = nc.createVariable('v','f',[tax,yax,xax])
    wvar = nc.createVariable('w','f',[tax,yax,xax])
    uvar[:,:] = uu[:].astype('float32')
    vvar[:,:] = vv[:].astype('float32')
    wvar[:,:] = ww[:].astype('float32')
    nc.close()


def Tank (deplevels,filedep,P0,unstick,outfile):
    ''' Tank below 0 + constante movers'''
    dep=GriddedTopo('depth',['dep'], file=filedep, zinvert=True)
    u = numpy.ones(len(deplevels))/3600.
    v = numpy.zeros(len(deplevels))
    current = ConstantMover('cur',['u','v'],u=u, v=v, is3d=True, depth=deplevels, topo=dep)
    bombs = PassiveTracer('bombs', nbuff=1000,geod=False,
                                movers=[current], stickers=[dep], unstick=unstick,
                                tstart=t1,tend=t1, tstep=0.,
                                reln=1,P0=P0,outfile=outfile)
    ercore=ERcore(geod=False)
    ercore.materials=[bombs]
    ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)

def TankCurrents (deplevels,filedep,filecur,P0,unstick,outfile):
    ''' Tank below 0 + constante movers'''
    dep=GriddedTopo('depth',['dep'], file=filedep, zinvert=True)
    u = numpy.ones(len(deplevels))/3600.
    v = numpy.zeros(len(deplevels))
    current = GriddedMover('cur',['u','v', 'w'], file=filecur)
    #, is3d=True, depth=deplevels, topo=dep)
    bombs = PassiveTracer('bombs', nbuff=1000,geod=False,
                                movers=[current], stickers=[dep], unstick=unstick,
                                tstart=t1,tend=t1, tstep=0.,
                                reln=1,P0=P0,outfile=outfile)
    ercore=ERcore(geod=False)
    ercore.materials=[bombs]
    ercore.run(t=datetime.datetime(2000,1,1),tend=datetime.datetime(2000,1,1,12),dt=3600)

###########
## Tank1 ##
###########
if False:
    deplevels=[2,1]
    filedep = 'dep1.nc'
    make_dep (deplevels=deplevels,filedep=filedep)

    # do not meet bottom - OK
    Tank(deplevels=deplevels,filedep=filedep,
        P0=[0,1,-0.5], unstick = 0, outfile = 'tank1_test0.out')

    # meet the bottom and stick - OK
# Time	id	x	y	z	state	age	mass
# 730121.041667	1	1.000000	1.000000	-1.499998	1	0.041667	1.000000
# 730121.083333	1	1.000000	1.000000	-1.499998	-1	0.041667	1.000000
    Tank(deplevels=deplevels,filedep=filedep,
        P0=[0,1,-1.5], unstick = 0, outfile = 'tank1_test1.out')

    # meet the bottom and do not stick - OK
# Time	id	x	y	z	state	age	mass
# 730121.041667	1	1.000000	1.000000	-1.499998	1	0.041667	1.000000
# 730121.083333	1	1.500004	1.000000	-1.499996	2	0.083333	1.000000
# 730121.125000	1	1.500004	1.000000	-1.499996	3	0.125000	1.000000
# 730121.166667	1	1.500004	1.000000	-1.499996	4	0.166667	1.000000
    Tank(deplevels=deplevels,filedep=filedep,
        P0=[0,1,-1.5], unstick = 1, outfile = 'tank1_test2.out')

    # start below de bottom - was nan, now is bottom
# Time	id	x	y	z	state	age	mass
# 730121.041667	1	0.000000	1.000000	-2.000000	2	0.041667	1.000000
# 730121.083333	1	1.000000	1.000000	-1.999998	2	0.083333	1.000000
# 730121.125000	1	1.000002	1.000000	-1.999998	3	0.125000	1.000000
    Tank(deplevels=deplevels,filedep=filedep,
        P0=[0,1,-3], unstick = 1, outfile = 'tank1_test3.out')


###########
## Tank2 ##
###########

if False:

    deplevels=[1,-1]
    filedep = 'dep2.nc'
    make_dep (deplevels=deplevels,filedep=filedep)

    # do not meet bottom - corrected - only had problem for constant movers when topo=0
    Tank(deplevels=deplevels,filedep=filedep,
        P0=[0,1,2], unstick = 0, outfile = 'tank2_test0.out')

    # meet the bottom and stick - OK
    Tank(deplevels=deplevels,filedep=filedep,
        P0=[0,1,-0.5], unstick = 0, outfile = 'tank2_test1.out')


    # meet the bottom and do not stick - corrected for topo = 0 and z > 0
    # import pdb; pdb.set_trace()
    Tank(deplevels=deplevels,filedep=filedep,
        P0=[0,1,-0.5], unstick = 1, outfile = 'tank2_test2.out')

    # start below de bottom - was nan, now is bottom
    Tank(deplevels=deplevels,filedep=filedep,
        P0=[0,1,-3], unstick = 1, outfile = 'tank2_test3.out')


## Tank2 with CURRENTS ##

make_currents (t1=datetime.datetime(2000,1,1),t2=datetime.datetime(2000,1,1,13),u=1./3600., fileout='currents.nc')
deplevels=[1,-1]
filedep = 'dep2.nc'
filecurrents='currents.nc'

# do not meet bottom
TankCurrents(deplevels=deplevels,filedep=filedep,filecur='currents.nc',
        P0=[0,1,2], unstick = 0, outfile = 'tank2_test4.out')

# meet the bottom and stick - OK
TankCurrents(deplevels=deplevels,filedep=filedep,filecur='currents.nc',
        P0=[0,1,-0.5], unstick = 0, outfile = 'tank2_test5.out')

# meet the bottom and do not stick - comeback to life?
TankCurrents(deplevels=deplevels,filedep=filedep,filecur='currents.nc',
        P0=[0,1,-0.5], unstick = 1, outfile = 'tank2_test6.out')
