#!/usr/bin/env python
import datetime,sys,os
#sys.path.insert(0,os.path.join(os.path.dirname(__file__),'..','..'))
sys.path.insert(1,'/home/rosa/git/ercore/')
from pylab import *
from ercore._flib_ercore import interp3d,interph,slipvel
from ercore import ERcore,ncep2dt,dt2ncep
from ercore.materials import PassiveTracer,BuoyantTracer
from ercore.fields import GriddedMover,ConstantDiffuser,VariableDiffuser,GriddedTopo,ConstantMover,GriddedDataGroup
from ercore.shoreline import Shoreline


def plot_particles(particles):
    tpast = 0    
    with open(particles.outfile) as f:     
        print f.readline()
        it = -1
        for line in f.readlines():
            tp,pid,px,py,pz,pstate,page,pmass = map(float, line.split('\t'))        
            if tp != tpast:            
                if it > 0:
                    plt.savefig('plt_%s_%03i.png' % (particles.outfile,it))
                    plt.close()
                tpast=tp
                it+=1
                print 'Plotting '+ncep2dt(tp).strftime('%Y-%m-%d %H:%M')
                figure()
                for i,j in enumerate(shoreline.polyi):
                    n=shoreline.polyn[i]
                    plot(shoreline.slx[j:j+n-1],shoreline.sly[j:j+n-1],'k')
                plot(P0[0],P0[1],'r+')
                title(ncep2dt(tp).strftime('%Y-%m-%d %H:%M'))
            plot(px,py, 'bo')
            #if it==2: break
    plt.savefig('plt_%s_%03i.png' % (particles.outfile,it))
    plt.close()


P0=[170.5,-46,-10]

dep=GriddedTopo('depth',['dep'],file='uds_gsb_test.nc',zinvert=True)
currents=GriddedMover('cur',['uo','vo'],is3d=True,file='uds_gsb_test.nc',topo=dep)
tide=GriddedMover('tide',['ut','vt'],is3d=False,file='uds_gsb_test.nc')

constantsurf=ConstantMover('constantsurf',['u','v'],u=1.0,v=0,topo=dep)
constant3D=ConstantMover('constant3D',['u','v'],is3d=True,levels=[0,-10,-30,-40],
                          u=[1.0,0.5,0.0,-0.5],v=[0,0,0,0],topo=dep)
constant3Dsub=ConstantMover('constant3Dsub',['u','v'],is3d=True,levels=[0,-10,-30,-40],
                          u=[0.0,-0.5,-1.0,-1.5],v=[0,0,0,0],surfsub=True,topo=dep)

diff=VariableDiffuser('diff',['diffx','diffy','diffz'],diffz=0.001,P0=P0)
shoreline=Shoreline('shoreline','shoreline2.bnd')
movergroup=GriddedDataGroup('group',['uo','vo'],[currents])

tstart = datetime.datetime(2009,1,1)
tend   = datetime.datetime(2009,1,3)
tstep  = 900.
tout   = 3600.

particles=BuoyantTracer('particlesA',10000,movers=[movergroup],stickers=[shoreline,dep],reln=10,P0=P0,w0=0.0,tstart=tstart, tend=tstart)

particles1=BuoyantTracer('particlesA2',10000,movers=[currents],stickers=[shoreline,dep],reln=10,P0=P0,w0=0.0,tstart=tstart, tend=tstart)

particles1t=BuoyantTracer('particlesA2t',10000,movers=[currents,tide],stickers=[shoreline,dep],reln=10,P0=P0,w0=0.0,tstart=tstart, tend=tstart)

particles2=BuoyantTracer('particlesB',10000,movers=[constantsurf,constant3Dsub],stickers=[shoreline,dep],reln=10,P0=P0,w0=0.0,tstart=tstart, tend=tstart)

particles3=BuoyantTracer('particlesC',10000,movers=[constant3D],stickers=[shoreline,dep],reln=10,P0=P0,w0=0.0,tstart=tstart, tend=tstart)

ercore=ERcore(tout=tout,rkorder=4)

###### start tests #####

ercore.materials=[particles]
ercore.run(tstart,tend,tstep)
plot_particles(particles)

ercore.materials=[particles1]
ercore.run(tstart,tend,tstep)
plot_particles(particles1) # same as before

ercore.materials=[particles1t]
ercore.run(tstart,tend,tstep)
plot_particles(particles1t)

## multiple releases of the same material ##

mdur = datetime.timedelta(seconds=7200) # release duration
mint = datetime.timedelta(hours=24)     # interval between releases

tstartp = tstart
mats = []
i=1
while tstartp < tend:
    tendp = tstartp + mdur
    print tstartp, tendp
    m1 = BuoyantTracer('p%i' % (i),10000,movers=[currents],stickers=[shoreline,dep],reln=10,P0=P0,w0=0.0,tstart=tstartp, tend=tendp)
    mats = mats + [m1]
    tstartp = tstartp + mint
    i+=1

ercore.materials=mats
ercore.run(tstart,tend,tstep)

## plot all pasrticles by timestamp ##
t1 = dt2ncep(tstart)
t2 = dt2ncep(tend)
dt = tout/86400.
for tc in numpy.arange(t1,t2,dt):
    tc1 = tc + dt
    tcs = ncep2dt(tc).strftime('%Y%m%d_%H%M')
    pxs = []
    pys = []
    for m in mats:
        with open(m.outfile) as f:     
            f.readline() # skip header            
            for line in f.readlines():
                tp,pid,px,py,pz,pstate,page,pmass = map(float, line.split('\t'))
                tps = ncep2dt(tp).strftime('%Y%m%d_%H%M')
                #if (tp >= tc) & (tp < tc1):
                if tps == tcs:
                    pxs = pxs + [px]
                    pys = pys + [py]    
    print ncep2dt(tc), len(pxs) # total number of particles at a given timestamp
    figure()
    for i,j in enumerate(shoreline.polyi):
        n=shoreline.polyn[i]
        plot(shoreline.slx[j:j+n-1],shoreline.sly[j:j+n-1],'k')
    plot(P0[0],P0[1],'r+')
    title(ncep2dt(tc).strftime('%Y-%m-%d %H:%M'))
    plot(pxs,pys, 'bo')
    plt.savefig(ncep2dt(tc).strftime('plt_multi_%Y%m%d_%H%M.png'))
    plt.close()

## 3D ##

ercore.rkorder=3
ercore.materials=[particles3]
ercore.run(tstart,tend,tstep)
plot_particles(particles3)  # function doesn't plot 3d particles
plt.plot(particles3.pos[:,0],particles3.pos[:,2])  # only final positions in the buffer




