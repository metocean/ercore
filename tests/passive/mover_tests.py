#!/usr/bin/env python
import datetime
from ercore import ERcore,ncep2dt,dt2ncep
from ercore.materials import PassiveTracer
from ercore.fields import GriddedTopo,GriddedMover,GriddedDataGroup,ConstantMover,VariableDiffuser
from ercore.shoreline import Shoreline

plot = False

def plot_particles(particles):
    import matplotlib.pyplot as plt
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
                plt.figure()
                for i,j in enumerate(shoreline.polyi):
                    n=shoreline.polyn[i]
                    plt.plot(shoreline.slx[j:j+n-1],shoreline.sly[j:j+n-1],'k')
                plt.plot(P0[0],P0[1],'r+')
                plt.title(ncep2dt(tp).strftime('%Y-%m-%d %H:%M'))
            plt.plot(px,py, 'bo')
            #if it==2: break
    plt.savefig('plt_%s_%03i.png' % (particles.outfile,it))
    plt.close()


P0=[170.5,-46,-10]

#import pdb;pdb.set_trace()

dep = GriddedTopo('depth',['dep'],file='uds_gsb_test.nc',zinvert=True)
shoreline  = Shoreline('shoreline','shoreline2.bnd')

currents = GriddedMover('cur',['uo','vo'],is3d=True,file='uds_gsb_test.nc',topo=dep)
tide     = GriddedMover('tide',['ut','vt'],is3d=False,file='uds_gsb_test.nc') #

movergroup = GriddedDataGroup('group',['uo','vo'],[currents])

constantsurf    = ConstantMover('constantsurf'  ,['u','v'],u=1.0,v=0,topo=dep)
constant3D      = ConstantMover('constant3D'    ,['u','v'],is3d=True,levels=[0,-10,-30,-40],u=[1.0,0.5,0.0,-0.5],v=[0,0,0,0],topo=dep)
constant3Dsub   = ConstantMover('constant3Dsub' ,['u','v'],is3d=True,levels=[0,-10,-30,-40],u=[0.0,-0.5,-1.0,-1.5],v=[0,0,0,0],surfsub=True,topo=dep)

diff       = VariableDiffuser('diff',['diffx','diffy','diffz'],diffz=0.001,P0=P0)


tstart = datetime.datetime(2009,1,1)
tend   = datetime.datetime(2009,1,3)
tstep  = 900.
tout   = 3600.

pa1 = PassiveTracer('particlesA1',10000,movers=[movergroup],stickers=[shoreline,dep],reln=10,P0=P0,tstart=tstart, tend=tstart)
pa2 = PassiveTracer('particlesA2',10000,movers=[currents],stickers=[shoreline,dep],reln=10,P0=P0,tstart=tstart, tend=tstart)
pa3 = PassiveTracer('particlesA3',10000,movers=[tide],stickers=[shoreline,dep],reln=10,P0=P0,tstart=tstart, tend=tstart)
pa4 = PassiveTracer('particlesA4',10000,movers=[currents,tide],stickers=[shoreline,dep],reln=10,P0=P0,tstart=tstart, tend=tstart)


pb1 = PassiveTracer('particlesB2',10000,movers=[constantsurf],stickers=[shoreline,dep],reln=10,P0=P0,tstart=tstart, tend=tstart)
pb2 = PassiveTracer('particlesB1',10000,movers=[constantsurf,constant3Dsub],stickers=[shoreline,dep],reln=10,P0=P0,tstart=tstart, tend=tstart)

pc1=PassiveTracer('particlesC1',10000,movers=[constant3D],stickers=[shoreline,dep],reln=10,P0=P0,tstart=tstart, tend=tstart)

ercore=ERcore(tout=tout,rkorder=4)

###################################
###     single release tests    ###
###################################

ercore.materials = [pa1]
ercore.run(tstart,tend,tstep)
if plot: plot_particles(pa1)


ercore.materials = [pa2]
ercore.run(tstart,tend,tstep)
if plot: plot_particles(pa2)

ercore.materials = [pa3]
ercore.run(tstart,tend,tstep)
if plot: plot_particles(pa3) # same as before


ercore.rkorder=3
ercore.materials=[pc1]
ercore.run(tstart,tend,tstep)
if plot: 
    plot_particles(pc1)  # function doesn't plot 3d particles
    plt.plot(pc1.pos[:,0],pc1.pos[:,2])  # only final positions in the buffer


###############################################
### multiple releases of the same material  ###
###############################################

mdur = datetime.timedelta(seconds=7200) # release duration
mint = datetime.timedelta(hours=24)     # interval between releases

tstartp = tstart
mats = []
i=1
while tstartp < tend:
    tendp = tstartp + mdur
    print tstartp, tendp
    m1 = PassiveTracer('p%i' % (i),10000,movers=[currents],stickers=[shoreline,dep],reln=10,P0=P0,tstart=tstartp, tend=tendp)
    mats = mats + [m1]
    tstartp = tstartp + mint
    i+=1

ercore.materials = mats
ercore.run(tstart,tend,tstep)

# plot all particles by timestamp
if plot:

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
            plt.plot(shoreline.slx[j:j+n-1],shoreline.sly[j:j+n-1],'k')
        plt.plot(P0[0],P0[1],'r+')
        plt.title(ncep2dt(tc).strftime('%Y-%m-%d %H:%M'))
        plt.plot(pxs,pys, 'bo')
        plt.savefig(ncep2dt(tc).strftime('plt_multi_%Y%m%d_%H%M.png'))
        plt.close()





