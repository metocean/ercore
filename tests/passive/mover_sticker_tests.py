#!/usr/bin/env python
import numpy,datetime,sys,os,time
from ercore import ERcore,ncep2dt,dt2ncep
from ercore.materials import PassiveTracer,BuoyantTracer
from ercore.fields import GriddedTopo,GriddedMover,GriddedDataGroup,ConstantMover,ConstantDiffuser,GriddedElevation
from ercore.shoreline import Shoreline

plot = False
test=int(sys.argv[1])


P0=[170.5,-46,-10]

#import pdb;pdb.set_trace()

dep = GriddedTopo('depth',['dep'],file='uds_gsb_test.nc',zinvert=True)
shoreline  = Shoreline('shoreline','otago_shoreline2.bnd')
currents = GriddedMover('cur',['uo','vo'],is3d=True,file='uds_gsb_test.nc',topo=dep)
constantcur    = ConstantMover('constantcur'  ,['u','v'],u=-1.0,v=1.0,topo=dep)
elevation=GriddedTopo
diff= ConstantDiffuser('diff',['diffx','diffy','diffz'],diffx=0.1,diffy=0.1,diffz=0.000)

tstart = datetime.datetime(2009,1,1)
tend   = datetime.datetime(2009,1,3)
tstart = 733409 
tend   = 733411
tstep  = 900.
tout   = 900.

ercore=ERcore(tout=tout,rkorder=4)


if test==1: ## shoreline sticker - material cannot unstick


    pa1 = PassiveTracer('shoreline_unstick',10000,movers=[constantcur],
                                            stickers=[shoreline],
                                            diffusers=[diff],
                                            unstick=1,
                                            reln=2*24*3600/900,
                                            P0=P0,
                                            tstart=tstart, 
                                            tend=tend)

if test==2:  ## shoreline sticker - material cannot unstick

    pa1 = PassiveTracer('shoreline_stick',10000,movers=[constantcur],
                                            stickers=[shoreline],
                                            diffusers=[diff],
                                            unstick=0,
                                            reln=2*24*3600/900,
                                            P0=P0,
                                            tstart=tstart, 
                                            tend=tend)
if test==3: ## griddedtopo sticker - material cannot unstick


    pa1 = BuoyantTracer('griddedtopo_unstick',10000,movers=[constantcur],
                                            stickers=[dep],
                                            diffusers=[diff],
                                            unstick=1,
                                            reln=2*24*3600/900,
                                            P0=P0,
                                            tstart=tstart, 
                                            tend=tend,
                                            w0=-10e-3)

if test==4:  ## griddedtopo sticker - material cannot unstick

    pa1 = BuoyantTracer('griddedtopo_stick',10000,movers=[constantcur],
                                            stickers=[dep],
                                            diffusers=[diff],
                                            unstick=0,
                                            reln=2*24*3600/900,
                                            P0=P0,
                                            tstart=tstart, 
                                            tend=tend,
                                            w0=-10e-3)

if test==5: ## griddedtopo sticker - material cannot unstick
# allow the unstick to be array sith same size as stickers so that the material can stick to one sticker
# but not the other

    pa1 = BuoyantTracer('griddedtopo_shore_unstick_both',10000,movers=[constantcur],
                                            stickers=[dep,shoreline],
                                            diffusers=[diff],
                                            unstick=1,
                                            reln=10*2*24*3600/900,
                                            P0=P0,
                                            tstart=tstart, 
                                            tend=tend+5,
                                            w0=-1e-3)
    ercore.materials = [pa1]
    ercore.run(tstart,tend,tstep)
    ## griddedtopo sticker - material cannot unstick

    pa1 = BuoyantTracer('griddedtopo_shore_stick_both',10000,movers=[constantcur],
                                            stickers=[dep,shoreline],
                                            diffusers=[diff],
                                            unstick=0,
                                            reln=10*2*24*3600/900,
                                            P0=P0,
                                            tstart=tstart, 
                                            tend=tend+5,
                                            w0=-1e-3)
    ercore.materials = [pa1]
    ercore.run(tstart,tend,tstep)

     ## griddedtopo sticker - material cannot unstick

    pa1 = BuoyantTracer('stick2dep_only',10000,movers=[constantcur],
                                            stickers=[dep,shoreline],
                                            diffusers=[diff],
                                            unstick=[0,1], 
                                            reln=10*2*24*3600/900,
                                            P0=P0,
                                            tstart=tstart, 
                                            tend=tend+5,
                                            w0=-3e-3)
    ercore.materials = [pa1]
    ercore.run(tstart,tend,tstep)

## griddedtopo sticker - material cannot unstick

    pa1 = BuoyantTracer('stick2shore_only',10000,movers=[constantcur],
                                            diffusers=[diff],
                                            stickers=[dep,shoreline],
                                            unstick=[1,0],
                                            reln=10*2*24*3600/900,
                                            P0=P0,
                                            tstart=tstart, 
                                            tend=tend+5,
                                            w0=-3e-3)
    ercore.materials = [pa1]
    ercore.run(tstart,tend,tstep)
    exit()

# > check that unstick work of for dep+shoreline
# > check that particles dont go above surface

if test==6: # positively buoyant particle reaching the surface
    pa1 = BuoyantTracer('positively_buoyant_particles',10000,movers=[constantcur],
                                            diffusers=[diff],
                                            stickers=[dep,shoreline],
                                            unstick=[1,0],
                                            reln=10*2*24*3600/900,
                                            P0=P0,
                                            tstart=tstart, 
                                            tend=tend+5,
                                            w0=+5e-3)


if test==7: # positively buoyant particle reaching the surface with ELEVATION FIELD - can unstick 
    elev = GriddedElevation('elev',['et'],file='uds_gsb_elevation.nc',zinvert=True)
    constantcur    = ConstantMover('constantcur'  ,['u','v'],u=-0.0,v=0.0,topo=dep)
    pa1 = BuoyantTracer('positivebuoyant_gridded_elevation_unstick',10000,movers=[constantcur],
                                            stickers=[dep,shoreline,elev],
                                            unstick=[1,1,1],
                                            reln=2*24*3600/900,
                                            P0=P0,
                                            tstart=tstart, 
                                            tend=tend+5,
                                            w0=+5e-3)
if test==8: # positively buoyant particle reaching the surface with ELEVATION FIELD - cannot unstick
    elev = GriddedElevation('elev',['et'],file='uds_gsb_elevation.nc',zinvert=True)
    constantcur    = ConstantMover('constantcur'  ,['u','v'],u=-0.0,v=0.0,topo=dep)
    pa1 = BuoyantTracer('positivebuoyant_gridded_elevation_stick',10000,movers=[constantcur],
                                            stickers=[dep,shoreline,elev],
                                            unstick=[1,1,0],
                                            reln=2*24*3600/900,
                                            P0=P0,
                                            tstart=tstart, 
                                            tend=tend+5,
                                            w0=+5e-3)

ercore.materials = [pa1]
ercore.run(tstart,tend,tstep)

if plot: plot_particles(pa1)




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

