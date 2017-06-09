#!/usr/bin/env python
import numpy,datetime,sys,os,time
from ercore import ERcore,ncep2dt,dt2ncep
from ercore.materials import PassiveTracer,BuoyantTracer
from ercore.fields import GriddedTopo,GriddedMover,GriddedDataGroup,ConstantMover,VariableDiffuser
from ercore.shoreline import Shoreline

plot = False

# 
test=int(sys.argv[1])

P0=[170.5,-46,-10]

#import pdb;pdb.set_trace()

dep = GriddedTopo('depth',['dep'],file='uds_gsb_test.nc',zinvert=True)
shoreline  = Shoreline('shoreline','shoreline2.bnd')
currents = GriddedMover('cur',['uo','vo'],is3d=True,file='uds_gsb_test.nc',topo=dep)
constant    = ConstantMover('currents'  ,['u','v'],u=0.0,v=0,topo=dep)


tstart = datetime.datetime(2009,1,1)
tend   = datetime.datetime(2009,1,3)
tstart = 733409 
tend   = 733411
tstep  = 900.
tout   = 900.

pa1 = PassiveTracer('particlesA2',10000,movers=[constant],stickers=[shoreline,dep],reln=2*24*3600/900,P0=P0,tstart=tstart, tend=tstart)

ercore=ERcore(tout=tout,rkorder=4)

if test==1:
#########################################
# TESTING RELEASE IN CIRCLE           ###
#########################################
    pa1 = PassiveTracer('particles_circle',10000,movers=[constant],stickers=[shoreline,dep],reln=10*2*24*3600/900,P0=P0,tstart=tstart,tend=tend,circular_radius=1000.0)
    ercore.materials = [pa1]
    ercore.run(tstart,tend,tstep)
    if plot: plot_particles(pa1)
if test==2:
#########################################
# TESTING RELEASE IN POLYGON          ###
#########################################
    #polygon: [ (x0,y0), (x1,y1), ..., (x0,y0) ] 
    pa1 = PassiveTracer('particles_polygon',10000,
                                                movers=[constant],
                                                stickers=[shoreline,dep],
                                                reln=10*2*24*3600/900,
                                                P0=P0,
                                                tstart=tstart,
                                                tend=tend,
                                                polygon=[[170.5,-46],[170.501,-46],[170.501,-46.001],[170.5,-46.001],[170.5,-46]])
    ercore.materials = [pa1]
    ercore.run(tstart,tend,tstep)
    if plot: plot_particles(pa1)
if test==3:
###############################################
# TESTING RELEASE IN TIME VARYING POLYGON   ###
###############################################
# 'variable_polys.txt' is a text file with X columns [time_in_days_since1-1-1   x1 y1 x2 y2 x3 y3 etc...]
    pa1 = PassiveTracer('particles_variable_poly',10000,movers=[constant],stickers=[shoreline,dep],reln=10*2*24*3600/900,P0=P0,tstart=tstart,tend=tend,variable_poly='variable_poly_file.txt')
    ercore.materials = [pa1]
    ercore.run(tstart,tend,tstep)
    if plot: plot_particles(pa1)
if test==4:
##############################################################
# TESTING RELEASE WITH TIME_VARYING NUMBER OF PARTICLES    ###
#############################################################3
# 'variable_release.txt' is a text file with 2 column [time_in_days_since1-1-1   number_of_particles_to_release]
#  The file must have the same timestep as the one used in the model
#  
    pa1 = PassiveTracer('particles_variable_reln',10000,movers=[currents],stickers=[shoreline,dep],reln=2*24*3600/900,P0=P0,tstart=tstart,tend=tend,variable_reln='variable_reln_file.txt')
    ercore.materials = [pa1]
    ercore.run(tstart,tend,tstep)
    if plot: plot_particles(pa1)
if test==5:
###################################################################################
# TESTING RELEASE WITH TIME_VARYING NUMBER OF PARTICLES AND TIME-VARYING POLYS##
###################################################################################
# 'variable_release.txt' is a text file with 2 column [time_in_days_since1-1-1   number_of_particles_to_release]
#  The file must have the same timestep as the one used in the model    
    pa1 = PassiveTracer('particles_variable_reln',10000,movers=[currents],stickers=[shoreline,dep],reln=2*24*3600/900,P0=P0,tstart=tstart,tend=tend,variable_reln='variable_reln_file.txt',variable_poly='variable_poly_file.txt')
    ercore.materials = [pa1]
    ercore.run(tstart,tend,tstep)
    if plot: plot_particles(pa1)

if test==6:
###################################################################################
# random depth within a specified range
###################################################################################
# 'variable_release.txt' is a text file with 2 column [time_in_days_since1-1-1   number_of_particles_to_release]
#  The file must have the same timestep as the one used in the model    
    pa1 = PassiveTracer('particles_random_depth',10000,movers=[currents],stickers=[shoreline,dep],reln=10*2*24*3600/900,P0=[170.5,-46,[-10,-12]],tstart=tstart,tend=tend)
    ercore.materials = [pa1]
    ercore.run(tstart,tend,tstep)
    if plot: plot_particles(pa1)

if test==7:
###################################################################################
# random depth within a circle and random depth within a specified range
###################################################################################
# 'variable_release.txt' is a text file with 2 column [time_in_days_since1-1-1   number_of_particles_to_release]
#  The file must have the same timestep as the one used in the model    
    pa1 = PassiveTracer('particles_random_depth_circle',10000,movers=[currents],stickers=[shoreline,dep],reln=10*2*24*3600/900,P0=[170.5,-46,[-10,-12]],tstart=tstart,tend=tend,circular_radius=1000.0)
    ercore.materials = [pa1]
    ercore.run(tstart,tend,tstep)
    if plot: plot_particles(pa1)

if test==8:
###################################################################################
# random depth within a circle and random depth within a specified range
###################################################################################
# 'variable_release.txt' is a text file with 2 column [time_in_days_since1-1-1   number_of_particles_to_release]
#  The file must have the same timestep as the one used in the model    
    pa1 = BuoyantTracer('particles_random_depth_circle_recycling',10000,movers=[currents],stickers=[shoreline,dep],reln=100*2*24*3600/900,P0=[170.5,-46,[-10,-12]],tstart=tstart,tend=tend,circular_radius=1000.0,w0=-0.01)
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

