#!/usr/bin python
import datetime, copy, os
from ercore import ERcore
from ercore.materials import BuoyantTracer
from ercore.fields import ConstantMover,VariableDiffuser

plot = True
if plot:
    import matplotlib.pyplot as plt 


tstart = datetime.datetime(2009,1,1)
tend   = datetime.datetime(2009,1,2)
tstep  = 3600.  # model time step
tout   = tstep  # output frequency
reln   = 24    # number of particles per release
tstep_release = 3  # hours

current = ConstantMover('cur',['uo','vo'],uo=1.0,vo=0.0)
diff    = VariableDiffuser('diff',['diffx','diffy','diffz'],diffz=0.001,P0=[0,0,0])

particles = BuoyantTracer(id = 'particles', nbuff = 10000, 
                          movers=[current], 
                          #diffusers = [diff],
                          tstart = tstart,tend = tend,
                          #tstep_release = tstep_release,
                          P0 = [0,0,0],
                          circular_radius = 1., 
                          reln = reln,
                          w0 = -0.1)  # positive upwards

ercore=ERcore(tout=tout, geod = False)

if plot:
    colors=['b+','r+','g+','m+']


for rk in range(4,0,-1):
    print 'Running rk%i' % rk
    part=copy.deepcopy(particles)
    ercore.materials = [part]
    if not os.path.isdir('rk'+str(rk)):os.mkdir('rk'+str(rk))
    ercore.rkorder = rk
    ercore.outpath = 'rk'+str(rk)
    #import pdb;pdb.set_trace()
    ercore.run(t = tstart ,tend = tend, dt = tstep)
    
    if plot: 
        # import pdb;pdb.set_trace()
        # last time stamp
        fig,axs = plt.subplots(1,2, figsize=(10,4))
        ax = axs[0]
        ax.plot(part.pos[:reln,0],part.pos[:reln,1],colors[rk-1], marker='o')
        ax.set_xlabel('x'); ax.set_ylabel('y')
        ax = axs[1]
        ax.plot(part.pos[:reln,0],part.pos[:reln,2],colors[rk-1], marker='o')
        ax.set_xlabel('x'); ax.set_ylabel('z')
        fig.suptitle('rk%i - hours since start = %.3f ' % (rk, part.tcum/3600.))
        if plot: plt.show()
        break

  

  
