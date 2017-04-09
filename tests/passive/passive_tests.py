#!/usr/bin python
import datetime, copy
# import matplotlib.pyplot as plt 
from ercore import ERcore
from ercore.materials import BuoyantTracer
from ercore.fields import ConstantMover,VariableDiffuser



# class _Material(object):
#   """Initialization:
#     <MaterialClass>(id,nbuff,movers,reactors,stickers,diffusers,tstart,tend,outfile,**properties)
#     Arguments:
#       id: Unique id
#       nbuff: Total number of particles in buffer
#       movers: List of mover id strings
#       reactors: List of reactor id strings
#       stickers: List of sticker id strings
#       diffusers: List of diffuser id strings
#       tstart: Starting time for release
#       tend: Ending time for release
#       tstep: Timestep of release
#       outfile: Filename of output file
#       P0: Initial position of release - Note the convention of particle vertical level Z is negative downards, where sea surface=0m i.e. -10 s 10 m below sea surface
#       spawn: Number of spawned particles (per day)
#       reln: Number of particles per release
#       R0: Total release of material
#       Q0: Flux of material (per day) 
#       **properties: Optional keyword arguments specifying additional properties
      
#     Properties:

# def __init__(self,id,nbuff,movers=[],reactors=[],stickers=[],diffusers=[],tstart=None,tend=None,tstep=0.,tstep_release=0.,outfile=None,P0=[0,0,0],spawn=1,reln=0,R0=1.,Q0=1.,unstick=0.,**prop):



tstart = datetime.datetime(2009,1,1)
tend   = datetime.datetime(2009,1,2)
tstep  = 900.
tout   = 3600.

current = ConstantMover('cur',['uo','vo'],uo=1.0,vo=0.0)
diff    = VariableDiffuser('diff',['diffx','diffy','diffz'],diffz=0.001,P0=[0,0,0])


particles = BuoyantTracer(id = 'particles', nbuff = 10000, 
                          movers=[current], diffusers = [],
                          tstart = tstart,tend = tend, tstep=0.,
                          P0 = [0,0,0], 
                          reln = 1000, 
                          w0 = -0.1)


ercore=ERcore(tout=900., geod = False)


colors=['b+','r+','g+','m+']
# plt.figure()
for rk in range(4,0,-1):
    print 'Running rk%i' % rk 
    part=copy.deepcopy(particles)
    ercore.materials = [part]
    if not os.path.isdir('rk'+str(rk)):os.mkdir('rk'+str(rk))
    ercore.rkorder = rk
    ercore.outpath = 'rk'+str(rk)
    ercore.run(t = tstart ,tend = tend ,dt = tstep)
  #plot(part.pos[:,0],part.pos[:,2],colors[rk-1], marker='o')

#show()
  

  

  
