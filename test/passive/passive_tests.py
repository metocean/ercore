#!/usr/bin/env python
import datetime,sys,os,copy
#sys.path.insert(0,os.path.join(os.path.dirname(__file__),'..','..'))
from pylab import figure,plot,show
from ercore import ERcore
from ercore.materials import PassiveTracer,BuoyantTracer
from ercore.fields import ConstantMover,VariableDiffuser

current=ConstantMover('cur',['uo','vo'],uo=1.0,vo=0.0)
diff=VariableDiffuser('diff',['diffx','diffy','diffz'],diffz=0.001,P0=[0,0,0])

particles=BuoyantTracer('particles',10000,geod=False,movers=[current],diffusers=[],stickers=[],reln=1000,P0=[0,0,0],w0=-0.1,tstart=datetime.datetime(2009,1,1),tend=datetime.datetime(2009,1,1))

colors=['b+','r+','g+','m+']

ercore=ERcore(tout=900.)
ercore.geod=False


#figure()
for rk in range(4,0,-1):
  part=copy.deepcopy(particles)
  ercore.materials=[part]
  if not os.path.isdir('rk'+str(rk)):os.mkdir('rk'+str(rk))
  ercore.rkorder=rk
  ercore.outpath='rk'+str(rk)
  ercore.run(datetime.datetime(2009,1,1),datetime.datetime(2009,1,2),900)
  #plot(part.pos[:,0],part.pos[:,2],colors[rk-1], marker='o')

#show()
  

  

  
