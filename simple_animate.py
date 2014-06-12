#!/usr/bin/env python
from pylab import *
from ercore import ncep2dt
import sys,time

ion()

f=open(sys.argv[1])
f.readline()

t0=-1
fig=figure()
show()

showshore=len(sys.argv)>2
if showshore:
  from ercore.shoreline import Shoreline
  shoreline=Shoreline('shoreline',sys.argv[2])
  for i,j in enumerate(shoreline.polyi):
    n=shoreline.polyn[i]
    plot(shoreline.slx[j:j+n-1],shoreline.sly[j:j+n-1],'k')

xstack=[]
ystack=[]

for line in f.readlines():
  bits=line.split()
  t,x,y,z,state,age,mass=map(float,bits)
  if (t!=t0) and state>0:
    p=plot(xstack,ystack,'.')
    if t0>700000:
      title(ncep2dt(t0).strftime('%Y-%m-%d %H:%M:%S'))
    else:
      title('T+%f' % (t))
    draw()
    p.pop(0).remove()
    t0=t
    time.sleep(0)
    xstack=[]
    ystack=[]
  else:
    xstack.append(x)
    ystack.append(y)
  


