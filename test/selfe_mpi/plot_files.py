## plot all pasrticles by timestamp ##
import sys; sys.path += ["/source/ercore"]

import numpy
import datetime 
import netCDF4
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

from ercore.shoreline import Shoreline
from ercore.fields import GriddedTopo,TidalMover
from ercore import dt2ncep,ncep2dt
from ercore.lib.tide import TideStr


# f = netCDF4.Dataset('tide_akl.nc')
# lon = f.variables['lon'][:]
# lat = f.variables['lat'][:]
# nv  = f.variables['nv'][:]-1
# dep = f.variables['dep'][:]
# triang = Triangulation(lon,lat,nv)

# depth_nz=GriddedTopo('depth_nz',['dep'], file='tide_nz.nc', zinvert=True)
# depth_tide_akl=GriddedTopo('depth_tide_akl',['dep'], file='tide_akl.nc', zinvert=True)
shoreline = Shoreline(id='shore_nz', file='NZ(WGS84).bnd') #shoreline_akl_eri.bnd')  #NZ(WGS84).bnd,

f = netCDF4.Dataset('tide_akl.nc')
lon = f.variables['lon'][:]
lat = f.variables['lat'][:]
nv  = f.variables['nv'][:]-1
dep = f.variables['dep'][:]
triang = Triangulation(lon,lat,nv)
cons = f.variables['cons'][:]
cons = ['M2','S2']
t0 = datetime.datetime.strptime(f.t0str, '%Y-%m-%d %H:%M:%S')
lat0 = f.lat0

lev = f.variables['lev'][:]
k=0 # surface
u_amp = f.variables['u_amp'][1:3,k,:]
u_pha = f.variables['u_pha'][1:3,k,:]
ustr=TideStr(u_amp,u_pha,cons,t0,lat0)

v_amp = f.variables['v_amp'][1:3,k,:]
v_pha = f.variables['v_pha'][1:3,k,:]
vstr=TideStr(v_amp,v_pha,cons,t0,lat0)

f.close()

f1 = netCDF4.Dataset('tide_nz.nc')
lon1 = f1.variables['lon'][:]
lat1 = f1.variables['lat'][:]
nv1  = f1.variables['nv'][:]-1
dep1 = f1.variables['dep'][:]
triang1 = Triangulation(lon1,lat1,nv1)
cons = f1.variables['cons'][:]
cons = ['M2','S2']
t01 = datetime.datetime.strptime(f1.t0str, '%Y-%m-%d %H:%M:%S')
lat01 = f1.lat0

lev1 = f1.variables['lev'][:]
k=0 # surface
u_amp1 = f1.variables['u_amp'][1:3,k,:]
u_pha1 = f1.variables['u_pha'][1:3,k,:]
ustr1=TideStr(u_amp1,u_pha1,cons,t01,lat01)

v_amp1 = f1.variables['v_amp'][1:3,k,:]
v_pha1 = f1.variables['v_pha'][1:3,k,:]
vstr1=TideStr(v_amp1,v_pha1,cons,t01,lat01)


# px = 174.705856
# py = -36.850966
# t = ncep2dt(731586.270833)
# # t = ncep2dt(731586.312500)

px,py = 174.705634, -36.851551
t = ncep2dt(731586.000000)
print t 

u = ustr.ts(t)[0]
v = vstr.ts(t)[0]
ws = numpy.sqrt(u*u + v*v)
print '%.3f, %.3f' % (u.min(), u.max())
print '%.3f, %.3f' % (v.min(), v.max())
print '%.3f, %.3f' % (ws.min(), ws.max())


u1 = ustr1.ts(t)[0]
v1 = vstr1.ts(t)[0]
ws1 = numpy.sqrt(u1*u1 + v1*v1)
print '%.3f, %.3f' % (u1.min(), u1.max())
print '%.3f, %.3f' % (v1.min(), v1.max())
print '%.3f, %.3f' % (ws1.min(), ws1.max())


## PLOT ##
levels = [0,0.001,0.15,0.3,0.45,0.60,0.75,1]# MaxNLocator(nbins=15).tick_values(vmin,vmax);
print levels
cmap = cm.get_cmap(name='jet')

plt.figure(figsize=(10,5))
plt.subplot(121)
plt.tricontourf(triang, ws, cmap=cmap, levels=levels)
plt.colorbar()
plt.plot(px,py, 'r^')
for i,j in enumerate(shoreline.polyi):
    n=shoreline.polyn[i]
    plt.plot(shoreline.slx[j:j+n-1],shoreline.sly[j:j+n-1],'grey')        
bnd = [174.70062597520379,174.71001584256734,-36.85465075572332, -36.849083065986633]
plt.axis(bnd)
plt.title('akl')
plt.subplot(122)
plt.tricontourf(triang1, ws1, cmap=cmap, levels=levels)
plt.colorbar()
plt.plot(px,py, 'r^')
plt.axis(bnd)
plt.title('nz')



plt.scatter(lon, lat, c=ws1, s=20, edgecolors='none')


