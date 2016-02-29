## plot all pasrticles by timestamp ##
import sys; sys.path += ["/source/ercore"]

import numpy
import matplotlib.pyplot as plt
from ercore.shoreline import Shoreline
from ercore.fields import GriddedTopo
from ercore import dt2ncep,ncep2dt

plt_acc = True
# P0 = [174.655, -36.818, 0] # too s hallow
# P0 = [174.686, -36.8655, 0] # near bridge
P0 =[174.6884,-36.8366,0] # andre
#bnd = [174.56873198401499,174.82807580347946,-36.927971300850707,-36.717879213428184]
filep = 'passive_akl.txt'
file_shore='shoreline_akl_eri.bnd'

def plot_shoreline():
    plt.figure()
    for i,j in enumerate(shoreline.polyi):
        n=shoreline.polyn[i]
        plt.plot(shoreline.slx[j:j+n-1],shoreline.sly[j:j+n-1],'k')        
    # plt.axis(bnd)
    plt.tricontourf(triang, dep) #, cmap=cmap, levels=levels)
    plt.plot(P0[0],P0[1],'r+')
    plt.colorbar()


#dep=GriddedTopo('depth',['dep'], file='Auckland_Cons_3D.nc', zinvert=True)
#d = dep.get(time=dt2ncep(datetime.datetime(2005,1,1)))[0]
import netCDF4
from matplotlib.tri import Triangulation
with netCDF4.Dataset('Auckland_Cons_3D.nc') as f:
    lon = f.variables['lon'][:]
    lat = f.variables['lat'][:]
    nv  = f.variables['nv'][:]-1
    dep = f.variables['dep'][:]
    triang = Triangulation(lon,lat,nv)



print 'Reading ', filep
with open(filep) as f: lines = f.readlines()
shoreline = Shoreline(id='shore_nz', file=file_shore)

first_time = True
times = []
told = None
for lin in lines[1:]:
    time,pid,px,py,pz,pstate,page,pmass,zbot,elev = map(float,lin.split('\t'))
    tc = ncep2dt(time)
    if tc != told:
        print 'Plotting %.3f %s' % (time, tc.strftime('%Y-%m-%d %H:%M'))        
        # import pdb; pdb.set_trace()
        if plt_acc and told is None: plot_shoreline()        
        plt.show()
        raw_input("Press Enter to continue...")        
        if not plt_acc: plt.close()
        told = tc
    print px,py
    plt.plot(px,py, 'ko')

plt.show()
plt.savefig('plt_part.png')
# plt.close()
