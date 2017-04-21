from plot_particles import *

# ts = read_release_txt('p1.dat')
# print ts.dtypes
# ts.to_pickle('p1.pickle')
# print ts.head()
# print ts.tail()

# plot_frames(filein='p1.pickle',
#         usemap = False,
#         lims = [177.15,177.35,-37.95, -37.85],
#         # dmer = 23
#         # dpar = 2,
#         itout = 744,
#         fileplotpref = 'plt_p1',
#         P0 = [177.2147235,-37.87237508,0],
#         polygon = [[177.214978,-37.908414],
#                     [177.214723,-37.872375],
#                     [177.322697,-37.871846],
#                     [177.322995,-37.907884],
#                     [177.214978,-37.908414]]
# )

# plot_all(filein='p1.pickle',
#         lims = [177.15,177.35,-37.95, -37.85],
#         fileplotpref = 'plt_p1',
#         P0 = [177.2147235,-37.87237508,0],
#         polygon = [[177.214978,-37.908414],
#                     [177.214723,-37.872375],
#                     [177.322697,-37.871846],
#                     [177.322995,-37.907884],
#                     [177.214978,-37.908414]]
# )

############################################
# ts = read_release_txt('p2.dat')
# print ts.dtypes
# ts.to_pickle('p2.pickle')
# print ts.head()
# print ts.tail()

# plot_all(filein='p2.pickle',
#         lims = [177.15,177.35,-37.95, -37.85],
#         fileplotpref = 'plt_p2',
#         P0 = [177.2147235,-37.87237508,0],
#         polygon = [[177.214978,-37.908414],
#                     [177.214723,-37.872375],
#                     [177.322697,-37.871846],
#                     [177.322995,-37.907884],
#                     [177.214978,-37.908414]]
# )

############################################
# ts = read_release_txt('p3.dat')
# print ts.dtypes
# ts.to_pickle('p3.pickle')
# print ts.head()
# print ts.tail()

# plot_all(filein='p3.pickle',
#         lims = [177.15,177.35,-37.95, -37.85],
#         fileplotpref = 'plt_p3',
#         P0 = [177.2147235,-37.87237508,0],
#         polygon = [[177.214978,-37.908414],
#                     [177.214723,-37.872375],
#                     [177.322697,-37.871846],
#                     [177.322995,-37.907884],
#                     [177.214978,-37.908414]],
#         shoreline = 'shoreline.bnd'
# )

############################################
# ts = read_release_txt('p4.dat')
# print ts.dtypes
# ts.to_pickle('p4.pickle')
# print ts.head()
# print ts.tail()
# import yaml
# config = yaml.load(open('bop4.yml'))

# plot_all(filein='p4.pickle',
#         lims = [177.15,177.35,-37.95, -37.85],
#         fileplotpref = 'plt_p4',
#         P0 = [177.2147235,-37.87237508,0],
#         polygon = [[177.214978,-37.908414],
#                     [177.214723,-37.872375],
#                     [177.322697,-37.871846],
#                     [177.322995,-37.907884],
#                     [177.214978,-37.908414]],
#         shoreline = 'shoreline.bnd'
# )

############################################

import sys

if len(sys.argv) < 2:
    print 'usage: %prog imp'
    sys.exit(-1)

imp = sys.argv[1] #'bop7_back'
print 'imp = ', imp
itout = int(sys.argv[2]) if len(sys.argv)>=2 else 6
print 'itout = ', itout

ts = read_release_txt('%s.dat' % imp)
print ts.dtypes
ts.to_pickle('%s.pickle' % imp)
print ts.head()
print ts.tail()
import yaml
config = yaml.load(open('%s.yml' % imp))

plot_all(filein='%s.pickle' % imp,
        #lims = [177.15,177.35,-37.95, -37.85],
	lims = [176.6,177.5,-38,-37.4],
        fileplotpref = 'plt_%s' % imp,
        P0 = [177.2147235,-37.87237508,0],
        polygon = [[177.214978,-37.908414],
                    [177.214723,-37.872375],
                    [177.322697,-37.871846],
                    [177.322995,-37.907884],
                    [177.214978,-37.908414]],
        shoreline = 'shoreline.bnd'
)

lims = [176.6,177.5,-38,-37.4]
#lims = [177.15,177.35,-37.95, -37.85]


import matplotlib.pyplot as plt 
import netCDF4
nc = netCDF4.Dataset('bop_currents.nc')
lon = nc.variables['lon'][:]
lat = nc.variables['lat'][:]
dep = nc.variables['dep'][:]

from ercore.shoreline import Shoreline
shore = Shoreline('shore', 'shoreline.bnd')

from shapely.geometry import Polygon
poly = Polygon(config['materials'][0]['polygon'])

plt.ion()
fig = plt.figure()
df = ts
times = df.index.get_level_values('time').unique()
releases = df.index.get_level_values('release').unique()

nt = len(times)
its = range(0,nt,itout)

for it in its:
    time = times[it]
    dft = df.loc[(df.index.get_level_values('time') == time)]
    fig.clear()
    plt.axis(lims)
    title_str = 'time = %s (%d/%d)' % (time, it, nt)
    print title_str
    plt.title(title_str)

    for i,rel in enumerate(releases):
        dftr= dft.loc[dft.index.get_level_values('release') == rel]
        plt.scatter(dftr['x'], dftr['y'], c='grey', edgecolor='', alpha=0.5)
        # plt.plot(ts['x'], ts['y'], c='grey', marker='.',ls='None', alpha=0.5)

    for i,j in enumerate(shore.polyi):
        n=shore.polyn[i]
        plt.plot(shore.slx[j:j+n-1],shore.sly[j:j+n-1],'k')

    x,y = poly.exterior.xy
    plt.plot(x,y, 'k--', lw=2)    

    cs = plt.contour(lon,lat,dep, range(0,100,5))
    # plt.clabel(CS)#, inline=1)#, fontsize=10)
    plt.colorbar()

    # plt.waitforbuttonpress()
    fileplot = 'plt_%s_it%03i.png' % (imp, it)
    plt.savefig(fileplot)

