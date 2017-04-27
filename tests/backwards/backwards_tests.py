# python -m unittest -v backwards_tests.TestCase.test_tide_point_multiple_files
import os
import numpy as np
import datetime
import unittest
from ercore import ERcore
from ercore.fields import ConstantMover,GriddedMover
from ercore.materials import PassiveTracer
from plot_particles import *

def create_gridded_depth (x,y,dep,filedep='dep.nc'):
    import netCDF4
    dep = dep.astype('float32')
    print 'Creating %s' % filedep
    nc = netCDF4.Dataset(filedep,'w')
    nc.createDimension('lon', len(x))
    nc.createDimension('lat', len(y))
    nc.createVariable('lon','float32',('lon'))
    nc.createVariable('lat','float32',('lat'))
    nc.createVariable('dep','float32',('lat','lon'))
    nc.variables['lon'][:] = x
    nc.variables['lat'][:] = y
    nc.variables['dep'][:] = dep
    nc.close()


def create_gridded_current (x,y,u,v,time,
                            units = 'days since 1-1-1', 
                            filecurr='cur.nc'):
    import netCDF4
    u = u.astype('float64')
    v = v.astype('float64')
    time = netCDF4.date2num(time, units)
    time = time.astype('float64')
    print 'Creating %s' % filecurr
    nc = netCDF4.Dataset(filecurr,'w')
    nc.createDimension('lon', len(x))
    nc.createDimension('lat', len(y))
    nc.createDimension('time', len(time))
    nc.createVariable('lon','float64',('lon'))
    nc.createVariable('lat','float64',('lat'))
    t = nc.createVariable('time','float64',('time'))
    t.units = units
    nc.createVariable('u','float32',('time', 'lat','lon'))
    nc.createVariable('v','float32',('time', 'lat','lon'))
    nc.variables['lon'][:] = x
    nc.variables['lat'][:] = y
    nc.variables['time'][:] = time
    nc.variables['u'][:] = u
    nc.variables['v'][:] = v
    nc.close()

def create_tide_current(tstart, tend, dt=3600,
                    units = 'days since 1-1-1', 
                    filecurr='tide_cur.nc'):

    if tstart > tend:
        raise Exception('tstart must be lower than tend in create_tide_current')

    x = np.arange(-100,100)
    y = np.arange(-100,100)
    lent = tend - tstart
    lent = lent.total_seconds()
    times = np.arange(0,lent+dt/2,dt)
    nt = len(times)
    u = np.zeros((nt, len(x),len(y)))
    v = np.zeros(u.shape)

    ut = 0.5/3600.*np.sin(np.linspace(0,2.*np.pi,nt))
    for it in range(nt): u[it,...] = ut[it]

    dtimes = [tstart+datetime.timedelta(seconds=t) for t in times]
    # import pdb;pdb.set_trace()
    create_gridded_current (x=x,y=x,
                            u=u,v=v,time=dtimes,
                            units = units,
                            filecurr=filecurr)
        

class TestCase(unittest.TestCase):

    def setUp(self):
        self.tstart = datetime.datetime(2000,1,3,12)
        self.tend   = datetime.datetime(2000,1,2,12)
        self.dt     = -3600 if self.tstart > self.tend else 3600

    def tearDown(self):
        pass

    def test_linear_point (self):
        imp = 'linear_point'
        if self.dt < 0: imp += '_back'
        print 'Running for imp %s' % imp

        current = ConstantMover('cur',['uo','vo'],uo=0.5/3600.,vo=0.0)
        p1 = PassiveTracer( id=imp,
                            outfile='%s.out' % imp,
                            movers=[current], 
                            reln=10000,
                            nbuff=10, 
                            P0=[0,0,0],
                            tstart = self.tstart,
                            tend =   self.tend)

        ercore=ERcore(geod=False)
        ercore.materials=[p1]
        ercore.run(t=self.tstart,tend=self.tend,dt=self.dt)

        df = read_release_txt('%s.out' % imp)
        print df.head()

        plot_all(df = df,
                colors=['b', 'r', 'g', 'k'],
                lims=None,
                fileplotpref='plt_%s' % imp,
                P0 = [0,0,0],
                polygon = None,
                shoreline = None)

    def test_linear_polygon (self):
        imp = 'linear_poly'
        if self.dt < 0: imp += '_back'
        print 'Running for imp %s' % imp

        current = ConstantMover('cur',['uo','vo'],uo=0.5/3600.,vo=0.0)
        p1 = PassiveTracer( id=imp,
                            outfile='%s.out' % imp,
                            movers=[current], 
                            reln=100,
                            nbuff=1000, 
                            polygon=[[0,0],[0,0.01],[0.01,1],[0,1]],
                            tstart = self.tstart,
                            tend =   self.tend)

        ercore=ERcore(geod=False)
        ercore.materials=[p1]
        ercore.run(t=self.tstart,tend=self.tend,dt=self.dt)

        df = read_release_txt('%s.out' % imp)
        print df.head()

        plot_all(df = df,
                colors=['b', 'r', 'g', 'k'],
                # lims=[-15,1,0,2],
                fileplotpref='plt_%s' % imp,
                polygon = [[0,0],[0,0.01],[0.01,1],[0,1]],
                dens = True)

    def test_linear_polygon2 (self):
        imp = 'linear_poly2'
        if self.dt < 0: imp += '_back'
        print 'Running for imp %s' % imp

        current = ConstantMover('cur',['uo','vo'],uo=0.5/3600.,vo=0.0)
        p1 = PassiveTracer( id=imp,
                            outfile='%s.out' % imp,
                            movers=[current], 
                            reln=1000,
                            nbuff=100, 
                            polygon=[[0,0],[0,0.01],[0.01,1],[0,1]],
                            tstart = self.tstart,
                            tend =   self.tend)

        ercore=ERcore(geod=False)
        ercore.materials=[p1]
        ercore.run(t=self.tstart,tend=self.tend,dt=self.dt)

        df = read_release_txt('%s.out' % imp)
        print df.head()

        plot_all(df = df,
                colors=['b', 'r', 'g', 'k'],
                # lims=[-15,1,0,2],
                fileplotpref='plt_%s' % imp,
                polygon = [[0,0],[0,0.01],[0.01,1],[0,1]],
                dens = True)


    def test_tide_point(self):
        imp = 'tide_point'
        if self.dt < 0: imp += '_back'
        print 'Running for imp %s' % imp
        filecurr = 'tide_cur.nc'
        if not os.path.exists(filecurr):
            print 'Creating tide current file'
            create_tide_current(tstart=min(self.tstart,self.tend)-datetime.timedelta(days=1), 
                                tend=max(self.tstart,self.tend)+datetime.timedelta(days=1), 
                                filecurr=filecurr)

        current=GriddedMover('cur',['u','v'],file=filecurr)
        p1 = PassiveTracer( id=imp,
                            outfile='%s.out' % imp,
                            movers=[current], 
                            reln=10,
                            nbuff=1000, 
                            P0=[0,0,0],
                            tstart = self.tstart,
                            tend =   self.tstart)

        ercore=ERcore(geod=False)
        ercore.materials=[p1]
        # import pdb;pdb.set_trace()
        ercore.run(t=self.tstart,tend=self.tend,dt=self.dt)

        df = read_release_txt('%s.out' % imp)
        print df.head()

        plot_all(df = df,
                lims=None,
                fileplotpref='plt_%s' % imp,
                P0 = [0,0,0],
                polygon = None,
                shoreline = None)

        plot_frames(df= df,
                    usemap=False,
                    lims=None,
                    itout=12,
                    fileplotpref='plt_%s' % imp,
                    P0 = [0,0,0])


    def test_tide_polygon (self):
        imp = 'tide_poly'
        if self.dt < 0: imp += '_back'
        print 'Running for imp %s' % imp

        filecurr = 'tide_cur.nc'
        if not os.path.exists(filecurr):
            print 'Creating tide current file'
            create_tide_current(tstart=min(self.tstart,self.tend)-datetime.timedelta(days=1), 
                                tend=max(self.tstart,self.tend)+datetime.timedelta(days=1), 
                                filecurr=filecurr)

        current=GriddedMover('cur',['u','v'],file=filecurr)

        p1 = PassiveTracer( id=imp,
                            outfile='%s.out' % imp,
                            movers=[current], 
                            reln=100,
                            nbuff=1000, 
                            polygon=[[0,0],[0,0.01],[0.01,1],[0,1]],
                            tstart = self.tstart,
                            tend =   self.tend)

        ercore=ERcore(geod=False)
        ercore.materials=[p1]
        ercore.run(t=self.tstart,tend=self.tend,dt=self.dt)

        df = read_release_txt('%s.out' % imp)
        print df.head()

        plot_all(df = df,
                colors=['b', 'r', 'g', 'k'],
                # lims=[-15,1,0,2],
                fileplotpref='plt_%s' % imp,
                polygon = [[0,0],[0,0.01],[0.01,1],[0,1]],
                dens = True)


        plot_frames(df= df,
                    usemap=False,
                    lims=None,
                    itout=12,
                    fileplotpref='plt_%s' % imp,
                    polygon = [[0,0],[0,0.01],[0.01,1],[0,1]])


    def test_tide_point_multiple_files(self):
        imp = 'tide_point_mfiles'
        if self.dt < 0: imp += '_back'
        print 'Running for imp %s' % imp

        current=GriddedMover('cur',['u','v'],file='tide_cur_%Y%m%d_%Hz.nc')

        p1 = PassiveTracer( id=imp,
                            outfile='%s.out' % imp,
                            movers=[current], 
                            reln=10,
                            nbuff=1000, 
                            P0=[0,0,0],
                            tstart = self.tstart,
                            tend =   self.tstart)

        ercore=ERcore(geod=False)
        ercore.materials=[p1]
        # import pdb;pdb.set_trace()
        ercore.run(t=self.tstart,tend=self.tend,dt=self.dt)

        df = read_release_txt('%s.out' % imp)
        print df.head()

        plot_all(df = df,
                lims=None,
                fileplotpref='plt_%s' % imp,
                P0 = [0,0,0],
                polygon = None,
                shoreline = None)

        plot_frames(df= df,
                    usemap=False,
                    lims=None,
                    itout=12,
                    fileplotpref='plt_%s' % imp,
                    P0 = [0,0,0])


if __name__ == '__main__':
    unittest.main()


