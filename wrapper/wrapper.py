import os
import sys
import copy
import datetime as dt
import logging
import zipfile
import yaml
import csv
import json
import time

from billiard import Process, Pool
from scipy.stats import gaussian_kde
import numpy as np
from matplotlib.figure import Figure


from msl_actions import read_xml, catch_exception
from msl_actions.download import UDSQuery

from ercore.materials import *
from ercore.fields import *
from ercore.shoreline import Shoreline
from ercore import ERcore

_DT0_=dt.datetime(2000,1,1)
_NCEPT0_=730121.0
_DTROUND=10
_NSECROUND=86400/_DTROUND

ncep2dt=lambda t:_DT0_+dt.timedelta(0,_DTROUND*round(_NSECROUND*(t-_NCEPT0_)))

class KernelDensity(object):
    """docstring for KernelDensity"""
    def __init__(self, release_id, outfile, outdir='', res=100, merge=False):
        super(KernelDensity, self).__init__()
        self.release_id = release_id
        self.res = res
        self.merge = merge
        self.outfile = outfile
        self.outdir = outdir
        self.read_output()
        if self.data:
            self.generate_grid()
            self.compute_density()
            #self.mask_result()
            #self.save_kde()

    def read_output(self):
        self.timelist = []
        self.data = {}
        with open(self.outfile, 'r') as openf:
            csvfile = csv.reader(openf, delimiter='\t')
            csvfile.next()
            for row in csvfile:
                timestep = float(row[0])
                if timestep not in self.timelist:
                    self.timelist.append(timestep)
                    self.data[timestep] = [map(float, (row[2], row[3]))]
                else:
                    self.data[timestep].append(map(float, (row[2], row[3]))) 

    def generate_grid(self):
        self.grids = {}
        self.points = {}
        self.gridsflat = {}
        for tstep in self.timelist:
            self.grids[tstep] = {}
            self.points[tstep] = np.array(self.data[tstep]).T
            buff = np.sqrt(self.points[tstep][0].std()**2 + self.points[tstep][1].std()**2 )    
            xflat = np.linspace(self.points[tstep][0].min()-buff,
                                self.points[tstep][0].max()+buff, self.res)
            yflat = np.linspace(self.points[tstep][1].min()-buff, 
                                self.points[tstep][1].max()+buff, self.res)
            X, Y = np.meshgrid(xflat, yflat)
            self.grids[tstep]['X'], self.grids[tstep]['Y'] = X,Y
            self.gridsflat[tstep] = np.vstack([X.ravel(), Y.ravel()])

    def compute_density(self):
        for tstep in self.timelist:
            kde = gaussian_kde(self.points[tstep])
            X = self.grids[tstep]['X']
            self.grids[tstep]['Z'] = np.reshape(kde(self.gridsflat[tstep]).T, X.shape)

    def save_json(self):
        fig = Figure()
        ax = fig.gca()
        output = {}
        for tstep in self.timelist:
            datestr = ncep2dt(tstep).strftime('%Y/%m/%d %H:%M:%S')
            output[datestr] = []
            x, y, z = self.grids[tstep]['X'], self.grids[tstep]['Y'],self.grids[tstep]['Z']
            contours = ax.contour(x,y,z)
            for c, line in enumerate(contours.collections):
                val = contours.cvalues[c]
                vertices = line.get_paths()[0].vertices
                x = vertices[:,0]
                y = vertices[:,1]
                output[datestr].append({'fac':z.max(), 
                                        'val':val,
                                        'np':len(x),
                                        'x':x.round(4).tolist(), 
                                        'y':y.round(4).tolist()})
        jsonfile = os.path.join(self.outdir, '%s.json' % self.release_id)
        with open(jsonfile, 'w') as jfile:
            json.dump(output, jfile)

def compute_kde(release_id, config_file, outdir):
    with open(config_file) as openfile:
        config = yaml.load(openfile)
    for material in config['materials']:        
        outfile = os.path.join(outdir, material['outfile'])
        kde = KernelDensity(release_id, outfile, outdir)
        kde.save_json()

def run_ercore(release_id, tstart, tstop, trun, tout, outdir, config_file):
    ercore_obj = ERcore(geod=True, release_id=release_id, tout=tout, outpath=outdir)
    ercore_obj.readYAML(config_file, globals())
    ercore_obj.run(t=tstart, tend=tstop, dt=trun)
    compute_kde(release_id, config_file, outdir)

class ERCoreWrapper(object):
    """Wraps the ercore model for one implementation,
    each implementation corresponds to area to be computed 
    and can have multiple sites and for each site is generated
    a kernel density estimation and saved as json format to
    be ingested by MOV.

    Arguments:
    id: short name id
    name: the implementation name
    rootdir: the final data root directory
    rundir: The run directory. A.k.a scratch directory
    static: Static directory with the static files
    boundary: list as [x0,y0,x1,y1]
    uds: dict
        host: the udshost to gather data from
        datasets: list of uds datasets
        vars: list of vars to get

    duration: Number of days the simulation will run for
    release_duration: Number of hours will release from
    delta_out: Number of seconds a simulation output will be written at
    delta_run: Simulation time step in seconds
    releases: dict with sites to be computed, Check Release class to see what values to set
    movers:  dict 
    reactors: dict
    diffuser: dict with diffusing parameters
    stickers: dict of stickers to use, defaul to ground and shoreline
    levels: list of levels to use in case of 3D simulation, default to surface only

    """
    def __init__(self, id, name, boundary, sitesfile, 
                 rootdir='/data/ercore',
                 staticdir= '/static/ercore',
                 rundir='/scratch/ercore',
                 uds={},
                 duration=4,
                 delta_out=1800,
                 delta_run=600,
                 levels=[0],
                 currents={},
                 wind={},
                 shoreline={},
                 diffuser={},
                 pipe=None, 
                 logger=logging, 
                 **kwargs):
        super(ERCoreWrapper, self).__init__()
        self.movers = { 'currents': {'id': 'currents',
                                        'class' : 'GriddedMover',
                                        'topo': 'bathy',
                                        'file'  : None,
                                        'vars'  : ['uso', 'vso']}}

        self.diffusers = {'diffuser': {'id' : 'diffuser',
                              'class':'ConstantDiffuser',
                              'diffx': 0.1,
                              'diffy': 0.1,
                              'diffz': 0.001,
                              'vars' : ['diffx', 'diffy', 'diffz']}}

        self.stickers = {'shoreline': { 'id':'shoreline',
                                           'class': 'Shoreline',
                                           'file':None}}

        self.reactors = { 'wind': {'id': 'wind',
                                   'class' : 'GriddedReactor',
                                   'topo' : 'bathy',
                                   'file'  : None,
                                   'vars'  : ['ugrd10m', 'vgrd10m']}}

        self.bathy = {'bathy': { 'id':'bathy',
                                 'class': 'ConstantTopo',
                                 'vars' : 'depth',
                                 'depth': 10}}

        self.default_material = {
            'class' : 'BuoyantTracer',
            'movers' : ['currents'],
            'reactors':['wind'],
            'stickers': ['shoreline'],
            'nbuff': 1000,
            'reln': 1000,
            'P0': None,
            'tstart': None,
            'tend': None,
            'outfile':None}

        self.uds = {
            'host': None,
            'datasets': ['gfs_sfc', 'rtofs_sfc'],
            'vars': ['ugrd10m', 'vgrd10m', 'uso','vso'],
            'outfile': None,
            'retry':1}

        self.uds.update(uds)

        self.outputs = {}
        self.id = id
        self.name = name
        self.sitesfile = sitesfile
        self.rootdir = rootdir
        self.staticdir = staticdir
        self.rundir = rundir
        self.boundary = boundary
        self.duration = duration
        self.delta_out = delta_out
        self.delta_run = delta_run
        self.levels = levels
        self.wind = wind
        self.currents = currents
        self.shoreline = shoreline
        self.diffuser = diffuser
        self.logger = logger
        self.pipe = pipe
        self.pool = None

    def set_cycle(self, cycle_dt):
        self.cycle_dt = cycle_dt
        self.cycle_str = cycle_dt.strftime('%Y%m%d_%Hz')
        self.rundir = os.path.join(self.rundir, self.id, self.cycle_str)
        self.rootdir = os.path.join(self.rootdir, self.id, self.cycle_str)
        self.start_time = self.cycle_dt
        self.end_time = self.cycle_dt + dt.timedelta(days=self.duration)
        self.outdir = os.path.join(self.rundir, 'out')
        self.indir = os.path.join(self.rundir, 'in')
        self.staticdir = os.path.join(self.rundir, 'static')
        
    def _makedirs(self):
        for path in [self.outdir, self.indir, self.rundir]:
            try:
                os.makedirs(path)
            except:
                pass

    def _get_movers_query(self):
        query = {
            'dset':self.uds['datasets'],
            'var':self.uds['vars'],
            'time': [self.start_time.strftime('%Y%m%d.%H%M%S'),
                     self.end_time.strftime('%Y%m%d.%H%M%S')],
            'type': 'fc,hc',
            'fmt': 'nc',
            'datum': 'msl',
            'sort': ['quality', '-res'],
            'spinup': '0.0',
            'stepback': 0,
            'nomissing':True,
            'slfmt':'bnd',
            'sl':'gshhs_h',
            'bnd': ["%.6f" % b for b in self.boundary],
        }
        return query

    def get_input_data(self):
        udsout = "uds_%s_%s.zip" % (self.id, self.cycle_str)
        query = self._get_movers_query()
        udsquery = UDSQuery(query, udshost=self.uds['host'], 
                            local_dir=self.indir, 
                            outfile=udsout,
                            retry=self.uds['retry'],
                            logger=self.logger)
        udsquery.run()
        zipfile.ZipFile(os.path.join(self.indir,udsout)).extractall(self.indir)
        self.inputfile = os.path.join(self.indir, 'uds_out.nc')
        self.shorelinefile = os.path.join(self.indir, 'shoreline.bnd')


    def get_ercore_config(self):
        self.movers['currents']['file'] = self.inputfile
        self.reactors['wind']['file'] = self.inputfile
        self.movers['currents'].update(self.currents)
        self.reactors['wind'].update(self.wind)
        self.stickers['shoreline']['file'] = self.shorelinefile
        self.stickers['shoreline'].update(self.shoreline)
        self.diffusers['diffuser'].update(self.diffuser)

        config = {
            'diffusers': self.diffusers.values(),
            'movers' : self.movers.values(),
            'stickers': self.stickers.values(),
            'reactors': self.reactors.values(),
            'topo' : self.bathy.values()
        }
    
        return config

    def get_material(self, site):
        this_material = copy.deepcopy(self.default_material)
        material = {
            'id' : site['id'],
            'tstart': self.cycle_dt.strftime('%Y-%m-%d %H:%M:%S'),
            'tend': (self.cycle_dt+dt.timedelta(hours=site['duration'])).strftime('%Y-%m-%d %H:%M:%S'),
            'P0': [site['lon'], site['lat'], site['depth']],
            'outfile': "%s_%s.txt" % (site['id'],self.cycle_str)}

        this_material.update(material)
        return this_material

    def get_sites(self):
        ingrid_sites = []
        x0,x1,y0,y1 = self.boundary
        sites = read_xml(self.sitesfile)
        for site in sites.getchildren():
            site_id = site.attrib['id']
            name = site.attrib['name']
            lon = float(site.attrib['x'])
            lat = float(site.attrib['y'])
            depth = float(site.attrib['z']) if site.attrib.has_key('z') else 0
            runlenght = float(site.attrib['z']) if site.attrib.has_key('z') else 0
            duration = (float(site.attrib['rend'])/60/60) if site.attrib.has_key('rend') else 6
            if x0 <= lon and lon <= x1 and\
               y0 <= lat and lat <= y1:
               ingrid_sites.append({'id' : site_id,
                                    'name': name,
                                    'lon' : lon,
                                    'lat' : lat,
                                    'depth': depth,
                                    'duration' : duration})
        return ingrid_sites
            
    def preprocess(self):
        self.get_input_data()
        sites = self.get_sites()
        self.config_files = {}
        self.configs = {}
        for site in sites:
            config = self.get_ercore_config()
            material = self.get_material(site)
            config['materials'] = [material]
            configfile = self._dump_config(config)
            self.config_files[site['id']] = configfile
            self.configs[site['id']] = config

        # Gather input data
        # Generate config files
        

    def _dump_config(self, config):
        site_id = config['materials'][0]['id']
        configpath = os.path.join(self.rundir, '%s_%s.yaml' % (site_id, self.cycle_str))
        with open(configpath, 'w') as configfile:
            yaml.dump(config, configfile)
        return configpath

    def set_mpiexec(self, ncores, hosts=None):
        self.pool = Pool(ncores)

    def run_model(self):
        if not self.pool: 
            self.pool = Pool(1)
        self.processes = {}
        for site_id, config_file in self.config_files.items():
            # run_ercore(site_id, 
            #           self.start_time,
            #           self.end_time, 
            #           self.delta_run,
            #           self.delta_out,
            #           self.outdir,
            #           config_file)
            self.processes[site_id] = self.pool.apply_async(run_ercore, 
                                                            args=[site_id, 
                                                                  self.start_time,
                                                                  self.end_time, 
                                                                  self.delta_run,
                                                                  self.delta_out,
                                                                  self.outdir,
                                                                  config_file])
        total = len(self.processes.keys())
        success = []
        while self.preprocess:
            for site_id, process in self.processes.items():
                if process.ready():
                    success.append(self.processes.pop(site_id))
                    result = process.get()
                    self.set_progress(len(success), total)
                else:
                    continue
            time.sleep(2)

    def set_progress(self, partial, total):
        if hasattr(self.pipe, 'send'):
            self.pipe.send((partial, total))

    def run(self):
        try:
            self._makedirs()
            self.preprocess()
            self.run_model()
        except Exception as exc:
            tb = catch_exception(exc)
            if hasattr(self.pipe, 'send'):
                self.pipe.send((tb))
            raise
