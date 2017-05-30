#!/usr/bin/env python
import datetime,os,numpy,re,copy
import base64
import gc 

_DT0_=datetime.datetime(2000,1,1)
_NCEPT0_=730120.99999
# conversion of fraction of days to NCEP/CF convention
ncep2dt=lambda t:_DT0_+datetime.timedelta(t-_NCEPT0_)
# conversion of datetime times to fraction of days
dt2ncep=lambda t: (1.+t.toordinal()+t.hour/24.+t.minute/1440.+t.second/86400.) if isinstance(t,datetime.datetime) else t

def get_summary_report(summary, dt, tout, outdir, materials):
  total_nrel = sum([m.reln for m in materials])
  materials_id = ','.join(['%d - %s (%s)' % (m.reln,m.id, m.__class__.__name__) for m in materials])
  report = """
===============================================
  ErCore Model Summary
--------------------------------------------
  Directory output :          %s
  Model time step:            %ss
  Output frequency:           %ss
  Materials:                  %s
  Total particles released:   %d
--------------------------------------------
  Time Summary
--------------------------------------------
  Start time:              %s
  End time:                %s
  Initialization time:     %0.5fs
  Avg. step time:          %0.5fs
  Avg. release time:       %0.5fs
  Avg. react time:         %0.5fs
  Avg. advect time:        %0.5fs
  Avg. diffuse time:       %0.5fs
  Avg. stick time:         %0.5fs
  Avg. ageing time:        %0.5fs
  Avg. dieing time:        %0.5fs
  Total modeling time:     %s
===============================================

+""" % (outdir, dt, tout, materials_id, total_nrel,
       summary['start_time'], summary['end_time'], summary['initialize'],
       summary['tstep'], summary['release'], summary['react'],
       summary['advect'], summary['diffuse'], summary['stick'],
       summary['ageing'], summary['die'],
       datetime.timedelta(seconds=summary['total_time']))
  return report

def parsetime(t):
    if isinstance(t,datetime.datetime):
        return dt2ncep(t)
    elif isinstance(t,(str, unicode)):
        try:
            return dt2ncep(datetime.datetime.strptime(t,'%Y%m%d_%Hz'))
        except:
            return dt2ncep(datetime.datetime.strptime(t,'%Y-%m-%d %H:%M:%S'))
    else:
        return t

#Decorator to allow inheritance of docstrings
def copydoc(fromfunc, sep="\n"):
    def _decorator(func):
        sourcedoc = fromfunc.__doc__
        if func.__doc__ == None:
            func.__doc__ = sourcedoc
        else:
            func.__doc__ = sep.join([sourcedoc, func.__doc__])
        return func
    return _decorator

class TerminateException(Exception): pass

#Enhanced list class to allow indexing by id
#All simulation objects (i.e. movers, materials etc. must have an id attribute)
class ObjectList(list):
    def __init__(self,items,copy=False):
        if not isinstance(items,list):items=[items]
        for i in items:
            if not hasattr(i,'id'):raise ERConfigException('Each item must have an id attribute')
        if copy:
            import copy
            items=copy.deepcopy(items)
        list.__init__(self,items)
        
    def __getitem__(self,key):
        if key is None:return None
        #if isinstance(key,str):
        if isinstance(key,(str, unicode)):
            idlist=self.idList()
            for idl in idlist:
                if key in idl:return self[idlist.index(idl)]
        elif isinstance(key,int):
            if key<len(self):return list.__getitem__(self,key)
        return None
            
    def __str__(self):
        return "\n".join([str(n) for n in self])
            
    def subset(self,keys):
        #if isinstance(keys,str):keys=[keys]
        if isinstance(keys,(str,unicode)):keys=[keys]
        idlist=self.idList()
        return [self[idlist.index(key)] for key in keys if key in idlist]
                        
    def sort(self,attr='id',reverse=False):
        sortkey=lambda o:getattr(o,attr)
        list.sort(self,key=sortkey,reverse=reverse)
    
    def idList(self):
        return [ds.id for ds in self]

class ERCoreException(Exception):
  pass
  
class ERConfigException(ERCoreException):
  pass

class ERRuntimeException(ERCoreException):
  pass


#Main class for ERcore simulation
class ERcore(object):
  geod=True
  tout=0.
  rkorder=4
  release_id = 0
  stopped = False
  def __init__(self,manager=None,**k):
    self.manager = manager
    self.fout={}
    self.outpath=k.get('outpath','.')
    self.save_summary = True
    self.__dict__.update(k)

    self.summary = {'initialize':None, 'release':None, 'advect':None, 
                    'stick':None, 'diffuse':None,'stick':None,'spawn':None,
                    'total_time':None,'end_time':None,'tstep':None}

  def __exit__(self,*args):
    if hasattr(self, 'materials'):
      del self.materials
    gc.collect()

  def __enter__(self):
    return self
 

    
  def readYAML(self,ctlfile, namespace=globals()):
    """Read a YAML configuration file for an ERcore simulation"""
    import yaml,inspect
    config = yaml.load(open(ctlfile).read())
    objects = ObjectList([])
    for key in config:
      conf=config[key]
      if not isinstance(conf,list):conf=[conf]
      for c in conf:
        if not c.has_key('id'):raise ERConfigException('Every configuration object must specify a unique id')
        if not c.has_key('class'):raise ERConfigException('Every configuration object must specify its class')
        if not namespace.has_key(c['class']):raise ERConfigException('Class %s not defined in program scope' % (c['class']))
        try:
          objects.append(namespace[c['class']](**c))
        except:
          import traceback
          raise ERConfigException('Configuration file error for object %s\n%s' % (c['id'],traceback.format_exc()))
    # intialize materials
    for obj in objects:  
        for prop in dir(obj):
            if prop in ['movers','reactors','diffusers','stickers','topo','members']:
                val=getattr(obj,prop)
                if isinstance(val,list):
                    val=ObjectList([objects[v] if isinstance(v,str) and objects[v] else v for v in val])
                    for v in val:
                        if isinstance(v,str):raise ERConfigException('Cannot find one of %s with id(s) %s specifed for %s' % (prop,v,obj.id))
                #elif isinstance(val,str) and objects[val]:
                elif isinstance(val,(str, unicode)) and objects[val]:
                    val=objects[val]
                # elif val is None: # allow for case when no topo are defined 
                #     val=objects[val]
                else:
                    raise ERConfigException('Cannot find one of %s with id(s) %s specifed for %s' % (prop,val,obj.id))
                setattr(objects[obj.id],prop,val)
    import materials
    import inspect
    self.materials=[o for o in objects if materials._Material in inspect.getmro(o.__class__)]

  def timestamp(self, section, start=None, avg=True):
    if self.save_summary:
      if not start:
        self.summary[section] = '%s' % datetime.datetime.now()
      else:
        start = start or datetime.datetime.now()
        delta = datetime.datetime.now()-start
        total_seconds=(delta.microseconds + (delta.seconds + delta.days * 24 * 3600) * 10**6) / 10**6 
        # see https://docs.python.org/2/library/datetime.html#datetime.timedelta.total_seconds
        if section not in self.summary or self.summary[section] is None:
          #self.summary[section] = delta.total_seconds() # delta.total_seconds() only works from Python 2.7 onwards
          self.summary[section] =total_seconds
        else:
          #self.summary[section] = numpy.average([self.summary[section],delta.total_seconds()]) # delta.total_seconds() only works from Python 2.7 onwards
          self.summary[section] = numpy.average([self.summary[section],total_seconds]) 
    return datetime.datetime.now() 


  def run(self,t,tend,dt,keep_sticked=False):
    """Run Ercore
    Arguments:
        t: Start time as datetime, ncep decimal time or string
        tend: Start time as datetime, ncep decimal time or string
        dt: Simulation timestep as seconds
    """
    start_time = self.timestamp('start_time')
    #t=parsetime(t)
    t=t0=parsetime(t)
    tend=parsetime(tend)
    #outofcompute = -1 if keep_sticked else 0
    outofcompute = -1 # this set the flag value for particle that needs to be removed from the computational pool
    iprint=int(self.tout/dt) if self.tout>0 else 1
    print iprint
    dt/=86400.
    etypes=[e.id for e in self.materials]
    i=0
    for e in self.materials:
      self.fout[e.id]=open(os.path.join(self.outpath,e.outfile),'w')
      self.fout[e.id].write(e.fheader())
      e.initialize(t,tend)
      print '%s: Initialized %s' % (ncep2dt(t).strftime('%Y%m%d %H:%M:%S'), e.id)
      if self.geod:e.geodcalc(True)
    print 'Running ERcore times %f to %f' % (t,tend)
    last_step = self.timestamp('initialize', start_time)
    while ((dt>0) and (t<tend)) or ((dt<0) and (t>tend)):
      if self.stopped:
        raise TerminateException('Model has been terminated')
      start_step = datetime.datetime.now() 
      i+=1
      t2=t+dt
      if self.stopped:
        raise TerminateException('Model has been terminated')
      start_step = datetime.datetime.now()
      if self.manager:
        r, rt = self.run_n
        total = (tend-t0)*rt
        partial = ((tend-t0)*(r-1))+(t-t0)
        self.manager.progress.emit(self.release_id, (partial/total)) 
      for e in self.materials:
        if e.tcum>=e.tstep:
          e.tcum-=e.tstep
          # release only if the material is not a child of another material
          # the release of spawned material is taken care of further below
          if not e.props['ischild']:
            e.release(t,t2)
          last_time = self.timestamp('release', start_step)
          if self.geod:e.geodcalc()
          e.react(t,t2)
          last_time = self.timestamp('react', last_time)
          if (e.state[:e.np]>0).any():
            e.advect(t,t2,self.rkorder)
            last_time = self.timestamp('advect', last_time)
            e.diffuse(t,t2)
            last_time = self.timestamp('diffuse', last_time)
            e.stick(t,t2)
            last_time = self.timestamp('stick', last_time)
          e.spawn(t,t2)
          last_time = self.timestamp('spawn', last_time)
        #e.pos[:e.np,:]=numpy.where(e.state[:e.np,numpy.newaxis]>0,e.post[:e.np,:],e.pos[:e.np,:])
        e.pos[:e.np,:] = e.post[:e.np,:]
        #print '%s: %s %d particles' % (t if t<700000 else ncep2dt(t).strftime('%Y%m%d %H:%M:%S'),e.id,e.np)
        print '[%s] %s: %s %d particles' % (self.release_id, t if t<700000 else ncep2dt(t).strftime('%Y%m%d %H:%M:%S'),e.id,e.np)
        #e.age[:e.np]+=abs(dt)*(e.state[:e.np]>0)
        e.age[:e.np]+=abs(dt)*(e.state[:e.np]>=outofcompute)
        last_time = self.timestamp('ageing', last_time)
        e.die(t,t2)
        last_time = self.timestamp('die', last_time)
        e.tcum+=86400.*abs(dt)
        # Output to file
        if i%iprint==0:
          self.fout[e.id].write(e.sfprint(t2))
          last_time = self.timestamp('write output', last_time)
          # remove particles flagged with state<=outofcompute (stuck to seabed,shoreline etc...)  outofcompute=-1
          #ind=(e.state<0)
          ind=(e.state<=outofcompute)  
          nind=ind.sum()
          if nind:e._reset(ind) #Recycle particles
          if nind>0:print '%d %s particles removed' % (nind,e.id)
      # check if the initial materials have "spawned" any children
      for e in self.materials:
        if dt>0:
          for spw in e.children:
            # import pdb;pdb.set_trace()
            spw0=spw.split('_')[0] # not sure why this is used ? maybe a convention to use same name with _plume, _sediment etc.. ?
            # if spw0 not in etypes:continue
            # self.materials[etypes.index(spw0)].release(t,t2,**e.children[spw])
            if spw not in etypes:continue
            self.materials[etypes.index(spw)].release(t,t2,**e.children[spw]) # **e.children[spw] allows passing the info from the parent material 
      t+=dt
      last_step = self.timestamp('tstep', start_step)
    for e in self.materials:
      self.fout[e.id].close()
    last_time = self.timestamp('total_time', start_time)
    end_time = self.timestamp('end_time')
    if self.save_summary:
      report = get_summary_report(self.summary, dt*86400, self.tout, self.outpath, self.materials)
      with open(os.path.join(self.outpath,'summary.txt'),'w') as sfile:
        sfile.write(report)
      print report 

  
  


