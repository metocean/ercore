#!/usr/bin/env python
import datetime,os,numpy

_DT0_=datetime.datetime(2000,1,1)
_NCEPT0_=730120.99999
ncep2dt=lambda t:_DT0_+datetime.timedelta(t-_NCEPT0_)
dt2ncep=lambda t: (1.+t.toordinal()+t.hour/24.+t.minute/1440.+t.second/86400.) if isinstance(t,datetime.datetime) else t

def parsetime(t):
    if isinstance(t,datetime.datetime):
        return dt2ncep(t)
    elif isinstance(t,str):
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
        if isinstance(key,str):
            idlist=self.idList()
            return self[idlist.index(key)] if key in idlist else None
        elif isinstance(key,int):
            return list.__getitem__(self,key) if key<len(self) else None
        else:
            return None
            
    def __str__(self):
        return "\n".join([str(n) for n in self])
            
    def subset(self,keys):
        if isinstance(keys,str):keys=[keys]
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

#Main class for ERcore simulation
class ERcore(object):
  geod=True
  tout=0.
  rkorder=4
  release_id = 0
  def __init__(self,**k):
    self.fout={}
    self.outpath=k.get('outpath','.')
    self.__dict__.update(k)
    
  def readYAML(self,ctlfile, namespace=globals()):
    """Read a YAML configuration file for an ERcore simulation"""
    import yaml,inspect
    try:
      config = yaml.load(open(ctlfile).read())
    except:
      raise ERConfigException('Cannot read '+ctlfile)
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
    for obj in objects:
        for prop in dir(obj):
            if prop in ['movers','reactors','diffusers','stickers','topo','members']:
                val=getattr(obj,prop)
                if isinstance(val,list):
                    val=[objects[v] if isinstance(v,str) and objects[v] else v for v in val]
                    for v in val:
                        if isinstance(v,str):raise ERConfigException('Cannot find one of %s with id(s) %s specifed for %s' % (prop,v,obj.id))
                elif isinstance(val,str) and objects[val]:
                    val=objects[val]
                else:
                    raise ERConfigException('Cannot find one of %s with id(s) %s specifed for %s' % (prop,val,obj.id))
                setattr(objects[obj.id],prop,val) 
    import materials
    import inspect
    self.materials=[o for o in objects if materials._Material in inspect.getmro(o.__class__)]

  def run(self,t,tend,dt):
    """Run Ercore
    Arguments:
        t: Start time as datetime, ncep decimal time or string
        tend: Start time as datetime, ncep decimal time or string
        dt: Simulation timestep as seconds
    """
    t=parsetime(t)
    tend=parsetime(tend)
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
    while ((dt>0) and (t<tend)) or ((dt<0) and (t>tend)):
      i+=1
      t2=t+dt
      for e in self.materials:
        if e.tcum>=e.tstep:
          e.tcum-=e.tstep
          e.release(t,t2)
          if self.geod:e.geodcalc()
          e.react(t,t2)
          if (e.state[:e.np]>0).any():
            e.advect(t,t2,self.rkorder)
            e.diffuse(t,t2)
            e.stick(t,t2)
          e.spawn(t,t2)
        e.pos[:e.np,:]=numpy.where(e.state[:e.np,numpy.newaxis]>0,e.post[:e.np,:],e.pos[:e.np,:])
        print '[%s] %s: %s %d particles' % (self.release_id, t if t<700000 else ncep2dt(t).strftime('%Y%m%d %H:%M:%S'),e.id,e.np)
        e.age[:e.np]+=abs(dt)*(e.state[:e.np]>0)
        e.die(t,t2)
        e.tcum+=86400.*abs(dt)
        if i%iprint==0:
          self.fout[e.id].write(e.sfprint(t2))
          ind=(e.state<0) 
          nind=ind.sum()
          if nind:e._reset(ind) #Recycle particles
          if nind>0:print '%d %s particles removed' % (nind,e.id)
      for e in self.materials:
        if dt>0:
          for spw in e.children:
            spw0=spw.split('_')[0]
            if spw0 not in etypes:continue
            self.materials[etypes.index(spw0)].release(t,t2,**e.children[spw])
      t+=dt
    for e in self.materials:
      self.fout[e.id].close()

  
  


