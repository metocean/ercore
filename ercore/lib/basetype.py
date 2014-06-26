#!/usr/bin/env python
import datetime,math
d2r=math.pi/180.
r2d=1./d2r


def add_month(d,x):
    """Add month to datetime and retain numerical day
    d: datetime
    x: number of months to add
    """
    newmonth = ((( d.month - 1) + x ) % 12 ) + 1
    newyear  = d.year + ((( d.month - 1) + x ) // 12 )
    newday=d.day
    while True:
        try:
            newdate=datetime.date( newyear, newmonth, newday,d.hour,d.minute,d.second)
            return newdate
        except ValueError:
            newday-=1

def add_year(d,x):
    """Add month to datetime and retain numerical month and day
    d: datetime
    x: number of years to add
    """
    return datetime.date(d.year+1,d.month,d.day,d.hour,d.minute,d.second)

def parsetime(t):
    """Convert string representation of time to datetime.
    Formats allowed:
    %Y%m%d.%H%M%S
    %Y%m%d.%H%M%Sz
    %Y%m%d_%Hz
    """
    if isinstance(t,str):
        if t[-1]=='z':
            if len(t)==16:
                t=datetime.datetime.strptime(t[0:15],'%Y%m%d.%H%M%S')
            else:
                t=datetime.datetime.strptime(t,'%Y%m%d_%Hz')
        else:
            t=datetime.datetime.strptime(t,'%Y%m%d.%H%M%S')
    return t

class Container(object):
    """Container class for non-specified entities in configuration files"""
    def __init__(self,**k):
        self.__dict__.update(k)
        if not k.has_key('id'):self.id='_undef_'

    def str(self,format='XML',name='Container'):
        if format=='XML':
            attrs=["%s=\"%s\"" % (k,self.__dict__[k]) for k in self.__dict__ if str(self.__dict__[k])[0]!='_']
            return "<%s %s />" % (name,' '.join(attrs))
        else:
            return name+':'+str(self.__dict__)
            
class MonthTime(datetime.datetime):
    """Class containing time intened for use with monthly timesteps"""
    def __add__(self,t):
        """Add t months to time"""
        newmonth = ((( self.month - 1) + t ) % 12 ) + 1
        newyear  = self.year + ((( self.month - 1) + t ) // 12 )
        newday=self.day
        while True:
            try:
                newdate=MonthTime( newyear, newmonth, newday)
                return newdate
            except ValueError:
                newday-=1
                
    def __gt__(self,other):
        return self.year>other.year or (self.year==other.year and self.month>other.month)
                
    def __cmp__(self,other):
        pass
        
    def __rcmp__(self,other):
        pass
        
    def __sub__(self,t):
        if isinstance(t,int):
            return self.__add__(-t)
        else:
            return 12*self.year+self.month-12*t.year-t.month
        
    def __rsub__(self,t):
        return -self.__sub__(t)

#Base time class
class Times(object):
    """Class representing a series of times with constant timestep""" 
    t2d={'DAY':datetime.timedelta(1,),'HR':datetime.timedelta(0,3600,),'MIN':datetime.timedelta(0,60,),'SEC':datetime.timedelta(0,1,),'WEEK':datetime.timedelta(7,0,0),'MONTH':datetime.timedelta(0,0,0)}
    def __init__(self,**k):
        """Argument options - keyword arguments:
        1.
        times: List of datetime objects
        2.
        comptimes: List of cdms component time objects
        3.
        t0: datetime of start
        t1: datetime of end
        dt: integer timestep
        dtunit: timestep unit -'DAY' | 'HR'(default) | 'MIN' | 'SEC' | 'WEEK' | 'MONTH'
        """
        self.dtunit='HR'
        if k.has_key('times'):
            self.t0=k['times'][0]
            self.t1=k['times'][-1]
            if len(k['times'])>1:
                self.dt=k['times'][1]-k['times'][0]
            else:
                self.dt=datetime.timedelta(0.)
        elif k.has_key('comptimes'):
            tax=k['comptimes']
            #timedelta hack below is ensure proper rounding to nearest second.
            self.t0=datetime.datetime(tax[0].year,tax[0].month,tax[0].day,tax[0].hour,tax[0].minute)+datetime.timedelta(0,10*round(0.1*tax[0].second))
            self.t1=datetime.datetime(tax[-1].year,tax[-1].month,tax[-1].day,tax[-1].hour,tax[-1].minute)+datetime.timedelta(0,10*round(0.1*tax[-1].second))
            self.dt=datetime.datetime(tax[1].year,tax[1].month,tax[1].day,tax[1].hour,tax[1].minute)+datetime.timedelta(0,round(tax[1].second))-self.t0 if len(tax)>1 else datetime.timedelta(0.0)
        else:
            self.t0=parsetime(k['t0'])
            self.t1=parsetime(k['t1'])
            if self.t1<self.t0:
                ttmp=self.t1
                self.t1=self.t0
                self.t0=ttmp
            if not k.has_key('dt'):
                self.dt=self.t1-self.t0
            else:
                self.dtunit=k.get('dtunit','HR').upper()
                if self.dtunit=='MONTH':
                    self.t0=MonthTime(self.t0.year,self.t0.month,self.t0.day)
                    self.t1=self.t0+(MonthTime(self.t1.year,self.t1.month,self.t1.day)-self.t0)
                    self.dt=int(k['dt'])
                elif (isinstance(k['dt'],int) or isinstance(k['dt'],float)):
                    self.dt=int(k['dt'])*self.t2d[self.dtunit]
                else:
                    self.dt=k['dt']
        self.gaps=[]
        
    def __add__(self,op):
        if isinstance(op,float) or isinstance(op,int):
            dt=datetime.timedelta(op)
        elif not isistance(op,datetime,timedelta):
            dt=datetime.timedelta(0.)
        return Times(t0=self.t0+dt,t1=self.t1+dt,dt=self.dt,dtunit=self.dtunit)
                    
    def __eq__(self,other):
        """Test for equality of 2 Times objects"""
        if self.t0<>other.t0:return False
        if self.t1<>other.t1:return False
        if self.t0==self.t1:return True
        if self.dt<>other.dt:return False
        return True
    
    def __cmp__(self,other):
        """Comparison test of 2 Times objects based on start time"""
        if self.t0>other.t0:
            return 1
        elif self.t0<other.t0:
            return -1
        return 0
    
    def __iter__(self):
        """Iterate over individual times. """
        t=self.t0
        igap=0
        dt=max(self.dt,datetime.timedelta(0,0,1))
        g=self.gaps[0] if len(self.gaps) else None
        while t<=self.t1:
            if not g or not g.intime(t):
                yield t
            if g and t>g.t1:
                igap+=1
                g=self.gaps[igap] if igap<len(self.gaps) else None
            if self.dtunit=='MONTH':
                t=add_month(t,self.dt)
            else:
                t+=dt
    
    def __len__(self):
        """Return number of indivdual times"""
        if self.t0 is None: return 0
        return int(self.ndt(self.t1))+1
    
    def __nonzero__(self):
        return self.t0 is not None
    
    def __getitem__(self,ind):
        """Return datetime or slice from indices"""
        if isinstance(ind,slice):
            return Times(t0=self[ind.start],t1=self[ind.stop-1],dt=self.dt)
        else:
            tz=self.t0+ind*self.dt
            if tz>self.t1:
                raise IndexError('Index greater than maximum time')
            else:
                return tz
    
    def ndt(self,t):
        """Return number of timesteps between specified time and start of Times object"""
        if self.t0 is None:return 0
        if self.dtunit=='MONTH':
            return (t-self.t0)/self.dt
        elif self.dt>datetime.timedelta(0.):
            tspan=t-self.t0
            return (86400.*tspan.days+tspan.seconds)/self.dtsecs()
        elif t==self.t0:
            return 0.
        else:
            return None
        
    def index(self,t):
        """Return nearest index for specified datetime"""
        ndt=self.ndt(t)
        return int(ndt) if ndt % 1.0==0. else None
    
        
    def sub(self,tsub):
        """Return overlapping subset of Times based on a second Times object"""
        if self.t0 is None:return False
        if tsub.t0 is None:return True
        if tsub.t0<=self.t0:
            if (tsub.t1>=self.t0):
                if (tsub.t1<self.t1):#Spans bottom half
                    self.t0+=(int(self.ndt(tsub.t1))+1)*self.dt
                else:#Spans whole - return Null time 
                    self.t0=None
                    self.t1=None
                    return True
            else:#Spans none below
                return False
        elif tsub.t1>=self.t1:
            if (tsub.t0<=self.t1):#Spans top half
                self.t1=self.t0+int(self.ndt(tsub.t0)-0.00000001)*self.dt
            else:#Spans none above
                return False
        else:#Is gap
            mergegap=0
            for g in self.gaps:#Try to merge existing gap
                if g.extend(tsub):mergegap+=1
            if mergegap==0:#Add new gap
                self.gaps.append(self.intersect(tsub)[0])
            self.gaps.sort()
            if mergegap>1:#Overlapping gaps
                i=1
                while i<len(self.gaps):
                    if self.gaps[i].t0<=self.gaps[i-1].t1:
                        self.gaps[i-1].t1=max(self.gaps[i].t1,self.gaps[i-1].t1)
                        del self.gaps[i]
        if len(self.gaps)>0:#Test no gaps overlap ends
            if self.t0>=self.gaps[0].t0:
                self.t0=self.gaps[0].t1+self.dt
                self.gaps=self.gaps[1:]
            elif self.t1<=self.gaps[-1].t1:
                self.t1=self.gaps[-1].t0-self.dt
                self.gaps=self.gaps[:-1]
        if self.t1<self.t0:
            self.t0=None
            self.t1=None
        return True


    def extend(self,tadd):
        """Extend based on a second Times object - return whether extended"""
        if (tadd.t0>self.t1+self.dt) or (tadd.t1<self.t0-self.dt):return False
        if (tadd.t0<self.t0):self.t0+=self.dt*int(self.ndt(tadd.t0))
        if (tadd.t1>self.t1):self.t1=self.t0+self.dt*int(self.ndt(tadd.t1))
        return True
        
    def dtstr(self,rounded=''):
        """String representation of timestep"""
        if not self.dt:
            return 'None'
        elif self.dtunit=='MONTH':
            return '%d MONTHS' % (self.dt)
        seconds=1*self.dt.seconds+86400*self.dt.days
        if seconds==0:
            return '0 SECS'
        dtstr=''
        if rounded[0:2]=='HR':
            seconds=int(3600.*round(seconds/3600.))
        elif rounded[0:3]=='MIN':
            seconds=int(60.*round(seconds/60.))
        if seconds>86400 and not rounded:
            dday=seconds/86400
            dtstr+='%d DAYS ' % (dday)
            seconds-=86400*dday
        if seconds>=3600 and not rounded[0:3]=='MIN':
            dthr=seconds/3600.
            dtstr+='%d HRS ' % (dthr)
            seconds-=dthr*3600
        if seconds>=60:
            dtmin=seconds/60.
            dtstr+='%d MINS ' % (dtmin)
            seconds-=60*dtmin
        if seconds>0 and not rounded:
            dtstr+='%d SECS ' % (seconds)
        return dtstr
    
    def dtsecs(self):
        """Return timestep as seconds"""
        return 86400*self.dt.days+self.dt.seconds
    
    def _toXMLobj(self):
        if self.dtunit=="MONTH":
            return "t0=\"%s\" t1=\"%s\" dt=\"%d\" dtunits=\"MONTH\"" % (self.t0.strftime('%Y%m%d.%H%M%S'),self.t1.strftime('%Y%m%d.%H%M%S'),self.dt)
        else:
            return "t0=\"%s\" t1=\"%s\" dt=\"%f\" dtunits=\"SEC\"" % (self.t0.strftime('%Y%m%d.%H%M%S'),self.t1.strftime('%Y%m%d.%H%M%S'),self.dtsecs())
                
    def str(self,format='txt'):
        """Return string representation
        format: 'CYCLE'   t0 as %Y%m%d_%Hz 
                'DOT'     t0 as %Y%m%d.%H%M%S
                'CDTIME'  t0 as %Y-%m-%d %H:%M:%S
                'XML'     Complete object as XML
                'TXT'     Complete object as readable string
        """
        if not self.t0:
            return 'Null time'
        elif format=='CYCLE':
            return self.t0.strftime('%Y%m%d_%Hz')
        elif format=='DOT':
            return self.t0.strftime('')
        elif format=='CDTIME':
            return self.t0.strftime('%Y-%m-%d %H:%M:%S')
        elif format=='XML':
            return "<Times %s />" % (self._toXMLobj())
        else:
            return 'tstart='+self.t0.strftime('%Y%m%d %H:%M:%S')+'  tend='+self.t1.strftime('%Y%m%d %H:%M:%S')+'  timestep='+self.dtstr()

        
    def __repr__(self):
        return 'MO Times: '+self.str('object')

    def get(self,fillgaps=False):
        """Return list of individual datetime times"""
        t=[self.t0]
        if self.t1>self.t0:
            if not fillgaps:
                for g in self.gaps:
                    while t[-1]<(g.t0-self.dt):
                        t.append(t[-1]+self.dt)
                    t.append(g.t1+self.dt)
            t2=t[-1]+self.dt
            while t2<=self.t1:
                t.append(t2)
                t2=t[-1]+self.dt
        return t

    def decdays(self,t0=(datetime.datetime(2000,1,1))):
        """Return list individual times as decimal days since t0""" 
        return [(t-t0).days+(t-t0).seconds/86400. for t in self.get()]
    
    def ncepdays(self):
        """Return list of individual times as NCEP decimal days (udunits) since 1-1-1"""
        return [1.+t.toordinal()+t.hour/24.+t.minute/1440.+t.second/86400. for t in self.get()]
    
    def intime(self,t1):
        """Return logical test for single timestamp within time"""
        if (t1>=self.t0) & (t1<=self.t1):
            return True
        else:
            return False
        
    def overlaps(self,tother):
        """Return logical test for other Times object overlapping"""
        if (tother.t0>=self.t0) and (tother.t0<=self.t1):
            return True
        elif (tother.t1>=self.t0) and (tother.t1<=self.t1):
            return True
        elif (tother.t0<=self.t0) and (tother.t1>=self.t1):
            return True
        elif (tother.t0>=self.t0) and (tother.t1<=self.t1):
            return True
        else:
            return False
        
    def intersect(self,tint,dt2=False):
        """Return intersection of this and other Times objects""" 
        if isinstance(tint,list):
            tstart=tint[0]
            tend=tint[1]
            dt=self.dt
        else:
            if not tint or tint==NullTime:return NullTime,[]
            tstart=tint.t0
            tend=tint.t1
            dt=tint.dt if (dt2 and tint.dt and type(tint.dt)==type(self.dt)) else self.dt
        if tstart is None or self.t0 is None or (self.t1<tstart) or (self.t0>tend):
            return NullTime,[]
        else:
            ind=[0,len(self)-1]
            if tstart<=self.t0:
                ts0=self.t0
            else:
                ndt=int(self.ndt(tstart))
                ind[0]=int(self.ndt(tstart))
                ts0=self.t0+ind[0]*self.dt
            if tend>=self.t1:
                ts1=self.t1
            else:
                ndt=self.ndt(tend)
                ind[1]=int(ndt+0.99999999)
                ts1=self.t0+ind[1]*self.dt
            tnew=Times(t0=ts0,t1=ts1,dt=dt,dtunit=self.dtunit)
            for gap in self.gaps:
                tnew.sub(gap)
            return tnew,ind
        
    def interpfac(self,tint,filled=False):
        """Return interpolation factors for interpolation of other Times object onto self"""
        t1=self.get(filled)
        ts=tint.get(filled)
        ind=[]
        fac=[]
        j1=0
        for it,t in enumerate(t1):
            if t<=ts[0]:
                ind.append(0)
                fac.append(1.0)
            elif t>=ts[-1]:
                ind.append(len(ts)-2)
                fac.append(0.0)
            else:
                while j1 < len(ts)-1:
                    if (t >= ts[j1] and t < ts[j1+1]):
                        ind.append(j1)
                        dt1=ts[j1+1]-t
                        dt2=ts[j1+1]-ts[j1]
                        fac.append((86400.*dt1.days+dt1.seconds)/(86400.*dt2.days+dt2.seconds))
                        break
                    else:
                        j1+=1
        return ind,fac
    
NullTime=Times(t0=None,t1=None,dt=0)
        
#Enhanced list class
class ObjectList(list):
    """Enhance list class for objects to allow access by object id"""
    def __init__(self,items,copy=False):
        """
        items: Object or list of objects - must have an 'id' attribute
        copy: Make copies of member items
        """
        if not isinstance(items,list):items=[items]
        for i in items:
            if not hasattr(i,'id'):raise ValueError('Member of an ObjectList must have an id attribute')
        if copy:
            import copy
            items=copy.deepcopy(items)
        list.__init__(self,items)
        
    def __getitem__(self,key):
        """Over-ridden list access by id attribute"""
        if key is None:return None
        if isinstance(key,str):
            idlist=self.idList()
            return self[idlist.index(key)] if key in idlist else None
        else:
            return list.__getitem__(self,key)
  
    def __str__(self):
        return "\n".join([str(n) for n in self])
            
    def subset(self,keys):
        """Return a sublist of member items from a specified key or list of keys"""
        if isinstance(keys,str):keys=[keys]
        idlist=self.idList()
        return ObjectList([self[idlist.index(key)] for key in keys if key in idlist])
            
    def pop(self,key=None):
        """Enhanced pop method to pop item based on key(optional)"""
        if isinstance(key,int):
            return list.pop(self,key)
        elif isinstance(key,str):
            ind=self.idList().index(key)
            return self.pop(ind)
        else:
            return None
            
    def sort(self,attr='id',reverse=False):
        """Sort list based on specified attribute"""
        sortkey=lambda o:getattr(o,attr)
        list.sort(self,key=sortkey,reverse=reverse)
    
    def idList(self):
        """Return list of member id attributes"""
        return [ds.id for ds in self]
        

#Base site class
class Site(object):
    def __init__(self,**k):
        self.x=k['x']
        self.y=k['y']
        self.id=str(k.get('id',str(self.x)+'_'+str(self.y)))
        self.name=k.get('name','')
        self.active=k.get('active',1)
        self.lev=k.get('lev',-1)
        self.prod=k.get('prod','')
        
    def str(self,format='txt'):
        if format=='XML':
            return '<Site x="%.8g" y="%.8g" id="%s" name="%s" />' % (self.x,self.y,self.id,self.name)
        elif format=='long':
            return 'name='+str(self.name)+' id='+str(self.id)+' Location: '+self.str()
        else:
            return 'x='+str(self.x)+' y='+str(self.y)
    
    def __repr__(self):
        return 'MO Site: '+self.str()    

#Base class representing an area
class Area(Site):
    def __init__(self,**k):
        if (not k.has_key('x')) & (k.has_key('bnd')):k['x']=k['bnd'][0]
        if (not k.has_key('y')) & (k.has_key('bnd')):k['y']=k['bnd'][2]
        Site.__init__(self,**k)
        if (not k.has_key('x2')) & (k.has_key('bnd')):k['x2']=k['bnd'][1]
        if (not k.has_key('y2')) & (k.has_key('bnd')):k['y2']=k['bnd'][3]
        self.ang=d2r*k['ang'] if k.has_key('ang') else 0.
        if (not k.has_key('x2') and k.has_key('lenx')):k['x2']=k['x']+k['lenx']*math.cos(self.ang)-k['leny']*math.sin(self.ang)
        if (not k.has_key('y2') and k.has_key('leny')):k['y2']=k['y']+k['lenx']*math.sin(self.ang)+k['leny']*math.cos(self.ang)
        self.x2=k['x2']
        self.y2=k['y2']
        self.lenx=abs((k['x2']-k['x'])*math.cos(self.ang)+(k['y2']-k['y'])*math.sin(self.ang))
        self.leny=abs((k['y2']-k['y'])*math.cos(self.ang)-(k['x2']-k['x'])*math.sin(self.ang))
        
    def str(self,format='txt'):
        if format=='XML':
            return '<Area x="%g" x2="%g" y="%g" y2="%g" angle="%g"/>' % (self.x,self.x2,self.y,self.y2,self.ang)
        elif format=='long':
            return 'name='+str(self.name)+' id='+str(self.id)+' Bounds: '+self.str('short')
        else:
            return 'BL=[%g %g] UR=[%g %g] angle=%g' % (self.x,self.y,self.x2,self.y2,r2d*self.ang)

    def __str__(self):
        return 'MO Area: '+self.str('long')
        
    def mod360(self,copy=False):
        if copy:
            import copy
            a=copy.copy(self)
        else:
            a=self
        if (a.x<0) and (a.x2<=0):
            a.x%=360.
            a.x2%=360.
        return a
    
    def enlarge(self,factor,copy=False):
        if copy:
            import copy
            area=copy.copy(self)
        else:
            area=self
        dx=0.5*(factor-1.)*area.lenx
        dy=0.5*(factor-1.)*area.leny
        area.x-=dx
        area.x2+=dx
        area.y-=dy
        area.y2+=dy
        return area

    
    def inside(self,pl):
        #Is a point p1 contained in area
        if isinstance(pl,tuple):
            ing=(pl[0]>=self.x) & (pl[0]<=self.x2) & (pl[1]>=self.y) & (pl[1]<=self.y2)
        elif isinstance(pl,list):
            ing=[self.inside(p) for p in p1]
        else:
            ing=(pl.x>=self.x) & (pl.x<=self.x2) & (pl.y>=self.y) & (pl.y<=self.y2)
        return ing
    
    def contains(self,a1,mod360=True):
        #Is another area a1 contained within area
        if isinstance(a1,list):
            ing=[self.contains(a,mod360) for a in a1]
        else:
            ing=(a1.x>=self.x) & (a1.x2<=self.x2) & (a1.y>=self.y) & (a1.y2<=self.y2)
        if mod360:
            ing=ing or self.contains(a1.mod360(True),False)
        return ing
    
    def overlap(self,a1,mod360=True):
        #Is another area a1 overlapping with area
        if isinstance(a1,list):
            return [self.overlap(a,mod360) for a in a1]
        else:
            sepx=abs(0.5*(a1.x+a1.x2)-0.5*(self.x+self.x2))
            sepy=abs(0.5*(a1.y+a1.y2)-0.5*(self.y+self.y2))
            ovlp=(sepx<0.5*(a1.lenx+self.lenx)) & (sepy<0.5*(a1.leny+self.leny))
        if mod360:
            ovlp=ovlp or self.overlap(a1.mod360(True),False)
        return ovlp
            
    def bnd(self):
        return [self.x,self.x2,self.y,self.y2]
        
    def merc_aspect(self):
        y1=math.log(math.tan(math.pi*self.y2/180.)+1./math.cos(math.pi*self.y2/180.))
        y0=math.log(math.tan(math.pi*self.y/180.)+1./math.cos(math.pi*self.y/180.))
        return 180*(y1-y0)/(self.x2-self.x)/math.pi
    
        
#Base class representing a computational or storage grid
class Grid(Area):
    def __init__(self,**k):
        if k.has_key('area'):
            self.__dict__.update(k['area'].__dict__)
        else:
            Area.__init__(self,**k)
        if (not k.has_key('dx')) & (k.has_key('res')):k['dx']=k['res'][0] if isinstance(k['res'],list) else k['res']
        if (not k.has_key('dy')) & (k.has_key('res')):k['dy']=k['res'][1] if isinstance(k['res'],list) else k['res']
        if (not k.has_key('nx')) & (k.has_key('dx')): k['nx']=int(self.lenx/k['dx'])+1
        if (not k.has_key('ny')) & (k.has_key('dy')): k['ny']=int(self.leny/k['dy'])+1
        if (not k.has_key('type')):k['type']='LL'
        self.type=k['type']
        self.nx=int(k['nx'])
        self.ny=int(k['ny'])
        self.dx=self.lenx/(self.nx-1) if self.nx>1 else 0. 
        self.dy=self.leny/(self.ny-1) if self.ny>1 else 0.
        self.buf=min(k['buf'],self.nx,self.ny) if k.has_key('buf') else 0
        
    def bufArea(self):
        if self.buf:
            return Area(x=self.x+self.buf*self.dx,x2=self.x2-self.buf*self.dx,y=self.y+self.buf*self.dy,y2=self.y2-self.buf*self.dy)
        else:
            return self
        
    def wrapped360(self):
        return (self.x+360.-self.x2-0.00001<=self.dx)
            
    def str(self,format='txt'):
        if format=='BOUNDS1':
            return str(self.x)+' '+str(self.x2)+' '+str(self.y)+' '+str(self.y2)
        elif format=='XML':
            return '<Grid x="%g" x2="%g" y="%g" y2="%g" dx="%g" dy="%g" nx="%i" ny="%i" ang="%g" buf="%i"/>' % (self.x,self.x2,self.y,self.y2,self.dx,self.dy,self.nx,self.ny,r2d*self.ang,self.buf)
        elif format=='long':
            return 'name='+str(self.name)+' id='+str(self.id)+' Bounds: '+self.str()
        elif format=='KML':
            desc=" Lon:%g to %g Lat:%g to %g dx:%g dy:%g nx:%d ny:%d" % (self.x,self.x2,self.y,self.y2,self.dx,self.dy,self.nx,self.ny)
            return "<Placemark><styleUrl>#%s</styleUrl><name>%s</name><description>%s</description><Polygon><outerBoundaryIs><LinearRing><coordinates>%f,%f,0 %f,%f,0 %f,%f,0 %f,%f,0</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>" \
                % (self.id,self.name,desc,self.x,self.y,self.x,self.y2,self.x2,self.y2,self.x2,self.y)
        else:
            return Area.str(self)+' %i by %i  dx=%g dy=%g' % (self.nx,self.ny,self.dx,self.dy)

        def __str__(self):
            return 'MO Grid: '+self.str('long')

    def __eq__(self,ogrid):
        tol=0.0001*(self.x2-self.x)
        if abs(self.x-ogrid.x)>tol:return False
        if abs(self.x2-ogrid.x2)>tol:return False
        if abs(self.y-ogrid.y)>tol:return False
        if abs(self.y2-ogrid.y2)>tol:return False
        if self.nx<>ogrid.nx:return False
        if self.ny<>ogrid.ny:return False
        return True
    
    def __ne__(self,ogrid):
        return not self.__eq__(ogrid)
        
   
    def ingrid(self,p1,mod360=False):
        if self.inside(p1):
            return True
        if mod360:
            import copy
            pnew=copy.copy(p1)
            pnew.x-=360.
            if self.inside(pnew):
                p1.x=pnew.x
                return True
            pnew.x+=720.
            if self.inside(pnew):
                p1.x=pnew.x
                return True
        return False
    
    def contains(self,a1,mod360=False):
        return Area.contains(self,a1,mod360) or (mod360 and self.wrapped360() and a1.y>=self.y and a1.y2<=self.y2)
        
    def overlap(self,a1,mod360=False):
        return Area.overlap(self,a1,mod360) or (mod360 and self.wrapped360() and abs(0.5*(a1.y+a1.y2)-0.5*(self.y+self.y2))<0.5*(a1.leny+self.leny))
    
    def extend(self,ng):
        import math
        extended=False
        if ng.x2<>self.x2:
            extra=math.floor((ng.x2-self.x2)/self.dx)
            self.nx+=int(extra)
            self.x2+=self.dx*extra
            extended=True
        if ng.x<>self.x:
            extra=math.floor((self.x-ng.x)/self.dx)
            self.nx+=int(extra)
            self.x-=self.dx*extra
            extended=True
        ymax=max(self.y,self.y2)
        if ng.y2<>ymax:
            extra=math.floor((ng.y2-ymax)/self.dy)
            self.ny+=int(extra)
            if self.y2==ymax:
                self.y2+=self.dy*extra
            else:
                self.y+=self.dy*extra
            extended=True
        ymin=min(self.y,self.y2)
        if ng.y<>ymin:
            extra=math.floor((ymin-ng.y)/self.dy)
            self.ny+=int(extra)
            if self.y==ymin:
                self.y-=self.dy*extra
            else:
                self.y2-=self.dy*extra
            extended=True
        return extended
        
    
    def boundary(self):
        import numpy
        cx=range(1,self.nx)+(self.ny-1)*[self.nx-1]+range(self.nx-2,-1,-1)+(self.ny-1)*[0]
        cy=(self.nx-1)*[0]+range(1,self.ny)+ (self.nx-2)*[self.ny-1] + range(self.ny-1,-1,-1)
        bpx=self.x+self.dx*numpy.array(cx)*math.cos(self.ang)-self.dy*numpy.array(cy)*math.sin(self.ang)
        bpy=self.y+self.dx*numpy.array(cx)*math.sin(self.ang)+self.dy*numpy.array(cy)*math.cos(self.ang)
        p=[Site(**{'x':bpx[ip],'y':bpy[ip],'tag':'bnd'+str(ip)}) for ip,p1 in enumerate(bpx)]
        coord=[(icx,cy[i]) for i,icx in enumerate(cx)]
        return p,coord
    
    def axes(self,ax):
        import numpy
        if ax=='x':
            return self.x+self.dx*numpy.arange(0,self.nx)
        elif ax=='y':
            return self.y+self.dy*numpy.arange(0,self.ny)
            
    def res(self):
        return min(self.dx,self.dy)
        
class LineString(list):
    def __init__(self,xy):
        if isinstance(xy[0],float):
            for ix,x in enumerate(xy[0::2]):
                self.append(Site(x=x,y=xy[2*ix+1]))
        elif isinstance(xy[0],Site):
            for s in xy:
                self.append(s)
                
    def list_xy(self):
        out=[]
        for s in self:
            out.extend([s.x,s.y])
        return out


if __name__=="__main__":
    mtime=Times(t0="20000101_00z",t1="20020101_00z",dtunit="MONTH",dt=6)
    mtime2=Times(t0="20000101_00z",t1="20100101_00z",dtunit="HR",dt=5)
    print mtime2+3./24
    print mtime.get()
    
    
