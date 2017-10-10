#!/usr/bin/env python
#Classes and functions to calculate tidal elevations and flows from constituents

import datetime,numpy
from pymo.core.basetype import *

DEFAULT_CONS=['Z0','M2','S2','N2','K2','K1','O1','P1','Q1']
##omg=2*math.pi*array([0,1.93227361,2.0,1.89598197,2.00547582,1.00273791,0.929535715,0.99726209,0.89324406])
#V0=array([0.,1.731557546,0.,6.050721243,3.487600001,0.173003674,1.558553872,6.110181633,5.877717569])
#tidet0=datetime.datetime(1992,1,1)
R2D=180./numpy.pi
D2R=1./R2D
PI2=2*numpy.pi
INV2PI=1./(2.*math.pi)
I2PI=1j*2.*math.pi

ASTRO_C=numpy.matrix([[ 270.434164,13.1763965268,-0.0000850, 0.000000039],
[ 279.696678, 0.9856473354, 0.00002267,0.000000000],
[ 334.329556, 0.1114040803,-0.0007739,-0.00000026],
[-259.183275, 0.0529539222,-0.0001557,-0.000000050],
[ 281.220844, 0.0000470684, 0.0000339, 0.000000070]])

S_CONS=['K1','K2','L2','M2','N2','O1','P1','S2']

#Array members: [freq(rad/day),doodson(6)+semi or shallow_factors]
CONS_STR={
'Z0':[0.0],
'K1':[6.300388096,1,1,0,0,0,0,-0.75,],
'K2':[12.60077619,2,2,0,0,0,0,0,],
'L2':[12.36886032,2,1,0,-1,0,0,-0.5,],
'M2':[12.1408332,2,0,0,0,0,0,0,],
'M4':[24.28166638,{'M2':2,}],
'MF':[0.4599430076,0,2,0,0,0,0,0,],
'MM':[0.2280271193,0,1,0,-1,0,0,0,],
'N2':[11.91280606,2,-1,0,1,0,0,0,],
'O1':[5.840445088,1,-1,0,0,0,0,-0.25,],
'P1':[6.265982514,1,1,-2,0,0,0,-0.25,],
'Q1':[5.612417969,1,-2,0,1,0,0,-0.25,],
'S2':[12.56637061,2,2,-2,0,0,0,0,],
}
CONS_STR={
'Z0':[ 0.000000000,0,0,0,0,0,0,0],
'SA':[ 0.017201969,0,0,1,0,0,-1,0],
'SSA':[ 0.034405582,0,0,2,0,0,0,0],
'MSM':[ 0.197510291,0,1,-2,1,0,0,0],
'MM':[ 0.228027119,0,1,0,-1,0,0,0],
'MSF':[ 0.425537426,0,2,-2,0,0,0,0],
'MF':[ 0.459943008,0,2,0,0,0,0,0],
'ALP1':[ 5.186880543,1,-4,2,1,0,0,-0.25],
'2Q1':[ 5.384390834,1,-3,0,2,0,0,-0.25],
'SIG1':[ 5.414907677,1,-3,2,0,0,0,-0.25],
'Q1':[ 5.612417969,1,-2,0,1,0,0,-0.25],
'RHO1':[ 5.642934796,1,-2,2,-1,0,0,-0.25],
'O1':[ 5.840445088,1,-1,0,0,0,0,-0.25],
'TAU1':[ 5.874850685,1,-1,2,0,0,0,-0.75],
'BET1':[ 6.037955394,1,0,-2,1,0,0,-0.75],
'NO1':[ 6.072360976,1,0,0,1,0,0,-0.75],
'CHI1':[ 6.102877804,1,0,2,-1,0,0,-0.75],
'PI1':[ 6.248780545,1,1,-3,0,0,1,-0.25],
'P1':[ 6.265982514,1,1,-2,0,0,0,-0.25],
'S1':[ 6.283186127,1,1,-1,0,0,1,-0.75],
'K1':[ 6.300388096,1,1,0,0,0,0,-0.75],
'PSI1':[ 6.317590065,1,1,1,0,0,-1,-0.75],
'PHI1':[ 6.334793677,1,1,2,0,0,0,-0.75],
'THE1':[ 6.497898387,1,2,-2,1,0,0,-0.75],
'J1':[ 6.528415230,1,2,0,-1,0,0,-0.75],
'2PO1':[ 6.691519940,{'P1':2,'O1':-1,}],
'SO1':[ 6.725925521,{'S2':1,'O1':-1,}],
'OO1':[ 6.760331103,1,3,0,0,0,0,-0.75],
'UPS1':[ 6.988358222,1,4,0,-1,0,0,-0.75],
'ST36':[11.061731227,{'M2':2,'N2':1,'S2':-2,}],
'2NS2':[11.259241519,{'N2':2,'S2':-1,}],
'ST37':[11.289758347,{'M2':3,'S2':-2,}],
'ST1':[11.293647101,{'N2':2,'K2':1,'S2':-2,}],
'OQ2':[11.456751810,2,-3,0,3,0,0,0],
'EPS2':[11.487268638,2,-3,2,1,0,0,0],
'ST2':[11.521674235,{'M2':1,'N2':1,'K2':1,'S2':-2,}],
'ST3':[11.646484609,{'M2':2,'S2':1,'K2':-2,}],
'O2':[11.680890191,{'O1':2,}],
'2N2':[11.684778945,2,-2,0,2,0,0,0],
'MU2':[11.715295773,2,-2,2,0,0,0,0],
'SNK2':[11.878400482,{'S2':1,'N2':1,'K2':-1,}],
'N2':[11.912806064,2,-1,0,1,0,0,0],
'NU2':[11.943322892,2,-1,2,-1,0,0,0],
'ST4':[11.981617228,{'K2':2,'N2':1,'S2':-2,}],
'OP2':[12.106427617,{'O1':1,'P1':1,}],
'GAM2':[12.110316356,2,0,-2,2,0,0,-0.5],
'H1':[12.123631230,2,0,-1,0,0,1,-0.5],
'M2':[12.140833199,2,0,0,0,0,0,0],
'H2':[12.158035168,2,0,1,0,0,-1,0],
'MKS2':[12.175238780,{'M2':1,'K2':1,'S2':-1,}],
'ST5':[12.209644362,{'M2':1,'K2':2,'S2':-2,}],
'ST6':[12.303937908,{'S2':2,'N2':1,'M2':-1,'K2':-1,}],
'LDA2':[12.338343490,2,1,-2,1,0,0,-0.5],
'L2':[12.368860318,2,1,0,-1,0,0,-0.5],
'2SK2':[12.531965028,{'S2':2,'K2':-1,}],
'T2':[12.549168640,2,2,-3,0,0,1,0],
'S2':[12.566370609,2,2,-2,0,0,0,0],
'R2':[12.583572578,2,2,-1,0,0,-1,-0.5],
'K2':[12.600776191,2,2,0,0,0,0,0],
'MSN2':[12.794397744,{'M2':1,'S2':1,'N2':-1,}],
'ETA2':[12.828803325,2,3,0,-1,0,0,0],
'ST7':[12.863208907,{'K2':2,'M2':1,'S2':-1,'N2':-1,}],
'2SM2':[12.991908035,{'S2':2,'M2':-1,}],
'ST38':[13.022424863,{'M2':2,'S2':1,'N2':-2,}],
'SKM2':[13.026313617,{'S2':1,'K2':1,'M2':-1,}],
'2SN2':[13.219935170,{'S2':2,'N2':-1,}],
'NO3':[17.753251167,{'N2':1,'O1':1,}],
'MO3':[17.981278286,{'M2':1,'O1':1,}],
'M3':[18.211249790,3,0,0,0,0,0,-0.5],
'NK3':[18.213194160,{'N2':1,'K1':1,}],
'SO3':[18.406815712,{'S2':1,'O1':1,}],
'MK3':[18.441221294,{'M2':1,'K1':1,}],
'SP3':[18.832353123,{'S2':1,'P1':1,}],
'SK3':[18.866758720,{'S2':1,'K1':1,}],
'ST8':[23.628101837,{'M2':2,'N2':1,'S2':-1,}],
'N4':[23.825612128,{'N2':2,}],
'3MS4':[23.856128971,{'M2':3,'S2':-1,}],
'ST39':[24.019233681,{'M2':1,'S2':1,'N2':1,'K2':-1,}],
'MN4':[24.053639263,{'M2':1,'N2':1,}],
'ST9':[24.088044844,{'M2':1,'N2':1,'K2':1,'S2':-1,}],
'ST40':[24.247260800,{'M2':2,'S2':1,'K2':-1,}],
'M4':[24.281666382,{'M2':2,}],
'ST10':[24.316071964,{'M2':2,'K2':1,'S2':-1,}],
'SN4':[24.479176673,{'S2':1,'N2':1,}],
'KN4':[24.513582270,{'K2':1,'N2':1,}],
'MS4':[24.707203808,{'M2':1,'S2':1,}],
'MK4':[24.741609390,{'M2':1,'K2':1,}],
'SL4':[24.935230927,{'S2':1,'L2':1,}],
'S4':[25.132741234,{'S2':2,}],
'SK4':[25.167146815,{'S2':1,'K2':1,}],
'MNO5':[29.894084351,{'M2':1,'N2':1,'O1':1,}],
'2MO5':[30.122111485,{'M2':2,'O1':1,}],
'3MP5':[30.156517067,{'M2':3,'P1':-1,}],
'MNK5':[30.354027358,{'M2':1,'N2':1,'K1':1,}],
'2MP5':[30.547648896,{'M2':2,'P1':1,}],
'2MK5':[30.582054478,{'M2':2,'K1':1,}],
'MSK5':[31.007591903,{'M2':1,'S2':1,'K1':1,}],
'3KM5':[31.041997485,{'K2':1,'K1':1,'M2':1,}],
'2SK5':[31.433129329,{'S2':2,'K1':1,}],
'ST11':[35.772823789,{'N2':3,'K2':1,'S2':-1,}],
'2NM6':[35.966445327,{'N2':2,'M2':1,}],
'ST12':[36.000850908,{'N2':2,'M2':1,'K2':1,'S2':-1,}],
'2MN6':[36.194472446,{'M2':2,'N2':1,}],
'ST13':[36.228878043,{'M2':2,'N2':1,'K2':1,'S2':-1,}],
'ST41':[36.388093999,{'M2':3,'S2':1,'K2':-1,}],
'M6':[36.422499581,{'M2':3,}],
'MSN6':[36.620009872,{'M2':1,'S2':1,'N2':1,}],
'MKN6':[36.654415454,{'M2':1,'K2':1,'N2':1,}],
'ST42':[36.813631425,{'M2':2,'S2':2,'K2':-1,}],
'2MS6':[36.848037006,{'M2':2,'S2':1,}],
'2MK6':[36.882442588,{'M2':2,'K2':1,}],
'NSK6':[37.079952880,{'N2':1,'S2':1,'K2':1,}],
'2SM6':[37.273574417,{'S2':2,'M2':1,}],
'MSK6':[37.307979999,{'M2':1,'S2':1,'K2':1,}],
'S6':[37.699111843,{'S2':3,}],
'ST14':[42.034917549,{'M2':2,'N2':1,'O1':1,}],
'ST15':[42.266833422,{'N2':2,'M2':1,'K1':1,}],
'M7':[42.492916172,{'M2':3.500000e+00,}],
'ST16':[42.688482094,{'M2':2,'S2':1,'O1':1,}],
'3MK7':[42.722887676,{'M2':3,'K1':1,}],
'ST17':[43.148425102,{'M2':1,'S2':1,'K2':1,'O1':1,}],
'ST18':[48.107278525,{'M2':2,'N2':2,}],
'3MN8':[48.335305645,{'M2':3,'N2':1,}],
'ST19':[48.369711226,{'M2':3,'N2':1,'K2':1,'S2':-1,}],
'M8':[48.563332779,{'M2':4,}],
'ST20':[48.760843071,{'M2':2,'S2':1,'N2':1,}],
'ST21':[48.795248652,{'M2':2,'N2':1,'K2':1,}],
'3MS8':[48.988870190,{'M2':3,'S2':1,}],
'3MK8':[49.023275772,{'M2':3,'K2':1,}],
'ST22':[49.220786078,{'M2':1,'S2':1,'N2':1,'K2':1,}],
'ST23':[49.414407616,{'M2':2,'S2':2,}],
'ST24':[49.448813197,{'M2':2,'S2':1,'K2':1,}],
'ST25':[54.407666621,{'M2':2,'N2':2,'K1':1,}],
'ST26':[54.635693740,{'M2':3,'N2':1,'K1':1,}],
'4MK9':[54.863720875,{'M2':4,'K1':1,}],
'ST27':[55.289258285,{'M2':3,'S2':1,'K1':1,}],
'ST28':[60.476138843,{'M2':4,'N2':1,}],
'M10':[60.704165962,{'M2':5,}],
'ST29':[60.901676254,{'M2':3,'N2':1,'S2':1,}],
'ST30':[61.129703388,{'M2':4,'S2':1,}],
'ST31':[61.361619262,{'M2':2,'N2':1,'S2':1,'K2':1,}],
'ST32':[61.555240814,{'M2':3,'S2':2,}],
'ST33':[67.430091484,{'M2':4,'S2':1,'K1':1,}],
'M12':[72.844999161,{'M2':6,}],
'ST34':[73.270536587,{'M2':5,'S2':1,}],
'ST35':[73.502452460,{'M2':3,'N2':1,'K2':1,'S2':1,}],
}
SAT_FAC={}
SAT_FAC['cons']=['O1','O1','O1','O1','O1','O1','O1','O1','P1','P1','P1','P1','P1','P1','K1','K1','K1','K1','K1','K1','K1','K1','K1','K1','N2','N2','N2','N2','M2','M2','M2','M2','M2','M2','M2','M2','M2','L2','L2','L2','L2','L2','S2','S2','S2','K2','K2','K2','K2','K2',]
SAT_FAC['amprat']=numpy.array([0.0003,0.0058,0.1885,0.0004,0.0029,0.0004,0.0064,0.001,0.0008,0.0112,0.0004,0.0004,0.0015,0.0003,0.0002,0.0001,0.0007,0.0001,0.0001,0.0198,0.1356,0.0029,0.0002,0.0001,0.0039,0.0008,0.0005,0.0373,0.0001,0.0004,0.0005,0.0373,0.0001,0.0009,0.0002,0.0006,0.0002,0.0366,0.0047,0.2505,0.1102,0.0156,0.0022,0.0001,0.0001,0.0024,0.0004,0.0128,0.298,0.0324,])
SAT_FAC['phcorr']=numpy.array([0.25,0.5,0,0.25,0.75,0.25,0.5,0.5,0,0.5,0.5,0.75,0.5,0.5,0,0.75,0.25,0.75,0,0.5,0,0.5,0.25,0.25,0.5,0,0,0.5,0.75,0.75,0,0.5,0.25,0.75,0.75,0,0,0.5,0,0.5,0.5,0.5,0,0.75,0,0.75,0.75,0.5,0,0,])
SAT_FAC['ilatfac']=numpy.array([1,0,0,1,1,1,0,0,0,0,0,1,0,0,0,1,1,1,0,0,0,0,1,1,0,0,0,0,2,2,0,0,2,2,2,0,0,0,0,0,0,0,0,2,0,2,2,0,0,0,])
SAT_FAC['deldood']=numpy.array([[-1,0,0,1,1,1,2,2,0,0,0,1,2,2,-2,-1,-1,-1,0,0,0,0,1,1,-2,-1,0,0,-1,-1,0,0,1,1,1,2,2,0,2,2,2,2,0,1,2,-1,-1,0,0,0,],
        [0,-2,-1,-1,0,1,0,1,-2,-1,0,0,0,1,-1,-1,0,1,-2,-1,1,2,0,1,-2,0,-2,-1,-1,0,-2,-1,-1,0,1,0,1,-1,-1,0,1,2,-1,0,0,0,1,-1,1,2,],
        [0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,]])


def astro_args(t):
    """list astro, list dastro = astro_args(datetime t)
    Computes the astronomical variables:
            astro=[tau,s,h,p,np,pp] (cycles) and time derivates dastro"""
    dz=t-datetime.datetime(1899,12,31,12,0,0)
    D=dz.days/10000.+dz.seconds/864000000.
    args=numpy.matrix([[1.],[10000.*D],[D*D],[D*D*D]])
    astro=(ASTRO_C*args)/360.0 % 1.
    tau=t.hour / 24. + t.minute / 1440. + t.second / 86400. + astro[1] - astro[0]
    args=numpy.matrix([[0.],[1.],[2.e-4*D],[3.e-4*D*D]])
    dastro=(ASTRO_C*args)/360.0 % 1.
    dtau=1.0 + dastro[1] - dastro[0]
    return numpy.hstack((tau.A1,astro.A1)), numpy.hstack((dtau.A1,dastro.A1))

def get_astro(icons,t0=datetime.datetime.now(),lat=None):
    icons=[ic.strip() for ic in icons]
    astro,ader=astro_args(t0)
    shallow=any([(len(CONS_STR[ic])==2 if CONS_STR.has_key(ic) else False)for ic in icons ])
    v=[]

    for ic in icons:
      cic=ic.strip()
      cic = ic.replace("\x00", "")  #added simon 10/11/2014
      cic = cic.replace("\x01", "")
      cic = cic.replace("\x80", "")
      # cic = ic.replace("\x00", "").replace("\x01", "").replace("\x80", "")
      if cic in CONS_STR:
        v.append(PI2* ((numpy.vdot(CONS_STR[cic][1:7],astro)+CONS_STR[cic][7]) % 1.) if len(CONS_STR[cic])>3 else 0.)
    
    if lat:
        if lat<5 and lat>=0:lat=5
        if lat>-5 and lat<0:lat=-5
        slat=math.sin(D2R*lat)
        rr=SAT_FAC['amprat']
        diufac=0.36309*(1.0-5.0*slat*slat)/slat
        rr=numpy.where(SAT_FAC['ilatfac']==1,diufac*rr,rr)
        rr=numpy.where(SAT_FAC['ilatfac']==2,2.59808*slat*rr,rr)
        uu=(numpy.dot(astro[3:],SAT_FAC['deldood'])+SAT_FAC['phcorr']) % 1.
        fsum=(1+0j)*numpy.ones(len(icons))
        if shallow:fssum=(1+0j)*numpy.ones(len(S_CONS))
        for isat,csat in enumerate(SAT_FAC['cons']):
            if csat in icons:
                fsum[icons.index(csat)]+=rr[isat]*numpy.exp(I2PI*uu[isat])
            if shallow and (csat in S_CONS):
                fssum[S_CONS.index(csat)]+=rr[isat]*numpy.exp(I2PI*uu[isat])
        f=abs(fsum)
        u=numpy.arctan2(fsum.imag,fsum.real)
        if shallow:
            fs=abs(fssum)
            us=numpy.arctan2(fssum.imag,fssum.real)
    else:
        u=numpy.zeros(len(icons))
        f=numpy.ones(len(icons))
    #Do shallow water v,u,f    
    if shallow:
        vs=[(numpy.vdot(CONS_STR[ic][1:7],astro)+CONS_STR[ic][7] % 1.) for ic in S_CONS]
        for ic,cc in enumerate(icons):
          if CONS_STR.has_key(cc):
            if len(CONS_STR[cc])==2:
              for isfac in CONS_STR[cc][1]:
                vsind=S_CONS.index(isfac)
                fac=CONS_STR[cc][1][isfac]
                v[ic]+=vs[vsind]*fac
                if lat:
                  u[ic]+=us[vsind]*fac
                  f[ic]*=math.pow(fs[vsind],fac)
                v[ic]=PI2*(v[ic] % 1.)
    return v,u,f
        

def get_freqs(icons):
    #return [CONS_STR.get(ic.strip(),[None])[0] for ic in icons]
    # return [CONS_STR.get(ic.replace("\x00", ""),[None])[0] for ic in icons] #Modif s.weppe 10/11/2014
    # convert coded string to "readable" format 
    freq = numpy.zeros(len(icons))
    for cnt,ic in enumerate(icons):
        cic = ic.replace("\x00", "")  #added simon 10/11/2014
        cic = cic.replace("\x01", "")
        cic = cic.replace("\x80", "")
       
        freq[cnt] =CONS_STR.get(cic,[None])[0]
    return freq

def ap2ep(u,v):
    u = u.amp*numpy.ma.exp(-1.j*u.pha)
    v = v.amp*numpy.ma.exp(-1.j*v.pha)
    
    wp=0.5*(u+1.j*v)
    wm=numpy.ma.conjugate(0.5*(u-1.j*v))
    thetap=numpy.ma.array(numpy.angle(wp),mask=wp.mask)
    thetam=numpy.ma.array(numpy.angle(wm),mask=wm.mask)
    wp=numpy.ma.abs(wp)
    wm=numpy.ma.abs(wm)
    
    sema=wp+wm
    semi=wp-wm
    pha=numpy.ma.mod(0.5*(thetam-thetap),2*numpy.pi)
    inc=numpy.ma.mod(0.5*(thetam+thetap),2*numpy.pi)
    k=numpy.fix(inc/numpy.pi)
    inc-=k*numpy.pi
    pha=numpy.ma.mod(pha*k*numpy.pi,2*numpy.pi)

    return sema,semi,inc,pha

#Base tidal data class
class TideStr:
    def __init__(self,amp,pha,icons=DEFAULT_CONS,t0=datetime.datetime.now(),lat=0.):
        if (amp.shape[0]<>pha.shape[0]) or (amp.shape[0]<>len(icons)):raise ValueError('Amplitude, phase and constituent arguments must have same number of members')
        self.cons=[ic.strip().upper() for ic in icons]
        self.amp=amp
        self.pha=pha
        self.freq=get_freqs(self.cons)
        self.t0=t0
        self.v0,self.u0,self.f=get_astro(self.cons,t0,lat)
        self.V=self.v0+self.u0
    
    #Return a timeseries of tidal quantity
    def ts(self,times,datum='msl'):
        if isinstance(times,datetime.datetime):times=[times]
        if isinstance(times,list):
            t1=numpy.array([(t-self.t0).days+(t-self.t0).seconds/86400. for t in times])
        elif isinstance(times,Times):
            t1=numpy.array(times.decdays(self.t0))
        else:
            return None
        ts1=numpy.ma.zeros((len(t1),)+self.amp.shape[1:])
        ndim=len(ts1.shape)
        one1=numpy.ones(ts1.shape)
        for i,f in enumerate(self.freq):
            ph1=self.pha[i]*one1
            f1=f*one1
            if ndim==3:
                ts1+=self.f[i]*self.amp[i]*numpy.ma.cos(f1*t1[:,numpy.newaxis,numpy.newaxis]-ph1+self.V[i])
            elif ndim==2:
                ts1+=self.f[i]*self.amp[i]*numpy.ma.cos(f1*t1[:,numpy.newaxis]-ph1+self.V[i])
            else:
                ts1+=self.f[i]*self.amp[i]*numpy.ma.cos(f1*t1-ph1+self.V[i])
            if (datum.lower()=='lat')and (f>0.):ts1+=self.amp[i]
        return ts1
    
    
