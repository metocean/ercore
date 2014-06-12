#Functions for sun rise and sun set
import numpy,datetime        

D2R=numpy.pi/180
M0 = 357.5291 * D2R
M1 = 0.98560028 * D2R
J0 = 0.0009
J1 = 0.0053
J2 = -0.0069
C1 = 1.9148 * D2R
C2 = 0.0200 * D2R
C3 = 0.0003 * D2R
P = 102.9372 * D2R
e = 23.45 * D2R
th0 = 280.1600 * D2R
th1 = 360.9856235 * D2R
h0 = -0.83 * D2R #sunset angle
d0 = 0.53 * D2R #sun diameter
h1 = -6 * D2R #nautical twilight angle
h2 = -12 * D2R #astronomical twilight angle
h3 = -18 * D2R #darkness angle
J2000 = 730121.0

#All functions below use radian lon,lat, time as NCEP time
def getJulianCycle( time, lon ): 
    return numpy.round(time - J2000 - J0 -0.5 + lon/(2 * numpy.pi)) 

def getApproxSolarTransit( Ht, lon, n ): 
    return J2000 + J0 + (Ht - lon)/(2 * numpy.pi) + n 

def getSolarMeanAnomaly( Js ): 
    return M0 + M1 * (Js - J2000 - 0.5)

def getEquationOfCenter( M ): 
    return C1 * numpy.sin(M) + C2 * numpy.sin(2 * M) + C3 * numpy.sin(3 * M) 

def getEclipticLongitude( M, C ): 
    return M + P + C + numpy.pi

def getSolarTransit( Js, M, Lsun ): 
    return Js + (J1 * numpy.sin(M)) + (J2 * numpy.sin(2 * Lsun)) 

def getSunDeclination( Lsun ): 
    return numpy.arcsin(numpy.sin(Lsun) * numpy.sin(e)) 

def getRightAscension( Lsun ):
    return numpy.arctan2(numpy.sin(Lsun) * numpy.cos(e), numpy.cos(Lsun))

def getSiderealTime( time, lon ):
    return th0 + th1 * (time -0.5 - J2000) + lon

def getAzimuth( th, a, lat, d ):
    H = th - a
    return numpy.arctan2(numpy.sin(H), numpy.cos(H) * numpy.sin(lat) - 
            numpy.tan(d) * numpy.cos(lat))

def getAltitude( th, a, lat, d ):
    H = th - a
    return numpy.arcsin(numpy.sin(lat) * numpy.sin(d) + 
            numpy.cos(lat) * numpy.cos(d) * numpy.cos(H))

def getHourAngle( h, lat, d ):
    return numpy.arccos((numpy.sin(h) - numpy.sin(lat) * numpy.sin(d)) / 
            (numpy.cos(lat) * numpy.cos(d))) 

def getSunsetJulianDate( w0, M, Lsun, lon, n ): 
    return getSolarTransit(getApproxSolarTransit(w0, lon, n), M, Lsun)+0.5 

def getSunriseJulianDate( Jtransit, Jset ): 
    return Jtransit - (Jset-0.5 - Jtransit)+0.5

   

class Sun(object):
    def __init__(self,lon,lat):
        self.lon=D2R*lon
        self.lat=D2R*lat
        
    def getTimes(self,date):
        n = getJulianCycle(date, self.lon)
        Js = getApproxSolarTransit(0, self.lon, n)
        M = getSolarMeanAnomaly(Js)
        C = getEquationOfCenter(M)
        Lsun = getEclipticLongitude(M, C)
        d = getSunDeclination(Lsun)
        Jtransit = getSolarTransit(Js, M, Lsun)
        w0 = getHourAngle(h0, self.lat, d)
        w1 = getHourAngle(h0 + d0, self.lat, d)
        Jset = getSunsetJulianDate(w0, M, Lsun, self.lon, n) #Set
        Jsetstart = getSunsetJulianDate(w1, M, Lsun, self.lon, n) #Set start 
        Jrise = getSunriseJulianDate(Jtransit, Jset) #Rise
        Jriseend = getSunriseJulianDate(Jtransit, Jsetstart) #Rise end
        w2 = getHourAngle(h1, self.lat, d)
        Jnau = getSunsetJulianDate(w2, M, Lsun, self.lon, n) #Dusk
        Jciv2 = getSunriseJulianDate(Jtransit, Jnau) #Dawn
        w3 = getHourAngle(h2, self.lat, d)
        w4 = getHourAngle(h3, self.lat, d)
        Jastro = getSunsetJulianDate(w3, M, Lsun, self.lon, n) #Nautical evening twilight end
        Jdark = getSunsetJulianDate(w4, M, Lsun, self.lon, n) #Astronomical evening twilight end
        Jnau2 = getSunriseJulianDate(Jtransit, Jastro) #Nautical morning twilight start
        Jastro2 = getSunriseJulianDate(Jtransit, Jdark) #Astronomical morning twilight start
        #AstroTwiStart->NautTwiStart->Dawn(CivilTwiStart)->Sunrise->SunriseEnd->SunsetStart->Sunset->Dusk(CivilTwiEnd)->NautTwiEnd->AstroTwiEnd

        return {
            'AstroTwiStart':Jastro2,
            'NautTwiStart':Jnau2,
            'Dawn':Jciv2,
            'Sunrise':Jrise,
            'SunriseEnd':Jriseend,
            'SunsetStart':Jsetstart,
            'Sunset':Jset,
            'Dusk':Jnau,
            'NautTwiEnd':Jastro,
            'AstroTwiEnd':Jdark
        }
        
    def getPosition(self,date):
        M = getSolarMeanAnomaly(date)
        C = getEquationOfCenter(M)
        Lsun = getEclipticLongitude(M, C)
        d = getSunDeclination(Lsun)
        a = getRightAscension(Lsun)
        th = getSiderealTime(date, self.lon)
        return {'azimuth':(90-getAzimuth( th, a, self.lat, d )/D2R) % 360.,'altitude':getAltitude( th, a, self.lat, d )/D2R}
            
    