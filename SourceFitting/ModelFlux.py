#ModelFlux
# Contain functions for modelling flux from known calibrator sources.
# See Weingartner 2001 (or 2011? Not sure)

import numpy as np

def CasAFlux(nu,mjd):
    '''
    Return expected flux of Cas A for a given frequency and time in Jy.

    Arguments
    nu -- Frequency in GHz
    mjd -- Time in modified Julian Date
    
    '''
    
    #Tobs = years from 2000
    Tobs = (mjd - 51543.)/365.25

    #Flux parameters for year 2000
    a =  2.204
    b = -0.682
    c = 0.038

    alpha = 5.3e-3 #per year

    return 10**(a + b*np.log10(nu/40.) + c*np.log10(nu/40.)**2 ) * (1 - alpha*Tobs)

def TauAFlux(nu,mjd):
    '''
    Return expected flux of Tau A for a given frequency and time in Jy.

    Arguments
    nu -- Frequency in GHz
    mjd -- Time in modified Julian Date
    
    '''
    
    #Tobs = years from 2005
    Tobs = (mjd - 51543.)/365.25 + 5.

    #Flux parameters for year 2000
    a =  2.502
    b = -0.35
    c = 0.0
    #a =  2.506
    #b = -0.302
    #c = 0.0

    alpha = 2.2e-3 #per year

    return 10**(a + b*np.log10(nu/40.) + c*np.log10(nu/40.)**2 ) * (1. - alpha*Tobs)


def SourceCoord(source,gal=False):
    '''
    Returns RA/DEC of source

    Arguments
    source -- Name of a known source
    
    '''    

    #Read in sources file
    dir =  __file__.split('/')
    srcs = np.loadtxt('/'.join(dir[:-1])+'/sources.lis',dtype='string')

    src_dict = {s[0]:[float(s[1]),float(s[2]),float(s[3]),float(s[4])] for s in srcs}

    
    if gal:
        return src_dict[source][2],src_dict[source][3]
    else:
        return src_dict[source][0],src_dict[source][1]

def EphemCoord(source,mjd,lat=28.300258*np.pi/180,lng=16.510112*np.pi/180.):
    '''
    Returns RA/DEC of ephemeris source (radians)

    Arguments
    source -- Name of a known source
    mjd -- MJD of observation

    Keyword Arguments
    lng -- Longitude of telescope (radians,+W) (Default: Izana Observatory)
    lat -- Latitude  of telescope (radians,+N) (Default: Izana Observatory)

    '''
    import Coordinates  

    Planets = {'Sun':0,
               'Mercury':1,
               'Venus':2,
               'Moon':3,
               'Mars':4,
               'Jupiter':5,
               'Saturn':6,
               'Uranus':7,
               'Neptune':8}

    cRa,cDec = Coordinates.Ephem.plpos(mjd,lng,lat,Planets[source])

    return cRa,cDec
    

    

def toJy(nu,beam):
    '''
    Return Jy/K for a given beam and frequency

    Arguments
    nu -- frequency in GHz
    beam -- beam in steradians.
    
    '''

    k = 1.3806488e-23
    c = 299792458.
    nu *= 1e9
    Jy = 1e26

    return 2.*k*nu**2/c**2 * beam * Jy
