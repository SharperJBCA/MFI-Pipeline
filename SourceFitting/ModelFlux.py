#ModelFlux
# Contain functions for modelling flux from known calibrator sources.
# See Weingartner 2001 (or 2011? Not sure)

import numpy as np

def CygAFlux(nu,mjd):
    '''
    Return expected flux of Cyg A for a given frequency and time in Jy.

    Arguments
    nu -- Frequency in GHz
    mjd -- Time in modified Julian Date (Not used)
    
    '''

    #Flux parameters:
    a = 1.482
    b = -1.2
    c = 0.

    if isinstance(mjd,type(list)) | isinstance(mjd,type(np.array)):
        return 10**(a + b*np.log10(nu/40.) + c*np.log10(nu/40.)**2 ) * np.ones(len(mjd))
    else:
        return 10**(a + b*np.log10(nu/40.) + c*np.log10(nu/40.)**2 ) 

    
def CasAFlux(nu,mjd):
    '''
    Return expected flux of Cas A for a given frequency and time in Jy.

    Arguments
    nu -- Frequency in GHz
    mjd -- Time in modified Julian Date
    
    '''
    
    #Flux parameters for year 2000
    a =  2.204
    b = -0.682
    c = 0.038

    alpha = (0.68 - 0.15 * np.log10(nu))/100.#5.3e-3 #per year

    return 10**(a + b*np.log10(nu/40.) + c*np.log10(nu/40.)**2 ) * CasASecular(nu,mjd)
    #return 179.9*(nu/33.)**-0.668 * np.exp( -alpha*Tobs)
    
def CasASecular(nu,mjd):
    '''
    Return expected flux of Tau A for a given frequency and time in Jy.

    Arguments
    mjd -- Time in modified Julian Date
    
    '''
    
    #Tobs = years from 2003
    Tobs = mjd/365.25 - ( 51543./365.25)
    alpha = (0.68 - 0.15 * np.log10(nu))/100.#5.3e-3 #per year

    return np.exp( -alpha*Tobs) #Macias-Perez 2010

def TauAFlux(nu,mjd):
    '''
    Return expected flux of Tau A for a given frequency and time in Jy.

    Arguments
    nu -- Frequency in GHz
    mjd -- Time in modified Julian Date
    
    '''
    
    #Tobs = years from 2003
    Tobs = mjd/365.25 -   ( 3. + 51543./365.25)

    #Flux parameters for year 2005
    a = 2.506
    b = -0.302
    c = 0.0
    #a =  2.506
    #b = -0.302
    #c = 0.0

    alpha = 0.22/100. #per year

    #print 'Percentage Flux Change:',(1. - alpha*Tobs)*100.,Tobs
    #print 'MEAN FLUX:', np.median( 10**(a + b*np.log10(nu/40.) ) * (1. - alpha*Tobs)), \
    #       np.median(10**(a + b*np.log10(nu/40.) )),np.mean(nu)

    return 973.*nu**(-0.296) * TauASecular(mjd)#* np.exp(-1.67e-3 * Tobs) #Macias-Perez 2010
    #return 10**(a + b*np.log10(nu/40.) ) * (1. - alpha*Tobs)


def TauASecular(mjd):
    '''
    Return expected flux of Tau A for a given frequency and time in Jy.

    Arguments
    mjd -- Time in modified Julian Date
    
    '''
    
    #Tobs = years from 2003
    Tobs = mjd/365.25 -   ( 3. + 51543./365.25)

    return np.exp(-1.67e-3 * Tobs) #Macias-Perez 2010
    #return 10**(a + b*np.log10(nu/40.) ) * (1. - alpha*Tobs)    


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
    inu = nu*1e9
    Jy = 1e26

    return 2.*k*inu**2/c**2 * beam * Jy
