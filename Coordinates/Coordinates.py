#Coordinate.py
# Goal:
# 1) Provide a safe wrapper function for fCoord_tpt function

import numpy as np
import fCoord_tpt

def Hor2Sky(az,el,mjd,lat=28.300258*np.pi/180,lng=16.510112*np.pi/180.,gal=0,TPoints=None,inType=type(dict())):
    '''
    Return sky coordinates and parallactic angle

    Arguments
    az -- Focal azimuth coordinates (radians)
    el -- Focal elevation coordinates (radians)
    mjd -- Time in modified Julian Date

    Keyword Arguments
    lat -- Latitude of telescope, positive towards North (radians)
    lng -- Longitude of telescope, postive towards West (radians)
    gal -- Set to 0 for equatorial coordinates or 1 for Galactic
    TPoints -- T-Point parameters for the telescope horn (radians)

    Notes: The default for lat and lng are the position of the QUIJOTE
    MFI, Izana Observatory, Tenerife.
    
    '''

    if (type(TPoints) != type(dict())) and (type(TPoints) != inType):
        TPoints = {'xpos':0.,
                   'ypos':0.,
                   'Pf':0.,
                   'Px':0.,
                   'Py':0.,
                   'Pc':0.,
                   'Pn':0.,
                   'Pa':0.,
                   'Pb':0.}

    if (np.max(az) > 2.*np.pi) | (np.min(az) < 0.):
        print 'WARNING: AZIMUTH OUT OF RANGE (0 < az < 2pi)'
        return None

    if (np.max(el) > np.pi) | (np.min(el) < -np.pi):
        print 'WARNING: ELEVATION OUT OF RANGE (-pi < el < pi)'
        return None

    if np.max(mjd) > 2400000:
        print 'WARNING: ARE YOU USING MODIFIED JULIAN DATE?'

    ra,dec,p = fCoord_tpt.get_h2e_tpoint(az,el,
                                         TPoints['xpos'],TPoints['ypos'],TPoints['Pf'],
                                         TPoints['Px']  ,TPoints['Py']  ,TPoints['Pc'],
                                         TPoints['Pn']  ,TPoints['Pa']  ,TPoints['Pb'],
                                         lat,lng,mjd,gal)

    return np.mod(ra,2.*np.pi),dec,p
                                         
def Sky2Hor(ra,dec,mjd,lat=28.300258*np.pi/180,lng=16.510112*np.pi/180.):
    '''
    Return horizon coordinates

    Arguments
    ra -- Focal right ascension coordinates (radians)
    dec -- Focal declination coordinates (radians)
    mjd -- Time in modified Julian Date

    Keyword Arguments
    lat -- Latitude of telescope, positive towards North (radians)
    lng -- Longitude of telescope, postive towards West (radians)

    Notes: The default for lat and lng are the position of the QUIJOTE
    MFI, Izana Observatory, Tenerife.
    
    '''
    if (np.max(ra) > 2.*np.pi) | (np.min(ra) < 0.):
        print 'WARNING: RIGHT ASCENSION OUT OF RANGE (0 < ra < 2pi)'
        return None

    if (np.max(dec) > np.pi) | (np.min(dec) < -np.pi):
        print 'WARNING: DECLINATION OUT OF RANGE (-pi < dec < pi)'
        return None

    if np.max(mjd) > 2400000:
        print 'WARNING: ARE YOU USING MODIFIED JULIAN DATE?'

    az,el = fCoord_tpt.get_e2h(ra, dec, mjd, lng, lat)

    return np.mod(az,2.*np.pi),el


def Horizon(az,el,TPoints=None,inType=type(dict())):
    '''
    Return sky coordinates and parallactic angle

    Arguments
    az -- Focal azimuth coordinates (radians)
    el -- Focal elevation coordinates (radians)
    mjd -- Time in modified Julian Date

    Keyword Arguments
    lat -- Latitude of telescope, positive towards North (radians)
    lng -- Longitude of telescope, postive towards West (radians)
    gal -- Set to 0 for equatorial coordinates or 1 for Galactic
    TPoints -- T-Point parameters for the telescope horn (radians)

    Notes: The default for lat and lng are the position of the QUIJOTE
    MFI, Izana Observatory, Tenerife.
    
    '''

    if (type(TPoints) != type(dict())) and (type(TPoints) != inType):
        TPoints = {'xpos':0.,
                   'ypos':0.,
                   'Pf':0.,
                   'Px':0.,
                   'Py':0.,
                   'Pc':0.,
                   'Pn':0.,
                   'Pa':0.,
                   'Pb':0.}

    if (np.max(az) > 2.*np.pi) | (np.min(az) < 0.):
        print 'WARNING: AZIMUTH OUT OF RANGE (0 < az < 2pi)'
        return None

    if (np.max(el) > np.pi) | (np.min(el) < -np.pi):
        print 'WARNING: ELEVATION OUT OF RANGE (-pi < el < pi)'
        return None


    iaz,iel= fCoord_tpt.get_horizon_tpoint(az,el,
                                     TPoints['xpos'],TPoints['ypos'],TPoints['Pf'],
                                     TPoints['Px']  ,TPoints['Py']  ,TPoints['Pc'],
                                     TPoints['Pn']  ,TPoints['Pa']  ,TPoints['Pb'])

    return np.mod(iaz,2.*np.pi),iel
