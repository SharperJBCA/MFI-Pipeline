#Coordinate.py
# Goal:
# 1) Provide a safe wrapper function for fCoord_tpt function

import numpy as np
import fCoord_tpt

def Hor2Sky(az,el,mjd,lat=28.3002*np.pi/180,lng=16.5090*np.pi/180.,gal=0,TPoints=None):
    '''
    Return sky coordinates and parallactic angle
    '''

    if type(TPoints) != type(dict()):
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

    if np.max(jd) > 2400000:
        print 'WARNING: ARE YOU USING MODIFIED JULIAN DATE?'


    ra,dec,p = fCoord_tpt.get_h2e_tpoint(az,el,
                                         TPoints['xpos'],TPoints['ypos'],TPoints['Pf'],
                                         TPoints['Px']  ,TPoints['Py']  ,TPoints['Pc'],
                                         TPoints['Pn']  ,TPoints['Pa']  ,TPoints['Pb'],
                                         lat,lng,mjd,gal)

    return np.mod(ra,2.*np.pi),dec,p
                                         
