#Ephem.py
# Wrapper for fEphem function.
# 1) Return Healpix coordinate for a planet

import numpy as np
import fEphem

def plpos(mjd,lon,lat,planet,gal=0):
    '''
    Return RA/DEC (or l/b) coordinates of a planet

    Arguments
    mjd   -- Modified julian date
    lon   -- Longitude of telescope (Positive West, radians)
    lat   -- Latitude of telescope (Positive North, radians)
    planet-- Planet number (see Notes)

    Keyword Arguments
    gal   -- Flag to return galactic coordinates (default = 0)


    Notes: <=0 - Sun
             1 - Mercury
             2 - Venus
             3 - Earth
             4 - Mars
             5 - Jupiter
             6 - Saturn
    '''

    rCen, dCen = fEphem.plpos(float(mjd),
                              float(lon),
                              float(lat),
                              int(planet),
                              int(gal))

    print rCen*180./np.pi,dCen*180./np.pi
    return np.mod(rCen,2.*np.pi), dCen
