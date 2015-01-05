#ModelFlux
# Contain functions for modelling flux from known calibrator sources.
# See Weingartner 2001 (or 2011? Not sure)

import numpy as np

def CasAFlux(nu,mjd):
    '''
    Return expected flux of Cas A for a given frequency and time.

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
