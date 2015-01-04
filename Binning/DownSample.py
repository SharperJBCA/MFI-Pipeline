#Write a save wrapper for the DSD.f90 routine

import numpy as np
import DSD

def DownSample(a,newlen,Errors=False):
    '''
    Return binned version of input array

    Arguments
    a -- Input array
    newlen -- Length of binned output array

    Keyword Arguments
    Errors -- Return errors on bin values (Default: False)


    Notes: The last datum of the output array may contain less samples
    than the rest of the data in the output array.

    '''

    if newlen > a.size:
        print 'WARNING: BIN SIZE GREATER THAN ARRAY SIZE'
        return None

    bins,errs = DSD.DownSample(a,newlen)

    if Errors:
        return bins, errs
    else:
        return bins
