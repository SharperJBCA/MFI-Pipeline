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

    bins,errs = DSD.downsample(a,newlen)

    if Errors:
        return bins, errs
    else:
        return bins

def bindata(X,Y,nbins):

    """
    A fast binning routine for regularly spaced data.

    Returns the bin values (Y) and the bin centers (X).

    """

    ybins = np.linspace(Y.min(),Y.max(),nbins)
    xbins = np.linspace(X.min(),X.max(),nbins)

    digitized = np.digitize(X,xbins)
    
    bins    = np.array([np.median(Y[digitized == i]) for i in range(1,nbins)])
    binmids = np.array([np.median(X[digitized == i]) for i in range(1,nbins)])
    binerr  = np.array([np.std(Y[digitized == i]) for i in range(1,nbins)])

    binmids = binmids[np.where(np.isnan(bins) == 0)[0]]
    binerr  = binerr[np.where(np.isnan(bins) == 0)[0]]
    bins    = bins[np.where(np.isnan(bins) == 0)[0]]

    return bins,binmids,binerr
