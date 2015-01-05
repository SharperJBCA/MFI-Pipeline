

import numpy as np
from lmfit import minimize, Parameters, Parameter,fit_report

def ModelFit(P,x,y):
    '''
    Return model

    Arguments
    P -- Parameter dictionary
    x -- X grid coordinates
    y -- Y grid coordinates
    
    '''
    
    U = ((x-P['cx'] )/(P['wx']))**2 + ((y-P['cy'])/P['wy'])**2
    return  (P['amp']*np.exp(-0.5*U)) + P['bkgd']


def residual(P,x,y,data,e):
    '''
    Return array of residuals.

    Arguments
    P -- Parameter dictionary
    x -- X grid coordinates
    y -- Y grid coordinates
    data -- Data array to compare model to
    e -- Error on data values
    
    '''

    d = ModelFit(P.valuesdict(),x,y)
    
    return (d - data)/e

def FitSource(data,x,y,beam=0.92/2.355,rvals=None,filename='Out.dat'):
    '''
    Return array of parameters defining a source

    Arguments
    data -- Input data to fit to
    x -- X grid coordinates
    y -- Y grid coordinates


    Keyword Arguments
    beam -- FWHM/2.355 of telescope beam
    rvals -- Dictionary of several fitting limits ({"r0":,"r1","r2","r3"})
    filename -- File to write parameters to (Default: Out.dat)


    Notes:
    Use of rvals:

    r0 -- Used to define whether source is within the input data.
    
    r1 -- Defines the extent of beam that can be accurately described
          as a Gaussian.

    r2 -- Defines inner radius of background.

    r3 -- Defines outer radius of background.

    '''

    #Define some default values for rvals:
    if type(rvals) != type(dict()):
        rvals = {'r0':beam*3.5,
                 'r1':beam*1.175,
                 'r2':beam*3.,
                 'r3':beam*6.}


    #Define limits:
    isSource     = (xi**2 + yi**2 < r0**2) 
    isGaussian   = (xi**2 + yi**2 < r1**2) 
    isBackground = (xi**2 + yi**2 < r3**2) & (xi**2 + yi**2 > r2**2)
    isNotSource  = (xi**2 + yi**2 > r2**2)    

    #Check source exists:
    if data[isSource].size > 10:

        params = Parameters()
        params.add('bkgd',np.median(data)  )
        params.add('amp' ,np.max(data-np.median(data))  ,min=0)
        params.add('wx'  ,beam,min=beam*0.95,max=beam*1.05,vary=True)
        params.add('wy'  ,beam,min=beam*0.95,max=beam*1.05,vary=True)
        params.add('cx'  ,0.,vary=False)
        params.add('cy'  ,0.,vary=False)

        #Index of the source peak:
        fwhm = np.where(isSource)[0]
        time = np.arange(data.size)
        maxarg = time[fwhm[(data[fwhm]).argmax()]]

        out = minimize(residual,params,args=(xi [(isBackground | isGaussian)],
                                             yi [(isBackground | isGaussian)],
                                             tod[(isBackground | isGaussian)],
                                             err[(isBackground | isGaussian)],),method='leastsq')
        vals = param.valuedict()

        return param.valuedict()
    else:
        return None
