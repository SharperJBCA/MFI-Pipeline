

import numpy as np
from lmfit import minimize, Parameters, Parameter,fit_report
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from matplotlib import pyplot
from scipy.interpolate import griddata

import Beams

def gModelFit(P,x,y):
    '''
    Return model

    Arguments
    P -- Parameter dictionary
    x -- X grid coordinates
    y -- Y grid coordinates
    
    '''
    
    U = ((x-P['cx'] )/(P['wx']))**2 + ((y-P['cy'])/P['wy'])**2
    return  (P['amp']*np.exp(-0.5*U)) + P['bkgd']


def gauss_residual(P,x,y,data,e):
    '''
    Return array of residuals.

    Arguments
    P -- Parameter dictionary
    x -- X grid coordinates
    y -- Y grid coordinates
    data -- Data array to compare model to
    e -- Error on data values
    
    '''

    d = gModelFit(P.valuesdict(),x,y)
    
    return (d - data)/e

def ModelFit(P,x,y,d2):
    '''
    Return model

    Arguments
    P -- Parameter dictionary
    x -- X grid coordinates
    y -- Y grid coordinates
    
    '''
    #dx = 180./float(d2.shape[1])
    #dy = 180./float(d2.shape[0])

    #r = np.sqrt(x**2 + y**2)
    #maxval = 10**(map_coordinates(d2,[[90./dx],[90./dx]])/10.)
    #data  = 10**(map_coordinates(d2,[((x+P['cx'])+90.)/dx,((y+P['cy'])+90.)/dy])/10.)/maxval
    #maxval = 10**(map_coordinates(d2,[[90./dx],[90./dx]])/10.)
    #data  = map_coordinates(d2,[((x+P['cx'])+90.)/dx,((y+P['cy'])+90.)/dy])

    #pyplot.plot(x/dx+d2.shape[1]/2,y/dy+d2.shape[0]/2)
    #pyplot.show()
    #data  = map_coordinates(d2.T,[(x-P['cx'])/dx+d2.shape[1]/2,(y-P['cy'])/dy+d2.shape[0]/2])


    #from matplotlib import pyplot
    #pyplot.plot(d2['X'],d2['Y'])
    #pyplot.plot(x,y)
    #pyplot.show()
    #print d2['X'].shape,d2['Y'].shape,d2['DATA'].shape,x.shape,y.shape
    #data  = 10**griddata((-d2['X'],-d2['Y']),np.log10(d2['DATA']),(x,y),fill_value=2)
    #data[data == 100] = 0.
    
    #print P['width']
    #data = 10**(data/10.)
    #data /= np.max(data)

    #pyplot.plot(np.reshape(d2,d2.shape[0]**2),'-')
    #pyplot.figure()
    #pyplot.plot(data,'-')    
    #pyplot.show()

    #data = 10**(ip(x,y)/10.)/maxval#np.zeros(len(x))
    #for i in range(len(x)):
   #     data[i] = 10**(ip(x[i],y[i])/10.)/maxval
    #return  P['amp']*np.reshape(data,data.shape[0]*data.shape[1]) + P['bkgd']

    #data = Beams.BeamModelTOD(d2,x,y)
    
    return  P['amp']*d2 + P['bkgd']


def residual(P,x,y,data,e,d2):
    '''
    Return array of residuals.

    Arguments
    P -- Parameter dictionary
    x -- X grid coordinates
    y -- Y grid coordinates
    data -- Data array to compare model to
    e -- Error on data values
    
    '''

    d = ModelFit(P.valuesdict(),x,y,d2)

    #pyplot.plot(data)
    #pyplot.plot(d)
    #pyplot.show()

    return (d - data)/e

def FitSource(data,x,y,ch,beamMap=None,beam={'X':[0.92/2.355,True],'Y':[0.92/2.355,True]},rvals=None,filename='Out.dat',err=None):
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
    #d2  = np.loadtxt('/nas/scratch/sharper/QUIJOTE/MFI-Pipeline/SourceFitting/beam_centred.dat')
    #d2 = np.reshape(d2,(721,721))



    #Define some default values for rvals:
    if type(rvals) != type(dict()):
        rvals = {'r0':beam['X'][0]*3.5,
                 'r1':beam['X'][0]*3.5,
                 'r2':beam['X'][0]*3.,
                 'r3':beam['X'][0]*6.}


    #Define limits:
    isSource     = (x**2 + y**2 < rvals['r0']**2) 
    isGaussian   = (x**2 + y**2 < rvals['r1']**2) 
    isBackground = (x**2 + y**2 < rvals['r3']**2) & (x**2 + y**2 > rvals['r2']**2)
    isNotSource  = (x**2 + y**2 > rvals['r2']**2)    


    #r_samp = d2[:,0]
    #ip = interp1d(r_samp,np.mean(d2[:,[ch]],axis=1),bounds_error=False,fill_value=0.)
    #ip = interp2d(np.linspace(-180,180,1441),np.linspace(-90,90,721),d2,bounds_error=False,fill_value=0.)
    #ip = interp2d(np.linspace(-90,90,721),np.linspace(-90,90,721),np.reshape(d2,d2.shape[0]*d2.shape[1]),bounds_error=False,fill_value=0.,kind='cubic')
    #print ip
    #Check source exists:
    if data[isSource].size > 10:


        params = Parameters()
        params.add('bkgd',0.,vary=True  )#np.median(data)
        params.add('amp' ,np.max(data-np.median(data))  ,min=0)
        params.add('wx'  ,beam['X'][0],min=beam['X'][0]*0.78,max=beam['X'][0]*1.22,vary=beam['X'][1])
        params.add('wy'  ,beam['Y'][0],min=beam['Y'][0]*0.78,max=beam['Y'][0]*1.22,vary=beam['Y'][1])
        params.add('cx'  ,0.,min=-beam['X'][0]*0.75,max=beam['X'][0]*0.75,vary=False)
        params.add('cy'  ,0.,min=-beam['Y'][0]*0.75,max=beam['Y'][0]*0.75,vary=False)
        params.add('width'  ,1.,vary=True)

        #Index of the source peak:
        fwhm = np.where(isSource)[0]
        time = np.arange(data.size)
        maxarg = time[fwhm[(data[fwhm]).argmax()]]

        if type(err) == type(None):
            err= data*0.+1.

        if type(beamMap) == type(None):
            out = minimize(gauss_residual,params,args=(x [(isBackground | isGaussian)],
                                                       y [(isBackground | isGaussian)],
                                                       data[(isBackground | isGaussian)],
                                                       err[(isBackground | isGaussian)]),method='leastsq')
        else:
            out = minimize(residual,params,args=(x [(isBackground | isGaussian)]*np.pi/180.,
                                                 y [(isBackground | isGaussian)]*np.pi/180.,
                                                 data[(isBackground | isGaussian)],
                                                 err[(isBackground | isGaussian)],
                                                 beamMap[(isBackground | isGaussian)]),method='leastsq')



        vals = params.valuesdict()

        return vals
    else:
        return None
