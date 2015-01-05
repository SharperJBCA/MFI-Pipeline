import numpy as np
import scipy.optimize as sp

import scipy.interpolate as ip
import scipy.fftpack as sfft

def FitSine(P,X):
    '''
    Return a sine function defined by P.

    Arguments
    P -- Parameters of sine wave
    X -- Frequencies

    (P[0] = background, P[1] = Amplitude, P[2] = Phase, P[3] = Wavelength)
    
    
    '''

    return P[0] + P[1]*np.cos( 2.*np.pi * (X)/P[3] + P[2]*2*np.pi )

def Error(P,X,Y):
    '''
    Return array of residuals
    '''

    P[3] = 360.
    if (P[3] > 0):
        P[2] = np.mod(P[2],1.)
        return (FitSine(P,X) - Y)
    else:
        return Y*0. + 1e24


def FFTMethod(x_in,y_in):
    '''
    Return best fit parameters to a discontinuous sinewave

    Arguments
    x_in -- time/angle/etc...
    y_in -- data
    '''

    from lmfit import minimize, Parameters, Parameter,fit_report


    #if y_in.size < 40000:
    #    return None

    sd = x_in.argsort()

    #Ensure data is ordered in time for FFT:
    X = x_in[sd]
    Y = y_in[sd]

    #Only use the middle 98% incase of end effects:

    start = len(X)*0.#.01
    end   = len(X)*1.#0.99

    X = X[start:end]
    Y = Y[start:end]

    pfit = np.poly1d(np.polyfit(X,Y,3))

    Xt = np.linspace(0,360,X.size)
    Yt = pfit(Xt)


    #Transform data to estimate wavelength and phase:
    fdata = sfft.fft(Yt-np.mean(Yt))
    freqs= sfft.fftfreq(fdata.size, d= (Xt[1] - Xt[0])/(2.*np.pi))

    fdata = fdata[1:len(freqs)/2]
    freqs = freqs[1:len(freqs)/2]
    freqs = freqs[::-1]

        
    fabs = np.abs(fdata)**2
    idx  = fabs.argmax() - 1



    wave  = np.abs((2.*np.pi)/freqs[idx])
    phase = np.mod(np.angle(fdata[idx])/(2.*np.pi),1)



    mid = (np.max(Y)+np.min(Y))/2.
    amp = np.max(Y)-mid

    P0 = [mid ,amp, phase ,wave]

    P1, s = sp.leastsq(Error,P0,args=(X,Y))

    if P1[1] < (np.max(Y)-np.mean(Y))/100.:
        P0 = [mid ,amp, np.mod(phase+0.5,1) ,wave]
        P1, s = sp.leastsq(Error,P0,args=(X,Y))


    return P1
    