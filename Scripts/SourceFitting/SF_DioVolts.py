#SourceFitting.py
# Goals:
# 1) Fit to the peak of source.
# 2) Save polarisation corrected diode signal.

import numpy as np
import pyfits
from DataAccess import ReadObs
from DataAccess import FileList
from DataAccess import TextFiles
from matplotlib import pyplot
import CalFitting
import WaveFitter
import Binning

import SourceFitting
import Coordinates

import scipy.interpolate as ip

import scipy.optimize as op
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

    if (P[3] > 0):
        P[2] = np.mod(P[2],1.)
        return (FitSine(P,X) - Y)
    else:
        return Y*0. + 1e24



def residual(P,X,Y,Pup):

    if (P[0] >0) & (P[0]< Pup):

        nsteps = 50
        stepsize = (X.size)/50

        #Sort the data based on time P0
        sd = (np.mod(X,P[0])).argsort()

        stds = 0.
        N = Y[sd[:nsteps*stepsize]]
        N = np.reshape(N,(nsteps,stepsize))
        N = np.std(N,axis=1)
        stds = np.sum(N**2)
        
        #for n in range(nsteps):
        #    stds += np.std(Y[sd[n*stepsize:(n+1)*stepsize]])**2

        return np.sqrt(stds)
    else:
        return 1e24

def FindWave(time,vals,WL=20.,Pup=None):
    '''
    Return best fit cosine wave for discontinous periodic data.

    Arguments
    time -- Array containing time coordinate of each vals sample
    vals -- Measured data array

    Keyword Arguments
    WL -- Initial guess at Wavelength (better to underestimate!)
    Pup -- Upper limit to wavelength guess (Default: 2*WL)
    '''

    P0 = [WL]

    if Pup == None:
        Pup = WL * 2.



    #Take off the first and last 10% of data due to end effects:
    start = int(time.size*0.1)
    end = time.size-start
    minjd = np.min(time)

    X = time[start:end]-minjd
    Y = vals[start:end]    

    #Loop the basin hopping until inital guess is good enough:
    while True:

        #Find initial guess at wavelength
        Pbasin = op.basinhopping(residual,P0,minimizer_kwargs={'args':(X,Y,Pup)},stepsize=4.,T = 1.)

        #Sort the data according to guess at wavelength
        sd = (np.mod(X,Pbasin['x'])).argsort()

        #Estimate the phase
        maxval_time = np.argmax(Y[sd])
        Phase = (Pbasin['x'] - np.mod(X[sd[maxval_time]],Pbasin['x']))/Pbasin['x']

        #Generate least-squares initial parameters
        Wavelength = Pbasin['x']
        Amplitude = np.max(Y) - np.median(Y)
        P0a = [Amplitude,np.median(Y),Phase,Wavelength]

        #Fit wavelength of original X and Y data.
        P1, s = op.leastsq(Error,P0a,args=(X,Y))

        #If fit is close to inital guess, accept fit.
        if P1[1] > Amplitude*0.75:
            break

    return P1



if __name__ == "__main__":

    #Define data files:
    datadir = '../DownSampData/data/'
    mfidir  = '/nas/scratch/sharper/QUIJOTE/MFI-Pipeline/MFI_INFO/'

    dioModel = np.loadtxt('/nas/scratch/sharper/QUIJOTE/MFI-Pipeline/MFI_INFO/DiodePolModel-Epoch2.dat')
    
    files = np.loadtxt('Cass',dtype='string',ndmin=1)

    #Read in T-Point parameters:
    tpoints = np.loadtxt(mfidir+'focal_plane.ini')
    tpoints *= np.pi/180.
    jd0 = 56244.

    chanPairs = np.array([[0,6],[1,7],[2,4],[3,5]])
    beams     = np.array([0.88,0.66,0.88,0.66])
    
    Params = np.zeros(12)
    rFacts = np.zeros(12)

    for f in files:
        print 'File',f

        #Read in data:
        filelist = FileList(f,dir=datadir)
        data = ReadObs.FullObs(filelist,['DATA','AZ','EL','JD','CAL','MOD','CAL_JD','ERR'],dir='')

        meanJD = np.mean(data['JD'][0,:,0])

        for horn in range(2,3):
            for ch,chPair in enumerate(chanPairs):

                Params[(horn-1)*4 + ch] = np.median(data['CAL'][0,:,chPair[0]+horn*8])
                rFacts[(horn-1)*4 + ch] = np.median(4.*np.mod(data['MOD'][0,:,horn],90.))
                       
        TextFiles.AppendFile('DIO_VOLTS',[f,
                                         Params[0],Params[1],Params[2],Params[3],
                                         Params[4],Params[5],Params[6],Params[7],
                                         Params[8],Params[9],Params[10],Params[11],
                                         rFacts[0],rFacts[1],rFacts[2],rFacts[3],
                                         rFacts[4],rFacts[5],rFacts[6],rFacts[7],
                                         rFacts[8],rFacts[9],rFacts[10],rFacts[11],meanJD])
