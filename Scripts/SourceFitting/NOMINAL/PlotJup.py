import numpy as np
import pyfits
from DataAccess import ReadObs
from DataAccess import FileList
from DataAccess import TextFiles
from matplotlib import pyplot
import CalFitting
import WaveFitter
import Binning
import DataAccess
from SourceFitting.ModelFlux import toJy

def ReadBandpasses(mfidir):

    #Read in bandpasses:
    bandpasses = [np.loadtxt(mfidir+'Pol1_bandpass.dat',skiprows=1),
                  np.loadtxt(mfidir+'Pol2_bandpass.dat',skiprows=1),
                  np.loadtxt(mfidir+'Pol3_bandpass.dat',skiprows=1),
                  np.loadtxt(mfidir+'Pol4_bandpass.dat',skiprows=1)]

    for i in range(len(bandpasses)):
        nu = bandpasses[i][:,0]
        for j in range(bandpasses[i].shape[1]-1):
            G = bandpasses[i][:,j+1]

            print 'FREQ: ', np.sum(nu*G)/np.sum(G) 
    
    #Read in avg-rfactors
    rfactors = np.loadtxt(mfidir+'AvgRFactors.dat')

    freqs = bandpasses[0][:,0]
    
    band1 = np.array([bandpasses[0][:,1] + rfactors[0  ]*bandpasses[0][:,8],
                      bandpasses[0][:,2] + rfactors[1  ]*bandpasses[0][:,5],
                      bandpasses[0][:,3] + rfactors[2  ]*bandpasses[0][:,6],
                      bandpasses[0][:,4] + rfactors[3  ]*bandpasses[0][:,7],
                      bandpasses[1][:,1] + rfactors[4  ]*bandpasses[1][:,7],
                      bandpasses[1][:,2] + rfactors[5  ]*bandpasses[1][:,8],
                      bandpasses[1][:,3] + rfactors[6  ]*bandpasses[1][:,5],
                      bandpasses[1][:,4] + rfactors[7  ]*bandpasses[1][:,6],
                      bandpasses[2][:,1] + rfactors[4+4]*bandpasses[2][:,7],
                      bandpasses[2][:,2] + rfactors[5+4]*bandpasses[2][:,8],
                      bandpasses[2][:,3] + rfactors[6+4]*bandpasses[2][:,5],
                      bandpasses[2][:,4] + rfactors[7+4]*bandpasses[2][:,6],
                      bandpasses[3][:,1] + rfactors[4+8]*bandpasses[3][:,7],
                      bandpasses[3][:,2] + rfactors[5+8]*bandpasses[3][:,8],
                      bandpasses[3][:,3] + rfactors[6+8]*bandpasses[3][:,5],
                      bandpasses[3][:,4] + rfactors[7+8]*bandpasses[3][:,6]])

    return [freqs,band1]


if __name__ == '__main__':

    mfidir = '/nas/scratch/sharper/QUIJOTE/MFI-Pipeline/MFI_INFO/'

    #Define telescope constants:
    jd0 = 56244.
    nside = 128.

    beams = np.loadtxt(mfidir+'BeamSizes.dat')    
    beams = [(beams[0,1]*np.pi/180.)**2,
             (beams[1,1]*np.pi/180.)**2,
             (beams[2,1]*np.pi/180.)**2,
             (beams[3,1]*np.pi/180.)**2] #Beams of H2 to H4
    pbeam = 4.*np.pi/(12.*256.**2)

    freqs = np.loadtxt(mfidir+'FreqCens.dat')

    TDiode = np.loadtxt(mfidir+'TDiodes.dat')[:,1]

    rfactors = np.loadtxt(mfidir+'AvgRFactors.dat')

    channels = np.arange(4,dtype='i')
    chPairs = [[0,6],[1,7],[2,4],[3,5]]
    horns = range(1,2)
    
    chan = ['H219c','H217c','H219u','H217u','H313c','H311c','H313u','H311u','H419c','H417c','H419u','H417u']

    #Read in bandpasses:
    bandpasses = ReadBandpasses(mfidir)


    #d = np.loadtxt('CasA_GaussianFits_PerChannel_All',dtype='string')
    d = np.loadtxt('H2JupSeptember2013_Good',dtype='string')

    mjd   = d[:,-1].astype('f') + jd0
    peaks = d[:,1:33].astype('f')
    srcs  = np.array([c[0:4] for c in d[:,0]])
    dios  = d[:,1+32*2:33+32*2].astype('f')
    sums  = d[:,1+32:33+32].astype('f')

    for horn in horns:
        for pair in chPairs:


            ratio = (peaks[:,pair[0]+horn*8] + rfactors[8+pair[0]]*peaks[:,pair[1]+horn*8]) / \
                    ( dios[:,pair[0]+horn*8] + rfactors[8+pair[0]]* dios[:,pair[1]+horn*8])
            ratio = (sums[:,pair[0]+horn*8] + rfactors[8+pair[0]]*sums[:,pair[1]+horn*8]) / \
                    ( dios[:,pair[0]+horn*8] + rfactors[8+pair[0]]* dios[:,pair[1]+horn*8])

            pixbeam = 4.*np.pi/(12.*nside**2)/beams[horn]
            ratio = TDiode[horn*4 + pair[0]]*ratio*toJy(freqs[horn*4 + pair[0]],beams[horn])*pixbeam

            gd = np.where((ratio < 27.) & (ratio > 5.6) & (mjd < 56612) & (mjd > 56508))[0]
            bd = np.where((ratio > 5.) & (ratio < 25.))[0]
            bd = np.where((ratio > 20.) & (ratio < 60.)& (mjd < 56612))[0]

            #for ft in d[bd,0]:
            #    print ft

            print np.mean(ratio[bd]),np.std(ratio[bd]),np.std(ratio[bd])/np.mean(ratio[bd])*100./np.sqrt(bd.size)
            pyplot.plot(mjd,ratio,'o')
            pyplot.show()
