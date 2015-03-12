#DiodeTemp.py
# Goals
# 1) Determine the Calibration Diode temperature from calibrator observations


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

import SourceFitting
from SourceFitting.ModelFlux import toJy

def SourceFlux(nu_bp,G_bp,mjd,srcs):
    '''
    Return list of fluxes for Tau A or Cas A
    '''
    #First determine central frequency:
    dnu = nu_bp[1] - nu_bp[0]
    nu_c = np.sum(nu_bp*G_bp)/np.sum(G_bp)
    BW = nu_bp*np.sum(G_bp)**2/np.sum(G_bp**2)



    #Central frequency of Tau A:
    spec_crab = SourceFitting.ModelFlux.TauAFlux(nu_bp,mjd[0])*10.
    nu_c_crab = np.sum(nu_bp*G_bp*spec_crab)/np.sum(G_bp*spec_crab)

    #Central frequency of Cas A:
    spec_casa = SourceFitting.ModelFlux.CasAFlux(nu_bp,mjd[0])*10.
    nu_c_casa = np.sum(nu_bp*G_bp*spec_casa)/np.sum(G_bp*spec_casa)     


    #Which vals are Cas/Tau:
    crabs = np.where((srcs == 'CRAB'))[0]
    casas = np.where((srcs == 'CASS'))[0]

    fluxs = np.zeros(len(srcs))
    corfact = np.zeros(len(srcs))

    fluxs[crabs] = SourceFitting.ModelFlux.TauAFlux(nu_c,mjd[crabs])
    fluxs[casas] = SourceFitting.ModelFlux.CasAFlux(nu_c,mjd[casas])

    
    print 'WMAP 22.8:',SourceFitting.ModelFlux.TauAFlux(22.8,mjd[0])
    print '13:',SourceFitting.ModelFlux.TauAFlux(13.,mjd[0])

    #Ratio between Source at nu_c and Source at nu_c_colour
    corfact[crabs] = SourceFitting.ModelFlux.TauAFlux(nu_c,mjd[crabs])/ \
                     SourceFitting.ModelFlux.TauAFlux(nu_c_crab,mjd[crabs]) 
    
    corfact[casas] = SourceFitting.ModelFlux.CasAFlux(nu_c,mjd[casas])/ \
                     SourceFitting.ModelFlux.CasAFlux(nu_c_casa,mjd[casas])


    return fluxs,corfact,nu_c

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

            #print 'FREQ: ', np.sum(nu*G)/np.sum(G) 
    
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


    beams = np.loadtxt(mfidir+'BeamSizes.dat')    
    beams = [(beams[0,1]*np.pi/180.)**2,
             (beams[1,1]*np.pi/180.)**2,
             (beams[2,1]*np.pi/180.)**2,
             (beams[3,1]*np.pi/180.)**2] #Beams of H2 to H4
    pbeam = 4.*np.pi/(12.*256.**2)
    TDiode = np.loadtxt(mfidir+'TDiodes.dat')[:,1]

    freqs = np.loadtxt(mfidir+'FreqCens.dat')


    rfactors = np.loadtxt(mfidir+'AvgRFactors.dat')

    channels = np.arange(4,dtype='i')
    chPairs = [[0,6],[1,7],[2,4],[3,5]]
    horns = range(2,3)
    
    chan = ['H219c','H217c','H219u','H217u','H313c','H311c','H313u','H311u','H419c','H417c','H419u','H417u']

    #Read in bandpasses:
    bandpasses = ReadBandpasses(mfidir)

    #Read in the data:
    taua = np.loadtxt('TauA_Nom60',dtype='string')

    mjd   = taua[:,-1].astype('f') + jd0
    peaks = taua[:,1:5].astype('f')
    srcs  = np.array(['CRAB' for c in taua[:,0]])
    sums = taua[:,1+4:5+4].astype('f')

    pbeam = 4.*np.pi/(12.*128.**2)
    
    for horn in horns:
        for pair in chPairs:
            Fluxs,Corfacts,nu_c = SourceFlux(bandpasses[0],bandpasses[1][pair[0]+horn*4,:],mjd,srcs)
            TDio = peaks[:,pair[0]]*toJy(freqs[horn*4 + pair[0]],beams[horn])#/0.86
            Sums = sums[:,pair[0]]*toJy(freqs[horn*4 + pair[0]],beams[horn])*pbeam/beams[horn]

            #print TDio, Fluxs
            gd = np.where(TDio > 0.)[0]
            #print horn,np.median(TDio/Fluxs)
            
            pyplot.plot(mjd[gd],Fluxs[gd]/Sums[gd],'o',color='gray',label='Diode Temperature',zorder=0,markersize=8,alpha=0.8)
            pyplot.plot(mjd[gd],Fluxs[gd]/TDio[gd],'o',color='blue',label='Diode Temperature',zorder=1,markersize=8,alpha=0.8)
            pyplot.ylim(0,3)

            pyplot.show()
    
    #for ch in channels:
    #    print chan[ch], '&', '%2.3f' % np.median(dios[:,ch]), '\\\\'
