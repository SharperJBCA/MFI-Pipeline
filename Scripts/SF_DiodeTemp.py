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

    #Ratio between Source at nu_c and Source at nu_c_colour
    corfact[crabs] = SourceFitting.ModelFlux.TauAFlux(nu_c,mjd[crabs])/ \
                     SourceFitting.ModelFlux.TauAFlux(nu_c_crab,mjd[crabs]) 
    
    corfact[casas] = SourceFitting.ModelFlux.CasAFlux(nu_c,mjd[casas])/ \
                     SourceFitting.ModelFlux.CasAFlux(nu_c_casa,mjd[casas])


    return fluxs,corfact,nu_c

def ReadBandpasses(mfidir):

    #Read in bandpasses:
    bandpasses = [np.loadtxt(mfidir+'Pol2_bandpass.dat',skiprows=1),
                  np.loadtxt(mfidir+'Pol3_bandpass.dat',skiprows=1),
                  np.loadtxt(mfidir+'Pol4_bandpass.dat',skiprows=1)]
    
    #Read in avg-rfactors
    rfactors = np.loadtxt(mfidir+'AvgRFactors.dat')

    freqs = bandpasses[0][:,0]
    band1 = np.array([bandpasses[0][:,1] + rfactors[4  ]*bandpasses[0][:,7],
                      bandpasses[0][:,2] + rfactors[5  ]*bandpasses[0][:,8],
                      bandpasses[0][:,3] + rfactors[6  ]*bandpasses[0][:,5],
                      bandpasses[0][:,4] + rfactors[7  ]*bandpasses[0][:,6],
                      bandpasses[1][:,1] + rfactors[4+4]*bandpasses[1][:,7],
                      bandpasses[1][:,2] + rfactors[5+4]*bandpasses[1][:,8],
                      bandpasses[1][:,3] + rfactors[6+4]*bandpasses[1][:,5],
                      bandpasses[1][:,4] + rfactors[7+4]*bandpasses[1][:,6],
                      bandpasses[2][:,1] + rfactors[4+8]*bandpasses[2][:,7],
                      bandpasses[2][:,2] + rfactors[5+8]*bandpasses[2][:,8],
                      bandpasses[2][:,3] + rfactors[6+8]*bandpasses[2][:,5],
                      bandpasses[2][:,4] + rfactors[7+8]*bandpasses[2][:,6]])

    return [freqs,band1]

    
if __name__ == '__main__':

    mfidir = '/nas/scratch/sharper/QUIJOTE/MFI-Pipeline/MFI_INFO/'

    #Define telescope constants:
    jd0 = 56244.
    
    beams = [0.99*(0.64*np.pi/180.)**2,
             0.92*(0.85*np.pi/180.)**2,
             0.99*(0.64*np.pi/180.)**2] #Beams of H2 to H4
    
    freqs = [[12.891,11.605,12.881,11.149],
             [12.891,11.605,12.881,11.149],
             [12.891,11.605,12.881,11.149]] #Freqs of H2 to H4

    channels = np.arange(12,dtype='i')
    chan = ['H219c','H217c','H219u','H217u','H313c','H311c','H313u','H311u','H419c','H417c','H419u','H417u']

    #Read in bandpasses:
    bandpasses = ReadBandpasses(mfidir)

    #Read in the data:
    #crab = np.loadtxt('TauA',dtype='string')
    crab = np.loadtxt('TauA',dtype='string')   
    casa = np.loadtxt('CasA',dtype='string')

    el    = np.append(crab[:,25].astype('f'),casa[:,25].astype('f')) + jd0
    mjd   = np.append(crab[:,26].astype('f'),casa[:,26].astype('f')) + jd0
    peaks = np.reshape(np.append(crab[:,1:13].astype('f'),casa[:,1:13].astype('f')),(len(mjd),12))
    srcs  = np.append(np.array([c[0:4] for c in crab[:,0]]),
                      np.array([c[0:4] for c in casa[:,0]]))
    dios  = np.reshape(np.append(crab[:,1+12:13+12].astype('f'),casa[:,1+12:13+12].astype('f')),(len(mjd),12))
    para  = np.append(crab[:,27].astype('f'),casa[:,27].astype('f'))



    gd    = np.ones(len(srcs),dtype='bool')

    #Fluxs,Corfacts = SourceFlux(bandpasses[0],bandpasses[1][2,:],mjd,srcs)
    #TDio = Fluxs/(peaks[:,2]*Corfacts)/1000.
    #gd    = gd & (TDio < 0.69) & (TDio > 100.52)
    
    #Fluxs,Corfacts = SourceFlux(bandpasses[0],bandpasses[1][5,:],mjd,srcs)
    #TDio = Fluxs/(peaks[:,5]*Corfacts)/1000.
    #gd    = gd & (TDio > 100)

    
    crabs = np.where((srcs[gd] == 'CRAB'))[0]

    for ch in channels:
        flatSpec = np.zeros(len(bandpasses[1][ch,:]))
        flatSpec[np.where(np.abs(bandpasses[0]-13.) < 1.)[0]] = 1.
        Fluxs,Corfacts,nu_c = SourceFlux(bandpasses[0],bandpasses[1][ch,:],mjd,srcs)
        
        TDio = Fluxs/(peaks[:,ch]*Corfacts)/toJy(nu_c,beams[ch/4])

        gd = (np.abs(TDio-np.median(TDio)) < 0.15) 
        pfit = np.poly1d(np.polyfit(np.arange(len(TDio[gd])), TDio[gd],1))
        gd = np.where(gd)[0]

        print ch+4,np.median(TDio[gd])#,np.std(TDio[gd])/np.median(TDio[gd]) * 100.,np.median(Corfacts)

        #print np.append(crab[:,0],casa[:,0])[(TDio[gd] > 1.31)]
        
        #pyplot.plot(mjd[gd],TDio[gd],'o')
        #pyplot.plot(mjd[gd[crabs]],TDio[gd[crabs]],'o')

        #pyplot.plot(para[gd],TDio[gd],'o')


        #pyplot.show()
    
    #for ch in channels:
    #    print chan[ch], '&', '%2.3f' % np.median(dios[:,ch]), '\\\\'
