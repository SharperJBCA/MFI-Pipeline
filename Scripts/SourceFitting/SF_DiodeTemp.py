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
    spec_crab = SourceFitting.ModelFlux.TauAFlux(nu_bp,mjd[0])
    nu_c_crab = np.sum(nu_bp*G_bp*spec_crab)/np.sum(G_bp*spec_crab)

    #Central frequency of Cas A:
    spec_casa = SourceFitting.ModelFlux.CasAFlux(nu_bp,mjd[0])
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


    beams = np.loadtxt(mfidir+'BeamSizes.dat')    
    beams = [1.133*(beams[0,1]*np.pi/180.)**2,
             1.133*(beams[1,1]*np.pi/180.)**2,
             1.133*(beams[2,1]*np.pi/180.)**2,
             1.133*(beams[3,1]*np.pi/180.)**2] #Beams of H2 to H4
    pbeam = 4.*np.pi/(12.*256.**2)
    TDiode = np.loadtxt(mfidir+'TDiodes.dat')[:,1]

    freqs = np.loadtxt(mfidir+'FreqCens.dat')


    rfactors = np.loadtxt(mfidir+'AvgRFactors.dat')

    channels = np.arange(4,dtype='i')
    chPairs = [[0,6],[1,7],[2,4],[3,5]]
    horns = range(0,4)
    
    chan = ['H219c','H217c','H219u','H217u','H313c','H311c','H313u','H311u','H419c','H417c','H419u','H417u']

    #Read in bandpasses:
    bandpasses = ReadBandpasses(mfidir)

    #Read in the data:
    casa = np.loadtxt('CasA_GaussianFits_PerChannel_All',dtype='string')
    taua = np.loadtxt('TauA_GaussianFits_PerChannel_All',dtype='string')

    mjd   = np.concatenate((taua[:,-1].astype('f'),casa[:,-1].astype('f'))) + jd0
    peaks = np.concatenate((taua[:,1:33].astype('f'),casa[:,1:33].astype('f')))
    srcs  = np.concatenate((np.array([c[0:4] for c in taua[:,0]]),np.array([c[0:4] for c in casa[:,0]])))
    dios  = np.concatenate((taua[:,1+32*2:33+32*2].astype('f'),casa[:,1+32*2:33+32*2].astype('f')))
    sums  = np.concatenate((taua[:,1+32:33+32].astype('f'),casa[:,1+32:33+32].astype('f')))


    crabs = np.where(srcs == 'CRAB')[0]
    casss = np.where(srcs == 'CASS')[0]

    for horn in horns:
        for pair in chPairs:
            Fluxs,Corfacts,nu_c = SourceFlux(bandpasses[0],bandpasses[1][pair[0]+horn*4,:],mjd,srcs)

            apRatio =(peaks[:,pair[0]+horn*8] + rfactors[8+pair[0]]*peaks[:,pair[1]+horn*8]) / \
                      (( sums[:,pair[0]+horn*8] + rfactors[8+pair[0]]* sums[:,pair[1]+horn*8])*pbeam/beams[horn])

            beamFactor = np.median(apRatio)/1.065

            inter = (sums[:,pair[0]+horn*8] + rfactors[8+pair[0]]*sums[:,pair[1]+horn*8]) / \
                    ( dios[:,pair[0]+horn*8] + rfactors[8+pair[0]]* dios[:,pair[1]+horn*8])# / beamFactor

            ratio = (peaks[:,pair[0]+horn*8] + rfactors[8+pair[0]]*peaks[:,pair[1]+horn*8]) / \
                    ( dios[:,pair[0]+horn*8] + rfactors[8+pair[0]]* dios[:,pair[1]+horn*8])# / beamFactor

            TDio = Fluxs/(ratio*Corfacts/beamFactor)/toJy(nu_c,beams[horn])#/beamFactor
            #TDio = TDiode[horn*4 + pair[0]]*ratio*toJy(freqs[horn*4 + pair[0]],beams[horn])

            gd = np.where((np.abs(TDio - np.median(TDio) ) < 0.19*(nu_c/13.)**2.7) & (np.abs(apRatio - np.median(apRatio) ) < 0.5) )[0]
            gd = np.where((np.abs(TDio - np.median(TDio) ) < 0.06) & (np.abs(apRatio - np.median(apRatio) ) < 0.5) )[0]             
            #gd = np.where(TDio > 0.)[0]
            #Random colours:
            c =  np.random.random(len(gd))

            #print  0.19*(nu_c/13.)**2.7

            print 'TDIO: ',horn,np.median(TDio),np.std(TDio[gd]),np.median(TDio[gd])*100.,nu_c, gd.size,beamFactor#rfactors[8+pair[0]]
            print np.median(np.sqrt(inter*pbeam/ratio/1.133)*180./np.pi)
            #pyplot.plot(np.abs(TDio[gd] - np.median(TDio) ),'o')
            #pyplot.show()

            '''
            pyplot.plot(mjd[gd],TDio[gd],'o',color='gray',label='Diode Temperature',zorder=0,markersize=8)
            pyplot.plot([np.min(mjd[gd]),np.max(mjd[gd])],[np.median(TDio[gd]),np.median(TDio[gd])],'--',lw=5,color='#F22000',label='Median')

            pyplot.legend(numpoints=1,frameon=False)
            
            pyplot.scatter(mjd[gd],TDio[gd],c=c,cmap='gray',s=120,lw=0,label='Diode Temperature',zorder=1)
            pyplot.plot([np.min(mjd[gd]),np.max(mjd[gd])],[np.median(TDio[gd]),np.median(TDio[gd])],'--',lw=5,color='#F22000',label='Median',zorder=2)
            pyplot.xlim(np.min(mjd[gd]),np.max(mjd[gd]))
            pyplot.xticks(size=14)
            pyplot.yticks(size=14)

            #cgd = np.where(np.abs(TDio[crabs] - np.median(TDio) ) < 0.19*(nu_c/13.)**2.7)[0]
            #pyplot.plot(mjd[crabs[cgd]],TDio[crabs[cgd]],'o',alpha=0.9,color='#0059F2',label='Tau A')
            pyplot.ylim(0,3)

            pyplot.title('Horn '+str(horn))
            pyplot.ylabel(r'T$_{\rm{dio}}$ (K)',size=14)
            pyplot.xlabel(r'MJD',size=14)
            pyplot.figure()

            pyplot.plot(np.arange(gd.size),apRatio[gd],'o',color='gray',label='Diode Temperature',zorder=0,markersize=8)
            pyplot.plot([0,gd.size],[np.median(apRatio[gd]),np.median(apRatio[gd])],'--',lw=5,color='#F22000',label='Median',zorder=0)

            pyplot.legend(numpoints=1,frameon=False)            
            pyplot.scatter(np.arange(gd.size),apRatio[gd],c=c,cmap='gray',s=120,lw=0,label='Fit-to-Photometry Ratio',zorder=1)
            pyplot.plot([0,gd.size],[np.median(apRatio[gd]),np.median(apRatio[gd])],'--',lw=5,color='#F22000',label='Median',zorder=2)
            pyplot.xlim(0,gd.size)
            pyplot.xticks(size=14)
            pyplot.yticks(size=14)
            pyplot.ylabel(r'$\frac{T_{fit}}{T_{phot} } \left(\frac{\Omega_{beam}}{\Omega_{pix}}\right)$',size=16)

            pyplot.show()
            '''
    #for ch in channels:
    #    print chan[ch], '&', '%2.3f' % np.median(dios[:,ch]), '\\\\'
