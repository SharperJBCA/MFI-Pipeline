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
import Coordinates

from scipy import optimize

def gradfit(P,X):

    return P[0]*X+ P[1]

def errfunc(P,X,Y):

    if (P[0] > -0.8) | (P[0] < -3.2):
        return 1e24 + np.zeros(len(X))
    else:
        return gradfit(P,X) - Y

def ReadBandpasses(mfidir):

    #Read in bandpasses:
    bandpasses = [np.loadtxt(mfidir+'Pol1_bandpass.dat',skiprows=1),
                  np.loadtxt(mfidir+'Pol2_bandpass.dat',skiprows=1),
                  np.loadtxt(mfidir+'Pol3_bandpass.dat',skiprows=1),
                  np.loadtxt(mfidir+'Pol4_bandpass.dat',skiprows=1)]

    nu_c = np.zeros(32)
    for i in range(len(bandpasses)):
        nu = bandpasses[i][:,0]
        for j in range(bandpasses[i].shape[1]-1):
            G = bandpasses[i][:,j+1]

            print 'FREQ: ', np.sum(nu*G)/np.sum(G)

            nu_c[i*8 + j] = np.sum(nu*G)/np.sum(G)

            print nu_c[i*8 + j],np.sum((nu/nu_c[i*8 + j])**(-0.3)*G)/np.sum(G*(nu/nu_c[i*8 + j])**(3.4) )
            
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

    return [freqs,band1],nu_c


if __name__ == '__main__':

    mfidir = '/nas/scratch/sharper/QUIJOTE/MFI-Pipeline/MFI_INFO/'

    #Define telescope constants:
    jd0 = 56244.
    nside = 128.

    beams = np.loadtxt(mfidir+'BeamSizes.dat')    
    beams = [1.133*(beams[0,1]*np.pi/180.)**2,
             1.133*(beams[1,1]*np.pi/180.)**2,
             1.133*(beams[2,1]*np.pi/180.)**2,
             1.133*(beams[3,1]*np.pi/180.)**2] #Beams of H2 to H4
    pbeam = 4.*np.pi/(12.*256.**2)
    pixbeam = 4.*np.pi/(12.*nside**2)#/beams[horn]

    freqs = np.loadtxt(mfidir+'FreqCens.dat')

    TDiode = np.loadtxt(mfidir+'TDiodes.dat')[:,1]

    rfactors = np.loadtxt(mfidir+'AvgRFactors.dat')

    channels = np.arange(4,dtype='i')
    chPairs = [[0,6],[1,7],[2,4],[3,5]]
    
    chan = ['H219c','H217c','H219u','H217u','H313c','H311c','H313u','H311u','H419c','H417c','H419u','H417u']

    #Read in bandpasses:
    bandpasses,nu_c = ReadBandpasses(mfidir)


    #d = np.loadtxt('CasA_GaussianFits_PerChannel_All',dtype='string')
    d = np.loadtxt('JupSeptember2013_Good',dtype='string')

    mjd   = d[:,-1].astype('f') + jd0
    peaks = d[:,1:33].astype('f')
    srcs  = np.array([c[0:4] for c in d[:,0]])
    dios  = d[:,1+32*2:33+32*2].astype('f')
    sums  = d[:,1+32:33+32].astype('f')

    means = np.zeros(8)
    errs = np.zeros(8)

    xjup = np.zeros(len(mjd))
    yjup = np.zeros(len(mjd))
    zjup = np.zeros(len(mjd))

    xear = np.zeros(len(mjd))
    year = np.zeros(len(mjd))
    zear = np.zeros(len(mjd))

    dist = np.zeros(len(mjd))
    for i in range(len(mjd)):
        xjup[i],yjup[i],zjup[i],a,b,c = Coordinates.Ephem.Planet(mjd[i],5)
        xear[i],year[i],zear[i],a,b,c = Coordinates.Ephem.Planet(mjd[i],3)

        dist[i] =  np.sqrt((xear[i]-xjup[i])**2 + (year[i]-yjup[i])**2 + (zear[i]-zjup[i])**2)

    #pyplot.plot(mjd,dist,'o')
    #pyplot.show()
    horns = range(1,3)
    for horn in horns:
        for pair in chPairs:


            ratio = (peaks[:,pair[0]+horn*8] + rfactors[8+pair[0]]*peaks[:,pair[1]+horn*8]) / \
                    ( dios[:,pair[0]+horn*8] + rfactors[8+pair[0]]* dios[:,pair[1]+horn*8])
            sumsr = (sums[:,pair[0]+horn*8] + rfactors[8+pair[0]]*sums[:,pair[1]+horn*8])*(pixbeam/beams[horn]) / \
                    ( dios[:,pair[0]+horn*8] + rfactors[8+pair[0]]* dios[:,pair[1]+horn*8])

            print np.median(ratio/sumsr), TDiode[horn*4 + pair[0]],horn

            r2 = TDiode[horn*4 + pair[0]]*ratio*toJy(freqs[horn*4 + pair[0]],beams[horn]) * (dist/4.04)**2#*pixbeam

            ratio = TDiode[horn*4 + pair[0]]*sumsr*1.065*toJy(freqs[horn*4 + pair[0]],beams[horn]) * (dist/4.04)**2#*pixbeam
            print np.median(ratio/r2), np.median(ratio),np.median(ratio),np.sqrt(beams[horn]/1.133)*180./np.pi


            #bd = np.where((ratio < 70.) & (ratio > 1.) & (np.isnan(ratio) == False))[0]# & (mjd < 56612) & (mjd > 56508))[0]
            #bd = np.where((ratio > 5.) & (ratio < 25.))[0]
            #bd = np.where((ratio > 20.) & (ratio < 60.)& (mjd < 56612))[0]
            bd = np.where((ratio > 10.))[0]

            #for ft in d[bd,0]:
            #    print ft
            if bd.size > 0:
                means[(horn - 1)*4 + pair[0]] = np.mean(ratio[bd])
                errs[(horn - 1)*4 + pair[0]]  = np.std(ratio[bd])/np.sqrt(bd.size)


                pfit = np.poly1d(np.polyfit(dist[bd],ratio[bd],2))
                pfit = np.poly1d(np.polyfit(np.log10(dist[bd]),np.log10(ratio[bd]),1))
                
                
                P0 = [-2.,np.log10(1e11/4./np.pi)]

                P1,s = optimize.leastsq(errfunc,P0,args=(np.log10(dist[bd]),np.log10(ratio[bd])))
            
                print pfit,P1,s
            
                plotdist = np.linspace(4,6.5,100)
                pyplot.plot(plotdist,10**gradfit(P1,np.log10(plotdist)),'--')                        
                pyplot.plot(dist[bd],ratio[bd],'o')

                pyplot.show()

                print np.median(ratio),np.std(ratio[bd])/np.sqrt(bd.size),np.std(ratio[bd])/np.mean(ratio[bd])*100./np.sqrt(bd.size)

    print means
    weightedMeans = (means[[0,1,4,5]]/errs[[0,1,4,5]]**2 + means[[2,3,6,7]]/errs[[2,3,6,7]]**2)/(1./errs[[0,1,4,5]]**2 + 1./errs[[2,3,6,7]]**2)
    weightedErrs  =  np.sqrt(1./(1./errs[[0,1,4,5]]**2 + 1./errs[[2,3,6,7]]**2))


    gd = np.where(np.isnan(weightedMeans) == False)[0]

    cents = np.array([8,9,16,17])
    print weightedMeans[gd],weightedErrs[gd]

    pyplot.plot(nu_c[cents[gd]],weightedMeans[gd],'o')
    pyplot.show()
