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

from MapMaker.Destriper import Control
import healpy as hp

if __name__ == "__main__":

    #Define data files (Change these to whatever you file structure is):
    datadir = '../../DownSampData/data/'
    mfidir  = '/nas/scratch/sharper/QUIJOTE/MFI-Pipeline/MFI_INFO/'
    foreground_mask_dir = '/nas/scratch/sharper/QUIJOTE/Pipeline/MAPS/' #Foreground mask map directory

    dioModel = np.loadtxt('/nas/scratch/sharper/QUIJOTE/MFI-Pipeline/MFI_INFO/DiodePolModel-Epoch2.dat')

    #Read in list of files (Just needs to be part of name e.g. CASS would read in all file containing substring CASS)
    #files = np.loadtxt('FileList.txt',dtype='string',ndmin=1)
    files = np.loadtxt('H2JUPFILES',dtype='string',ndmin=1)

    #Read in T-Point parameters:
    tpoints = np.loadtxt(mfidir+'focal_plane.ini')
    tpoints *= np.pi/180.
    jd0 = 56244.

    chanPairs = np.array([[0,6],[1,7],[2,4],[3,5]])

    beams = np.loadtxt(mfidir+'BeamSizes.dat')
    beams     = [{'X':[beams[0,1]/2.355,True],'Y':[beams[0,1]/2.355,True]},
                 {'X':[beams[1,1]/2.355,True],'Y':[beams[1,1]/2.355,True]},
                 {'X':[beams[2,1]/2.355,True],'Y':[beams[2,1]/2.355,True]},
                 {'X':[beams[3,1]/2.355,True],'Y':[beams[3,1]/2.355,True]}]
    
    #beams     = [{'X':[54.8/60./2.355,True],'Y':[54.8/60./2.355,True]},
    #             {'X':[0.65/2.355,True],'Y':[0.87*0.65/2.355,True]},
    #             {'X':[54.8/60./2.355,True],'Y':[54.8/60./2.355,True]},
    #             {'X':[0.64/2.355,True],'Y':[0.99*0.64/2.355,True]}]
    
    #beams     = [{'X':[54.8/60./2.355,True],'Y':[54.8/60./2.355,True]},
    #             {'X':[0.65/2.355,True],'Y':[0.87*0.65/2.355,True]},
    #             {'X':[0.92/2.355,True],'Y':[0.90/2.355,True]},
    #             {'X':[0.64/2.355,True],'Y':[0.99*0.64/2.355,True]}]

    #Noise limits:
    noiseLims = np.array([0.0764244215806,
                          0.0514327622598,
                          0.0889688890257,
                          0.054032925059,
                          0.00664899914188,
                          0.014294220537,
                          0.00656608744065,
                          0.0119769612796,
                          0.0599282051666,
                          0.124438793376,
                          0.0575965942848,
                          0.0426464567593])*2.

    

    #Foreground mask:
    foremap = hp.read_map(foreground_mask_dir+'maskmap_N512.fits')
    foremap = (foremap == 0)

    
    Params = np.zeros(32)
    rFacts = np.zeros(32)
    Sums = np.zeros(32)


    cals = np.zeros(len(files))
    mdls = np.zeros(len(files))
    cals2 = np.zeros(len(files))
    mdls2 = np.zeros(len(files))
    mods = np.zeros(len(files))
    jds  = np.zeros(len(files))

    sallmaps = np.zeros(12*512**2)
    wallmaps = np.zeros(12*512**2)

    for fi,f in enumerate(files):
        print 'File',f

        #Read in data:
        filelist = FileList(f,dir=datadir)
        data = ReadObs.FullObs(filelist,['DATA','AZ','EL','JD','CAL','MOD','CAL_JD','ERR','MASK'],dir='')
        #data = ReadObs.FullObs(filelist,['DATA','AZ','EL','JD','CAL','MOD','CAL_JD','ERR'],dir='')

        
        #Get some mean values:
        meanEl = np.median(data['EL'][0,:,0])
        meanJd = np.median(data['JD'][0,:,0])
        mask = np.where((data['MASK'][0,:,0] == 0))[0]

        ajd = data['JD'][0,mask,0]
        #mask = np.ones(data['EL'][0,:,0].size,dtype='bool')


        #Find the AZ scan starts/ends
        seconds = 24.*60.**2        
        P1 = WaveFitter.FindWave(data['JD'][0,mask,0]*seconds,data['AZ'][0,mask,0])
        scanLength = P1[3]/seconds/2.

        phaseStart = P1[2]*P1[3]/seconds*2
        nScans = int( ((np.max(data['JD'][0,mask,0])-np.min(data['JD'][0,mask,0])) - phaseStart) /scanLength )
        
        X = data['JD'][0,mask,0]-np.min(data['JD'][0,mask,0])    
        scanStarts = np.zeros(nScans+2)
        maskHits  = np.ones(nScans+2,dtype='bool')
        for i in range(1,nScans+1):
            hits = np.where((X > phaseStart + (i-1)*scanLength) & (X < phaseStart + (i)*scanLength))[0]
            if hits.size > 0:
                scanStarts[i]= hits[0]
            else:
                maskHits[i] = False

        scanStarts = scanStarts[maskHits]
        scanStarts[-1] = len(X)


        #sra,sdec = SourceFitting.ModelFlux.SourceCoord('CASA')
        sra,sdec = SourceFitting.ModelFlux.EphemCoord('Jupiter',meanJd+jd0)
        sra *= 180./np.pi
        sdec *= 180./np.pi
        
        c = ['#F29900','#0059F2','#F22000','#59F200']
        l = ['Horn 1','Horn 2','Horn 3','Horn 4']
        for horn in range(1,2): #Cannot calibrate horn 1 (0)


            #Get the RA/DEC coordinates:
            ra,dec,p = Coordinates.Hor2Sky(data['AZ'][0,mask,0]*np.pi/180.,
                                           data['EL'][0,mask,0]*np.pi/180.,
                                           data['JD'][0,mask,0]+jd0,TPoints={'xpos':tpoints[horn,1],
                                                                             'ypos':tpoints[horn,2],
                                                                             'Pf':  tpoints[horn,3],
                                                                             'Px':  tpoints[horn,4],
                                                                             'Py':  tpoints[horn,5],
                                                                             'Pc':  tpoints[horn,6],
                                                                             'Pn':  tpoints[horn,7],
                                                                             'Pa':  tpoints[horn,8],
                                                                             'Pb':  tpoints[horn,9]})


            #Convert RA/DEC to telescope frame:
            #xi = (np.mod(ra*180./np.pi+180.,360.) - (sra-180.))*np.cos(dec) #!!! Remove the 180 if not fitting for Cas A!!!
            xi = (np.mod(ra*180./np.pi,360.) - (sra))*np.cos(dec) #!!! Remove the 180 if not fitting for Cas A!!!

            yi = (dec*180./np.pi-sdec) #!!! Remove the 180 if not fitting for Cas A!!!

            #Rotate sky coordinates to telescope frame
            x =  xi*np.cos(p) + yi*np.sin(p) 
            y = -xi*np.sin(p) + yi*np.cos(p)
                    
            beam = beams[horn]
            
            isNotSource  = (x**2 + y**2 > (beam['X'][0]*1.25*2.355)**2)
            
            stds = np.zeros(len(scanStarts))

            meanMod = np.mean(np.mod(np.median(data['MOD'][0,mask,horn]),90)*4.)

            xpos = np.array([])
            ypos = np.array([])
    
            for ch,chPair in enumerate(chanPairs):

                for chP in chPair:

                    
                    volts = data['DATA'][0,mask,chP+horn*8]
                    volts_bad = data['DATA'][0,mask,chP+horn*8]*1.
                    chErr = data['ERR'][0,mask,chP+horn*8]

                    isNotSource  = (x**2 + y**2 > (beam['X'][0]*2.)**2)

                    
                    for i in range(len(scanStarts)-1):
                        dat = volts[scanStarts[i]:scanStarts[i+1]]*1.
                        jd = ajd[scanStarts[i]:scanStarts[i+1]]
                        xi = x[scanStarts[i]:scanStarts[i+1]]
                        yi = y[scanStarts[i]:scanStarts[i+1]]

                        isNotSource  = (xi**2 + yi**2 > (beam['X'][0]*2.355)**2)

                        #pyplot.plot(xi,yi,'.b',alpha=0.5)
                        #pyplot.plot(xi[isNotSource],yi[isNotSource],'.r',alpha=0.5)

                        if len(xi[isNotSource]) > 4:
                            bkgdfit = np.poly1d(np.polyfit(xi[isNotSource],dat[isNotSource],3))
                            lastfit = bkgdfit
                        else:
                            bkgdfit = lastfit
                        volts[scanStarts[i]:scanStarts[i+1]] -=bkgdfit(x[scanStarts[i]:scanStarts[i+1]])
    
                        #volts[scanStarts[i]:scanStarts[i+1]] -= 1.
                        
                        stds[i] = np.std(dat[isNotSource])
                    #pyplot.show()

                    
                    if (np.sqrt(np.median(stds**2)) < noiseLims[ch+(horn-1)*4]):
                        P1 = SourceFitting.lm_SourceFit.FitSource(volts,x,y,beam=beam,err=chErr)

                        if (P1 != None):
                            #pyplot.plot(volts,',')
                            ##pyplot.plot(SourceFitting.lm_SourceFit.ModelFit(P1,x,y),'-')
                            #pyplot.show()
                            nside = 128
                            npix  = 12*nside**2
                            pix   = hp.ang2pix(nside,(np.pi/2. - y*np.pi/180.),x*np.pi/180.)
                            #pix   = hp.ang2pix(nside,(np.pi/2. - dec),ra)

                            gd = np.where(volts < 100.001)[0]
                            Maps  = Control.Destriper(volts[gd].astype('Float64'),gd.size,pix[gd],npix,Medians=True)

                            ipix = hp.query_disc(nside,hp.ang2vec((90.-0.)*np.pi/180.,0.*np.pi/180.),2.355*beam['X'][0]*np.pi/180.)
                            bpix = hp.query_disc(nside,hp.ang2vec((90.-(0.-1.))*np.pi/180.,(0.-1.)*np.pi/180.),0.5*np.pi/180.)

                            Sums[chP+horn*8] = np.sum(Maps.m[ipix])# - float(ipix.size)*np.median((Maps.m[bpix])[Maps.m[bpix] != 0.])

                            #hp.gnomview(Maps.m)
                            #pyplot.show()
                            #Maps.m[ipix] = 0.
                            #amaxp = np.argmax(Maps.m)
                            #c2,c1 = hp.pix2ang(nside,amaxp)
                            #c2 = (np.pi/2. - c2)*180./np.pi
                            #c1 = c1 * 180./np.pi

                            argmax = np.argmax(SourceFitting.lm_SourceFit.ModelFit(P1,x,y))
                            maxjd = data['JD'][0,argmax,0]

                            print maxjd

                            gd = np.where((data['CAL'][0,:,ch] != 0) & (np.abs(data['CAL_JD'][0,:,0]-maxjd) < 1./60./24.  )  )[0]
                            #r-factor
                            r = np.median(data['CAL'][0,gd,chP+horn*8])

                            Params[chP+horn*8] = P1['amp']
                            rFacts[chP+horn*8] = r                             
                            print 'P AMP',P1['amp'],r,chP+horn*8


                    if type(P1) != type(None):
                        xpos = np.append(xpos,P1['cx'])
                        ypos = np.append(ypos,P1['cy'])
                        
                        wxpos = np.append(xpos,P1['wx'])
                        wypos = np.append(ypos,P1['wy'])

                        if chP+horn*8 == 17:
                            bpix = hp.query_disc(nside,hp.ang2vec((90.-(0.))*np.pi/180.,0.*np.pi/180.),1.*np.pi/180.)


                            #hp.gnomview(Maps.m)
                            #pyplot.show()

                            bdpix = np.where(np.abs(Maps.m) > 5e-4)[0]
                            Maps.m[bpix] = 0.
                            print 'NOISE!!!: ',np.std(Maps.m[Maps.m != 0])

                            #if np.std(Maps.m[Maps.m != 0]) < 5e-4:
                            Maps.sw[bdpix] = 0.
                            Maps.hw[bdpix] = 0.

                            sallmaps[Maps.hw != 0] += Maps.sw[Maps.hw != 0]
                            wallmaps[Maps.hw != 0] += Maps.hw[Maps.hw != 0]

                            #allmaps = sallmaps*0.
                            #allmaps[wallmaps != 0] = sallmaps[wallmaps != 0]/wallmaps[wallmaps != 0]
                            #allmaps[allmaps == 0.] = hp.UNSEEN


        print Params[0],Params[1],Params[2],Params[3],Params[4],Params[5],Params[6],Params[7]
        print Sums[0],Sums[1],Sums[2],Sums[3],Sums[4],Sums[5],Sums[6],Sums[7]            
        print rFacts[0],rFacts[1],rFacts[2],rFacts[3],rFacts[4],rFacts[5],rFacts[6],rFacts[7]
        print 'POSITION:',np.nanmean(xpos),np.nanmean(ypos)
        #if type(P1) != type(None):
        #    TextFiles.AppendFile('H2JupSeptember2013_Good',np.concatenate((np.array([f]),Params,Sums,rFacts,np.array([meanEl]),np.array([meanJd]))))#,np.mean(xpos),np.mean(ypos),np.mean(wxpos),np.mean(wypos))))

        print '--'

    allmaps = sallmaps*0.
    allmaps[wallmaps != 0] = sallmaps[wallmaps != 0]/wallmaps[wallmaps != 0]
    allmaps[allmaps == 0.] = hp.UNSEEN
    hp.write_map('Jupiter_H311c_Healpix.fits', allmaps)
