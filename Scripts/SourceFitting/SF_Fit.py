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

import healpy as hp

if __name__ == "__main__":

    #Define data files (Change these to whatever you file structure is):
    datadir = '../DownSampData/data/'
    mfidir  = '/nas/scratch/sharper/QUIJOTE/MFI-Pipeline/MFI_INFO/'
    foreground_mask_dir = '/nas/scratch/sharper/QUIJOTE/Pipeline/MAPS/' #Foreground mask map directory

    dioModel = np.loadtxt('/nas/scratch/sharper/QUIJOTE/MFI-Pipeline/MFI_INFO/DiodePolModel-Epoch2.dat')

    #Read in list of files (Just needs to be part of name e.g. CASS would read in all file containing substring CASS)
    #files = np.loadtxt('FileList.txt',dtype='string',ndmin=1)
    files = np.loadtxt('CASSFILES',dtype='string',ndmin=1)

    #Read in T-Point parameters:
    tpoints = np.loadtxt(mfidir+'focal_plane.ini')
    tpoints *= np.pi/180.
    jd0 = 56244.

    chanPairs = np.array([[0,6],[1,7],[2,4],[3,5]])
    beams     = [{'X':[0.89/2.355,False],'Y':[0.99*0.89/2.355,False]},
                 {'X':[0.65/2.355,False],'Y':[0.87*0.65/2.355,False]},
                 {'X':[0.85/2.355,False],'Y':[0.92*0.85/2.355,False]},
                 {'X':[0.64/2.355,False],'Y':[0.99*0.64/2.355,False]}]

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

    
    Params = np.zeros(12)
    rFacts = np.zeros(12)

    for f in files:
        print 'File',f

        #Read in data:
        filelist = FileList(f,dir=datadir)
        data = ReadObs.FullObs(filelist,['DATA','AZ','EL','JD','CAL','MOD','CAL_JD','ERR'],dir='')

        #Get some mean values:
        meanEl = np.median(data['EL'][0,:,0])
        meanJd = np.median(data['JD'][0,:,0])

        for i in range(data['DATA'].shape[2]):
            data['CAL'][0,:,i] /= (1. + dioModel[i,1]*np.cos( np.pi/180.*np.mod(np.median(data['MOD'][0,:,i/8]),90)*4. + dioModel[i,2]*2.*np.pi))



        #Find the AZ scan starts/ends
        seconds = 24.*60.**2        
        P1 = WaveFitter.FindWave(data['JD'][0,:,0]*seconds,data['AZ'][0,:,0])
        scanLength = P1[3]/seconds/2.
        #nScans = (Phase Start time) -> Range(jd) / scanLength
        phaseStart = P1[2]*P1[3]/seconds*2
        nScans = int( ((np.max(data['JD'][0,:,0])-np.min(data['JD'][0,:,0])) - phaseStart) /scanLength )
        
        X = data['JD'][0,:,0]-np.min(data['JD'][0,:,0])    
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


        for horn in range(2,3): #Cannot calibrate horn 1 (0)
            for ch,chPair in enumerate(chanPairs):

                gd = np.where(data['CAL'][0,:,ch] != 0)[0]
                #r-factor
                r = np.median(data['CAL'][0,gd,chPair[0]+horn*8])/np.median(data['CAL'][0,gd,chPair[1]+horn*8])
                
                steplen = 1440
                nsteps = gd.size/steplen
                dioVolt = data['CAL'][0,gd,chPair[0]+horn*8] + r*data['CAL'][0,gd,chPair[1]+horn*8]

                if nsteps > 0:
                    newlen = nsteps * steplen
                    #Produce a model of the gain drifts with JD
                    bca    = np.median(np.reshape(data['CAL'][0,gd[:newlen],chPair[0]+horn*8]+ \
                                                  r*data['CAL'][0,gd[:newlen],chPair[1]+horn*8],(nsteps,steplen)),axis=1)
        
                    bjd    = np.median(np.reshape(data['CAL_JD'][0,gd[:newlen],0],(nsteps,steplen)),axis=1)
                else:
                    bca = np.array([np.median(data['CAL'][0,gd,chPair[0]+horn*8]+ \
                                              r*data['CAL'][0,gd,chPair[1]+horn*8])])
                    
                if bca.size > 4:
                    bjd[0]  = np.min([np.min(data['JD'][0,:,0]),np.min(data['CAL_JD'][0,gd,0])])
                    bjd[-1] = np.max([np.max(data['JD'][0,:,0]),np.max(data['CAL_JD'][0,gd,0])])        
                    pfit    = ip.interp1d(bjd,bca)
                    pfit = pfit(data['JD'][0,:,0])
                else:
                    pfit = np.median(bca) + np.zeros(data['DATA'][0,:,chPair[0]+horn*8].size)
                    

                chVoltage = (data['DATA'][0,:,chPair[0]+horn*8] + r*data['DATA'][0,:,chPair[1]+horn*8])/pfit


                ratio = np.median(data['DATA'][0,:,chPair[0]+horn*8]/(data['DATA'][0,:,chPair[1]+horn*8]))

                pyplot.plot(      data['DATA'][0,:,chPair[0]+horn*8]-np.median(      data['DATA'][0,:,chPair[0]+horn*8]),',')
                pyplot.plot(ratio*data['DATA'][0,:,chPair[1]+horn*8]-np.median(ratio*data['DATA'][0,:,chPair[1]+horn*8]),',')

                pyplot.figure()
                pyplot.plot(data['DATA'][0,:,chPair[0]+horn*8]-(ratio*data['DATA'][0,:,chPair[1]+horn*8]),',')
                print ratio/r
                print r
                pyplot.show()
                
                chErr = data['ERR'][0,:,chPair[0]+horn*8]

                #Get the RA/DEC coordinates:
                ra,dec,p = Coordinates.Hor2Sky(data['AZ'][0,:,0]*np.pi/180.,
                                               data['EL'][0,:,0]*np.pi/180.,
                                               data['JD'][0,:,0]+jd0,TPoints={'xpos':tpoints[(chPair[0]+horn*8)/8,1],
                                                                              'ypos':tpoints[(chPair[0]+horn*8)/8,2],
                                                                              'Pf':tpoints[(chPair[0]+horn*8)/8,3],
                                                                              'Px':tpoints[(chPair[0]+horn*8)/8,4],
                                                                              'Py':tpoints[(chPair[0]+horn*8)/8,5],
                                                                              'Pc':tpoints[(chPair[0]+horn*8)/8,6],
                                                                              'Pn':tpoints[(chPair[0]+horn*8)/8,7],
                                                                              'Pa':tpoints[(chPair[0]+horn*8)/8,8],
                                                                              'Pb':tpoints[(chPair[0]+horn*8)/8,9]})

                #Convert RA/DEC to telescope frame:
                sra,sdec = SourceFitting.ModelFlux.SourceCoord('CASA')
                xi = (np.mod(ra*180./np.pi+180.,360.) - (sra-180.))*np.cos(dec) #!!! Remove the 180 if not fitting for Cas A!!!
                #xi = (np.mod(ra*180./np.pi,360.) - (sra))*np.cos(dec) #!!! Remove the 180 if not fitting for Cas A!!!

                yi = (dec*180./np.pi-sdec) #!!! Remove the 180 if not fitting for Cas A!!!

                #Rotate sky coordinates to telescope frame
                x = xi*np.cos(p) - yi*np.sin(p) 
                y = xi*np.sin(p) + yi*np.cos(p)

                beam = beams[horn]
                
                isNotSource  = (x**2 + y**2 > (beam['X'][0]*1.25*2.355)**2)

                stds = np.zeros(len(scanStarts))
                
                for i in range(len(scanStarts)-1):
                    dat = chVoltage[scanStarts[i]:scanStarts[i+1]]*1.
                    jd = data['JD'][0,scanStarts[i]:scanStarts[i+1],0]
                    xi = x[scanStarts[i]:scanStarts[i+1]]
                    yi = y[scanStarts[i]:scanStarts[i+1]]

                    isNotSource  = (xi**2 + yi**2 > (beam['X'][0]*1.25*2.355)**2)

                    bkgdfit = np.poly1d(np.polyfit(xi[isNotSource],dat[isNotSource],3))
                
                    chVoltage[scanStarts[i]:scanStarts[i+1]] -=bkgdfit(xi)

                    stds[i] = np.std(dat[isNotSource])

                #print 'NOISE:',np.sqrt(np.median(stds**2))
                #pyplot.plot(stds,'o')
                #pyplot.show()

                if (np.sqrt(np.median(stds**2)) < noiseLims[ch+(horn-1)*4]):
                    P = SourceFitting.lm_SourceFit.FitSource(chVoltage,x,y,beam=beam,err=chErr)

                    if (P != None):
                        Params[ch+(horn-1)*4] = P['amp']
                        rFacts[ch+(horn-1)*4] = r                             
                        print 'P AMP',P['amp']
                else:
                    Params[ch+(horn-1)*4] = 0.
                    rFacts[ch+(horn-1)*4] = 0.

        '''
        TextFiles.AppendFile('CasA',[f,
                                         Params[0],Params[1],Params[2],Params[3],
                                         Params[4],Params[5],Params[6],Params[7],
                                         Params[8],Params[9],Params[10],Params[11],
                                         rFacts[0],rFacts[1],rFacts[2],rFacts[3],
                                         rFacts[4],rFacts[5],rFacts[6],rFacts[7],
                                         rFacts[8],rFacts[9],rFacts[10],rFacts[11],
                                         meanEl,meanJd,np.mean(p)])
        '''
        print '--'

