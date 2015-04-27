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
from MapMaker.Destriper import Control
from scipy.interpolate import interp1d

import scipy.optimize as op
import scipy.fftpack as sfft

import healpy as hp

import Masking

import Beams

def Rotate(c1,c2,a):
    Mat = np.array([[ np.cos(a), np.sin(a)],
                    [-np.sin(a), np.cos(a)]])

    vec = np.array([c1,c2])

    return Mat.dot(vec)

def RotateCoords(ra,dec,cra,cdec,cp):

    print cra*180./np.pi
    #Convert into cartisian coordinates
    x = np.cos(ra) * np.cos(dec)
    y = np.sin(ra) * np.cos(dec)
    z = np.sin(dec)

    #First Rotate by RA
    xr1,yr1 = Rotate(x,y,cra)
    zr1 = z
    #Then Rotate by Dec
    xr1,zr1 = Rotate(xr1,z,cdec)

    #Then Rotate by Para Ang
    yr1,zr1 = Rotate(yr1,zr1,-cp)
    

    new_ra  = np.arctan2(yr1,xr1)
    new_dec = np.arcsin(zr1)

    return new_ra,new_dec

if __name__ == "__main__":

    #Define data files (Change these to whatever you file structure is):
    datadir = 'CasA/'
    mfidir  = '/nas/scratch/sharper/QUIJOTE/MFI-Pipeline/MFI_INFO/'
    foreground_mask_dir = '/nas/scratch/sharper/QUIJOTE/Pipeline/MAPS/' #Foreground mask map directory

    dioModel = np.loadtxt('/nas/scratch/sharper/QUIJOTE/MFI-Pipeline/MFI_INFO/DiodePolModel-Epoch2.dat')

    #Read in list of files (Just needs to be part of name e.g. CASS would read in all file containing substring CASS)
    #files = np.loadtxt('FileList.txt',dtype='string',ndmin=1)
    files = np.loadtxt('Cass2Files.lis',dtype='string',ndmin=1)

    #Read in T-Point parameters:
    tpoints = np.loadtxt(mfidir+'focal_plane.ini')
    tpoints *= np.pi/180.
    jd0 = 56244.

    chanPairs = [0,1,2,3]#np.array([[1,7],[3,5]])#,[1,7],[3,5]])[0,6],[2,4]

    beams = np.loadtxt(mfidir+'BeamSizes.dat')
    beams     = [{'X':[beams[0,1]/2.355,True],'Y':[beams[0,1]/2.355,True]},
                 {'X':[beams[2,1]/2.355,True],'Y':[beams[2,1]/2.355,True]},
                 {'X':[beams[4,1]/2.355,True],'Y':[beams[4,1]/2.355,True]},
                 {'X':[beams[6,1]/2.355,True],'Y':[beams[6,1]/2.355,True]}]

    

    #Foreground mask:
    foremap = hp.read_map(foreground_mask_dir+'maskmap_N512.fits')
    foremap = (foremap == 0)

    
    Params = np.zeros(4)
    Sums = np.zeros(4)
    FullSums = np.zeros(4)


    cals = np.zeros(len(files))
    mdls = np.zeros(len(files))
    cals2 = np.zeros(len(files))
    mdls2 = np.zeros(len(files))
    mods = np.zeros(len(files))
    jds  = np.zeros(len(files))

    nside = 256
    sallmaps = np.zeros(12*nside**2)
    wallmaps = np.zeros(12*nside**2)


    import pyfits
    #Import the Beammaps
    beamMaps = {'217':pyfits.open('/nas/scratch/sharper/QUIJOTE/MFI-Pipeline/Beams/Horn2_17GHz_CTS_Beam.fits'),
                '219':pyfits.open('/nas/scratch/sharper/QUIJOTE/MFI-Pipeline/Beams/Horn2_19GHz_CTS_Beam.fits'),
                '311':pyfits.open('/nas/scratch/sharper/QUIJOTE/MFI-Pipeline/Beams/Horn3_11GHz_CTS_Beam.fits'),
                '313':pyfits.open('/nas/scratch/sharper/QUIJOTE/MFI-Pipeline/Beams/Horn3_13GHz_CTS_Beam.fits')}
    channels = [0,0,0,0,'219','217','219','217','313','311','313','311',0,0,0,0]
    
    for fi,f in enumerate(files):
        print 'File',f

        #Read in data:
        filelist = FileList(f,dir=datadir)
        print filelist
        data = ReadObs.FullObs(filelist,['DATA','AZ','EL','JD','MASK'],dir='')

        #Get some mean values:
        meanEl = np.median(data['EL'][0,:,0])

        jd = np.roll(data['JD'][0,:,0],1)-data['JD'][0,:,0]
        stepmax = np.argmax(np.abs(jd[1:])) + 1
        steptime = np.max(np.abs(jd[1:]))

        if steptime > 0.02:
            steps = [[0,stepmax],[stepmax,len(data['JD'][0,:,0])]]
        else:
            steps = [[0,len(data['JD'][0,:,0])]]
            
        for step in steps:
        
            meanJd = np.median(data['JD'][0,step[0]:step[1],0])


            sra,sdec = SourceFitting.ModelFlux.SourceCoord('CASA')
            #sra,sdec = SourceFitting.ModelFlux.EphemCoord('Jupiter',meanJd+jd0)
            #sra  *= 180./np.pi
            #sdec *= 180./np.pi


            c = ['#F29900','#0059F2','#F22000','#59F200']
            l = ['Horn 1','Horn 2','Horn 3','Horn 4']
            for horn in range(2,3): #Cannot calibrate horn 1 (0)


                #Get the RA/DEC coordinates:
                ra,dec,p = Coordinates.Hor2Sky(data['AZ'][0,step[0]:step[1],0]*np.pi/180.,
                                               data['EL'][0,step[0]:step[1],0]*np.pi/180.,
                                               data['JD'][0,step[0]:step[1],0]+jd0,TPoints={'xpos':tpoints[horn,1],
                                                                                            'ypos':tpoints[horn,2],
                                                                                            'Pf':  tpoints[horn,3],
                                                                                            'Px':  tpoints[horn,4],
                                                                                            'Py':  tpoints[horn,5],
                                                                                            'Pc':  tpoints[horn,6],
                                                                                            'Pn':  tpoints[horn,7],
                                                                                            'Pa':  tpoints[horn,8],
                                                                                            'Pb':  tpoints[horn,9]})
                x,y = RotateCoords(ra,dec,sra*np.pi/180.,sdec*np.pi/180.,np.mean(p))
                x *= 180./np.pi
                y *= 180./np.pi

                #Generate the beam model for this horn and data:
                
                xpos = 0.#83.6
                ypos = 0.#22.
                xo = np.mod(ra*180./np.pi,360.)
                yo = dec*180./np.pi
                
                print 'tt',x.shape,y.shape
                beam = beams[horn]
            
                isNotSource  = ((x)**2 + (y)**2 > (beam['X'][0]*1.25*2.355)**2)
                

                mask = data['MASK'][0,step[0]:step[1],(horn-1)*4].astype('bool')
                x = x[mask]
                y = y[mask]
                xo = xo[mask]
                yo = yo[mask]

                for chi,ch in enumerate(chanPairs):

                    mdlTOD = Beams.BeamModelTOD(beamMaps[channels[horn*4+ch]],x*np.pi/180.,y*np.pi/180.)

                    print 'CHANNEL:',ch
                    temps = data['DATA'][0,step[0]:step[1],ch+(horn-1)*4]
                    temps = temps[mask]


                    chsel = ch + (horn-1)*4.
                    P1 = SourceFitting.lm_SourceFit.FitSource(temps,
                                                              x,
                                                              y,
                                                              chsel,beam=beam,beamMap=mdlTOD)

                    if temps.size < 300:
                        P1 = None
                        
                    if (P1 != None):
                        print 'HELLo',P1

                        #pyplot.plot(temps,'-')
                        #pyplot.plot(temps-SourceFitting.lm_SourceFit.ModelFit(P1,x*np.pi/180.,y*np.pi/180.,mdlTOD),'-')
                        #pyplot.plot(SourceFitting.lm_SourceFit.ModelFit(P1,x*np.pi/180.,y*np.pi/180.,mdlTOD),'-')
                        #pyplot.show()
                        pbeam = 4.*np.pi/(12.*nside**2)
                        npix  = 12*nside**2 
                        pix   = hp.ang2pix(nside,(np.pi/2. - y*np.pi/180.),x*np.pi/180.)

                        Maps  = Control.Destriper(temps.astype('Float64'),temps.size,pix,npix,Medians=True)

                        ipix1 = hp.query_disc(nside,hp.ang2vec((90.-ypos)*np.pi/180.,xpos*np.pi/180.),   2.355*beam['X'][0]*np.pi/180.)
                        ipix2 = hp.query_disc(nside,hp.ang2vec((90.-ypos)*np.pi/180.,xpos*np.pi/180.),2.*2.355*beam['X'][0]*np.pi/180.)
                        ipix3 = hp.query_disc(nside,hp.ang2vec((90.-ypos)*np.pi/180.,xpos*np.pi/180.),3.*2.355*beam['X'][0]*np.pi/180.)

                        
                        mapmask = Maps.m*0.
                        mapmask[ipix2] +=1
                        mapmask[ipix3] +=1

                        bkgd = np.median(Maps.m[(mapmask == 1) & (Maps.m != 0)])
                        npixels = np.where((mapmask == 2) & (Maps.m != 0))[0].size
                        Sums[ch] = np.sum(Maps.m[ipix2]) -bkgd*npixels #- bkgd*ipix2.size#

                        cpix = hp.ang2pix(nside,(90.-sdec)*np.pi/180,sra*np.pi/180.)

                        #Maps.m[Maps.m == 0] = hp.UNSEEN
                        #Maps2.m[Maps2.m == 0] = hp.UNSEEN

                        #hp.gnomview(Maps.m,reso=2)
                        #hp.gnomview(Maps2.m,reso=2)                        
                        #print bkgd
                        #pyplot.show()
                        FullSums[ch] = bkgd*ipix2.size#(np.sum(Maps.m[ipix2])- bkgd*ipix1.size)/P1['amp']

                        #Maps.m[Maps.m == 0] = np.nan

                        if chsel == 4:
                            from SourceFitting.ModelFlux import toJy
                            if np.where(Maps.m[ipix2] != 0)[0].size > 140:
                                print 'FLUX DETERMINATION:',Sums[ch]*toJy(13.,pbeam),np.sum(Maps.m[ipix2])*toJy(13.,pbeam) , bkgd*toJy(13.,pbeam)*190.,np.where(Maps.m[ipix2] != 0)[0].size
                            sallmaps += Maps.sw#/P1['amp']
                            wallmaps += Maps.hw

                                       
                        Params[ch] = P1['amp'] 


                print Params[0],Params[1],Params[2],Params[3]
                if (type(P1) != type(None)):# & (Maps.m[cpix] != 0):
                    TextFiles.AppendFile('CasA_Nom_Fits_H3-All-FullBeam-Bkgd',np.concatenate((np.array([f]),Params,Sums,FullSums,np.array([npixels]),np.array([meanEl]),np.array([meanJd]))))

        print '--'
        #if np.where(Maps.m[ipix2] != 0)[0].size > 140:
        #    allmaps = sallmaps*0.
        #    allmaps[wallmaps != 0] = sallmaps[wallmaps != 0]/wallmaps[wallmaps != 0]
        #    allmaps[allmaps == 0.] = hp.UNSEEN
        #    hp.write_map('TAUA_MAP_STAMP.fits', allmaps)
        #    #print stop
