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

import MyMaths

if __name__ == "__main__":

    #Define data files (Change these to whatever you file structure is):
    datadir = '/nas/scratch/sharper/QUIJOTE/datos/caltod/'#'../../DownSampData/data/'
    datadir = '/nas/scratch/sharper/QUIJOTE/Observations/SimTele/TELE/SIMTOD/'#'../../DownSampData/data/'    
    mfidir  = '/nas/scratch/sharper/QUIJOTE/MFI-Pipeline/MFI_INFO/'
    foreground_mask_dir = '/nas/scratch/sharper/QUIJOTE/Pipeline/MAPS/' #Foreground mask map directory

    dioModel = np.loadtxt('/nas/scratch/sharper/QUIJOTE/MFI-Pipeline/MFI_INFO/DiodePolModel-Epoch2.dat')

    #Read in list of files (Just needs to be part of name e.g. CASS would read in all file containing substring CASS)
    #files = np.loadtxt('FileList.txt',dtype='string',ndmin=1)
    files = np.loadtxt('Jup.lis',dtype='string',ndmin=1)
    print files.shape
    files = np.array(['SIM_TAUA_HORN3_NONOISE_0-85x1-5.fits'])
    print files.shape
    files = np.array(['SIM_TAUA_HORN3_NONOISE.fits'])
    files = np.array(['DensePack.fits'])
    files = np.array(['VeryDensePack.fits'])

    beams = np.loadtxt(mfidir+'BeamSizes.dat')
    beam ={'X':[beams[2,1]/2.355,True],'Y':[beams[2,1]/2.355,True]}
    beam ={'X':[0.85/2.355,True],'Y':[1.5/2.355,True]}


    
    Params = np.zeros(4)
    rFacts = np.zeros(4)
    Sums = np.zeros(4)

    sallmaps = np.zeros(12*512**2)
    wallmaps = np.zeros(12*512**2)

    #for fi,f in enumerate(files):

    widths = np.array([2.,1.75,1.5,1.25,1.,0.75])#,0.5,0.25
    widths = np.linspace(0.25,2.1,10.)
    nloops = len(widths)

    meanw = np.zeros(nloops)

    gauss = lambda P,x,y: P[0]*np.exp(-0.5*(((x-P[1])/P[3])**2 + ((y-P[2])/P[4])**2) )
    
    for fi in range(nloops):#enumerate(files):

        fwhmconst = (2.*np.sqrt(2.*np.log(2.)))
        beam ={'X':[widths[fi]/fwhmconst,True],'Y':[widths[fi]/fwhmconst,True]}

        f = files[0]
        print 'File',f

        #Read in data:
        filelist = FileList(f,dir=datadir)
        data = ReadObs.FullObs(filelist,['DATA','RA','DEC'],dir='')

        data['RA'][0,:,0] -= 83.
        data['DEC'][0,:,0]-= 22.
        sra,sdec = 0.,0.

        data['DATA'][0,:,0] = gauss([1.,sra,sdec,widths[fi]/fwhmconst,widths[fi]/fwhmconst],data['RA'][0,:,0],data['DEC'][0,:,0])

        data['DATA'][0,:,0] = data['DATA'][0,:,0] * 300. / MyMaths.toJy(12.85,1.133*(widths[fi]*np.pi/180.)**2)

        nside2 = 512
        pix   = hp.ang2pix(nside2,(np.pi/2. - data['DEC'][0,:,0]*np.pi/180.),data['RA'][0,:,0]*np.pi/180.)
        npix  = 12*nside2**2
        temps = data['DATA'][0,:,0] #+ np.random.normal(loc=0.,scale=1e-8,size=len(data['DATA'][0,:,0]))
        Maps  = Control.Destriper(temps,temps.size,pix,npix,Medians=True,cn=np.ones(temps.size))
        #data['DATA'][0,:,0] = Maps.m[pix]

        
        

        #data['DATA'][0,:,0] -= np.min(data['DATA'][0,:,0])


        ra,dec = data['RA'][0,:,0],data['DEC'][0,:,0]

        #Convert RA/DEC to telescope frame:
        xi = (ra- sra)*np.cos(dec*np.pi/180.) #!!! Remove the 180 if not fitting for Cas A!!!
        yi = (dec-sdec) #!!! Remove the 180 if not fitting for Cas A!!!

        #Rotate sky coordinates to telescope frame
        x =  xi
        y =  yi

        #pyplot.plot(x,y,',')
        #pyplot.show()

        #Estimate widths:
        #xc = np.sum(x*data['DATA'][0,:,0])/np.sum(data['DATA'][0,:,0])
        lims = np.linspace(-20,20,2000)
        r = np.max(lims)-np.min(lims)
        xs = gauss([1.,sra,sdec,widths[fi]/fwhmconst,widths[fi]/fwhmconst],lims,np.zeros(2000))
        xc = np.sum(xs*lims)/np.sum(xs)

        width_X = np.sum(np.sqrt(np.abs(lims)**2*xs))/np.sum(xs)*1.476
        print 'WIDTH ESTIMATE',width_X, widths[fi],widths[fi]/width_X
                
        isNotSource  = (x**2 + y**2 > (beam['X'][0]*1.25*2.355)**2)

        temps = data['DATA'][0,:,0]#/MyMaths.toJy(12.85,0.85*1.5*(np.pi/180.)**2)
        #pyplot.plot(temps)
        #pyplot.show()
            
        P1 = SourceFitting.lm_SourceFit.FitSource(temps,x,y,beam=beam)#,err=chErr)
        if (P1 != None):
            nside = 512
            npix  = 12*nside**2
            pix   = hp.ang2pix(nside,(np.pi/2. - y*np.pi/180.),x*np.pi/180.)

            Maps  = Control.Destriper(temps.astype('Float64'),temps.size,pix,npix,Medians=True,cn=np.ones(temps.size))

            ipix = hp.query_disc(nside,hp.ang2vec((90.-0.)*np.pi/180.,0.*np.pi/180.),1.*2.355*beam['X'][0]*np.pi/180.)
            bpix = hp.query_disc(nside,hp.ang2vec((90.-(0.-1.))*np.pi/180.,(0.-1.)*np.pi/180.),0.5*np.pi/180.)
            
            Sums = np.sum(Maps.m[ipix])# - float(ipix.size)*np.median((Maps.m[bpix])[Maps.m[bpix] != 0.])

            #pyplot.plot(x,temps)
            #pyplot.plot(x,SourceFitting.lm_SourceFit.ModelFit(P1,x,y))
            #pyplot.show()

            #hp.gnomview(Maps.m,norm='hist',reso=3)
            #pyplot.show()

            #argmax = np.argmax(SourceFitting.lm_SourceFit.ModelFit(P1,x,y))
            #maxjd = data['JD'][0,argmax,0]

            fit = P1['amp']
            #aper = Sums * (4.*np.pi/(12.*float(nside)**2))/(0.85*np.pi/180.)**2

            
            notFittedWidths = fit*MyMaths.toJy(12.85,1.133*widths[fi]**2*(np.pi/180.)**2)
            fittedWidths    = fit*MyMaths.toJy(12.85,1.133*((P1['wx']+P1['wy'])/2.*2.355 * np.pi/180.)**2)

            meanwidth = (P1['wx']+P1['wy'])/2.*2.355
            
            #Params[ch] = P1['amp']
            print meanwidth, widths[fi]
            print 'P AMP', notFittedWidths
            print 'P AMP', fittedWidths    

            print 'Ratio', fittedWidths/notFittedWidths        

            meanw[fi] = fittedWidths/notFittedWidths       #(P1['wx']+P1['wy'])/2.*2.355

            print MyMaths.toJy(12.85,np.sqrt( (widths[fi]*np.pi/180.)**2 - (4.*np.pi/(12.*float(nside)**2)) )  )
            print 'P AMP',Sums*MyMaths.toJy(12.85,(4.*np.pi/(12.*float(nside)**2)))#*(meanwidth*np.pi/180.)**2/(widths[fi]*np.pi/180.)**2

        #if type(P1) != type(None):
        #    TextFiles.AppendFile('Jupiter_Cal',np.concatenate((np.array([f]),Params,Sums,np.array([meanEl]),np.array([meanJd]))))#,np.mean(xpos),np.mean(ypos),np.mean(wxpos),np.mean(wypos))))

        print '--'

    pyplot.plot(widths,meanw,'o')
    pyplot.show()
    #allmaps = sallmaps*0.
    #allmaps[wallmaps != 0] = sallmaps[wallmaps != 0]/wallmaps[wallmaps != 0]
    #allmaps[allmaps == 0.] = hp.UNSEEN
    #hp.write_map('Jupiter_H311c_Healpix.fits', allmaps)
