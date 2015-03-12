#CalData.py
# Goals:
# 1) Fit a polynomial to G_ch*T_diode voltages
# 2) Divide V_ch by G_ch*T_diode
# 3) Multiply by 1 + P sin(4x + O) diode model
# 4) Record amplitude of source

import DataAccess
import numpy as np
from matplotlib import pyplot
import scipy.interpolate as ip
import Binning

import scipy.signal as sig
   
import WaveFitter
import Coordinates
import healpy as hp
import SourceFitting

def MedianFilter(data,jd,jdlen,mask,foremask,cutoff=0.03):
    '''
    '''


    
    RFIMask = np.ones(mask.size,dtype='bool')

    #Subtract off median of rings:
    ijd = jd- np.min(jd)
    
    jdmax = np.max(ijd)
    nrings = int(jdmax)/int(jdlen)

    stds = np.zeros(nrings)
    maxvals = np.zeros(nrings)

    for i in range(nrings):
        if i == nrings-1:
            hi = jdmax+1.
        else:
            hi = (i+1)*jdlen
            
        thisRing = ((ijd >= i*jdlen) & (ijd < hi))
        thisBkgd = data[(thisRing & foremask & mask)]

        if thisBkgd.size > 0:
            data[thisRing] -= np.median(thisBkgd)
            maxvals[i] = np.max(np.abs(thisBkgd-np.median(thisBkgd)))
            stds[i] = np.std(thisBkgd)
        if (stds[i] > cutoff) | (stds[i] == 0) | (maxvals[i] > 0.2):
            RFIMask[thisRing] = False

    #pyplot.plot(maxvals,'o')
    #pyplot.figure()
    #pyplot.plot(data,',')    
    #pyplot.plot(jd[(mask & RFIMask)],data[(mask & RFIMask)],',')
    #pyplot.show()

    return RFIMask
    
if __name__ == '__main__':

    #Data access parameters:
    mfidir  = '/nas/scratch/sharper/QUIJOTE/MFI-Pipeline/MFI_INFO/'    
    dioModel = np.loadtxt(mfidir+'DiodePolModel-Epoch1.dat')

    files = np.loadtxt('Filelist.txt',dtype='string')
    dir = '/nas/scratch/sharper/QUIJOTE/datos/caltod/'

    cols = ['DATA','MASK','JD','EL','AZ']

    channels = np.array([0,1,2,3])

    #Read in T-Point parameters:
    tpoints = np.loadtxt(mfidir+'focal_plane.ini')
    tpoints *= np.pi/180.    

    #Constants and extras:
    nHorns = 4
    nside = 512
    maskmap = np.ones(12*nside**2,dtype='bool')
    jd0 = 56244.

    set = 0.8
    rise= 0.25
    
    #Foreground mask:
    foremap = hp.read_map('/nas/scratch/sharper/QUIJOTE/Pipeline/MAPS/maskmap_N512.fits')
    foremap = (foremap == 0)

    #Output columns:
    calcols = ['DATA','MASK','JD','EL','AZ']
    caldir = '/nas/scratch/sharper/QUIJOTE/datos/caltod_noground/'


    #Open all NOMINALXX files:
    filelist = DataAccess.FileList('NOMINAL60?-13',dir=dir)
    #filelist = filelist[0:len(filelist)/4]
    #print filelist

    #m = np.ones(len(filelist),dtype='bool')
    #m[24:28] = False
    #filelist = filelist[(m)]
    
    data = DataAccess.FullObs(filelist,cols)
        
    #for horn in range(2,3):
    #Get the RA/DEC coordinates:
    horn = 2
    ra,dec,p = Coordinates.Hor2Sky(data['AZ'][0,:,0]*np.pi/180.,
                                   data['EL'][0,:,0]*np.pi/180.,
                                   data['JD'][0,:,0]+jd0,TPoints={'xpos':tpoints[horn,1],
                                                                  'ypos':tpoints[horn,2],
                                                                  'Pf':  tpoints[horn,3],
                                                                  'Px':  tpoints[horn,4],
                                                                  'Py':  tpoints[horn,5],
                                                                  'Pc':  tpoints[horn,6],
                                                                  'Pn':  tpoints[horn,7],
                                                                  'Pa':  tpoints[horn,8],
                                                                  'Pb':  tpoints[horn,9]})
        
    pix = hp.ang2pix(nside,np.pi/2.-dec,ra)
        
    #Mask any bright Galactic sources:
    foremask = foremap[pix]

    #Sort AZ
    saz = np.argsort(data['AZ'][0,:,0])
    az  = data['AZ'][0,saz,0]
    jd  = np.mod(data['JD'][0,saz,0],1.)
    tod = data['DATA'][0,saz,:]
    mask = (data['MASK'][0,saz,:].astype('bool') )

    del data

    day   = np.where((jd > rise) & (jd < set))[0]
    night = np.where((jd < rise) | (jd > set))[0]
    
    formask = foremask[saz]
    #Loop through channel pairs:
    for ch in channels:
        m = mask[:,ch] & foremask & ((jd > rise) & (jd < set))
        #Sort data by azimuth:

        
        nbins = 360/5
        btod,berr = Binning.DownSample(tod[m,ch],nbins,Errors=True)
        baz = Binning.DownSample(az[m],nbins)

        gd = np.where(np.abs(baz[0:-2]-baz[1:-1]) < 5)[0]

        #ax = pyplot.subplot(111,polar=True)
        #ax.plot(baz[gd]*np.pi/180.,btod[gd])#-np.min(btod[gd]))
        
        DataAccess.AppendFile('GroundSpill_H3I60-Epoch1-day.dat',np.append(baz[gd],btod[gd]))

        m = mask[:,ch] & foremask & ((jd < rise) | (jd > set))

        nbins = 360/5
        btod,berr = Binning.DownSample(tod[m,ch],nbins,Errors=True)
        baz = Binning.DownSample(az[m],nbins)

        gd = np.where(np.abs(baz[0:-2]-baz[1:-1]) < 5)[0]
        #ax.plot(baz[gd]*np.pi/180.,btod[gd])#-np.min(btod[gd]))
        #pyplot.show()
        
        DataAccess.AppendFile('GroundSpill_H3I60-Epoch1-night.dat',np.append(baz[gd],btod[gd]))
        
        

'''
        
       bl = 60. #seconds
                RFImask = MedianFilter(chVoltage,data['JD'][0,:,0]*24.*3600.,bl,mask,foremask)

                caldata[0,:,ch] = chVoltage[:maxlen]
                calmask[0,:,ch] = mask[:maxlen] & RFImask[:maxlen]


        #Write the downsampled data to new files.
        nFiles = len(filelist)
        nData = maxlen/nFiles

        #Generate temporary storage containers
        Out_caldata = np.zeros((1,nData,4))
        Out_calmask = np.ones((1,nData,4),dtype='bool')
        Out_caljd = np.zeros((1,nData,1))
        Out_calel = np.zeros((1,nData,1))
        Out_calaz = np.zeros((1,nData,1))
        
        for ifile,file in enumerate(filelist):
            ft = (file.split('/'))[-1]

            lo = nData*ifile

            if ifile < nFiles - 1:
                hi = nData*(ifile + 1)
            else:
                hi = np.min([maxlen-lo,nData]) + lo
                

            Out_caldata[0,0:hi-lo,:] = caldata[0,lo:hi,:]
            Out_calmask[0,0:hi-lo,:] = calmask[0,lo:hi,:]
            Out_caljd[0,0:hi-lo,0] = data['JD'][0,lo:hi,0]
            Out_calel[0,0:hi-lo,0] = data['EL'][0,lo:hi,0]
            Out_calaz[0,0:hi-lo,0] = data['AZ'][0,lo:hi,0]

            DataAccess.WriteFits.WriteTable(caldir+ft,[Out_caldata[:,0:hi-lo,:],
                                                           Out_calmask[:,0:hi-lo,:],
                                                           Out_caljd[:,0:hi-lo,:],
                                                           Out_calel[:,0:hi-lo,:],
                                                           Out_calaz[:,0:hi-lo,:]],calcols,clobber=True)
'''
