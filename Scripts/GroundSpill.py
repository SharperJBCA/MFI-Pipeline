#CalData.py
# Goals:
# 1) Fit a polynomial to G_ch*T_diode voltages
# 2) Divide V_ch by G_ch*T_diode
# 3) Multiply by 1 + P sin(4x + O) diode model
# 4) Record amplitude of source
#
#NOTE: CHANGE HORN MANUALLY!!! (For now)

import DataAccess
import numpy as np
from matplotlib import pyplot
import scipy.interpolate as ip
import Binning
   
import Coordinates
import healpy as hp
from scipy.interpolate import interp1d as ip
    
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
    caldir = '/nas/scratch/sharper/QUIJOTE/datos/caltod_spill/'


    #Open all NOMINALXX files:
    files = np.loadtxt('Nom40.lis',dtype='string')

    for f in files:

        filelist = DataAccess.FileList(f,dir=dir)
        data = DataAccess.FullObs(filelist,cols)
        
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
        tod = data['DATA'][0,saz,:]
        mask = (data['MASK'][0,saz,:].astype('bool') )
    
        formask = foremask[saz]
        #Loop through channels:
        for ch in channels:
            m = mask[:,ch] & foremask

            #Bin data and generate template:
            nbins = 360/5
            btod,berr = Binning.DownSample(tod[m,ch],nbins,Errors=True)
            baz = Binning.DownSample(az[m],nbins)

            gd = np.where(np.abs(baz[0:-2]-baz[1:-1]) < 6)[0]

            
            spillfit = ip(np.concatenate(( np.array(baz[gd]-360.)
                                           ,baz[gd],
                                           np.array(baz[gd]+360.) )),
                          np.concatenate(( np.array(btod[gd]),
                                           btod[gd],
                                           np.array(btod[gd]) )))

            #Subtract template from data
            data['DATA'][0,:,ch] -= spillfit(data['AZ'][0,:,0])        
                

        #Write the downsampled data to new files.
        nFiles = len(filelist)
        maxlen = (data['JD'][0,:,0]).size
        nData = (data['JD'][0,:,0]).size/nFiles
        
        for ifile,file in enumerate(filelist):
            ft = (file.split('/'))[-1]

            lo = nData*ifile

            if ifile < nFiles - 1:
                hi = nData*(ifile + 1)
            else:
                hi = np.min([maxlen-lo,nData]) + lo
                
            DataAccess.WriteFits.WriteTable(caldir+ft,[data['DATA'][:,lo:hi,:],
                                                       data['MASK'][:,lo:hi,:],
                                                       data['JD'][:,lo:hi,:],
                                                       data['EL'][:,lo:hi,:],
                                                       data['AZ'][:,lo:hi,:]],cols,clobber=True)
