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
import sys

import Masking

import jdcal

def BrightnessModel(stellaData):
    tmdl = np.loadtxt(stellaData,dtype='string')
    temp = np.array(tmdl[:,2],dtype='float')
    date = np.array([day.split('.') for day in tmdl[:,0]],dtype='float')
    time = np.array([(t.split(':'))[0:2] for t in tmdl[:,1]],dtype='float')
    
    jd = np.array([jdcal.gcal2jd(d[2],d[1],d[0]) for d in date])
    mjd = jd[:,1] + time[:,0]/24. + time[:,1]/(60.*24.)
    
    return ip(mjd,temp)


if __name__ == '__main__':

    filelistname = sys.argv[1]

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
    
    #Cutoffs
    cutoffs = [0.056,0.029,0.059,0.032,
               0.033,0.026,0.021,0.023]
    peaks   = [0.28 ,0.13 ,0.28 ,0.13 ,
               0.09 ,0.09 ,0.09 ,0.09 ]

    #Constants and extras:
    nHorns = 4
    nside = 512
    maskmap = np.ones(12*nside**2,dtype='bool')
    jd0 = 56244.

    set = 0.8
    rise= 0.25


    #Sky Brightness Model:
    stellaData = 'stella-brightness-2013-2014.txt'
    skymdl = BrightnessModel(stellaData)
    
    #Foreground mask:
    foremap = hp.read_map('/nas/scratch/sharper/QUIJOTE/Pipeline/MAPS/maskmap_N512.fits')
    foremap = (foremap == 0)

    #Output columns:
    calcols = ['DATA','MASK','JD','EL','AZ']
    caldir = '/nas/scratch/sharper/QUIJOTE/datos/caltod_spill/'


    #Open all NOMINALXX files:
    files = np.loadtxt(filelistname,dtype='string')

    for f in files:

        filelist = DataAccess.FileList(f,dir=dir)
        data = DataAccess.FullObs(filelist,cols)
        
        for horn in range(2):
            ra,dec,p = Coordinates.Hor2Sky(data['AZ'][0,:,0]*np.pi/180.,
                                           data['EL'][0,:,0]*np.pi/180.,
                                           data['JD'][0,:,0]+jd0,TPoints={'xpos':tpoints[horn+1,1],
                                                                          'ypos':tpoints[horn+1,2],
                                                                          'Pf':  tpoints[horn+1,3],
                                                                          'Px':  tpoints[horn+1,4],
                                                                          'Py':  tpoints[horn+1,5],
                                                                          'Pc':  tpoints[horn+1,6],
                                                                          'Pn':  tpoints[horn+1,7],
                                                                          'Pa':  tpoints[horn+1,8],
                                                                          'Pb':  tpoints[horn+1,9]})
        
            pix = hp.ang2pix(nside,np.pi/2.-dec,ra)
        
            #Mask any bright Galactic sources:
            foremask = foremap[pix]

            #Sort AZ
            saz = np.argsort(data['AZ'][0,:,0])
            az  = data['AZ'][0,saz,0]
            tod = data['DATA'][0,saz,:]
            mjd = data['JD'][0,saz,0] + jd0       
            mask = (data['MASK'][0,saz,:].astype('bool') )
    
            formask = foremask[saz]
            #Loop through channels:
            for ch in channels:
                m = mask[:,ch+horn*4] & foremask

                #Loops twice, day and night:
                for period in ['DAY','NIGHT']:
                    if period == 'DAY':
                        dayNight = (skymdl(mjd) > 30)
                    else:
                        dayNight = (skymdl(mjd) <= 30)

                    print dayNight.shape,m.shape
                    #Bin data and generate template:
                    nbins = 360/5
                    try:
                        btod,berr = Binning.DownSample(tod[(m & dayNight),ch+horn*4],nbins,Errors=True)
                        baz = Binning.DownSample(az[(m & dayNight)],nbins)
                
                        gd = np.where(np.abs(baz[0:-2]-baz[1:-1]) < 6)[0]

                
                        spillfit = ip(np.concatenate(( np.array(baz[gd]-360.)
                                                       ,baz[gd],
                                                       np.array(baz[gd]+360.) )),
                                      np.concatenate(( np.array(btod[gd]),
                                                       btod[gd],
                                                       np.array(btod[gd]) )))

                        #Subtract template from data
                        data['DATA'][0,(dayNight),ch+horn*4] -= spillfit(data['AZ'][0,(dayNight),0])        
                    except TypeError:
                        continue
                    

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
