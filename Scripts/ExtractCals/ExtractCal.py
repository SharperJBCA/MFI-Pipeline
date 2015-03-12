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
import SourceFitting

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


    
    #Foreground mask:
    foremap = hp.read_map('/nas/scratch/sharper/QUIJOTE/Pipeline/MAPS/maskmap_N512.fits')
    foremap = (foremap == 0)

    #Output columns:
    calcols = ['DATA','MASK','JD','EL','AZ']
    caldir = '/nas/scratch/sharper/QUIJOTE/datos/caltod_spill/'


    #Open all NOMINALXX files:
    files = np.loadtxt(filelistname,dtype='string')

    sra_cas,sdec_cas= SourceFitting.ModelFlux.SourceCoord('CASA')
    #sra_tau,sdec_tau = SourceFitting.ModelFlux.SourceCoord('CRAB')
    #sra_tau,sdec_tau = SourceFitting.ModelFlux.SourceCoord('CYGA')
    #sra_tau,sdec_tau = SourceFitting.ModelFlux.SourceCoord('3C274')


    nsampring = 30*50

    for f in files:

        filelist = DataAccess.FileList(f,dir=dir)
        data = DataAccess.FullObs(filelist,cols)
        
        for horn in range(1,2):
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

            nrings = pix.size/nsampring
            set = np.arange(nsampring)

            for ch in channels:
                foremask = (foremap[pix] == 1) & (data['MASK'][0,:,horn*4 + ch]==1)

                for i in range(nrings):
                    dat = data['DATA'][0,i*nsampring:(i+1)*nsampring,ch+horn*4]
                    mask = foremask[i*nsampring:(i+1)*nsampring]
                    adec = dec[i*nsampring:(i+1)*nsampring]*180./np.pi
                    #pyplot.plot(adec,dat,'-')
                    #pyplot.plot(adec[mask],dat[mask])
                    #pyplot.show()
                    setsize = len(set[mask])
                    if setsize > 4:
                        bkgdfit = np.poly1d(np.polyfit(set[mask],dat[mask],3))
                        data['DATA'][0,i*nsampring:(i+1)*nsampring,ch+horn*4] -= bkgdfit(set)
                    else:
                        data['DATA'][0,i*nsampring:(i+1)*nsampring,ch+horn*4] -= np.median(dat)



            radius = 3.
            cas = ((np.mod(ra*180./np.pi-180.,360) - (sra_cas-180.))**2 + (dec*180./np.pi - sdec_cas)**2 < radius**2)
            #tau = ((ra*180./np.pi - sra_tau)**2 + (dec*180./np.pi - sdec_tau)**2 < radius**2)

            #pyplot.plot(ra*180./np.pi,dec*180./np.pi,'.')
            #pyplot.plot(ra[tau]*180./np.pi,dec[tau]*180./np.pi,'.')
            #pyplot.show()


            calibs = np.where(cas & (data['MASK'][0,:,horn]==1))[0]

            if len(calibs) > 1:
                DataAccess.WriteFits.WriteTable('CasA/'+f,
                                                [data['DATA'][:,calibs,4:],
                                                 data['MASK'][:,calibs,4:],
                                                 data['JD'][:,calibs,:],
                                                 data['EL'][:,calibs,:],
                                                 data['AZ'][:,calibs,:]],cols,clobber=True)
