#CalData.py
# Goals:
# 1) Fit a polynomial to G_ch*T_diode voltages
# 2) Divide V_ch by G_ch*T_diode
# 3) Multiply by 1 + P sin(4x + O) diode model
# 4) Record amplitude of source

import DataAccess
import numpy as np
import scipy.interpolate as ip
   
import Coordinates
import healpy as hp

import Masking
    
if __name__ == '__main__':

    #Directories:
    mfidir  = '/nas/scratch/sharper/QUIJOTE/MFI-Pipeline/MFI_INFO/'    #Where MFI parameter files are
    dir = '/nas/scratch/sharper/QUIJOTE/Pipeline/Preprocessing/DownSampData/data/' #Where downsampled data are
    caldir = '/nas/scratch/sharper/QUIJOTE/datos/caltod/' #Where to write calibrated data too
    foreground_mask_dir = '/nas/scratch/sharper/QUIJOTE/Pipeline/MAPS/' #Foreground mask map directory


    #Data access parameters:
    dioModel = np.loadtxt(mfidir+'DiodePolModel-Epoch2.dat')

    files = np.loadtxt('Nom75.lis',dtype='string')

    cols = ['DATA','MOD','AZ','EL','JD','CAL','CAL_JD','ERR']

    chanPairs = np.array([[0,6],[1,7],[2,4],[3,5]])
    TDiode = [1.,1.,1.,1.,0.38,0.41,0.41,0.41,2.12,1.36,2.12,1.55,0.49,0.20,0.54,0.50]
    freq   = [13.,11.,13.,11.,19.,17.,19.,17.,13.,11.,13.,11.,19.,17.,19.,17.]

    #Read in T-Point parameters:
    tpoints = np.loadtxt(mfidir+'focal_plane.ini')
    tpoints *= np.pi/180.    

    #Constants and extras:
    nHorns = 4
    nside = 512
    maskmap = np.ones(12*nside**2,dtype='bool')
    jd0 = 56244.
    
    #Foreground mask:
    foremap = hp.read_map(foreground_mask_dir+'maskmap_N512.fits')
    foremap = (foremap == 0)

    #Output columns:
    calcols = ['DATA','MASK','JD','EL','AZ']

    for f in files:
        #Read in the data
        print 'OPENING:',f
        filelist = DataAccess.FileList(f,dir=dir)
        data = DataAccess.FullObs(filelist,cols)

        #Output Containers:
        maxlen = (data['JD'].shape[1]/1000)*1000
        caldata = np.zeros((1,maxlen,4))
        calmask = np.ones((1,maxlen,4),dtype='bool')
        
        #Correct data for Diode Polarisation
        for i in range(data['DATA'].shape[2]):
            data['CAL'][0,:,i] /= (1. + dioModel[i,1]*np.cos( np.pi/180.*np.mod(np.median(data['MOD'][0,:,i/8]),90)*4. + dioModel[i,2]*2.*np.pi))

        for horn in range(2,3):
            #Get the RA/DEC coordinates:
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
            foremask = foremap[pix]
            
            #Generate Planet Mask:
            mask = Masking.GetPlanetMask(maskmap,data['JD'][0,:,0]+jd0,pix)

            #Loop through channel pairs:
            for ch,chPair in enumerate(chanPairs):
                
                gd      = np.where((data['CAL'][0,:,chPair[0]+horn*8] != 0))[0]
                #Median value of channel ratios ~ equal the ratio of the channels gains :
                r  = np.median(data['CAL'][0,gd,chPair[0]+horn*8])/np.median(data['CAL'][0,gd,chPair[1]+horn*8])
                         
        
                steplen = 1440
                nsteps = gd.size/steplen

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
                    pfit    = ip.interp1d(bjd,bca,kind='cubic')
                    pfit = pfit(data['JD'][0,:,0])
                else:
                    pfit = np.median(bca) + np.zeros(data['DATA'][0,:,chPair[0]+horn*8].size)

                #Calibrated data:
                chErr = np.sqrt(1./(1./data['ERR'][0,:,chPair[0]+horn*8]**2+1./(r*data['ERR'][0,:,chPair[1]+horn*8])**2))
                chVoltage = (data['DATA'][0,:,chPair[0]+horn*8]+ \
                             r*data['DATA'][0,:,chPair[1]+horn*8])/pfit #Divide by diode voltage
                chVoltage *= TDiode[horn*4 + ch] #Multiply by diode temperature


                
                bl = 60. #seconds
                RFImask = Masking.GetRingNoiseMask(chVoltage,data['JD'][0,:,0]*24.*3600.,bl,mask,foremask)

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
