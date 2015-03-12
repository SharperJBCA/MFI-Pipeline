import DataAccess
import numpy as np
import scipy.interpolate as ip

from matplotlib import pyplot
import Coordinates
import healpy as hp

import Masking

import sys
   
if __name__ == '__main__':

    filelistname = sys.argv[1]

    #Directories:
    mfidir  = '/nas/scratch/sharper/QUIJOTE/MFI-Pipeline/MFI_INFO/'    #Where MFI parameter files are
    dir = '/nas/scratch/sharper/QUIJOTE/Pipeline/Preprocessing/DownSampData/data/' #Where downsampled data are
    caldir = '/nas/scratch/sharper/QUIJOTE/datos/caltod/' #Where to write calibrated data too
    foreground_mask_dir = '/nas/scratch/sharper/QUIJOTE/Pipeline/MAPS/' #Foreground mask map directory


    #Data access parameters:
    dioModel = np.loadtxt(mfidir+'DiodePolModel-Epoch2.dat')
    rfactors = np.loadtxt(mfidir+'AvgRFactors.dat')

    files = np.loadtxt(filelistname,dtype='string')

    cols = ['DATA','MOD','AZ','EL','JD','CAL','CAL_JD','ERR','MASK']

    chanPairs = np.array([[0,6],[1,7],[2,4],[3,5]])
    TDiode = np.loadtxt(mfidir+'TDiodes.dat')[:,1]
    freq   = np.loadtxt(mfidir+'FreqCens.dat')
    cutoffs = [1.,1.,1.,1.,0.045 ,0.045,0.045 ,0.045,0.023 ,0.024 ,0.023 ,0.024 ,1.,1.,1.,1.]
    peaks   = [1.,1.,1.,1.,0.15  ,0.15,0.15  ,0.15  ,0.1   ,0.1   ,0.1   ,0.1   ,1.,1.,1.,1.]
    
    #Read in T-Point parameters:
    tpoints = np.loadtxt(mfidir+'focal_plane.ini')
    tpoints *= np.pi/180.    

    #Constants and extras:
    nHorns = 4
    nside = 512
    maskmap = np.ones(12*nside**2,dtype='bool')
    jd0 = 56244.
    Tatm = 274.
    
    #Foreground mask:
    foremap = hp.read_map(foreground_mask_dir+'maskmap_N512.fits')
    foremap = (foremap == 0)

    #Output columns:
    calcols = ['DATA','MASK','JD','EL','AZ','SPIKEMASK']

    #Absorbs
    absorbs = np.zeros((17,len(files)))

    for ifile,f in enumerate(files):
        #Read in the data
        print 'OPENING:',f
        filelist = DataAccess.FileList(f,dir=dir)
        data = DataAccess.FullObs(filelist,cols)

        #Output Containers:
        maxlen = (data['JD'].shape[1]/1000)*1000
        caldata = np.zeros((1,maxlen,8))
        calmask = np.ones((1,maxlen,8),dtype='bool')
        calspikes = np.zeros((1,maxlen,8),dtype='bool')
        
        #Correct data for Diode Polarisation
        #for i in range(data['DATA'].shape[2]):
        #    data['CAL'][0,:,i] /= (1. + dioModel[i,1]*np.cos( np.pi/180.*np.mod(np.median(data['MOD'][0,:,i/8]),90)*4. + dioModel[i,2]*2.*np.pi))

        for horn in range(1,3):
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
            mask = (mask) & (data['MASK'][0,:,0] == 0)

            #Loop through channel pairs:
            for ch,chPair in enumerate(chanPairs):
                
                gd      = np.where((data['CAL'][0,:,chPair[0]+horn*8] != 0))[0]                         
        
                steplen = 1440
                nsteps = gd.size/steplen

                #Bin the G*T_dio signal and timestamps:
                if nsteps > 0:
                    #If Data file is large enough to have main gain step measurements:
                    newlen = nsteps * steplen
                    bca    = np.median(np.reshape(data['CAL'][0,gd[:newlen],chPair[0]+horn*8]+ \
                                                  rfactors[chPair[0]+horn*4]*data['CAL'][0,gd[:newlen],chPair[1]+horn*8],(nsteps,steplen)),axis=1)

                    undervoltage = np.where(bca < 0)[0]
                    if undervoltage.size > 0:
                        from matplotlib import pyplot
                        print f
                        pyplot.plot(data['CAL'][0,gd[:newlen],chPair[0]+horn*8],'o')
                        pyplot.show()
        
                    bjd    = np.median(np.reshape(data['CAL_JD'][0,gd[:newlen],0],(nsteps,steplen)),axis=1)
                else:
                    #Else just measure the median gain for the whole file:
                    bca = np.array([np.median(data['CAL'][0,gd,chPair[0]+horn*8]+ \
                                              rfactors[chPair[0]+horn*4]*data['CAL'][0,gd,chPair[1]+horn*8])])

                
                #Produce a model of the gain drifts with JD
                if bca.size > 4:
                    #If the file is big enough to allow for a cubic spline interpolation:
                    bjd[0]  = np.min([np.min(data['JD'][0,:,0]),np.min(data['CAL_JD'][0,gd,0])])
                    bjd[-1] = np.max([np.max(data['JD'][0,:,0]),np.max(data['CAL_JD'][0,gd,0])])        
                    pfit    = ip.interp1d(bjd,bca,kind='cubic')
                    pfit = pfit(data['JD'][0,:,0])
                else:
                    #Else just use a single median for the whole file:
                    pfit = np.median(bca) + np.zeros(data['DATA'][0,:,chPair[0]+horn*8].size)


                #Calibrated channel error:
                chErr = np.sqrt(data['ERR'][0,:,chPair[0]+horn*8]**2 + (rfactors[chPair[0]+horn*4]*data['ERR'][0,:,chPair[1]+horn*8])**2)
         
                #Calibrated data:

                sw1 =                            data['DATA'][0,:,chPair[0]+horn*8]#/data['ERR'][0,:,chPair[0]+horn*8]**2
                sw2 = rfactors[chPair[0]+horn*4]*data['DATA'][0,:,chPair[1]+horn*8]#/data['ERR'][0,:,chPair[1]+horn*8]**2
                
                
                chVoltage = (sw1 + sw2)/pfit#/(1./chErr**2) #Divide by diode voltage and weighted average
                scaleFactor = np.median((sw1 + sw2)/pfit)

                #from matplotlib import pyplot
                #pyplot.plot(chVoltage[mask]*TDiode[horn*4 + ch] ,',')
                #pyplot.show()


                print horn, ch,TDiode[horn*4 + ch] ,scaleFactor
                chVoltage *= TDiode[horn*4 + ch] #Multiply by diode temperature

                pfit = np.poly1d(np.polyfit(1./np.sin(data['EL'][0,mask,0]*np.pi/180.),chVoltage[mask],1))

                absorbs[(horn-1)*4 + ch,ifile] = 100.-np.exp(-pfit[1]/Tatm)*100.
                absorbs[(horn-1)*4 + ch+8,ifile] = pfit[0]
                
                print 'Horn:',horn,100.-np.exp(-pfit[1]/Tatm)*100.
                vals = np.linspace(1,2,100)

                colors = np.random.uniform(size=chVoltage.size,low=0.4)
                #pyplot.scatter(1./np.sin(data['EL'][0,:,0]*np.pi/180.),chVoltage,c=colors,cmap='gray',s=80,lw=0)
                pyplot.plot(1./np.sin(data['EL'][0,mask,0]*np.pi/180.),chVoltage[mask],',',color='#0059F2')
                pyplot.plot(vals,pfit(vals),'--',lw=3,color='#F29900')     
                pyplot.ylabel('Voltage',size=16)
                pyplot.xlabel('Airmass',size=16)
                #pyplot.ylim(0.5,4.)
                #pyplot.xlim(0,100)
                pyplot.xticks(size=14)
                pyplot.yticks(size=14)
                pyplot.show()

                
                #pyplot.plot(1./np.sin(data['EL'][0,:,0]*np.pi/180.),chVoltage,',')
                #pyplot.plot(vals,pfit(vals),'-',lw=3)                
                #pyplot.show()
            
        absorbs[-1,ifile] = np.median(data['JD'][0,:,0])
    #np.savetxt('DipAbsorption.txt',absorbs)
