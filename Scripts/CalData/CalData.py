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
    
    #Foreground mask:
    foremap = hp.read_map(foreground_mask_dir+'maskmap_N512.fits')
    foremap = (foremap == 0)

    #Output columns:
    calcols = ['DATA','MASK','JD','EL','AZ','SPIKEMASK']

    for f in files:
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
            emask = Masking.GetPlanetMask(maskmap,data['JD'][0,:,0]+jd0,pix)
            mask = (emask) & (data['MASK'][0,:,0] == 0)

            #Loop through channel pairs:
            for ch,chPair in enumerate(chanPairs):



                gd = np.where((data['CAL_JD'][0,:,0] > np.min(data['JD'][0,:,0])) & (data['CAL_JD'][0,:,0] < np.max(data['JD'][0,:,0])))[0]

                cals1 =data['CAL'][0,gd,chPair[0]+horn*8]
                cals2 =data['CAL'][0,gd,chPair[1]+horn*8]
                caljd =data['CAL_JD'][0,gd,0]

                imaskcal = ip.interp1d(data['JD'][0,:,0],emask.astype('i'),kind='nearest')

                maskcal = imaskcal(caljd)

                
                gd      = np.where((cals1 > 0) & (maskcal == 1))[0]                         
        
                steplen =1440
                nsteps = gd.size/steplen

                #Bin the G*T_dio signal and timestamps:
                if nsteps > 1:
                    #If Data file is large enough to have main gain step measurements:
                    newlen = nsteps * steplen
                    bca    = np.median(np.reshape(cals1[gd[:newlen]]+ \
                                                  rfactors[chPair[0]+horn*4]*cals2[gd[:newlen]],(nsteps,steplen)),axis=1)

                    undervoltage = np.where(bca < 0)[0]
                    if undervoltage.size > 0:
                        from matplotlib import pyplot
                        print f
                        pyplot.plot(data['CAL'][0,gd[:newlen],chPair[0]+horn*8],'o')
                        pyplot.show()
        
                    bjd    = np.median(np.reshape(caljd,(nsteps,steplen)),axis=1)
                else:
                    #Else just measure the median gain for the whole file:
                    bca = np.array([np.median(cals1[gd]+ \
                                              rfactors[chPair[0]+horn*4]*cals2[gd])])

                
                #Produce a model of the gain drifts with JD
                if bca.size > 4:
                    #If the file is big enough to allow for a cubic spline interpolation:
                    bjd[0]  = np.min([np.min(data['JD'][0,:,0]),np.min(caljd[gd])])
                    bjd[-1] = np.max([np.max(data['JD'][0,:,0]),np.max(caljd[gd])])        
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

#                 scaleFactor = np.median((sw1 + sw2)/pfit)

#                 baseline = 1400
#                 a = np.zeros(chVoltage.size/baseline)
#                 temp = sw1 + sw2
#                 for i in range(a.size):
#                     a[i] = np.median(temp[i*baseline:(i+1)*baseline])

#                 a = np.repeat(a,baseline)
#                 from matplotlib import pyplot
#                 day = data['JD'][0,:,0]
#                 print np.mean(data['JD'][0,:,0])+jd0

#                 day -= int(np.min(day))
#                 day -= 1.
#                 #np.mod(data['JD'][0,:,0],1.)
#                 day = day[:a.size]
#                 sday = np.argsort(day)
#                 a = a/a[0]
#                 t = pfit*scaleFactor/a[0] - ( pfit[0]*scaleFactor/a[0] - 1.)
#                 t = t[:a.size]
#                 pyplot.plot(day[sday],a[sday] ,'-',color='#ADADAD',lw=2,label='Channel Voltage')
#                 pyplot.plot(day[sday],t[sday],'--',color='#0059F2',lw=2,label='Diode Voltage')

#                 ylims = pyplot.gca().axes.get_ylim()
#                 pyplot.plot([0.26,0.26],[ylims[0],ylims[1]],'-',color='#F29900',lw=2,label='Sunrise',alpha=0.8)
#                 pyplot.plot([0.81,0.81],[ylims[0],ylims[1]],'-',color='#F22000',lw=2,label='Sunset',alpha=0.8)
#                 pyplot.ylim(ylims[0],ylims[1])
#                 pyplot.legend(frameon=False,loc='upper left')
#                 pyplot.xlabel('Fraction of Day')
#                 pyplot.show()
#                 print stop
                #print horn, ch,TDiode[horn*4 + ch] ,scaleFactor
                chVoltage *= TDiode[horn*4 + ch] #Multiply by diode temperature


                
                bl = 60. #seconds
                print cutoffs[horn*4 + ch],peaks[horn*4+ch]
                spikemask,RFImask = Masking.GetRingNoiseMask(chVoltage,data['JD'][0,:,0]*24.*3600.,bl,mask,foremask,
                                                             std_cutoff=cutoffs[horn*4 + ch],
                                                             peak_cutoff=peaks[horn*4+ch])

                #from matplotlib import pyplot
                #pyplot.plot(dec[mask & RFImask]*180./np.pi,chVoltage[mask & RFImask],',')
                #pyplot.show()

                #
                #stopFact = -300000
                #pyplot.plot(np.mod(data['JD'][0,(mask[:stopFact] & RFImask[:stopFact]),0],1)*24.,
                #            sw1[mask[:stopFact] & RFImask[:stopFact]]+sw2[mask[:stopFact] & RFImask[:stopFact]],',',color='#ADADAD',label='System')
                #pyplot.plot(np.mod(data['JD'][0,(mask[:stopFact] & RFImask[:stopFact]),0],1)*24.,
                #            pfit[mask[:stopFact] & RFImask[:stopFact]]*scaleFactor,',',color='#0059F2',label='Diode')

                #pyplot.legend(frameon=False)
                #pyplot.xlabel('Time (Hours)')
                #pyplot.ylabel('Temperature (K)')

                #pyplot.plot(np.mod(data['JD'][0,(mask & RFImask),0],1),
                #            chVoltage[mask & RFImask],',',color='#ADADAD')
                #pyplot.show()

                #Horn-1 because we arent using horn 1 data.
                caldata[0,:,(horn-1)*4 + ch] = chVoltage[:maxlen]
                calmask[0,:,(horn-1)*4 + ch] = mask[:maxlen] & RFImask[:maxlen]
                calspikes[0,:,(horn-1)*4 + ch] =  spikemask[:maxlen]#(data['MASK'][0,:maxlen,0] == 0) &

        #Write the downsampled data to new files.
        nFiles = len(filelist)
        nData = maxlen/nFiles

        #Generate temporary storage containers
        Out_caldata = np.zeros((1,nData,8))
        Out_calmask = np.ones((1,nData,8),dtype='bool')
        Out_caljd = np.zeros((1,nData,1))
        Out_calel = np.zeros((1,nData,1))
        Out_calaz = np.zeros((1,nData,1))
        Out_calspikemask = np.zeros((1,nData,8),dtype='bool')

        for ifile,file in enumerate(filelist):
            ft = (file.split('/'))[-1]

            lo = nData*ifile

            if ifile < nFiles - 1:
                hi = nData*(ifile + 1)
            else:
                hi = np.min([maxlen-lo,nData]) + lo
                

            Out_caldata[0,0:hi-lo,:] = caldata[0,lo:hi,:]
            Out_calmask[0,0:hi-lo,:] = calmask[0,lo:hi,:]
            Out_calspikemask[0,0:hi-lo,:] = calspikes[0,lo:hi,:]


            
            Out_caljd[0,0:hi-lo,0] = data['JD'][0,lo:hi,0]
            Out_calel[0,0:hi-lo,0] = data['EL'][0,lo:hi,0]
            Out_calaz[0,0:hi-lo,0] = data['AZ'][0,lo:hi,0]

            DataAccess.WriteFits.WriteTable(caldir+ft,[Out_caldata[:,0:hi-lo,:],
                                                       Out_calmask[:,0:hi-lo,:],
                                                       Out_caljd[:,0:hi-lo,:],
                                                       Out_calel[:,0:hi-lo,:],
                                                       Out_calaz[:,0:hi-lo,:],
                                                       Out_calspikemask[:,0:hi-lo,:]],calcols,clobber=True)

