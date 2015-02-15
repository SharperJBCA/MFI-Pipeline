#DownSampleData.py
#Goals: -Generate QUIJOTE data at X Hz instead of 1000 Hz
#       -Save an array of weights, just an array of std(TOD[chunk])
#       -Save mean cal signal from every 1 second pulse for each channel

import numpy as np
from matplotlib import pyplot
from DataAccess import FullObs
from DataAccess import WriteFits

import glob
from DataAccess import Count
import CalFitting
import Binning


files = np.loadtxt('FileList.txt',dtype='string')
dir   = '/nas/scratch/sharper/QUIJOTE/datos/tod/'
outdir= 'data/'
filesize_limit = 1000000

SR_in = 1000
SR_out= 50

cols = ['DATA','MOD','AZ','EL','JD','ERR','MASK','CAL','CAL_JD']
colsin = ['DATA','MOD_ANGLE','AZ','EL','JD','CAL']

for f in files:

    filelist = glob.glob(dir+f+'*')
    filelist = np.sort(filelist)
    dShapes,nFiles = Count(filelist,colsin,dir='')
    nSamples = dShapes[0][1]
    
    newlen   = nSamples/SR_in*SR_out
    nFilesOut= int(np.ceil(nFiles/16.))#int(np.ceil(newlen/filesize_limit))

    nLoops = 16#int(np.ceil(nFiles/nFilesOut))

    
    for iloop in range(nFilesOut):    
        #First read in the full file. e.g. 000.tod -> XXX.tod
        if iloop == nFilesOut - 1:
            data = FullObs(filelist[iloop*nLoops:],colsin,dir='')
        else:
            data = FullObs(filelist[iloop*nLoops:(iloop+1)*nLoops],colsin,dir='')

        newlen   = data['DATA'].shape[1]/SR_in*SR_out

        #from matplotlib import pyplot
        #pyplot.plot(data['DATA'][0,:,0],',')
        #pyplot.plot(data['DATA'][0,:,6],',')        
        #pyplot.show()
        #Record each calsig for each channel.
        calamps,calbkgd,caljd = CalFitting.AvgCalSig(data['DATA'],data['CAL'],jd=data['JD'])

        #Remove CalSignal:
        #c    = (data['CAL'][0,:,0] == 0)
        #data['DATA']  = data['DATA'][:,c,:]
        #data['MOD_ANGLE']  = data['MOD_ANGLE'][:,c,:]
        #data['AZ']  = data['AZ'][:,c,:]
        #data['EL']  = data['EL'][:,c,:]
        #data['JD']  = data['JD'][:,c,:]
        #c = None

        #Downsample the data to whatever using DSD.
        ds_data = np.zeros((1,newlen,data['DATA'].shape[2]))
        ds_derr = np.zeros((1,newlen,data['DATA'].shape[2]))
        ds_mod  = np.zeros((1,newlen,data['MOD_ANGLE'].shape[2]))
        ds_az   = np.zeros((1,newlen,1))
        ds_el   = np.zeros((1,newlen,1))
        ds_jd   = np.zeros((1,newlen,1))
        ds_cal  = np.zeros((1,newlen,1))

        for i in range(data['DATA'].shape[2]):
            ds_data[0,:,i],ds_derr[0,:,i] = Binning.DownSample(data['DATA'][0,:,i],newlen,Errors=True)

        for i in range(data['MOD_ANGLE'].shape[2]):
            ds_mod[0,:,i] = Binning.DownSample(data['MOD_ANGLE'][0,:,i],newlen)

        ds_az[0,:,0] = Binning.DownSample(data['AZ'][0,:,0],newlen)

        ds_el[0,:,0] = Binning.DownSample(data['EL'][0,:,0],newlen)

        ds_jd[0,:,0]  = Binning.DownSample(data['JD'][0,:,0],newlen)
        ds_cal[0,:,0] = Binning.DownSample(data['CAL'][0,:,0],newlen)

        top = np.roll(ds_cal[0,:,0],5)
        top[0:5] = 0.
        bot =  np.roll(ds_cal[0,:,0],-5)
        bot[-6:] = 0.
        ds_cal[0,:,0] += top + bot
        ds_cal[0,(ds_cal[0,:,0] > 0),0] = 1

        data = None

        #nrings = int(np.max(ds_jd-np.min(ds_jd))*24.*60.*60./60)
        #print nrings

        #ds_rings = np.zeros((1,ds_jd.size,1))
        #ds_rings[0,0:nrings*60*SR_out,0] = np.repeat(np.arange(nrings),60*SR_out)

        c = np.where((ds_cal == 0))[0]

        #Write the downsampled data to a new file.
        WriteFits.WriteTable(outdir+f+'-'+str(iloop).zfill(3)+'-ds.fits',[ds_data,
                                                                          ds_mod,
                                                                          ds_az,
                                                                          ds_el,
                                                                          ds_jd,
                                                                          ds_derr,
                                                                          ds_cal,
                                                                          calamps,
                                                                          caljd],cols)
        ds_data = None
        ds_mod  = None
        ds_az   = None
        ds_el   = None
        ds_jd   = None
        ds_derr = None
        calamps = None
        caljd   = None
