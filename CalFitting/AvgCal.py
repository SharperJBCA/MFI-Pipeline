#AvgCal.py
# Provide method for averaging cal diode signal in QUIJOTE TOD.

import numpy as np
from CalFitting import DriftCal 

def AvgCalSig(data,c,jd=[None],mod=None,DriftModel=None):
    '''
    Return array of cal voltages

    Arguments
    data -- Input TOD in the format (1,nSamples,nChannels)
    c -- Cal signal array, 0 = off, 1 = on, format: (1,nSamples,1)

    Keyword Arguments
    jd = Return the julian date of each diode signal
    DriftModel = Model describing how cal signal timing drifts
    '''


    ijd = jd[0]
    imod = mod[0]

    #Get the number of channels in data:
    nChan = data.shape[2]

    #Make a boolean array of calsig on == True
    cal = (c[0,:,0] == 1)

    #Check first and last block is a full block
    if (np.sum(c[0,0:1000,0]) < 1000)  & (np.sum(c[0,1000:2000,0]) == 0):
        cal[0:1000] = False

    if (np.sum(c[0,-1000:,0]) < 1000)  & (np.sum(c[0,-2000:-1000,0]) == 0):
        cal[-1000:] = False

    #Put the calsigs into a new container

    calSig = data[0,cal,:]
    pulseTime = 50
    nPulses = len(c[0,cal,0])/pulseTime

    if (type(ijd) != type(None)):
        calJd = jd[0,cal,0]
        cjd = np.zeros((1,nPulses,1))

    if (type(imod) != type(None)):
        calMod = mod[0,cal,:]        
        cmod = np.zeros((1,nPulses,mod.shape[2]))

    
    cal = None

    #Calculate the drift model for this observation:
    pfit = DriftCal.CalDriftModel(data,c)

    #Calculate average signal - Calsignal = 50 samples
    pulses = np.zeros((pulseTime,nChan))
    amps  = np.zeros((1,nPulses,nChan))
    bkgds = np.zeros((1,nPulses,nChan))

    for i in range(nPulses):
        hi = (i+1)*pulseTime+int(pfit(i))
        lo = i*pulseTime+int(pfit(i))

        if (hi-lo > 0) & (lo > 0) & (hi < calSig.shape[0]):
            pulses[:,:] = calSig[i*pulseTime+int(pfit(i)):(i+1)*pulseTime+int(pfit(i)),:]


            upper = np.mean(pulses[11 :21,:],axis=0)
            lower = np.mean(pulses[35:45,:],axis=0)
            amps[0,i,:]  = upper - lower
            bkgds[0,i,:] = lower


            if (type(ijd) != type(None)):
                cjd[0,i,0] = np.mean(calJd[i*pulseTime+int(pfit(i)):(i+1)*pulseTime+int(pfit(i))])

            if (type(imod) != type(None)):
                cmod[0,i,:] = np.mean(calMod[i*pulseTime+int(pfit(i)):(i+1)*pulseTime+int(pfit(i)),:],axis=0)

    outlist = [amps,bkgds]
    if (type(ijd) != type(None)):
        outlist = outlist + [cjd]
    if (type(imod) != type(None)):
        outlist = outlist + [cmod]

    return outlist
