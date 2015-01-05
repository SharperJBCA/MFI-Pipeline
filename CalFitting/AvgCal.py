#AvgCal.py
# Provide method for averaging cal diode signal in QUIJOTE TOD.

import numpy as np

def AvgCalSig(data,c,jd=[None],DriftModel=[None]):
    '''
    Return array of cal voltages

    Arguments
    data -- Input TOD in the format (1,nSamples,nChannels)
    c -- Cal signal array, 0 = off, 1 = on, format: (1,nSamples,1)

    Keyword Arguments
    jd = Return the julian date of each diode signal
    DriftModel = Model describing how cal signal timing drifts
    '''

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
    if jd[0]:
        calJd = jd[0,cal,0]

    calSig = data[0,cal,:]
    pulseTime = 50
    nPulses = len(c[0,cal,0])/pulseTime
    cal = None

    if DriftModel[0]: 
        pfit = np.poly1d(DriftModel)
    else:
        pfit = np.poly1d([0.,0.])

    #Calculate average signal - Calsignal = 50 samples
    pulses = np.zeros((pulseTime,nPulses,nChan))
    amps = np.zeros((1,nPulses,nChan))
    cjd = np.zeros((1,nPulses,1))

    for j in range(nChan):
        for i in range(nPulses):

            if (i+1)*pulseTime+int(pfit(i))-i*pulseTime+int(pfit(i)) > 0:
                pulses[:,i,j] = calSig[i*pulseTime+int(pfit(i)):(i+1)*pulseTime+int(pfit(i)),j]
                upper = np.mean(pulses[8 :21,i,j])
                lower = np.mean(pulses[33:43,i,j])
                amps[0,i,j] = upper - lower

            if (jd[0]) & (j == 0):
                cjd[0,i,0] = np.mean(calJd[i*pulseTime+int(pfit(i)):(i+1)*pulseTime+int(pfit(i))])


    if (jd[0]):
        return amps,cjd
    else:
        return amps

