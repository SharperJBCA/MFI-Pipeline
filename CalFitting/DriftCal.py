#DriftCal.py
# Calculate model for drift of calibration diode

import numpy as np
import scipy.interpolate as ip

def CalDriftModel(data,c):
    '''
    Return model of drift in cal signal

    Arguments
    data -- Input TOD in the format (1,nSamples,nChannels)
    c -- Cal signal array, 0 = off, 1 = on, format: (1,nSamples,1)

    Notes: Locks onto the voltage spike causes when the cal diode
    switches off. Only accurately measured in low frequency channels
    where the diode voltage >> background noise.
    
    '''


    #Make a boolean array of calsig on == True
    cal = (c[0,:,0] == 1)
    calSig = data[0,cal,:]

    #First Estimate the drift model
    pulseTime = 50
    nPulses = len(c[0,cal,0])/pulseTime

    spikeVoltage = 26
    amax = np.zeros(nPulses)
    
    time = np.arange(nPulses)
    for i in range(nPulses):
        amax[i] = (calSig[i*pulseTime:(i+1)*pulseTime,0]).argmax() - spikeVoltage

    
    pfit = ip.interp1d(time,amax)
    
    return pfit
