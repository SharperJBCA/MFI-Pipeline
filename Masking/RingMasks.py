import numpy as np

def GetRingNoiseMask(data,jd,jdlen,mask,foremask,std_cutoff=0.03,peak_cutoff=0.2):
    '''
    Return mask for spikes and noisey data in the TOD. Data is modified in place.

    Arguments
    data     -- Array containing TOD
    jd       -- Array of times associated with each TOD sample
    jdlen    -- Length of time to inspect TOD over
    mask     -- Current data mask (E.g. Mask of ephemeris sources)
    foremask -- Mask of bright Galactic foregrounds

    Keyword argument

    std_cutoff  -- Maximum standard-deviation of TOD over a period of jdlen.
    peak_cutoff -- Maximum value the TOD over a period of jdlen can contain.


    '''


    
    RFIMask = np.ones(mask.size,dtype='bool')

    #Make time run from 0 -> T:
    ijd = jd- np.min(jd)

    #Calculate the number of rings/loops:
    jdmax = np.max(ijd)
    nrings = int(jdmax)/int(jdlen)

    stds = np.zeros(nrings)
    maxvals = np.zeros(nrings)

    
    for i in range(nrings):

        #Define upper step in data. Catch left-over data on last loop
        if i == nrings-1:
            hi = jdmax+1.
        else:
            hi = (i+1)*jdlen

        #Get the data in this loop that is not masked:
        thisRing = ((ijd >= i*jdlen) & (ijd < hi))
        thisBkgd = data[(thisRing & foremask & mask)]


        #If there is data, subtract median, check noise and spikes. 
        if thisBkgd.size > 0:
            data[thisRing] -= np.median(thisBkgd)
            maxvals[i] = np.max(np.abs(thisBkgd-np.median(thisBkgd)))
            stds[i] = np.std(thisBkgd)

            #Filter ring? 
            if (stds[i] > std_cutoff) | (stds[i] == 0) | (maxvals[i] > peak_cutoff):
                RFIMask[thisRing] = False

    return RFIMask
