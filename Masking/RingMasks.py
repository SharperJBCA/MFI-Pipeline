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

    notmask = np.where((foremask & mask))[0]
    notmaskjd = ijd[notmask]
    notmaskdata = data[notmask]

    BigStep = 180*50 #ms
    
    for i in range(nrings):

        #Define upper step in data. Catch left-over data on last loop
        if i == nrings-1:
            hi = jdmax+1.
        else:
            hi = (i+1)*jdlen

        mid = i*jdlen + jdlen/2
        mhi = int(np.min([len(notmaskjd),mid+BigStep]))
        mlo = int(np.max([0,mid-BigStep]))
     
        #Get the data in this loop that is not masked:
        thisRing = np.where((notmaskjd[mlo:mhi] >= i*jdlen) & (notmaskjd[mlo:mhi] < hi))[0]
        thisBkgd = notmaskdata[mlo+thisRing]


        #If there is data, subtract median, check noise and spikes. 
        if thisBkgd.size > 0:
            data[notmask[mlo+thisRing]] -= np.median(thisBkgd)
            stds[i] = np.std(thisBkgd)
            maxvals[i] = np.max( np.abs(thisBkgd-np.median(thisBkgd)) )

            #Filter ring? 
            if (stds[i] > std_cutoff) | (stds[i] == 0) | (maxvals[i] > peak_cutoff):
                RFIMask[notmask[mlo+thisRing]] = False

    return RFIMask
