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


    notmask = (foremask & mask)

    BigStep = 180*50*10 #ms

    lastStep = 0


    j = 0 
    for i in range(nrings):

        #Define upper step in data. Catch left-over data on last loop
        if i == nrings-1:
            hi = jdmax+1.
        else:
            hi = (i+1)*jdlen

        #if lastRingSize > 0:
        mid = j*jdlen*50 + jdlen*50/2
        mhi = int(np.min([len(ijd),mid+BigStep]))
        mlo = int(np.max([0,mid-BigStep]))

        #Get the data in this loop that is not masked:
        #thisRing = np.where((notmaskjd[mlo:mhi] >= i*jdlen) & (notmaskjd[mlo:mhi] < hi))[0]
        thisRing = np.where((ijd[mlo:mhi] >= i*jdlen) & (ijd[mlo:mhi] < hi))[0]


        #If there is data, subtract median, check noise and spikes.

        lastRingSize = thisRing.size


        if lastRingSize > 0:
            j = j + 1

            if i == nrings-1:
                dhi = len(data)
            else:
                dhi = mlo+thisRing[-1]

            thisMask = notmask[lastStep:dhi]
            thisBkgd = (data[lastStep:dhi])[thisMask]
            data[lastStep:dhi] -= np.median(thisBkgd)

            
            if thisBkgd.size > 0:
                RFI = RFIMask[lastStep:dhi]

                #Filter spikes:
                maxvals[i] = np.max( np.abs(thisBkgd - np.median(thisBkgd) ))
                if (maxvals[i] > peak_cutoff):

                    RFI_short = RFI[thisMask]
                    spikes = np.where( np.abs(thisBkgd - np.median(thisBkgd) ) > peak_cutoff)[0]
                
                    #Get position of each spike:
                    spike_diff = spikes[1:-1]-spikes[0:-2]
                    spike_pos = np.where(spike_diff > 1)[0] + 1

                    RFI_short[np.max([spikes[0]-60,0.]):np.min([spikes[0]+60,len(thisBkgd)])] = False

                    if spike_pos.size > 0:
                        for spike in spikes[spike_pos]:
                            RFI_short[np.max([spike-60,0.]):np.min([spike+60,len(thisBkgd)])] = False
                
                    RFI[thisMask] = RFI_short
                    RFIMask[lastStep:dhi] = RFI
                
                #Filter ring?
                stds[i] = np.std(thisBkgd[RFI[[thisMask]]])

                if (stds[i] > std_cutoff) | (stds[i] == 0):
                    RFIMask[lastStep:dhi] = False

            lastStep = dhi
            

    #print 'SIZES:',data.size,RFIMask.size,std_cutoff,peak_cutoff
    #from matplotlib import pyplot
    #pyplot.plot(stds,'o')
    #pyplot.figure()
    #pyplot.plot(maxvals,'o')
    #pyplot.figure()
    #pyplot.plot(ijd[mask & foremask],data[mask & foremask],',')
    #pyplot.plot(ijd[RFIMask & mask & foremask],data[RFIMask & mask & foremask],',')
    #pyplot.plot(ijd[mask],data[mask],',')
    #pyplot.plot(ijd[mask & foremask],RFIMask[mask & foremask],',')

    #pyplot.show()

    return RFIMask
