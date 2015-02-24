import numpy as np
import scipy.fftpack as fftpack

#----Smooth Func-------#

def smoothListGaussian(list,strippedXs=False,degree=5):

    listframe = np.concatenate((list[::-1],list,list[::-1]))

    padlen = 2.**np.ceil(np.log(listframe.size)/np.log(2.))

    fft = fftpack.rfft(listframe,n=padlen)
    
    window = (np.cos(np.arange(degree)/float(degree) * np.pi) *0.5 + 0.5)
    filter = np.append(window,np.zeros(padlen-degree))

    smoothed = fftpack.irfft(fft*filter)
    smoothed = smoothed[len(list):len(list)*2]
    
    return smoothed


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
    SpikeMask = np.zeros(mask.size,dtype='bool')

    #notforemask = (foremask & mask)
    notmask = (foremask & mask)
    #notmask = np.ones(notmask.size,dtype='bool')

    nSampRing = 6000 #1 ring #ms
    nrings = data.size/nSampRing
    
    stds = np.zeros(nrings)
    maxvals = np.zeros(nrings)

    lastStep = 0





    j = 0 
    for i in range(nrings):

        #Define upper step in data. Catch left-over data on last loop
        if i == nrings-1:
            hi = data.size
            lo = i*nSampRing
        else:
            hi = (i+1)*nSampRing
            lo = i*nSampRing
            

        thisMask = np.where(notmask[lo:hi])[0]
        thisBkgd = (data[lo:hi])[thisMask]
        data[lo:hi] -= np.median(thisBkgd)
        thisBkgd -= np.median(thisBkgd)

        #from matplotlib import pyplot
        #print peak_cutoff
        #pyplot.plot(data,',')
        #pyplot.figure()
        #pyplot.plot(data[lo:hi],',')
        

        if thisBkgd.size > 0:
            RFI = RFIMask[lo:hi]

            #Filter spikes:
            smthBkgd = thisBkgd-smoothListGaussian(thisBkgd,degree=50.)
            #from matplotlib import pyplot
            #print i,nrings
            #pyplot.plot(smthBkgd,'-')
            #pyplot.show()

            maxvals[i] = np.max( np.abs(smthBkgd ))
            if (maxvals[i] > peak_cutoff):
                #print 'FLAG SPIKE'

                #print 'DATA SIZE:',len(data[lo:hi])

                
                #RFI_short = RFI[thisMask]
                spikes = np.where( np.abs(smthBkgd ) > peak_cutoff)[0]
                
                #Get position of each spike:
                spike_diff = spikes[1:-1]-spikes[0:-2]
                spike_pos = np.where(spike_diff > 1)[0] + 1

                #RFI_short[np.max([spikes[0]-60,0.]):np.min([spikes[0]+60,len(thisBkgd)])] = False
                RFI[np.max([thisMask[spikes[0]]-60,0]):np.min([thisMask[spikes[0]]+60,RFI.size])] = False
                #print np.max([thisMask[spikes[0]]-60,0]),np.min([thisMask[spikes[0]]+60,RFI.size])

                if spike_pos.size > 0:
                    for spike in spikes[spike_pos]:
                        #RFI_short[np.max([spike-60,0.]):np.min([spike+60,len(thisBkgd)])] = False
                        #print np.max([thisMask[spike]-60,0]),np.min([thisMask[spike]+60,RFI.size])
                        RFI[np.max([thisMask[spike]-60,0]):np.min([thisMask[spike]+60,RFI.size])] = False
                #RFI[thisMask] = RFI_short
                RFIMask[lo:hi] = RFI

                #pyplot.plot(np.abs(smthBkgd),'-')
                #pyplot.plot(RFI_short,'-')
                                
                #pyplot.show()


                #print 'YOY OY YO'
                #from matplotlib import pyplot
                #print peak_cutoff
                #pyplot.plot(np.abs(smthBkgd),'-')
                #pyplot.plot(RFI_short,'-')
                #pyplot.figure()
                #pyplot.plot(data[lo:hi])
                #pyplot.plot(RFI)#RFIMask[lo:hi])
                #pyplot.show()

                SpikeMask[lo:hi] = (RFI == False)


            hispikes = np.where(data[lo:hi] >= 1.)[0]
            RFIMask[lo+hispikes] = False

            #thisMask = notforemask[lastStep:dhi]
            #thisBkgd = (data[lastStep:dhi])[thisMask]                
            #Filter ring?

            stds[i] = np.std(thisBkgd[RFI[thisMask]])

            if (stds[i] > std_cutoff) | (stds[i] == 0):
                RFIMask[lo:hi] = False
                #SpikeMask[lo:hi] = False

        
    #print 'SIZES:',data.size,RFIMask.size,std_cutoff,peak_cutoff
    #from matplotlib import pyplot
    #pyplot.plot(stds,'o')
    #pyplot.figure()

    #pyplot.plot(maxvals,'o')
    #pyplot.figure()
    #pyplot.plot(ijd[mask & foremask],data[mask & foremask],',')
    #pyplot.plot(ijd[RFIMask & mask & foremask],data[RFIMask & mask & foremask],',')
    #pyplot.plot(jd,data,',')#jd[mask & RFIMask],    
    #pyplot.plot(jd[mask & RFIMask],data[mask & RFIMask],',')#jd[mask & RFIMask],
    #pyplot.plot(jd[(foremask == False)],data[(foremask == False)],'o',alpha=0.7)    

    #pyplot.plot(ijd[mask & foremask],RFIMask[mask & foremask],',')

    #pyplot.show()

    return SpikeMask,RFIMask


'''
    #from matplotlib import pyplot
    #pyplot.plot(SpikeMask)
'''
