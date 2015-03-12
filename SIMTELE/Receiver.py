#Define the receiver noise of a telescope
import numpy as np

import scipy.fftpack as sfft
from matplotlib import pyplot
import matplotlib.ticker as ticker

def Noise(obslist,recinfo):
    '''
    Generate simulated 1/f noise for receivers
    '''

    NSamplesPerScan = int(recinfo['SampRate']/obslist['Speed'] * obslist['DSlew'])
    NSamples = int(NSamplesPerScan * obslist['NStep'])
    bit2 = int(2**np.ceil(np.log10(NSamples)/np.log10(2)))

    
    sig = recinfo['sigma']
    fknee = recinfo['fknee']
    alpha = recinfo['Alpha']

    #FFT of white noise signal:
    inputsig = np.random.normal(size=NSamples,scale=sig)
    whitenoise = sfft.fft(inputsig,n=bit2)

    #Generate temporal frequencies:
    nu = (np.arange(bit2)+1.) * (recinfo['SampRate']/float(NSamples))
    #nu = (np.arange(bit2/2)+1.) * (recinfo['SampRate']/float(NSamples/2.))

    #1/f Power spectrum:
    P = 1. + (fknee/nu)**alpha
    P2 = 1. + (1./nu)**-2

    #P = np.concatenate((P[-5:], (P[::-1])[:P.size-5],P))

    #pyplot.plot(P,',')
    #pyplot.show()
    #Generate 1/f noise:

    print P.shape,whitenoise.shape
    #noise = np.real( sfft.ifft( whitenoise * np.sqrt(P) ) )
    noise = np.real( sfft.ifft( whitenoise * np.sqrt(P) ,n=bit2) )
    noise = noise[:NSamples]

    
    noise2 = np.real( sfft.ifft( whitenoise * np.sqrt(P2) ,n=bit2) )
    noise2 = noise2[:NSamples]


    gd = np.where((nu < 50.) & (nu > 1./280.))[0]
    gd2 = np.where((nu < 50.) & (nu > 1./280.))[0]
    
    print np.sum(P2[gd2]-1.) ,np.sum(P[gd]-1.) ,np.sum(P[gd]-1.)/np.sum(P2[gd2]-1.) 

    print np.std(noise2[0:14000])/sig,np.std(noise[:14000])/sig,(np.std(noise[:14000])/np.std(noise2[0:14000]))**2
    print stop
    fig = pyplot.figure()



    fig.add_subplot(211)
    print noise.shape,inputsig.shape
    pyplot.plot(noise[0:NSamples],color='#0059F2')
    pyplot.plot(noise[:NSamples]-inputsig[:NSamples],color='#F29900')
    
    pyplot.xlim(0,14000)
    pyplot.xticks([])
    pyplot.yticks([-30,0,30],size=14)
    pyplot.ylim(-50,50)
    pyplot.ylabel('mK',size=16)

    fig.add_subplot(212)
    fig.subplots_adjust(hspace=0)
    freqs = nu#np.append(-nu[::-1],nu)#sfft.fftfreq(bit2,d=1./recinfo['SampRate']/2.)
    minfreq = recinfo['SampRate']/float(bit2)
    gd = np.where(freqs > minfreq)[0]
    pyplot.plot(freqs[gd], np.abs(whitenoise[gd])**2 * P[gd],color='#0059F2')
    pyplot.ylim(1,180000000000000.0)

    ylims = pyplot.gca().get_ylim()
    line, = pyplot.plot([minfreq,recinfo['SampRate']],[np.median(np.abs(whitenoise[gd]))**2,np.median(np.abs(whitenoise[gd]))**2],'-',color='#ADADAD',lw=2)
    line.set_dashes([8,4,2,4,2,4])

    
    pyplot.plot(freqs[gd],np.median(np.abs(whitenoise))**2 *(fknee/freqs[gd])**alpha,'--',color='#ADADAD',lw=2)
    
    pyplot.plot([recinfo['fknee'],recinfo['fknee']],[ylims[0],ylims[1]],':',color='#ADADAD',lw=2)

    print ylims
    pyplot.ylim(ylims[0],ylims[1])

    pyplot.plot(freqs[gd], np.median(np.abs(whitenoise[gd]))**2 * P[gd],'--',lw=3,color='#F29900')
    
    pyplot.semilogx()
    pyplot.semilogy()

    #pyplot.yticks([1e0,1e2,1e4,1e6,1e8,1e10,1e12,1e14])
    #pyplot.gca().xaxis.set_major_formatter(ticker.ScalarFormatter())
    #pyplot.gca().yaxis.set_major_formatter(ticker.ScalarFormatter())    
    pyplot.xlim(np.min(freqs[gd]),recinfo['SampRate'])

    pyplot.xlabel(r'$\nu$ (Hz)',size=16)
    pyplot.ylabel(r'P($\nu$) (mK$^2$)',size=16)
    pyplot.show()

    return noise


class Receiver:

    def __init__(self,obslist,recinfo):
        '''
        Contains the receiver information/noise
        '''

        #Generate noise:
        self.noise = Noise(obslist,recinfo)        
        self.noisespikes = NoiseWithSpikes(obslist,recinfo)
        self.NSampls = self.noise.size



def NoiseWithSpikes(obslist,recinfo):
    '''
    Generate simulated 1/f noise for receivers
    '''

    NSamplesPerScan = int(recinfo['SampRate']/obslist['Speed'] * obslist['DSlew'])
    NSamples = int(NSamplesPerScan * obslist['NStep'])
    bit2 = int(2**np.ceil(np.log10(NSamples)/np.log10(2)))

    
    sig = recinfo['sigma']
    fknee = recinfo['fknee']
    alpha = recinfo['Alpha']

    #FFT of white noise signal:
    white = np.random.normal(size=NSamples,scale=sig)
    whitenoise = sfft.fft(white,n=bit2)

    #Generate temporal frequencies:
    nu = (np.arange(bit2/2)+1.) * (recinfo['SampRate']/float(bit2))

    #1/f Power spectrum:
    P = 1. + (fknee/nu)**alpha

    #Add some noise spikes:
    spike2Hz = np.where((nu > 1.95) & (nu < 2.1))[0]
    #P[spike2Hz[0]] += np.mean(whitenoise)**2*200000.    
    spike2Hz = np.where((nu > 1.9) & (nu < 1.95))[0]
    #P[spike2Hz[0]] += np.mean(whitenoise)**2*10000.    
    spike2Hz = np.where((nu > 2.13) & (nu < 2.3))[0]
    #P[spike2Hz[0]] += np.mean(whitenoise)**2*10000.

    spike2Hz = np.where((nu > 1.8) & (nu < 1.85))[0]
    #P[spike2Hz[0]] += np.mean(whitenoise)**2*10000.    
    spike2Hz = np.where((nu > 2.2) & (nu < 2.3))[0]
    #P[spike2Hz[0]] += np.mean(whitenoise)**2*10000.    

    
    spike4Hz = np.where((nu > 3.95) & (nu < 4.1))[0]
    #P[spike4Hz[0]] += np.mean(whitenoise)**2*20000.


    #log-normal
    lnorm = lambda P,X: P[0]*np.exp(- ((np.log10(X)-P[1])/P[2])**2)

    Spike4Hz  = lnorm([np.mean(np.abs(whitenoise)**2)*200./1e7,np.log10(4.),2e-6],nu)
    Spike4Hz += lnorm([np.mean(np.abs(whitenoise)**2)*100./1e8,np.log10(4.),5e-5],nu)    
    Spike4Hz += lnorm([np.mean(np.abs(whitenoise)**2)*100./1e7,np.log10(4.05),2e-6],nu) 
    Spike4Hz += lnorm([np.mean(np.abs(whitenoise)**2)*100./1e7,np.log10(3.95),2e-6],nu) 
    Spike4Hz += lnorm([np.mean(np.abs(whitenoise)**2)*10./1e7,np.log10(4.1),2e-6],nu) 
    Spike4Hz += lnorm([np.mean(np.abs(whitenoise)**2)*10./1e7,np.log10(3.9),2e-6],nu) 

    P += Spike4Hz

    Spike7Hz  = lnorm([np.mean(np.abs(whitenoise)**2)*200./1e7,np.log10(7.),2e-6],nu)
    Spike7Hz += lnorm([np.mean(np.abs(whitenoise)**2)*100./1e8,np.log10(7.),5e-5],nu)    
    Spike7Hz += lnorm([np.mean(np.abs(whitenoise)**2)*100./1e7,np.log10(7.05),2e-6],nu) 
    Spike7Hz += lnorm([np.mean(np.abs(whitenoise)**2)*100./1e7,np.log10(6.95),2e-6],nu) 
    Spike7Hz += lnorm([np.mean(np.abs(whitenoise)**2)*10./1e7,np.log10(7.1),2e-6],nu) 
    Spike7Hz += lnorm([np.mean(np.abs(whitenoise)**2)*10./1e7,np.log10(6.9),2e-6],nu) 
    P += Spike7Hz
    


    Spike11Hz = lnorm([np.mean(np.abs(whitenoise)**2)/2000000.,np.log10(11.),1e-2],nu) 
    P += Spike11Hz

    P = np.append(P[::-1],P)
    nu = np.append(-nu[::-1],nu)


    #Generate 1/f noise:
    noise = np.real( sfft.ifft(  whitenoise*np.sqrt(P) ) ) 
    noise = noise[:NSamples]- white

    #freqs = nu
    #minfreq = recinfo['SampRate']/float(NSamples)
    #gd = np.where(freqs > minfreq)[0]
    #from matplotlib import pyplot
    #pyplot.plot(freqs[gd], np.abs(whitenoise[gd])**2*P[gd],color='#0059F2')#
    #pyplot.ylim(1,180000000000000.0)
    #pyplot.semilogx()
    #pyplot.semilogy()
    #pyplot.show()

    return noise
