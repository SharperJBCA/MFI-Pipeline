#DiodePol.py
# Goals:
# 1) Read in Cal Diode signal at many different polarisations.
# 2) Determine polarisation fraction of each MFI channel.

import numpy as np
import pyfits
from DataAccess import ReadObs
from DataAccess import TextFiles
from matplotlib import pyplot
import CalFitting
import WaveFitter
import Binning

def SelectNumbers(a,n):
    '''
    Return a selection of random numbers from a -bag- of numbers

    Arguments
    a -- input array
    n -- how many numbers to select
    '''

    bag = np.arange(len(a),dtype='i')
    num = np.zeros(n,dtype='i')

    mask = np.ones(len(a),dtype='bool')

    for i in range(n):

        ind = int(np.random.uniform(low=0,high=len(a)-i ))

        num[i] = (bag[mask])[ind]

        mask[num[i]] = False

    return num
        
        
    

def JackKnifeErrors(mod,amps,niter=10):
    '''
    '''

    PolFracs = np.zeros(niter)
    Phases   = np.zeros(niter)

    for i in range(niter):

        #Randomly select indices:
        indexs = SelectNumbers(amps,len(amps)/2)
        
        P = WaveFitter.FFTMethod(mod[indexs],amps[indexs])

        PolFracs[i] = P[1]/P[0]
        Phases[i] = P[2]
        
    return np.std(PolFracs)*np.sqrt(len(amps)),np.std(Phases)*np.sqrt(len(amps))

if __name__ == '__main__':

    filelist = np.loadtxt('TESTABCD05.lis',dtype='string')
    data = ReadObs.FullObs(filelist,['DATA','CAL','JD','MOD_ANGLE'],dir='/nas/scratch/sharper/QUIJOTE/datos/tod/')

    #model = CalFitting.CalDriftModel(data['DATA'],data['CAL'])
    model = np.array([-0.001542, 0.6685])
    
    
    amps,bkgds,mjd,mod = CalFitting.AvgCalSig(data['DATA'],data['CAL'],jd=data['JD'],mod=data['MOD_ANGLE'])

    pairs = [[0,6],[1,7],[2,4],[3,5]]

    for i in range(16,amps.shape[2]):
    #for ch in pairs:

        #Remove any bad diode voltages
        gd = (amps[0,:,i] != 0)

        #Define the modulator angles:
        m = 4.*np.mod(mod[0,gd,(i)/8],90)
        a = amps[0,gd,i]#/amps[0,gd,16+ch[1]] #- 1.)/2. + 1.
        jd= mjd[0,gd,0]

        ranges = [[355.,5.],[85.,95.],[175.,185.],[265.,275.]]
        grad = 0.
        for rad in ranges:
            #print rad
            if rad[0] > rad[1]:
                gd = np.where((m > rad[0]) | (m < rad[1]))[0]
            else:
                gd = np.where((m > rad[0]) & (m < rad[1]))[0]
                
            pfit = np.poly1d(np.polyfit(jd[gd],a[gd],1.))
            grad += pfit[1]/4.
            #print pfit[1]

        pfit = np.poly1d([grad,0.])
        pyplot.plot(a-np.mean(a),',')
        a = a/pfit(jd)
        pyplot.plot(a-np.mean(a),',')
        pyplot.show()
        
        #print grad
        
        #b = bkgds[0,gd,i]

        #Fit to the wave:
        P = WaveFitter.FFTMethod(m,a,wave=360.)
        P2 = WaveFitter.FFTMethod(m,a-WaveFitter.FitSine(P,m),wave=180.)

        
        polerr, phserr = JackKnifeErrors(m,a)
        n1 = '%1.3f' % float(P[1]/P[0]*100.)
        e1 = '%1.5f' % float(polerr/2.*100.)
                                 
        n2 = '%1.3f' % float(P[2])
        e2 = '%1.5f' % float(phserr)

        n1_2 = '%1.3f' % float(P2[1]/2.*100.)                                 
        n2_2 = '%1.3f' % float(P2[2])

        print n1,'&',e1 , '&',n2,'&',e2,'&',n1_2,'&',n2_2, '\\\\'
        #float(P[1]/P[0]),float(P[2])
        #'&',n1 , '&',n2, '\\\\'

                

        #pyplot.plot(m,a-WaveFitter.FitSine(P,m),'o')#-WaveFitter.FitSine(P2,m)
        #pyplot.plot(m,a2/np.mean(a2),'o')        
        #pyplot.plot(m,WaveFitter.FitSine(P2,m),'or',linewidth=3,label='Best Fit')        
        #pyplot.ylabel(r'$\frac{V}{\left< V \right> }$',size='large')
        #pyplot.xlabel(r'4$\phi$',size='large')
        #pyplot.show()


        #TextFiles.AppendFile('DiodePolParams_Epoch1.dat',[i,P[1]/P[0],polerr,P[2],phserr])

        #Plotting:
        #md = np.linspace(0,360,500)
        #pyplot.plot(mbig,1000.*dbig ,',',label=' Voltages')
        #pyplot.plot(md,1000.*WaveFitter.FitSine(Pd,md),'--r',linewidth=3,label='Best Fit')
        #pyplot.plot(m,1000.*a ,'o',label='Diode Voltages')
        #pyplot.plot(md,1000.*WaveFitter.FitSine(P,md),'--r',linewidth=3,label='Best Fit')        
        #pyplot.xlabel(r'4$\phi$')
        #pyplot.ylabel(r'Voltage (mV)')
        #pyplot.legend(numpoints=1)
        #pyplot.show()
        #pyplot.plot(m,1000.*a ,'o',label='Diode Voltages')
        #pyplot.plot(md,1000.*WaveFitter.FitSine(P,md),'--r',linewidth=3,label='Best Fit')

'''
        pyplot.plot(m,bkgds[0,gd,i]/np.median(bkgds[0,gd,i]),'o')
        pyplot.ylabel(r'$\frac{V}{\left< V \right> }$',size='large')
        pyplot.xlabel(r'4$\phi$',size='large')
        pyplot.figure()
        pyplot.plot(m,bkgds[0,gd,i+2]/np.median(bkgds[0,gd,i+2]),'o')
        pyplot.ylabel(r'$\frac{V}{\left< V \right> }$',size='large')
        pyplot.xlabel(r'4$\phi$'  ,size='large')                    
        pyplot.figure()
        pyplot.plot(m,(bkgds[0,gd,i+2]+bkgds[0,gd,i+4])/np.median(bkgds[0,gd,i+2]+bkgds[0,gd,i+4]),'o')
        pyplot.ylabel(r'$\frac{V}{\left< V \right> }$',size='large')
        pyplot.xlabel(r'4$\phi$',size='large')
        pyplot.figure()
        pyplot.plot(m,(bkgds[0,gd,i]+bkgds[0,gd,i+6])/np.median(bkgds[0,gd,i]+bkgds[0,gd,i+6]),'o')        
        pyplot.ylabel(r'$\frac{V}{\left< V \right> }$',size='large')
        pyplot.xlabel(r'4$\phi$',size='large')
        pyplot.show()
'''
