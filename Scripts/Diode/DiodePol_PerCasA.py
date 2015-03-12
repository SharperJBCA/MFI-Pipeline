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
        
        
    

def JackKnifeErrors(mod,amps,niter=300):
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

    gd = np.where(PolFracs < 1.)[0]
    #pyplot.hist(PolFracs,bins=np.arange(0.,1.+0.01,0.01))
    #pyplot.show()
    return np.std(PolFracs[gd])*np.sqrt(len(amps)),np.std(Phases[gd])*np.sqrt(len(amps))

if __name__ == '__main__':

    data = np.loadtxt('../SourceFitting/RFactors-NoFix',dtype='string')
    

    jd = data[:,-3].astype('f')

    pyplot.plot(jd,data[:,9+2].astype('f')/data[:,9+4].astype('f'),'o')
    pyplot.plot(jd,data[:,1+2].astype('f')/data[:,1+4].astype('f'),'o')    
    pyplot.show()

    gd = np.where(jd > 517.)[0]
    m = data[gd,-1].astype('f')


    nchans = 4
    pairs = [[0,6],[1,7],[2,4],[3,5]]
    #for i in range(nchans):
    for i in range(8):

        #Define the modulator angles:
        
        #a =data[gd,9+pairs[i][0]].astype('f')/data[gd,9+pairs[i][1]].astype('f')
        #a2=data[gd,9+pairs[i][0]].astype('f')/data[gd,1+pairs[i][0]].astype('f')# 
        a =data[gd,9+i].astype('f')/data[gd,1+i].astype('f')

        mask = np.ones(a.size,dtype='bool')
        for j in range(2):
            #Fit to the wave:
            P = WaveFitter.FFTMethod(m[mask],a[mask])
            a2 = a - WaveFitter.FitSine(P,m)

            pyplot.plot(np.abs(a2[mask]/np.mean(a[mask])),'o')
            pyplot.show()
            bd = np.where((np.abs(a2/np.mean(a[mask])) > 0.036))[0]
            mask[bd]= False
            

            
        polerr, phserr = JackKnifeErrors(m[mask],a[mask])
            
        n1 = '%1.3f' % float(P[1]/P[0]*100.)
        e1 = '%1.5f' % float(polerr*100.)
            
        n2 = '%1.3f' % float(P[2])
        e2 = '%1.5f' % float(phserr)
        

        print n1 ,'&',e1 , '&',n2,'&',e2 , '\\\\'


        pyplot.plot(m[mask],a[mask],'o')
        #pyplot.plot(m,a2/np.mean(a2),'o')
        md = np.linspace(0,360,500)
        pyplot.plot(md,WaveFitter.FitSine(P,md),'--r',linewidth=3,label='Best Fit')        
        pyplot.ylabel(r'$\frac{V}{\left< V \right> }$',size='large')
        pyplot.xlabel(r'4$\phi$',size='large')
        pyplot.show()

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
