import numpy as np
from matplotlib import pyplot
import MyMaths

import matplotlib as mpl

from scipy import optimize

def gradfit(P,X):

    return P[0]*X+ P[1]

def errfunc(P,X,Y,e):

    return (gradfit(P,X) - Y)/e



#Klein Gulkis 1978 fluxs
freqs = np.array([20.05,20.295,20.735,20.970,21.225,21.55,22.055,22.22,22.37,22.945,23.41,23.84,24.1,
                  4.62,7.9 , 9.9,15.4,15.8,16.3,19.,
                  9.53,9.1,8.0 ,6.65,5.,3.])
fluxs = np.array([72.27,74.33,76.65,77.82,81.35,81.58,85.31,83.85,84.6,89.92,93.88,98.99,100.56,
                  6.2 ,15.0,28.8,54.1,54.6,61.0,64.,
                  21.9,20.,18.5,15.5,10.7,7.77])


errs = np.sqrt((np.repeat(np.array([0.5]),len(freqs)))**2 + (0.01*fluxs)**2)
errs[-13:] = fluxs[-13:]*0.05

print len(fluxs),len(freqs)
#Jupiter solid angle
solid = 2.481e-8 #sr At 5.2 AU




mfidir = '/nas/scratch/sharper/QUIJOTE/MFI-Pipeline/MFI_INFO/'

beams = np.loadtxt(mfidir+'BeamSizes.dat')    
beams = np.array([1.133*(beams[0,1]*np.pi/180.)**2,
                  1.133*(beams[1,1]*np.pi/180.)**2,
                  1.133*(beams[2,1]*np.pi/180.)**2,
                  1.133*(beams[3,1]*np.pi/180.)**2]) #Beams of H1 to H4
                 

qfluxs = np.array([49.89706784,42.71858291,32.4,26.5]) #Jy
qerrs  = np.array([2.25110345,1.19874893,0.8250022,0.7425473])*np.sqrt(41.) #Jy
qfreqs = np.array([18.7,16.8,12.88,11.16]) #GHz
qerrs  = np.sqrt(qerrs**2 + (qfluxs*0.04)**2)


P0 = [0.,0.]
P1,s = optimize.leastsq(errfunc,P0,args=(np.log10(freqs),np.log10(fluxs),errs ))
print np.log10(errs**2) 
print P1
plotfreqs = np.linspace(1,30,100)

pyplot.plot(plotfreqs,10**gradfit(P1,np.log10(plotfreqs)),'--k',lw=2)
pyplot.errorbar(qfreqs,qfluxs,yerr=qerrs,fmt='s',color='#F22000',label='QUIJOTE',markersize=10)
pyplot.errorbar(freqs,fluxs,yerr=errs,fmt='o',color='#ADADAD',label='Ancillary')
pyplot.semilogx()
pyplot.semilogy()
pyplot.gca().xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
pyplot.gca().yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
pyplot.xlim(1,30)
pyplot.ylim(1,100)
pyplot.legend(frameon=False,numpoints=1,loc='upper left')
pyplot.xticks(size=14)
pyplot.yticks(size=14)
pyplot.ylabel('Flux Density (Jy)',size=16)
pyplot.xlabel('Frequency (GHz)',size=16)
pyplot.show()

qdist  = 4.04 #AU

#*MyMaths.toJy(qfreqs*1.,np.array([beams[1],beams[1],beams[2],beams[2]]))/1000.
print beams
#qfluxs = qfluxs * (qdist/4.04)**2# (5.2/qdist)**2
print freqs,qfreqs

fluxs  *= (4.04/5.2)**2
qfluxs *= (4.04/5.2)**2

Tb     = fluxs/MyMaths.toJy(freqs*1.,np.repeat(np.array([solid]),len(freqs)))
qTb    = qfluxs/MyMaths.toJy(qfreqs*1.,np.array([solid,solid,solid,solid]))  
qerrTb = qerrs/MyMaths.toJy(qfreqs*1.,np.array([solid,solid,solid,solid]))

gd = np.where(freqs < 22.5)[0]

pfit = np.poly1d(np.polyfit(np.log10(freqs[gd ]),np.log10(fluxs[gd ]),1))#[[0,1,2,3,4,6,7]] 



plotfreqs = np.linspace(1,30,100)

fitTb = (10**pfit(np.log10(plotfreqs)))/MyMaths.toJy(plotfreqs*1.,np.repeat(np.array([solid]),len(plotfreqs)))

nu = qfluxs/10**pfit(np.log10(qfreqs))
print nu,pfit

#pyplot.plot(qfreqs,nu,'o')
#pyplot.show()

#pyplot.plot(freqs,Tb* (4.04/5.2)**2,'o')
#pyplot.plot(qfreqs,qTb* (4.04/5.2)**2/nu,'o')

pyplot.plot(plotfreqs,fitTb,'--')
pyplot.errorbar(qfreqs,qTb*1.133,fmt='o',yerr=qerrs)
pyplot.plot(freqs,Tb,'o')

pyplot.show()
