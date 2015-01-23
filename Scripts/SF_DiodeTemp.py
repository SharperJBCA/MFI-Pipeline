#DiodeTemp.py
# Goals
# 1) Determine the Calibration Diode temperature from calibrator observations


import numpy as np
import pyfits
from DataAccess import ReadObs
from DataAccess import FileList
from DataAccess import TextFiles
from matplotlib import pyplot
import CalFitting
import WaveFitter
import Binning

import SourceFitting
from SourceFitting.ModelFlux import toJy


if __name__ == '__main__':

    n = np.loadtxt('TauA',dtype='string')
    jd0 = 56244.

    print n.shape

    beam13 = (0.66*np.pi/180.)**2#(0.88*np.pi/180.)**2#
    beam11 = (0.66*np.pi/180.)**2#(0.86*np.pi/180.)**2#
    f1 = 19.
    f2 = 17.
    h = 2

    gd = ((n[:,5].astype('f') != 0.))

    CasFlux13GHz = SourceFitting.ModelFlux.TauAFlux(f1,n[gd,26].astype('f')+jd0)
    CasFlux11GHz = SourceFitting.ModelFlux.TauAFlux(f2,n[gd,26].astype('f')+jd0)

    #pyplot.hist(CasFlux13GHz/n[:,1].astype('f')/ toJy(f1,beam13))
    #pyplot.hist(CasFlux13GHz/n[:,3].astype('f')/ toJy(f2,beam13),alpha=0.6)

    #Calculate stability of correlated channels
    cTd11 = np.median(CasFlux11GHz/n[gd,2+h*4].astype('f')/ toJy(f2,beam11))
    stdcTd11 = np.std(CasFlux11GHz/n[gd,2+h*4].astype('f')/ toJy(f2,beam11))    
    cTd13 = np.median(CasFlux13GHz/n[gd,1+h*4].astype('f')/ toJy(f1,beam13))
    stdcTd13 = np.std(CasFlux13GHz/n[gd,1+h*4].astype('f')/ toJy(f1,beam13))

    print '11GHz: ',cTd11,stdcTd11,stdcTd11/cTd11 * 100.
    print '13GHz: ',cTd13,stdcTd13,stdcTd13/cTd13 * 100.
    print np.log10(cTd11/cTd13)/np.log10(f2/f1)

    #Calculate stability of uncorrelated channels
    uTd11 = np.median(CasFlux11GHz/n[gd,4+h*4].astype('f')/ toJy(f2,beam11))
    stduTd11 = np.std(CasFlux11GHz/n[gd,4+h*4].astype('f')/ toJy(f2,beam11))    
    uTd13 = np.median(CasFlux13GHz/n[gd,3+h*4].astype('f')/ toJy(f1,beam13))
    stduTd13 = np.std(CasFlux13GHz/n[gd,3+h*4].astype('f')/ toJy(f1,beam13))

    print '11GHz: ',uTd11,stduTd11,stduTd11/uTd11 * 100.
    print '13GHz: ',uTd13,stduTd13,stduTd13/uTd13 * 100.
    print np.log10(uTd11/uTd13)/np.log10(f2/f1)

    #b = CasFlux13GHz/n[gd,3+h*4].astype('f')/ toJy(f1,beam13)
    #bd = ( (b > 2))
    #print n[bd,0]
    #pyplot.plot(n[gd,26].astype('f')+jd0,n[gd,3+h*4].astype('f'),'o')    
    #pyplot.plot(CasFlux13GHz/n[gd,7].astype('f')/ toJy(f1,beam13),'o')
    #pyplot.show()

    pyplot.plot(n[gd,26].astype('f')+jd0,CasFlux13GHz/n[gd,1+h*4].astype('f')/ toJy(f1,beam13),'o',label='Correlated 19GHz',alpha=0.8)
    pyplot.plot(n[gd,26].astype('f')+jd0,CasFlux11GHz/n[gd,4+h*4].astype('f')/ toJy(f2,beam11),'o',label='Uncorrelated 17GHz',alpha=0.8)
    pyplot.plot(n[gd,26].astype('f')+jd0,CasFlux13GHz/n[gd,3+h*4].astype('f')/ toJy(f1,beam13),'o',label='Uncorrelated 19GHz',alpha=0.8)
    pyplot.plot(n[gd,26].astype('f')+jd0,CasFlux11GHz/n[gd,2+h*4].astype('f')/ toJy(f2,beam11),'o',label='Correlated 17GHz',alpha=0.8)

    ax = pyplot.gca()
    ylims = ax.get_ylim()
    pyplot.plot([56757,56757],[0.,ylims[1]*2.],'--k',linewidth=3)
    pyplot.ylim(ylims[0],ylims[1])

    pyplot.xlabel(r'MJD')
    pyplot.ylabel(r'Antenna Temperature (K)')
    pyplot.legend(numpoints=1)

    pyplot.figure()
    pyplot.plot(n[gd,26].astype('f')+jd0,n[gd,1+h*4+12].astype('f'),'o',label='Correlated 19GHz')
    pyplot.plot(n[gd,26].astype('f')+jd0,n[gd,4+h*4+12].astype('f'),'o',label='Uncorrelated 17GHz')
    pyplot.plot(n[gd,26].astype('f')+jd0,n[gd,3+h*4+12].astype('f'),'o',label='Uncorrelated 19GHz')
    pyplot.plot(n[gd,26].astype('f')+jd0,n[gd,2+h*4+12].astype('f'),'o',label='Correlated 17GHz')

    ax = pyplot.gca()
    ylims = ax.get_ylim()
    pyplot.plot([56757,56757],[0.,ylims[1]*2.],'--k',linewidth=3)
    pyplot.ylim(ylims[0],ylims[1])
    
    pyplot.xlabel(r'MJD')
    pyplot.ylabel(r'r-factor',style='italic')
    pyplot.legend(numpoints=1)    
    pyplot.show()
