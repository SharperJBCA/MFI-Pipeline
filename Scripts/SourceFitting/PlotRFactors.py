import numpy as np
import pyfits
from DataAccess import ReadObs
from DataAccess import FileList
from DataAccess import TextFiles
from matplotlib import pyplot
import CalFitting
import WaveFitter
import Binning
import DataAccess

import SourceFitting
from SourceFitting.ModelFlux import toJy

if __name__ == '__main__':

    casa = np.loadtxt('CasA_GaussianFits_PerChannel_All',dtype='string')

    casamps = casa[:,1:33].astype('f')
    dates = casa[:,-1].astype('f')


    chanpairs = [[0,6],[1,7],[2,4],[3,5]]

    horns = range(4)

    for h in horns:

        for chp in chanpairs:

            #yplot.plot(casamps[:,h*8 + chp[0]]/casamps[:,h*8 + chp[1]],'o')

            apRatio = casamps[:,h*8 + chp[0]]/casamps[:,h*8 + chp[1]]
            gd = np.where(np.abs(apRatio-np.median(apRatio)) < np.median(apRatio)*0.09)[0]
            c =  np.random.random(len(gd))

            print gd.size
            print h, chp,np.std(apRatio[gd])/np.mean(apRatio[gd])


            pyplot.scatter(dates[gd],apRatio[gd],c=c,cmap='gray',s=120,lw=0,label='Fit-to-Photometry Ratio',zorder=1)
            pyplot.xticks(size=14)
            pyplot.yticks(size=14)
            pyplot.ylabel(r'$\frac{G_x}{G_y}$',size=16)
            pyplot.xlabel(r'Observing Day',size=16)

            #pyplot.ylim(0.5,1.5)
            pyplot.show()
