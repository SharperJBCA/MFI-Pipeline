import numpy as np
import pyfits
from matplotlib import pyplot
import sys

fname = sys.argv[1]

hdu = pyfits.open(fname)

data = hdu[1].data

#pyplot.plot(data['DATA'][0,500:,3]*1000.,'-')#,color='#0059F2')
pyplot.plot(data['DATA'][0,500:2500,3]*1000.,'-',color='#0059F2')
#pyplot.plot(data['DATA'][0,:,3]*1000.,'-',color='#0059F2')
pyplot.xlabel('Sample',size=16)
pyplot.ylabel('Antenna Temperature (mK)',size=16)
pyplot.xticks(size=14)
pyplot.yticks(size=14)
pyplot.show()
