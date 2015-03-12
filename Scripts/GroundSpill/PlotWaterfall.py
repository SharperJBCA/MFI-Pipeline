import numpy as np
from matplotlib import pyplot
import jdcal
from skimage import exposure


d = np.loadtxt('H313uNom30_GSWaterfall.dat')

data = d[:,1:]
time = d[:,0]
angs = np.arange(360)

stime = time.argsort()

data = data[stime,:]
itime = np.round(time[stime]).astype('i')
time  = time[stime]

times = np.unique(itime)

for i in range(times.size):
    hits = np.where(itime == times[i])[0]

    day = np.mod(time[hits],1.)

    calday = jdcal.jd2jcal(2400000.,np.mean(time[hits]))

    dat = np.transpose(data[hits,:])
    datshape = dat.shape
    dat = np.reshape(dat,dat.shape[0]*dat.shape[1])
    bd = np.where(dat > 1)[0]
    gd = np.where(dat < 1)[0]
    
    dat -= np.min(dat)
    hist_dat = dat*1.
    hist_dat[gd] = exposure.equalize_hist(dat[gd],nbins=100)#np.reshape(dat,dat.shape[0]*dat.shape[1]),nbins=10000)
    hist_dat[bd] = np.nan
    dat[bd] = np.nan

    r = np.max(dat[gd]) - np.min(dat[gd])
    mind = np.min(dat[gd])


    hist_dat = np.reshape(hist_dat,datshape)
    dat = np.reshape(dat,datshape)

    
    print np.max(hist_dat),np.min(hist_dat)
    pyplot.imshow(dat)#,vmax=1.,vmin=-0.5,aspect='auto')
    cb = pyplot.colorbar(orientation='horizontal',pad=0.1)    
    pyplot.figure()
    pyplot.imshow(hist_dat*r + mind)#,vmax=1.,vmin=-0.5,aspect='auto')
    pyplot.title('%i/%i/%i' % (calday[::-1])[1:] )
    pyplot.ylabel('Azimuth')
    #pyplot.gca().xaxis_date()
    cb = pyplot.colorbar(orientation='horizontal',pad=0.1)
    cb.set_label('K')
    #pyplot.xticks([])
    pyplot.show()
