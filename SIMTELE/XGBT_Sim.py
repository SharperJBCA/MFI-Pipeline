import numpy as np
import pyfits
import mlmapper
import nBinning as Binning
from matplotlib import pyplot

import DataAccess
import Control

from Telescope import Telescope
import Sky

Observation = {'RA':202.46841,
               'DEC':47.194814,
               'OBSMODE':'RALONMAP',
               'MJD': 51545.3,
               'DSlew':0.15,
               'DStep':0.25/60.,
               'NStep':36,
               'Speed':10./60./60.,
               'LAT':38. + 25./60. + 59./60.**2,
               'LNG':79. + 50./60. + 23./60.**2.}

horns = [[0.,0.]]

FocalPlane = [{'xpos':horns[i][0],
               'ypos':horns[i][1],
               'Pf':0,
               'Px':0,
               'Py':0,
               'Pc':0,
               'Pn':0,
               'Pa':0,
               'Pb':0} for i in range(len(horns))]

#3.8,
Receiver = [{'sigma':2.2,
             'fknee':1.,
             'SampRate':5.,
             'Alpha':1.25} for i in range(len(horns))]

import MyMaths
hdu = pyfits.open('m51_xband_ds.fits')
hdu[0].data = hdu[0].data/MyMaths.toJy(10.,(1.4/60.*np.pi/180.)**2)*0.

print np.max(hdu[0].data),MyMaths.toJy(10.,(1.4/60.*np.pi/180.)**2)

#Setup data containers:
NObservations = 2
NSamplesPerScan = int(Receiver[0]['SampRate']/Observation['Speed'] * Observation['DSlew'])
print 'NSAMPLES PER SCAN:', NSamplesPerScan
NSamples = int(NSamplesPerScan * Observation['NStep'])

tod = np.zeros(NSamples*NObservations*len(Receiver))
pix = np.zeros(NSamples*NObservations*len(Receiver))
ra  = np.zeros((NSamples*NObservations,len(Receiver)))
dec = np.zeros((NSamples*NObservations,len(Receiver)))

print 'TIME:' ,NSamples /Receiver[0]['SampRate'] / 60.,NSamples*NObservations


Tele = Telescope(FocalPlaneInfo=FocalPlane, ReceiverInfo=Receiver,SkyMap=hdu)

for i in range(NObservations):
    if np.mod(i,2) == 0:
        Observation['OBSMODE'] = 'RALONMAP'
    else:
        Observation['OBSMODE'] = 'DECLATMAP'
            
    Tele.Observation(Observation)

    for j in range(len(Tele.Receivers)):

        np.random.seed(seed=102313 + i*len(Tele.Receivers) + j)
        
        tod[i*NSamples*len(Tele.Receivers) + j*NSamples:i*NSamples*len(Tele.Receivers) + (j+1)*NSamples] = Tele.Skies[j].signal + Tele.Receivers[j].noise
        
        if i == 0:
            ra[ i*NSamples:(i+1)*NSamples,j] = Tele.Scans[j].ra
            dec[i*NSamples:(i+1)*NSamples,j] = Tele.Scans[j].dec
            pix[i*NSamples*len(Tele.Receivers) + j*NSamples:i*NSamples*len(Tele.Receivers) + (j+1)*NSamples] = Tele.Skies[j].pix            
        else:
            rot = 0. * np.pi/180.
            rt = Tele.Scans[j].ra  - Observation['RA']
            dt = Tele.Scans[j].dec - Observation['DEC']
            ra[ i*NSamples:(i+1)*NSamples,j] =  rt*np.cos(rot) + dt*np.sin(rot) + Observation['RA']
            dec[i*NSamples:(i+1)*NSamples,j] = -rt*np.sin(rot) + dt*np.cos(rot) + Observation['DEC']

            pt,x,y = Sky.GetPixlistSqr(Tele.Skies[j].WCS,dec[ i*NSamples:(i+1)*NSamples,j],ra[i*NSamples:(i+1)*NSamples,j],hdu)
            
            pix[i*NSamples*len(Tele.Receivers) + j*NSamples:i*NSamples*len(Tele.Receivers) + (j+1)*NSamples] = pt


        
            

    Observation['MJD'] +=  0.045 #Move the sky a little bit

#Remove any pixels outside range:
gd = np.where((pix >= 0) & (pix < hdu[0].data.shape[0]*hdu[0].data.shape[1]))[0]

print 'hj',tod.size-gd.size
pyplot.plot(ra[:,0],dec[:,0],'-',color='#0059F2')
pyplot.xlabel('Right Ascension')
pyplot.ylabel('Declination')
pyplot.show()

tod = tod[gd]
pix = pix[gd]


s,h,sig,jk = mlmapper.mlmapper(tod,pix,hdu[0].data.shape[0]*hdu[0].data.shape[1],1,100)
sw,hw,h= Control.Destriper(tod,135,pix,hdu[0].data.shape[0]*hdu[0].data.shape[1],bl_long=270,Verbose=True,maxiter=800)

sd = sw*0.
sd[hw != 0] = sw[hw != 0]/hw[hw != 0]

s2 = np.zeros(hdu[0].data.shape[0]*hdu[0].data.shape[1])
h = np.zeros(hdu[0].data.shape[0]*hdu[0].data.shape[1])
cn = tod*0. + 1.
Binning.BinMap(tod,1,pix,cn,s2,hits = h)

bad = (s == 0)
good = (s != 0)
bot = np.median(s[good])
s -= bot
s[bad] = np.nan
h[bad] = np.nan
sd[bad] = np.nan
s = np.reshape(s,hdu[0].data.shape)
sd = np.reshape(sd,hdu[0].data.shape)
h = np.reshape(h,hdu[0].data.shape)
m = s

print hdu[0].data.shape
print m.shape
hdu[0].data = s
#hdu.verify('fix')
hdu.writeto('RASTERSIMULATIONS/XBand-Knee1Hz-2Maps-Rot0.fits',clobber=True)
hdu[0].data = sd
#hdu.verify('fix')
hdu.writeto('RASTERSIMULATIONS/XBand-Knee1Hz-2Maps-Destriper-Rot0.fits',clobber=True)


#for i in range(len(Tele.Skies)):
#    pyplot.plot(ra[:,i],dec[:,i],'-',alpha=0.7,label='Horn %i' % i)

#pyplot.xlabel('Right Ascension')
#pyplot.ylabel('Declination')
#pyplot.legend()

#pyplot.show()

