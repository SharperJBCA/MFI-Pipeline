import numpy as np
import pyfits
#import mlmapper
from Polarisation.Tools import nBinning as Binning
import Binning as Bin
from matplotlib import pyplot
import healpy as hp

import DataAccess
#import Control

from Polarisation.Destriper import Control
#from MapMaker.MLMapper import Control

from Telescope import Telescope

files = ['Mod0.fits','Mod1.fits','Mod2.fits']
count = 0
for f in files:
    hdu = pyfits.open(f)
    count += hdu[1].data['QTOD'].shape[1]
    hdu.close()
    del hdu

qtod = np.zeros(count)
utod = np.zeros(count)
mod  = np.zeros(count)
pix  = np.zeros(count)
phi  = np.zeros(count)
noise= np.zeros(count)

lastLen = 0
for f in files:
    hdu = pyfits.open(f)
    fLen = hdu[1].data['QTOD'].shape[1] + lastLen

    qtod[lastLen:fLen] = hdu[1].data['QTOD'][0,:,0]*1.
    utod[lastLen:fLen] = hdu[1].data['UTOD'][0,:,0]*1.
    mod[lastLen:fLen]  = hdu[1].data['MOD'][0,:,0]*1.
    pix[lastLen:fLen]  = hdu[1].data['PIX'][0,:,0]*1.
    phi[lastLen:fLen]  = hdu[1].data['PHI'][0,:,0]*1.
    noise[lastLen:fLen]  = hdu[1].data['NOISE'][0,:,0] *1. 
    lastLen = fLen
    hdu.close()
    del hdu


ang = (2.*phi + 4.*mod)*np.pi/180.
tod = qtod*np.sin(ang) + utod*np.cos(ang) + noise
nside = 256

del phi, qtod, mod, utod, noise

Maps1    = Control.Destriper(tod,1500,pix,12*nside**2,ang,bl_long=3000,Verbose=True,maxiter=800)

class Maps:
    def __init__(self,npix):
        self.m = np.zeros(npix)
        self.q = np.zeros(npix)
        self.u = np.zeros(npix)

        self.c2 = np.zeros(npix)
        self.s2 = np.zeros(npix)
        self.sc = np.zeros(npix)

        self.vs = np.zeros(npix)
        self.vc = np.zeros(npix)

        self.DetA = np.zeros(npix)

        self.sw = np.zeros(npix)
        self.hw = np.zeros(npix)

nside = 256
Maps2 = Maps(12*nside**2)

cn = tod*0. +1.

Binning.BinMapPol_Angs(tod,1,pix.astype('i'),ang,cn,Maps2)
Binning.BinMapPol(tod,1,pix.astype('i'),ang,cn,Maps2)

Binning.BinMap(tod,1,pix.astype('i'),cn,Maps2.m)

#hp.write_map('DES_QNOM.fits',Maps1.q)
#hp.write_map('DES_UNOM.fits',Maps1.u)
#hp.write_map('QNOM.fits',qmap)
#hp.write_map('UNOM.fits',umap)

hp.mollview(np.sqrt(Maps1.q**2+Maps1.u**2)-np.sqrt(Maps2.q**2+Maps2.u**2),max=1,min=-1)
hp.mollview(np.sqrt(Maps1.q**2+Maps1.u**2),max=1,min=-1)
hp.mollview(np.sqrt(Maps2.q**2+Maps2.u**2),max=1,min=-1)
pyplot.show()

'''
bad = (h == 0)np.sqrt(Maps2.q**2+Maps2.u**2)
good = (h != 0)
bot = np.median(s[good])
s -= bot
s[bad] = hp.UNSEEN
s2[bad] = hp.UNSEEN
h[bad] = hp.UNSEEN
'''
#hp.write_map('NOMINALSIMULATIONS/NOM_EL60_EL80_48h_fk10_sig6_Signal.fits',s2)
#hp.write_map('NOMINALSIMULATIONS/NOM_EL60_EL80_48h_fk10_sig6-hits.fits',h)



'''
#m[h != 0] = s[h != 0]/[h != 0]
pyplot.imshow(m,interpolation='nearest',origin='lower')
pyplot.figure()

print hdu[0].data.shape
print m.shape
hdu[0].data = m
hdu.verify('fix')
hdu.writeto('XBand-Knee1Hz-BW1250-2Maps.fits',clobber=True)
hdu[0].data = h
hdu.verify('fix')
hdu.writeto('XBand-Knee1Hz-BW1250-2Maps-hits.fits',clobber=True)


for i in range(len(Tele.Skies)):
    pyplot.plot(ra[:,i],dec[:,i],'-',alpha=0.7,label='Horn %i' % i)

pyplot.xlabel('Right Ascension')
pyplot.ylabel('Declination')
pyplot.legend()

pyplot.show()

'''
