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

Observation = {'RA':83.,
               'DEC':70.,
               'OBSMODE':'NOMINAL',
               'MJD': 51545.0+0.775,
               'DSlew':360.,
               'DStep':0.,
               'NStep':360,
               'Speed':6.,
               'LAT':28. + 18./60. + 00./60.**2,
               'LNG':16. + 30./60. + 35./60.**2.}

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
Receiver = [{'sigma':0.1,
             'fknee':10.,
             'SampRate':50.,
             'Alpha':1.25} for i in range(len(horns))]

import MyMaths
#hdu = pyfits.open('m51_xband_ds.fits')
#hdu[0].data = hdu[0].data/MyMaths.toJy(10.,(1.4/60.*np.pi/180.)**2)
imap = hp.ud_grade(hp.read_map('../../../Pipeline/MAPS/wmap_K_iqu_RADEC.fits',0)*(23./11.)**2.15,256)
qmap = hp.ud_grade(hp.read_map('../../../Pipeline/MAPS/wmap_K_iqu_RADEC.fits',1)*(23./11.)**2.15,256)
umap = hp.ud_grade(hp.read_map('../../../Pipeline/MAPS/wmap_K_iqu_RADEC.fits',2)*(23./11.)**2.15,256)

maps = [imap,qmap,umap]
mods = [0.,22.5,45.,67.5]

#Setup data containers:
NObservations = 4
NSamplesPerScan = int(Receiver[0]['SampRate']/Observation['Speed'] * Observation['DSlew'])
NSamples = int(NSamplesPerScan * Observation['NStep'])

tod  = np.zeros( NSamples*NObservations*len(Receiver) )
qtod = np.zeros( NSamples*NObservations*len(Receiver) )
utod = np.zeros( NSamples*NObservations*len(Receiver) )
noise= np.zeros( NSamples*NObservations*len(Receiver) )

pix = np.zeros( NSamples*NObservations*len(Receiver) )
ra  = np.zeros((NSamples*NObservations,len(Receiver)))
dec = np.zeros((NSamples*NObservations,len(Receiver)))
p   = np.zeros((NSamples*NObservations,len(Receiver)))

print 'TIME:' ,NSamples /Receiver[0]['SampRate'] / 60./60.,NSamples*NObservations
print 'Samps Per Scan:',NSamplesPerScan

Tele = Telescope(FocalPlaneInfo=FocalPlane, ReceiverInfo=Receiver,SkyMap=maps,Healpix=True)


for i in range(NObservations):

    #f np.mod(i,2) == 0:
    #    Observation['OBSMODE']= 'RALONMAP'
    #else:
    #    Observation['OBSMODE']= 'DECLATMAP'
    
        
    Tele.Observation(Observation)
    


    for j in range(len(Tele.Receivers)):
            
        #np.random.seed(seed=1021310 + i*len(Tele.Receivers) + j ) #Set a seed so every map is the same
        
        tod[i*NSamples*len(Tele.Receivers) + j*NSamples:i*NSamples*len(Tele.Receivers) + (j+1)*NSamples] =  Tele.Skies[j].signal 
        tod[i*NSamples*len(Tele.Receivers) + j*NSamples:i*NSamples*len(Tele.Receivers) + (j+1)*NSamples] -= \
                                           np.median(tod[i*NSamples*len(Tele.Receivers) + j*NSamples:i*NSamples*len(Tele.Receivers) + (j+1)*NSamples])
        
        noise[i*NSamples*len(Tele.Receivers) + j*NSamples:i*NSamples*len(Tele.Receivers) + (j+1)*NSamples] =  Tele.Receivers[j].noise
        
        qtod[i*NSamples*len(Tele.Receivers) + j*NSamples:i*NSamples*len(Tele.Receivers) + (j+1)*NSamples] =  Tele.Skies[j].qsignal
        utod[i*NSamples*len(Tele.Receivers) + j*NSamples:i*NSamples*len(Tele.Receivers) + (j+1)*NSamples] =  Tele.Skies[j].usignal
        
        pix[i*NSamples*len(Tele.Receivers) + j*NSamples:i*NSamples*len(Tele.Receivers) + (j+1)*NSamples] = Tele.Skies[j].pix
        ra[ i*NSamples:(i+1)*NSamples,j] = Tele.Scans[j].ra
        dec[i*NSamples:(i+1)*NSamples,j] = Tele.Scans[j].dec
        p[i*NSamples:(i+1)*NSamples,j]   = Tele.Scans[j].p

    Observation['MJD'] += (NSamples/Receiver[0]['SampRate'])/24./60.**2#Move the sky a little bit




P = qtod*np.sin(2.*p[:,0]*np.pi/180.) + utod*np.cos(2.*p[:,0]*np.pi/180.) + noise

mod = np.zeros(tod.size) + 45.
DataAccess.WriteTable('Mod2.fits',[np.reshape(qtod,(1,tod.size,1)),
                                   np.reshape(utod,(1,tod.size,1)),
                                   np.reshape(noise,(1,tod.size,1)),
                                   np.reshape(pix,(1,pix.size,1)),
                                   np.reshape(mod,(1,mod.size,1)),
                                   np.reshape(p[:,0],(1,mod.size,1))],cols=['QTOD','UTOD','NOISE','PIX','MOD','PHI'])


print stop

tod = P
cn = tod*0. +1.

#Remove any pixels outside range:
gd = np.where((pix >= 0) & (pix < maps[0].size))[0]
#tod = tod[gd]
#pix = pix[gd]

print tod.size,pix.size,p[:,0].size

#s,h,sig,jk = mlmapper.mlmapper(tod,pix,map.size,1,1000)
#m = Control.MLMapper(tod,pix,map.size,1,1000)
#s = Control.MLMapper(tod,pix,map.size,bl_long=1000,Verbose=True)
Maps1    = Control.Destriper(tod,3000,pix,maps[0].size,p[:,0]*np.pi/180.,bl_long=3000,Verbose=True,maxiter=800)
#s = Maps.m
#s = sw*0.
#s[hw != 0] = sw[hw != 0]/hw[hw != 0]

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

Binning.BinMapPol_Angs(tod,1,pix.astype('i'),p[:,0]*np.pi/180.,cn,Maps2)
Binning.BinMapPol(tod,1,pix.astype('i'),p[:,0]*np.pi/180.,cn,Maps2)

Binning.BinMap(tod,1,pix.astype('i'),cn,Maps2.m)

#hp.write_map('DES_QNOM.fits',Maps1.q)
#hp.write_map('DES_UNOM.fits',Maps1.u)
#hp.write_map('QNOM.fits',qmap)
#hp.write_map('UNOM.fits',umap)

hp.mollview(Maps1.q,rot=[83,22])
hp.gnomview(Maps1.q,rot=[83,22])
hp.gnomview(Maps1.u,rot=[83,22])
hp.mollview(Maps2.m)
hp.mollview(Maps2.q,rot=[83,22])
hp.gnomview(Maps2.u,rot=[83,22])

#hp.gnomview(Maps.sw,rot=[83,22])
hp.gnomview(qmap,rot=[83,22])
hp.gnomview(umap,rot=[83,22])

pyplot.show()

'''
bad = (h == 0)
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
