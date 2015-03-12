import numpy as np
import pyfits
import mlmapper
import nBinning as Binning
from matplotlib import pyplot
import healpy as hp

import DataAccess
import Control

from Telescope import Telescope

Observation = {'RA':88.558,
               'DEC':1.82,
               'OBSMODE':'RALONMAP',
               'MJD': 51545.05,
               'DSlew':12./60.,
               'DStep':12./60./20.,
               'NStep':20,
               'Speed':1./60.,
               'LAT':38. + 25./60. + 59./60.**2,
               'LNG':79. + 50./60. + 23./60.**2.}

horns = [[0.,0.],
         [300./60.**2*np.pi/180.,0.000001*np.pi/180.]]

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
Receiver = [{'sigma':1.,
             'fknee':0.5,
             'SampRate':50.,
             'Alpha':1.25} for i in range(len(horns))]

import MyMaths
#hdu = pyfits.open('m51_xband_ds.fits')
#hdu[0].data = hdu[0].data/MyMaths.toJy(10.,(1.4/60.*np.pi/180.)**2)
map = hp.read_map('../../../Pipeline/MAPS/wmap_K_iqu_RADEC.fits')*(23./11.)**2.15


#Setup data containers:
NObservations = 6
NSamplesPerScan = int(Receiver[0]['SampRate']/Observation['Speed'] * Observation['DSlew'])
NSamples = int(NSamplesPerScan * Observation['NStep'])

tod = np.zeros( NSamples*NObservations*len(Receiver) )
pix = np.zeros( NSamples*NObservations*len(Receiver) )
ra  = np.zeros((NSamples*NObservations,len(Receiver)))
dec = np.zeros((NSamples*NObservations,len(Receiver)))

print 'TIME:' ,NSamples /Receiver[0]['SampRate'] / 60./60.,NSamples*NObservations


Tele = Telescope(FocalPlaneInfo=FocalPlane, ReceiverInfo=Receiver,SkyMap=map,Healpix=True)


for i in range(NObservations):

    if np.mod(i,2) == 0:
        Observation['OBSMODE'] = 'RALONMAP'
    else:
        Observation['OBSMODE'] = 'DECLATMAP'
        
            
    Tele.Observation(Observation)

    for j in range(len(Tele.Receivers)):

        np.random.seed(seed=1021310 + i*len(Tele.Receivers) + j ) #Set a seed so every map is the same
        
        tod[i*NSamples*len(Tele.Receivers) + j*NSamples:i*NSamples*len(Tele.Receivers) + (j+1)*NSamples] =  Tele.Receivers[j].noise + Tele.Skies[j].signal
        tod[i*NSamples*len(Tele.Receivers) + j*NSamples:i*NSamples*len(Tele.Receivers) + (j+1)*NSamples] -= \
                                           np.median(tod[i*NSamples*len(Tele.Receivers) + j*NSamples:i*NSamples*len(Tele.Receivers) + (j+1)*NSamples])
        pix[i*NSamples*len(Tele.Receivers) + j*NSamples:i*NSamples*len(Tele.Receivers) + (j+1)*NSamples] = Tele.Skies[j].pix
        ra[ i*NSamples:(i+1)*NSamples,j] = Tele.Scans[j].ra
        dec[i*NSamples:(i+1)*NSamples,j] = Tele.Scans[j].dec

    Observation['MJD'] +=  1./24. + (NSamples/Receiver[0]['SampRate'])/24./60.**2#Move the sky a little bit


fig = pyplot.figure()
ax  = fig.add_subplot(111)
ax.plot(ra[:,0],dec[:,0],'-',alpha=0.8,label='Feed 1', color='#0059F2')
ax.plot(ra[:,1],dec[:,1],'-',alpha=0.8,label='Feed 2', color='#F22000')


ax.set_xlabel('R.A.')
ax.set_ylabel('Decl.')
ax.set_title('Simulated GBT Ku-Band Feed Scan Tracks')
pyplot.legend()
ax.set_xlim(ax.get_xlim()[::-1])

pyplot.show()

'''
#Remove any pixels outside range:
gd = np.where((pix >= 0) & (pix < map.size))[0]
tod = tod[gd]
pix = pix[gd]

#s,h,sig,jk = mlmapper.mlmapper(tod,pix,map.size,1,1000)
sw,hw,h= Control.Destriper(tod,2500,pix,map.size,bl_long=10000,Verbose=True,maxiter=800)

s = sw*0.
s[hw != 0] = sw[hw != 0]/hw[hw != 0]
    
s2 = np.zeros(map.size)
h = np.zeros(map.size)
cn = tod*0. + 1.
Binning.BinMap(tod,1,pix,cn,s2,hits = h)

bad = (h == 0)
good = (h != 0)
bot = np.median(s[good])
s -= bot
s[bad] = hp.UNSEEN
s2[bad] = hp.UNSEEN
h[bad] = hp.UNSEEN

hp.write_map('NOMINALSIMULATIONS/NOM_EL60_24h_fk10_sig6_BL50sec.fits',s)
#hp.write_map('NOMINALSIMULATIONS/NOM_EL60_24h_fk10_sig6_BL10sec-hits.fits',h)
'''


#hp.mollview(s2-s,min=-10,max=20)
#hp.mollview(s2,min=-10,max=20)
#hp.mollview(s,min=-10,max=20)
#hp.mollview(s2-s)
#pyplot.show()

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
