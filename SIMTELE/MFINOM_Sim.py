import numpy as np
import pyfits
#import mlmapper
from MapMakerTest.Tools import nBinning as Binning
import Binning as Bin
from matplotlib import pyplot
import healpy as hp

import DataAccess
#import Control

from MapMaker.Destriper import Control as Control
#from MapMakerTest.Destriper import Control
#from MapMaker.Destriper import Control
#from MapMaker.MLMapper import Control

from Telescope import Telescope

Observation = {'RA':90,
               'DEC':30.,
               'OBSMODE':'NOMINAL',
               'MJD': 51545.0+0.85,
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
Receiver = [{'sigma':6.,
             'fknee':1.,
             'SampRate':50.,
             'Alpha':1.} for i in range(len(horns))]

import MyMaths
#hdu = pyfits.open('m51_xband_ds.fits')
#hdu[0].data = hdu[0].data/MyMaths.toJy(10.,(1.4/60.*np.pi/180.)**2)
nside = 256
map = hp.ud_grade(hp.read_map('../../../Pipeline/MAPS/wmap_K_iqu_RADEC.fits')*(23./11.)**2.15,nside)


#Setup data containers:
NObservations = 16
NSamplesPerScan = int(Receiver[0]['SampRate']/Observation['Speed'] * Observation['DSlew'])
NSamples = int(NSamplesPerScan * Observation['NStep'])

tod = np.zeros( NSamples*NObservations*len(Receiver) )
tod2= np.zeros( NSamples*NObservations*len(Receiver) )
pix = np.zeros( NSamples*NObservations*len(Receiver) )
ra  = np.zeros((NSamples*NObservations,len(Receiver)))
dec = np.zeros((NSamples*NObservations,len(Receiver)))
pang = np.zeros((NSamples*NObservations,len(Receiver)))

print 'TIME:' ,NSamples /Receiver[0]['SampRate'] / 60./60.,NSamples*NObservations


Tele = Telescope(FocalPlaneInfo=FocalPlane, ReceiverInfo=Receiver,SkyMap=map,Healpix=True)


for i in range(NObservations):
        
    Tele.Observation(Observation)
    if np.mod(i+1,4) == 0:
        Observation['MJD'] += 13./60./60./24.*51. #Jump some random number ahead
    #    Observation['DEC'] = 80.


    for j in range(len(Tele.Receivers)):

        np.random.seed(seed=1021310 + i*len(Tele.Receivers) + j ) #Set a seed so every map is the same
        
        tod[i*NSamples*len(Tele.Receivers) + j*NSamples:i*NSamples*len(Tele.Receivers) + (j+1)*NSamples] =  Tele.Receivers[j].noise+Tele.Receivers[j].noisespikes#Tele.Skies[j].signal +  
        tod[i*NSamples*len(Tele.Receivers) + j*NSamples:i*NSamples*len(Tele.Receivers) + (j+1)*NSamples] -= \
                                           np.median(tod[i*NSamples*len(Tele.Receivers) + j*NSamples:i*NSamples*len(Tele.Receivers) + (j+1)*NSamples])


        tod2[i*NSamples*len(Tele.Receivers) + j*NSamples:i*NSamples*len(Tele.Receivers) + (j+1)*NSamples] =  Tele.Receivers[j].noise#Tele.Skies[j].signal +  
        tod2[i*NSamples*len(Tele.Receivers) + j*NSamples:i*NSamples*len(Tele.Receivers) + (j+1)*NSamples] -= \
                                            np.median(tod2[i*NSamples*len(Tele.Receivers) + j*NSamples:i*NSamples*len(Tele.Receivers) + (j+1)*NSamples])

        
        pix[i*NSamples*len(Tele.Receivers) + j*NSamples:i*NSamples*len(Tele.Receivers) + (j+1)*NSamples] = Tele.Skies[j].pix
        ra[ i*NSamples:(i+1)*NSamples,j] = Tele.Scans[j].ra
        dec[i*NSamples:(i+1)*NSamples,j] = Tele.Scans[j].dec
        pang[i*NSamples:(i+1)*NSamples,j] = Tele.Scans[j].p

    Observation['MJD'] +=  (NSamples/Receiver[0]['SampRate'])/24./60.**2#Move the sky a little bit


#fig = pyplot.figure(figsize=(8,5))
# matplotlib is doing the mollveide projection
#ax = fig.add_subplot(111,projection='mollweide')

print np.mean(dec[:,0]),np.max(dec[:,0]),np.min(dec[:,0])

#ax.plot((180.-ra[:,0])*np.pi/180.,dec[:,0]*np.pi/180.,'-',color='#0059F2',lw=2)
#ax.set_longitude_grid(60)

#ax.set_latitude_grid(45)
#ax.set_longitude_grid_ends(90)
#ax.xaxis.set_ticklabels([])
#pyplot.grid(True,linewidth=2,color='#ADADAD',linestyle=':',alpha=0.8)

#pyplot.show()

#DataAccess.WriteTable('TestData.fits',[np.reshape(tod,(1,tod.size,1)),
#                                       np.reshape(pix,(1,pix.size,1))],cols=['TOD','PIX'])

#Remove any pixels outside range:

gd = np.where((pix >= 0) & (pix < map.size))[0]
tod = tod[gd]
tod2= tod2[gd]
pix = pix[gd]

deta = np.zeros(12*nside**2)
Binning.BinMapPol_Angs(tod,250,pix,2.*pang*np.pi/180.,tod*0.+1.,deta)

#hp.mollview(deta)
#pyplot.show()
#s,h,sig,jk = mlmapper.mlmapper(tod,pix,map.size,1,1000)
#m = Control.MLMapper(tod,pix,map.size,1,1000)
#s = Control.MLMapper(tod,pix,map.size,bl_long=1000,Verbose=True)


mask = np.ones(250,dtype='bool')
#mask[0:5] = False
mask = np.tile(mask,tod.size/250)

mask[deta[pix.astype('i')] < 1e-4] = False

#spikes
gd = np.where(mask)[0]

#np.random.seed(seed=100)
#spikes = np.random.uniform(low=0,high=len(gd)-1,size=300).astype('i')
#tod[gd[spikes]] = 10000.
#mask[gd[spikes]] = False

pyplot.plot(tod,',')
pyplot.show()

#Maps1    = Control1.Destriper(tod,250,pix,map.size,bl_long=10000,maxiter=800,Verbose=True)
#Maps2    = Control.Destriper(tod,250,pix,map.size,bl_long=10000,maxiter=800,Verbose=True,mask=mask)
Maps1    = Control.Destriper(tod,tod.size,pix,map.size,maxiter=800,Verbose=True,Medians=True)
Maps2    = Control.Destriper(tod2,tod.size,pix,map.size,maxiter=800,Verbose=True,Medians=True)

ipx = hp.query_disc(nside,hp.ang2vec((90.-28.)*np.pi/180.,0.),5.*np.pi/180.)

print 'NOISE:', np.std(Maps1.m[ipx]),6. * np.sqrt(1./50.)#*np.sqrt(np.mean(Maps1.hits[ipx])/50.)
print 'NOISE:', np.std(Maps2.m[ipx]),6. * np.sqrt(1./50.)

Maps1.m *= np.sqrt(Maps1.hits)
Maps1.m[Maps1.m == 0] = hp.UNSEEN
Maps2.m *= np.sqrt(Maps2.hits)
Maps2.m[Maps2.m == 0] = hp.UNSEEN

print 'NOISE:', np.std(Maps1.m[ipx]),6. * np.sqrt(1./50.)#*np.sqrt(np.mean(Maps1.hits[ipx])/50.)
print 'NOISE:', np.std(Maps2.m[ipx]),6. * np.sqrt(1./50.)
hp.mollview(Maps1.m-Maps2.m,coord=['C','C'])
hp.mollview(Maps1.m,coord=['C','C'])#,min=-10,max=10)
hp.mollview(Maps2.m,coord=['C','C'])

pyplot.show()

#Maps    = Control.Destriper(tod,250,pix,map.size,bl_long=10000,Verbose=True,maxiter=800)
s = Maps.m
#s = sw*0.
#s[hw != 0] = sw[hw != 0]/hw[hw != 0]
    
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

#hp.write_map('NOMINALSIMULATIONS/NOM_EL60_EL80_48h_fk10_sig6_Signal.fits',s2)
#hp.write_map('NOMINALSIMULATIONS/NOM_EL60_EL80_48h_fk10_sig6-hits.fits',h)



#hp.mollview(s2-s,min=-10,max=20)
hp.mollview(s2,min=-10,max=20)
hp.mollview(s,min=-10,max=20)
hp.mollview(s2-s)
pyplot.show()

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
