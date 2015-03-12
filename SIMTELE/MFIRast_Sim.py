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

Observation = {'RA':83.,
               'DEC':22.,
               'OBSMODE':'RALONMAP',
               'MJD': 51545.8 - 19.25/24.,
               'DSlew':10.,
               'DStep':0.025,
               'NStep':300,
               'Speed':0.6125,
               'LAT':28. + 18./60. + 00./60.**2,
               'LNG':16. + 30./60. + 35./60.**2.}

horns = [[-1.39,2.41],[2.41,1.39],[1.39,-2.41],[-2.41,-1.39]]
horns = [[0.,0.]]#,[2.41,1.39],[1.39,-2.41],[-2.41,-1.39]]

FocalPlane = [{'xpos':horns[i][0]*np.pi/180.,
               'ypos':horns[i][1]*np.pi/180.,
               'Pf':0,
               'Px':0,
               'Py':0,
               'Pc':0,
               'Pn':0,
               'Pa':0,
               'Pb':0} for i in range(len(horns))]

#3.8,
Receiver = [{'sigma':0.0006,
             'fknee':0.00001,
             'SampRate':50.,
             'Alpha':1.25} for i in range(len(horns))]

import MyMaths
#hdu = pyfits.open('m51_xband_ds.fits')
#hdu[0].data = hdu[0].data/MyMaths.toJy(10.,(1.4/60.*np.pi/180.)**2)
nside = 1024
map = hp.ud_grade(hp.read_map('../../../Pipeline/MAPS/wmap_K_iqu_RADEC.fits')*(23./11.)**2.15,nside)
map = hp.read_map('Horn3_Simulated_TauA_NoNoise_0-85x1-5.fits')

#Setup data containers:
NObservations = 2
NSamplesPerScan = int(Receiver[0]['SampRate']/Observation['Speed'] * Observation['DSlew'])
NSamples = int(NSamplesPerScan * Observation['NStep'])

tod = np.zeros(( NSamples*NObservations,len(Receiver)) )
pix = np.zeros(( NSamples*NObservations,len(Receiver)) )
ra  = np.zeros((NSamples*NObservations,len(Receiver)))
dec = np.zeros((NSamples*NObservations,len(Receiver)))
pang = np.zeros((NSamples*NObservations,len(Receiver)))

print 'TIME:' ,NSamples /Receiver[0]['SampRate'] / 60./60.,NSamples*NObservations


Tele = Telescope(FocalPlaneInfo=FocalPlane, ReceiverInfo=Receiver,SkyMap=map,Healpix=True)


for i in range(NObservations):
        
    Tele.Observation(Observation)
    print np.mod(i,2)
    if np.mod(i,2) == 0:
        Observation['RA'] = 80.
        Observation['MJD'] = 51545.8 - 0.2/24.
    else:
        Observation['RA'] = 285.
        Observation['MJD'] = 51545.8- 18.25/24.
        
    print Observation['RA'] ,Observation['MJD']

    for j in range(len(Tele.Receivers)):

        np.random.seed(seed=1021310 + i*len(Tele.Receivers) + j ) #Set a seed so every map is the same
        
        tod[ i*NSamples:(i+1)*NSamples,j] =  Tele.Skies[j].signal + Tele.Receivers[j].noise #+ 
        tod[ i*NSamples:(i+1)*NSamples,j] -= \
             np.median(tod[ i*NSamples:(i+1)*NSamples,j])
        pix[ i*NSamples:(i+1)*NSamples,j] = Tele.Skies[j].pix
        ra[ i*NSamples:(i+1)*NSamples,j] = Tele.Scans[j].ra
        dec[i*NSamples:(i+1)*NSamples,j] = Tele.Scans[j].dec
        pang[i*NSamples:(i+1)*NSamples,j] = Tele.Scans[j].p
        
    Observation['MJD'] +=  (NSamples/Receiver[0]['SampRate'])/24./60.**2#Move the sky a little bit




#Remove any pixels outside range:

gd = np.where((pix[:,0] >= 0) & (pix[:,0] < map.size))[0]



mask = np.ones(250,dtype='bool')
hits = np.ones(12*nside**2)

cols =['DATA','RA','DEC']
DataAccess.WriteTable('SIM_TAUA_HORN3_NONOISE.fits',[np.reshape(tod[:,0],(1,len(tod[:,0]),1)),
                                                     np.reshape(ra[:,0],(1,len(tod[:,0]),1)),
                                                     np.reshape(dec[:,0],(1,len(tod[:,0]),1))],cols)
print stop

for i in [0]:
    tod2 = tod[:,i]
    pix2 = pix[:,i]
    pyplot.plot(tod2)
    print tod2.size
    Maps2    = Control.Destriper(tod2,1000,pix2,map.size,maxiter=800,Verbose=True,mask=mask,Medians=True)
    #hp.gnomview(Maps2.m,rot=[54,31],reso=5)
    hits *= Maps2.hits
pyplot.show()

#hp.mollview(Maps1.m-Maps2.m)
#hp.mollview(Maps1.m,min=-10,max=20)

nhits = np.where(hits > 0.)[0]
nobs  = np.where(Maps2.hits > 0.)[0]
print 'Percentage area:', float(nhits.size)/float(nobs.size) * 100.

hp.gnomview(Maps2.m,rot=[83,22],reso=5)
hp.gnomview(Maps2.hits,rot=[83,22],reso=5)
hp.gnomview(hits,rot=[83,22],reso=5)

#hp.mollview(Maps2.m,min=-10,max=20)

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
