import numpy as np
import pyfits
import mlmapper
import nBinning as Binning
from matplotlib import pyplot

import DataAccess
import Control

from Telescope import Telescope

Observation = {'RA':202.46841,
               'DEC':47.194814,
               'OBSMODE':'RALONMAP',
               'MJD': 51545.3,
               'DSlew':0.15,
               'DStep':11.7/60.**2,
               'NStep':47,
               'Speed':10./60./60.,
               'LAT':38. + 25./60. + 59./60.**2,
               'LNG':79. + 50./60. + 23./60.**2.}

beamSep = 94.88/60.**2
d2r = np.pi/180.
horns = [[0.,0.],
         [0.000001*d2r, beamSep*d2r],
         [0.000001*d2r,-beamSep*d2r],
         [-82.16/60.**2*d2r, 47.44/60.**2*d2r ],
         [ 82.16/60.**2*d2r,-47.44/60.**2*d2r ],
         [-82.16/60.**2*d2r,-47.44/60.**2*d2r ],
         [ 82.16/60.**2*d2r, 47.44/60.**2*d2r ]]

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
Receiver = [{'sigma':1.18,
             'fknee':1.,
             'SampRate':5.,
             'Alpha':1.25} for i in range(len(horns))]

import MyMaths
hdu = pyfits.open('m51_kband_ds.fits')
hdu[0].data = hdu[0].data/MyMaths.toJy(26.,(32./60./60.*np.pi/180.)**2)/2.*0. #HalfPower

print np.max(hdu[0].data),MyMaths.toJy(26.,(32./60./60.*np.pi/180.)**2)

#Setup data containers:
NObservations = 2
NSamplesPerScan = int(Receiver[0]['SampRate']/Observation['Speed'] * Observation['DSlew'])
NSamples = int(NSamplesPerScan * Observation['NStep'])

tod = np.zeros(NSamples*NObservations*len(Receiver))
pix = np.zeros(NSamples*NObservations*len(Receiver))
ra  = np.zeros((NSamples*NObservations,len(Receiver)))
dec = np.zeros((NSamples*NObservations,len(Receiver)))

print 'TIME:' ,NSamples /Receiver[0]['SampRate'] / 60./ 60.,NSamples*NObservations


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
        pix[i*NSamples*len(Tele.Receivers) + j*NSamples:i*NSamples*len(Tele.Receivers) + (j+1)*NSamples] = Tele.Skies[j].pix
        ra[ i*NSamples:(i+1)*NSamples,j] = Tele.Scans[j].ra
        dec[i*NSamples:(i+1)*NSamples,j] = Tele.Scans[j].dec

    Observation['MJD'] += 0.0225 #Move the sky a little bit


#Remove any pixels outside range:
gd = np.where((pix >= 0) & (pix < hdu[0].data.shape[0]*hdu[0].data.shape[1]))[0]

print 'hj',tod.size-gd.size,NSamplesPerScan

tod = tod[gd]
pix = pix[gd]


#s,h,sig,jk = mlmapper.mlmapper(tod,pix,hdu[0].data.shape[0]*hdu[0].data.shape[1],1,100)
sw,hw,h= Control.Destriper(tod,270,pix,hdu[0].data.shape[0]*hdu[0].data.shape[1],bl_long=270,Verbose=True,maxiter=800)

sd = sw*0.
sd[hw != 0] = sw[hw != 0]/hw[hw != 0]
s = sd*1.



s2 = np.zeros(hdu[0].data.shape[0]*hdu[0].data.shape[1])
h = np.zeros(hdu[0].data.shape[0]*hdu[0].data.shape[1])
cn = tod*0. + 1.
Binning.BinMap(tod,1,pix,cn,s2,hits = h)


bad = (s == 0)
good = (s != 0)
bot = np.median(s[good])
s -= bot
s2-= bot
s[bad] = np.nan
sd[bad]= np.nan
h[bad] = np.nan
s = np.reshape(s,hdu[0].data.shape)
h = np.reshape(h,hdu[0].data.shape)
s2= np.reshape(s2,hdu[0].data.shape)
sd= np.reshape(sd,hdu[0].data.shape)

m = s


#m[h != 0] = s[h != 0]/[h != 0]
#pyplot.imshow(m,interpolation='nearest',origin='lower')
#pyplot.figure()

print hdu[0].data.shape
print m.shape

#del hdu[0].header['NAXIS3']
#del hdu[0].header['NAXIS4']

#hdu[0].data = s
#hdu.verify('fix')
#hdu[0].header.set('UNITS','mK')
#hdu.writeto('RASTERSIMULATIONS/KBand-Knee1Hz-2Maps-7Rec-ML.fits',clobber=True)
#pyplot.imshow(sd)
#pyplot.show()
hdu[0].data = sd
hdu.verify('fix')
hdu.writeto('RASTERSIMULATIONS/KBand-Knee1Hz-2Maps-7Rec-Des-BL54sec.fits',clobber=True)

#hdu[0].data = s2
#hdu.verify('fix')
#hdu.writeto('KBand-Knee1Hz-BW1250-10Maps-HalfPower-NOISE.fits',clobber=True)

'''
c = ['#ADADAD','#F29900','#0059F2','#9900F2','#F22000','#59F200','#D2F200']
for i in range(len(Tele.Skies)):
    minra = np.min(ra[:,i])
    maxra = np.max(ra[:,i])
    minde = np.min(dec[:,i])
    maxde = np.max(dec[:,i])

    mra = (maxra + minra)/2.
    mde = (maxde + minde)/2.
    rra = (maxra - minra)
    rde = (maxde - minde)
    
    pyplot.plot([minra,minra,maxra,maxra,minra] ,[minde,maxde,maxde,minde,minde],'-',alpha=0.9,label='Horn %i' % i,lw=2,color=c[i])
    
    pyplot.fill_between([mra - rra/2. ,mra + rra/2.] ,[mde - rde/2.,mde - rde/2.],[mde + rde/2.,mde + rde/2.],alpha=0.2,color=c[i])

    
    #pyplot.plot(ra[:,i],dec[:,i],'-',alpha=0.7,label='Horn %i' % i)

pyplot.xlabel('Right Ascension')
pyplot.ylabel('Declination')
pyplot.legend()

i = 0
minra = np.min(ra[:,i])
maxra = np.max(ra[:,i])
minde = np.min(dec[:,i])
maxde = np.max(dec[:,i])

mra = (maxra + minra)/2.
mde = (maxde + minde)/2.
rra = (maxra - minra)
rde = (maxde - minde)
pyplot.plot([minra,minra,maxra,maxra,minra] ,[minde,maxde,maxde,minde,minde],'-',alpha=1.,label='Horn %i' % i,lw=2,color='#3D3D3D')
pyplot.fill_between([mra - rra/2. ,mra + rra/2.] ,[mde - rde/2.,mde - rde/2.],[mde + rde/2.,mde + rde/2.],alpha=0.5,color=c[i])

pyplot.show()

'''
