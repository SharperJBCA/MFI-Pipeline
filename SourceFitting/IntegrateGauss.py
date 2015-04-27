import numpy as np
from matplotlib import pyplot
from scipy.interpolate import interp1d

import healpy as hp
from MapMaker.Destriper import Control

gauss = lambda P,x,y: P[0]*np.exp(-0.5*( ((x-P[1])/P[3])**2 + ((y-P[2])/P[4])**2))

p0 = [1.,0.,0.,1.,1.]

fullIntegral = 2.*np.pi * p0[0] * p0[3]*p0[4]

#x,y = np.meshgrid(np.linspace(-2.*np.sqrt(2.*np.log(2.)),2.*np.sqrt(2.*np.log(2.)),300),np.linspace(-2.*np.sqrt(2.*np.log(2.)),2.*np.sqrt(2.*np.log(2.)),300))

x,y = np.meshgrid(np.linspace(-10,10,3200),np.linspace(-10,10,3200))


diff = (x[0,0] - x[0,1])**2

xshape = x.shape
x = np.reshape(x,x.shape[0]*x.shape[1])
y = np.reshape(y,y.shape[0]*y.shape[1])

gd = np.where((x**2 + y**2 < (2.*np.sqrt(2.*np.log(2.)))**2))[0]#2.*
#gd = np.where((x**2 + y**2 < 1.))[0]#2.*

x = x[gd]
y = y[gd]

print np.sum(gauss(p0,x,y))*diff
print 'Ratio of FWHM to Full Gauss:',np.sum(gauss(p0,x,y))*diff/fullIntegral
print fullIntegral/(np.sum(gauss(p0,x,y))*diff)


#Compare two different beams:

pa = [1.,0.,0.,0.90/2.355,0.90/2.355]
pb = [1.,0.,0.,0.85/2.355,0.85/2.355]

x,y = np.meshgrid(np.linspace(-10,10,3200),np.linspace(-10,10,3200))

sa = np.sum(gauss(pa,x,y))*diff
sb = np.sum(gauss(pb,x,y))*diff

print 'Ratio of beams: ', sb/sa


#FullBeam

d2  = np.loadtxt('MFIBeams/SATB-130222-0111.txt')

x = np.reshape(x,x.shape[0]*x.shape[1])
y = np.reshape(y,y.shape[0]*y.shape[1])

r = np.sqrt(x**2 + y**2)
out = np.where(r > 1.7154)[0]
out = np.where(r > 0.85*3.)[0]

r_samp = d2[:,0]
b113 = np.mean(d2[:,[7]],axis=1)

ip = interp1d(r_samp,np.mean(d2[:,[7]],axis=1),bounds_error=False,fill_value=0.)#,kind='cubic')
#print x.shape,x.shape[1]/2
#r = (x[x.shape[1]/2,:])
#print np.min(r),np.max(r)


#pyplot.plot(r,gauss(pb,r,np.zeros(r.size)))
#pyplot.plot(r,ip(r))
#pyplot.show()

modelvals = ip(r)
gaussvals = gauss(pb,x,y)
modelvals /= np.max(modelvals)
gaussvals /= np.max(gaussvals)

modelvals[out] =0.
gaussvals[out] =0.

print np.sum(modelvals)/np.sum(gaussvals)
print np.sum(gaussvals)/np.sum(modelvals)
print np.sum(modelvals)*diff/(2.*np.pi*pb[3]*pb[4])

modelvals = np.reshape(modelvals,xshape)
gaussvals = np.reshape(gaussvals,xshape)
x = np.reshape(x,xshape)
y = np.reshape(y,xshape)

print r.shape
cs = pyplot.contourf(x,y,np.log10(modelvals))
pyplot.contour(x,y,np.log10(gaussvals),levels=cs.levels,colors='k')
pyplot.xlim(-8,8)
pyplot.ylim(-8,8)
pyplot.show()

modelvals = ip(r)/np.max(modelvals) + np.random.normal(size=r.size,scale=0.000001)
nside = 256
npix = 12*nside**2
pbeam = 4.*np.pi/(12.*nside**2)

x = np.reshape(x,x.shape[0]*x.shape[1])
y = np.reshape(y,y.shape[0]*y.shape[1])
gaussvals = gauss(pb,x,y) + np.random.normal(size=r.size,scale=0.000001)

pix = hp.ang2pix(nside,(np.pi/2. - y*np.pi/180.),x*np.pi/180.)

print modelvals.shape,pix.shape
Maps = Control.Destriper(modelvals.astype('f'),modelvals.size,pix.astype('i'),npix,Medians=True)
Maps2 = Control.Destriper(gaussvals.astype('f'),modelvals.size,pix.astype('i'),npix,Medians=True)


ipix = hp.query_disc(nside,hp.ang2vec(np.pi/2.,0.),0.85*np.pi/180.)

print np.sum(Maps.m[ipix])*pbeam/(1.133*(0.85*np.pi/180.)**2) *0.943
#[ipix])
print np.sum(Maps2.m[ipix])*pbeam/(1.133*(0.85*np.pi/180.)**2)

print stop
print np.sum(Maps2.m[ipix])*pbeam/(1.133*(0.85*np.pi/180.)**2)

print np.sum(Maps2.m[ipix])/np.sum(Maps.m[ipix])#MB/Full Beam efficiency?
print np.sum(Maps2.m)/np.sum(Maps.m)#MB/Full Beam efficiency?

hp.gnomview(10.*np.log10(Maps.m/np.max(Maps.m)),reso=3,min=-80)
hp.gnomview(10.*np.log10(Maps2.m/np.max(Maps2.m)),reso=3,min=-80)
pyplot.show()
