#PlanetMask.py
# This module will generate a TOD-level mask of planets and a band of coordinates 

import numpy as np
import healpy as hp
import Coordinates

def GetPlanetMask(data,maskmap,mjd,pix,lon=16.5090*np.pi/180.,lat=28.3002*np.pi/180,Planets=[-1],FWHM=[4.],bandlims={'up':-0.5,'down':-12.5},satMask=None):
    '''
    Return mask for planets and satellites

    Arguments
    maskmap -- Healpix map to contain mask
    mjd     -- Modified julian date
    pix     -- Healpix pixel coordinate for each sample

    Keyword Arguments
    lon  -- Longitude of telescope (Positive West, radians)
    lat  -- Latitude of telescope (Positive North, radians)
    
    FWHM -- List of radii to mask
    Planets -- List of planet id-codes
    bandlims -- Dictionary of shape {"up":up,"down":down} where up is the upper limit of band and lower is the lower limit.

    Note: Planets id-codes: -1 = Sun
                             1 = Mercury
                             2 = Venus
                             3 = Moon
                             4 = Mars
                             5 = Jupiter
                             6 = Saturn
                             7 = Uranus
                             8 = Neptune
    
    
    '''
    
    #Clear mask map
    maskmap[:] = True
    nside = int(np.sqrt(maskmap.size/12))

    #Loop through planets
    for i in range(len(Planets)):
        cRa,cDec = Coordinates.Ephem.plpos(np.mean(mjd),lon,lat,Planets[i])
        cPix   = hp.ang2pix(nside,(np.pi/2.-cDec),cRa)
        vec    = hp.pixelfunc.pix2vec(nside,cPix)
        ipix   = hp.query_disc(nside,vec,FWHM[i]*np.pi/180.)
        maskmap[ipix] = False

    #The Moon moves a bit fast so do this separately:
    for t in np.linspace(mjd[0],mjd[-1],6):
        cRa,cDec = Coordinates.Ephem.plpos(t,lon,lat,3)
        cPix   = hp.ang2pix(nside,(np.pi/2.-cDec),cRa)
        vec    = hp.pixelfunc.pix2vec(nside,cPix)
        ipix   = hp.query_disc(nside,vec,2.5*np.pi/180.)
        maskmap[ipix] = False
        
    #Pixels containing satellite band
    up     = bandlims['up']
    down   = bandlims['down']

    if isinstance(satMask,type(None)):
        satpix = hp.query_strip(nside, (90.-up)*np.pi/180., (90.-down)*np.pi/180.)
        maskmap[satpix] = False
  
        return maskmap[pix]
 
    else:
        #from matplotlib import pyplot
        #from MapMaker.Destriper import Control

        #nside = 512
        #npix = 12*nside**2
        #pyplot.plot(satMask[1],data,',')
        #pyplot.show()
        #m = maskmap[pix].astype('bool')
        #s = satMask[1][m]
        #d = data[m]
        #print s.size,d.size
        #Maps = Control.Destriper(d,d.size,s.astype('i'),npix,Medians=True,cn=np.ones(d.size))
        #e,a = hp.pix2ang(nside,np.arange(npix))
        #pyplot.plot(a*180./np.pi,Maps.m)
        #pyplot.show()
        #Maps.m[Maps.m == 0] = hp.UNSEEN        
        #hp.mollview(satMask[0],rot=[180,0])#*np.abs(satMask[0]-1.)
        #hp.cartview(Maps.m,rot=[180,0],lonra=[-180,180],latra=[42,52],norm='hist')        
        #pyplot.show()
        return maskmap[pix].astype('bool') & (satMask[0][satMask[1]]==0).astype('bool')
