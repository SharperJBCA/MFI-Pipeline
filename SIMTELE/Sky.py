#Generate a simulated sky TOD from a FITS image
import numpy as np
import pyfits
from astropy import wcs
from matplotlib import pyplot
import healpy as hp

def GetPixlistSqr(WCS,theta,phi,hdu):

    dcrds = np.array([phi,theta]).T
    pixcrds = WCS.wcs_world2pix(dcrds,0)
    pixcrds = pixcrds.T
    xpix = np.round(pixcrds[0,:])
    ypix = np.round(pixcrds[1,:])

    

    return np.array(xpix + hdu[0].data.shape[1]*ypix,dtype='Int64'),xpix,ypix


class Sky:

    def __init__(self,hdu,Scan,obslist,recinfo,Healpix=False):
        '''
        Hold Sky signal TOD
        '''


        if Healpix:
            if isinstance(hdu,list):
                map  = hdu[0]
                qmap = hdu[1]
                umap = hdu[2]
                Pol = True
            else:
                map = hdu
                Pol = False
                
            nside = int(np.sqrt(float(map.size)/12.))
            self.pix = hp.ang2pix(nside,(90.-Scan.dec)*np.pi/180.,Scan.ra*np.pi/180.)
        else:
            #Generate pixel coordinates
            self.WCS = wcs.WCS(naxis=2)
            self.WCS.wcs.crval = [hdu[0].header['CRVAL1'],hdu[0].header['CRVAL2']]
            self.WCS.wcs.cdelt = [hdu[0].header['CDELT1'],hdu[0].header['CDELT2']]
            self.WCS.wcs.crpix = [hdu[0].header['CRPIX1'],hdu[0].header['CRPIX2']]
            self.WCS.wcs.ctype = ['RA---TAN','DEC--TAN']
            
            self.pix,self.xpix,self.ypix = GetPixlistSqr(self.WCS,Scan.dec,Scan.ra,hdu)
            map = np.reshape(hdu[0].data,hdu[0].data.shape[0]*hdu[0].data.shape[1])

        gd = np.where((self.pix >= 0) & (self.pix < map.size))[0]

        NSamplesPerScan = int(recinfo['SampRate']/obslist['Speed'] * obslist['DSlew'])
        NSamples = int(NSamplesPerScan * obslist['NStep'])


        
        self.signal = np.zeros(NSamples)
        self.signal[gd] = map[self.pix[gd]] - np.median(map[self.pix[gd]])

        bd = np.where(np.isnan(self.signal))[0]
        self.signal[bd] = 0.

        #Generate polarised tods
        if Pol:
            self.qsignal = np.zeros(NSamples)
            self.usignal = np.zeros(NSamples)

            self.qsignal[gd] = qmap[self.pix[gd]] - np.median(qmap[self.pix[gd]])
            self.qsignal[bd] = 0.

            self.usignal[gd] = umap[self.pix[gd]] - np.median(umap[self.pix[gd]])
            self.usignal[bd] = 0.
            
