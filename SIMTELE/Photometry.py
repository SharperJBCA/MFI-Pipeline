import sys
import numpy as np
import pyfits
from matplotlib import pyplot

from scipy.ndimage.filters import gaussian_filter

import MyMaths


def Gauss2D(P,x,y):
    k = 1.3806488e-23
    c = 299792458.
    Jy = 1e26
    nu = 13.7e9

    #P[0] is in Jy.
    A = 2.*np.pi*(P[1]*P[2])*(np.pi/180.)**2 * 2. * k * nu**2 / c**2 * Jy  #2.*k*(10e9)**2/c**2*(P[1]*P[2])*(np.pi/180.)**2 * Jy #10GHz

    X = (x-P[4])*np.cos(P[3]*np.pi/180.) + (y-P[5])*np.sin(P[3]*np.pi/180.)
    Y =-(x-P[4])*np.sin(P[3]*np.pi/180.) + (y-P[5])*np.cos(P[3]*np.pi/180.)

    return (P[0]/A)*np.exp(-0.5*((X/P[1])**2 + (Y/P[2])**2))


def SimSource(hdr,map,P,sig):
    k = 1.3806488e-23
    c = 299792458.
    Jy = 1e26
    nu = 13.7e9
    
    crval1 = hdr['CRVAL1']
    crval2 = hdr['CRVAL2']
    cdelt1 = hdr['CDELT1']
    cdelt2 = hdr['CDELT2']  

    crmax1 = crval1 - map.shape[1]/2. * cdelt1 #Zero x
    crmin2 = crval2 - map.shape[0]/2. * cdelt2 #Zero y

    xvals = np.abs(cdelt1)*np.arange(-map.shape[1]/2,map.shape[1]/2) + crval1
    yvals = np.abs(cdelt2)*np.arange(-map.shape[0]/2,map.shape[0]/2) + crval2
    
    xpix = np.tile(xvals,map.shape[0])
    ypix = np.repeat(yvals,map.shape[1])

    tod = Gauss2D(P,xpix,ypix)

    sig_K = sig  /MyMaths.toJy(13.7,(np.pi/180. * 1.2/60.)**2)
    m = np.reshape(tod*1. + np.random.normal(scale=sig_K,size = tod.size),map.shape)

    print 'SUM OF SIMULATED SOURCE: ', np.sum(m)*np.abs(cdelt1*cdelt2)*(np.pi/180.)**2 * 2. * k * nu**2 / c**2 * Jy   
    #print 'TEST2: ', np.sum(tod2) *np.abs(cdelt1*cdelt2)*(np.pi/180.)**2 * 2. * k * nu**2 / c**2 * Jy 
    #print m[m.shape[1]/2,m.shape[0]/2]*1000.

    return m#np.reshape(m,m.shape[1]*m.shape[0])

def Query_Disc(xw,x,y,r,yw=None):
    '''
    Return pixels within disc of radius r.

    Arguments:
    xw -- Width of image in pixels.
    x -- X-Coordinate of disc centre in pixels.
    y -- Y-Coordinate of disc centre in pixels.
    r -- Radius of disc in pixels.

    Keyword Arguments:
    yw -- Height of image in pixels. If not set, equals xw.
    '''

    if yw == None:
        yw = xw
    
    xpix = np.tile(np.arange(xw),yw)
    ypix = np.reshape(np.tile(np.array([np.arange(yw)]).T,xw),yw*xw)

    return np.where(((xpix-x)**2 + (ypix-y)**2) < r**2)[0]

def Query_Ellipse(xw,x,y,a,b,PA=0.,yw=None):
    '''
    Return pixels within ellipse of axes a and b.

    Arguments:
    xw -- Width of image in pixels.
    x -- X-Coordinate of disc centre in pixels.
    y -- Y-Coordinate of disc centre in pixels.
    a -- Semi-major axis in pixels.
    b -- Semi-minor axis in pixels.
    
    Keyword Arguments:
    PA -- Position angle in radians. Measured clockwise.
    yw -- Height of image in pixels. If not set, equals xw.
    '''

    if yw == None:
        yw = xw
    
    xpix = np.tile(np.arange(xw),yw)
    ypix = np.reshape(np.tile(np.array([np.arange(yw)]).T,xw),yw*xw)

    X =  (xpix-x)*np.cos(PA) + (ypix-y)*np.sin(PA)
    Y = -(xpix-x)*np.sin(PA) + (ypix-y)*np.cos(PA)

    return np.where(( ((X)/a)**2 + ((Y)/b)**2) <= 1**2)[0]


def Query_Aperture(pixwidth,x,y,a1,b1,a2,b2,PA=0.):
    xpix = np.tile(np.arange(pixwidth),pixwidth)
    ypix = np.reshape(np.tile(np.array([np.arange(pixwidth)]).T,pixwidth),pixwidth**2)

    din  = Query_Ellipse(pixwidth,x,y,a1,b1,PA=PA)
    dout = Query_Ellipse(pixwidth,x,y,a2,b2,PA=PA)

    X =  (xpix-x)*np.cos(PA) + (ypix-y)*np.sin(PA)
    Y = -(xpix-x)*np.sin(PA) + (ypix-y)*np.cos(PA)

    appix = np.where(((X[dout]/a1)**2 + (Y[dout]/b1)**2) > 1**2)[0]

    return dout[appix]

def toJy(nu,beam):
    '''
    toJy(nu,beam)

    nu: frequency in GHz
    beam: beam in steradians.
    
    '''

    k = 1.3806488e-23
    c = 299792458.
    nu *= 1e9
    Jy = 1e26

    return 2.*k*nu**2/c**2 * beam * Jy


if __name__ == "__main__":

    sourcefile = sys.argv[1]
    srcs = np.loadtxt(sourcefile,skiprows=1,dtype='string',ndmin=2)


    #Have X sources I want to analyse.
    for src in srcs:
        #Read in the fits data:
        print src
        hdu = pyfits.open(src[0])
        crval1 = hdu[0].header['CRVAL1']
        crval2 = hdu[0].header['CRVAL2']
        cdelt1 = hdu[0].header['CDELT1']
        cdelt2 = hdu[0].header['CDELT2']
        map = hdu[0].data

        

        #Simulate data to test things work:
        #sig = 0.05/1000.#Jy/beam
        #P = [5./1000.,(1.2/60.)/2.355,(1.2/60.)/2.355,0.,crval1,crval2]
        #map = SimSource(hdu[0].header,map,P,sig)

        #Calculate central pixel coordinate of aperture and annulus:
        crmax1 = crval1 - map.shape[1]/2. * cdelt1
        crmin2 = crval2 - map.shape[0]/2. * cdelt2

        print cdelt1, cdelt2

        x = np.round(map.shape[1]  + (np.abs(crval1 - crmax1)/cdelt1 -  (crval1 - float(src[1]))/cdelt1)/np.cos(crval2*np.pi/180.) )
        y = np.round(np.abs(crval2 - crmin2)/cdelt2 -  (crval2 - float(src[2]))/cdelt2)

        a  = float(src[3])/60./cdelt2
        b  = float(src[4])/60./cdelt2
        PA = float(src[5])*np.pi/180.


        #Calculate the pixels within aperture and annulus:
        pix  = Query_Ellipse(map.shape[0],x,y,a,b,PA=PA)
        ape  = Query_Aperture(map.shape[0],x,y,a*1.4,b*1.4,a*2.2,b*2.2,PA=PA)
        beam = Query_Ellipse(map.shape[0],x,y,float(src[6])/60./cdelt2/2.,float(src[6])/60./cdelt2/2.,PA=PA)
        beam2= Query_Ellipse(map.shape[0],x,y,float(src[6])/60./cdelt2/2.355,float(src[6])/60./cdelt2/2.355,PA=PA)

        Nbeams = a*b/(float(src[6])/60./cdelt2/2.)**2
        
        #Extract the pixels from the map:
        mapflat = np.reshape(map,map.shape[0]*map.shape[1])
        disc = np.sum(mapflat[pix])
        aper = 0.*np.median(mapflat[ape])*float(len(pix))

        #Calculate flux + error (error uses method from haperflux, therefore expects uncorrelated white noise e.g. underestimates).
        flux = (disc - aper)*MyMaths.toJy(float(src[7]), (np.pi/180. )**2 * np.abs(cdelt1*cdelt2) )
        err1  = (np.std(mapflat[ape]))*MyMaths.toJy(float(src[7]),(np.pi/180. )**2 * np.abs(cdelt1*cdelt2)) * float(pix.size)/np.sqrt(Nbeams)/2.

        print 'FLUX: ', flux
        print 'ERROR:', err1
        print 'SNR: :', flux/err1#,np.array([a,b,a*1.3,b*1.3,a*1.8,b*1.8])*cdelt2*60. 

        #Check that the correct pixels were selected:
        mapflat[pix] = 0.
        mapflat[ape] = 0.
        
        map = np.reshape(mapflat,map.shape)

        pyplot.imshow(map,interpolation='nearest')
        pyplot.show()
        


       #err2  = np.std(mapflat[ape])*MyMaths.toJy(float(src[7]),(np.pi/180. )**2 * np.abs(cdelt1*cdelt2)) * np.sqrt(float(pix.size) +  float(pix.size)**2/float(ape.size) * np.pi/2. )
