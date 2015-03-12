#Class defining the scans of the telescope
import numpy as np
from matplotlib import pyplot

import Coordinates

def Time(obslist,recinfo):
    '''
    Return mjd array for a particular observation
    '''

    NSamplesPerScan = int(recinfo['SampRate']/obslist['Speed'] * obslist['DSlew'])
    NSamples = int(NSamplesPerScan * obslist['NStep'])

    Period = NSamples/recinfo['SampRate']/60./60./24. #In fractions of days

    mjd = np.linspace(obslist['MJD'],obslist['MJD']+Period,NSamples)
    return mjd

def Nominal(obslist,recinfo):
    '''
    Return nominal AZEL tracks for a survey telescope
    '''

    NSamplesPerScan = int(recinfo['SampRate']/obslist['Speed'] * obslist['DSlew'])
    
    #Define slew scan:
    slew = np.linspace(0,obslist['DSlew'],NSamplesPerScan)
    slew = np.tile(slew,obslist['NStep'])

    #Define plane latitude:
    plane = np.zeros(obslist['NStep']*NSamplesPerScan) + obslist['DEC']

    return slew, plane


def Tracks(obslist,recinfo):
    '''
    Return slew and step tracks for a particular observation
    '''
    
    NSamplesPerScan = int(recinfo['SampRate']/obslist['Speed'] * obslist['DSlew'])
    
    #Define slew scan:
    if obslist['OBSMODE'] == 'AZLONMAP':
        slews = np.linspace(0,obslist['DSlew'],NSamplesPerScan)
        slews = np.append(slews,slews[::-1])
        slews = np.tile(slews,obslist['NStep']/2)

        if np.mod(obslist['NStep'],2) != 0:
            slews = np.append(slews,np.linspace(0,obslist['DSlew'],NSamplesPerScan))
    else:
        slews = np.linspace(0,obslist['DSlew'],NSamplesPerScan)
        slews = np.tile(slews,obslist['NStep'])

    #Define step scan:
    if obslist['DStep'] != 0:
        steps = np.arange(0,obslist['NStep']*obslist['DStep'],obslist['DStep'])
        steps = np.repeat(steps,NSamplesPerScan)
    else:
        steps = np.zeros(slews.size)

    #Return the tracks:
    if (obslist['OBSMODE'] == 'RALONMAP') | (obslist['OBSMODE'] == 'AZLONMAP'):
        slews = slews + obslist['RA']  - obslist['DSlew']/2.
        steps = steps + obslist['DEC'] - obslist['NStep']*obslist['DStep']/2.
        return slews,steps
    else:
        slews = slews + obslist['DEC']  - obslist['DSlew']/2.
        steps = steps + obslist['RA']   - obslist['NStep']*obslist['DStep']/2.
        return steps,slews


class Scan:

    def __init__(self,obslist,recinfo,focalplane):
        '''
        Contains the scan information of the telescope

        
        '''
        #Generate time vector:
        self.mjd = Time(obslist,recinfo)
        

        if (obslist['OBSMODE'] == 'RALONMAP') | (obslist['OBSMODE'] == 'DECLATMAP'):
            #Generate Focal Centre RA/DEC:
            ra,dec = Tracks(obslist,recinfo)


            #Generate Focal Centre AZ/EL:
            self.az,self.el = Coordinates.Sky2Hor(ra*np.pi/180.,dec*np.pi/180.,self.mjd,
                                                  lng=obslist['LNG']*np.pi/180.,lat=obslist['LAT']*np.pi/180.)

            print 'AVERAGE EL: ', np.mean(self.el*180./np.pi)

        if (obslist['OBSMODE'] == 'AZLONMAP'):
            #Generate Focal Centre AZ/EL:
            self.az,self.el = Tracks(obslist,recinfo)            
            self.az *= np.pi/180.
            self.el *= np.pi/180.


        if (obslist['OBSMODE'] == 'NOMINAL'):

            print 'heloo'
            obslist['DSlew'] = 360. #degrees, ring distance
            self.az,self.el = Nominal(obslist,recinfo)

            self.az *= np.pi/180.
            self.el *= np.pi/180.
            
        #Generate the Sky position of scan:
        self.ra,self.dec,self.p = Coordinates.Hor2Sky(self.az,self.el,self.mjd,TPoints=focalplane,
                                                 lng=obslist['LNG']*np.pi/180.,lat=obslist['LAT']*np.pi/180.)

        self.ra *= 180./np.pi
        self.dec *= 180./np.pi
        self.p *= 180./np.pi
