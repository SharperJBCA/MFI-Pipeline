#Telescope class that will handle all the simulated data.
import numpy as np

from Scan import Scan
from Receiver import Receiver
from Sky import Sky

class Telescope:

    def __init__(self, FocalPlaneInfo=None, ReceiverInfo=None,SkyMap=None,Healpix=False):
        '''
        Simulated Telescope

        Arguments

        ObsInfo-- List of user defined observations
        

        KeyWord Arguments
        
        FocalPlaneInfo -- Information on focal positions of horns and TPoint values
        ReceiverInfo   -- Information of receiver noise and 1/f characteristics
        SkyMap         -- Fits file containing image to be sampled
        '''

        if type(FocalPlaneInfo) == type(None):
            self.FocalPlaneInfo = [{'xpos':0.,
                                    'ypos':0.,
                                    'Pf':0,
                                    'Px':0,
                                    'Py':0,
                                    'Pc':0,
                                    'Pn':0,
                                    'Pa':0,
                                    'Pb':0}]
        else:
            self.FocalPlaneInfo = FocalPlaneInfo

        if type(ReceiverInfo) == type(None):
            self.ReceiverInfo = [{'sigma':1.,
                                  'fknee':1.,
                                  'SampRate':100.,
                                  'Alpha':1.}]
        else:
            self.ReceiverInfo = ReceiverInfo

        if type(SkyMap) == type(None):
            self.SkyMap = None
        else:
            self.SkyMap = SkyMap

        self.Healpix = Healpix

    def Observation(self,ObsList):
        '''
        Prime telescope for a particular observation

        Arguments

        ObsList -- A dict of variables defining the observation

                   Format: {RA,DEC:,MJD:,OBSMODE:,DSlew:,DStep:,NStep:,Speed:,LAT:,LNG:}
        '''

        self.ObsList = ObsList

        #Define the scans
        self.Scans = [Scan(self.ObsList,self.ReceiverInfo[i],focalplane) for i,focalplane in enumerate(self.FocalPlaneInfo)]

        #Define the receivers
        self.Receivers = [Receiver(self.ObsList,recinfo) for recinfo in self.ReceiverInfo]

        #Define the Sky signal
        if type(self.SkyMap) != type(None):
            self.Skies = [Sky(self.SkyMap,self.Scans[i],self.ObsList,self.ReceiverInfo[i],Healpix=self.Healpix) for i in range(len(self.Scans))]
