#ReadObs.py
# Goal:
# 1) A simple suite of functions for reading Tabular fits files.

import numpy as np
import pyfits
import glob

def FileList(substring,dir=''):
    '''
    Return list of filenames that contain the input substring

    Arguments
    substring -- string part of files to search for
    '''

    filelist = glob.glob(dir+'*'+substring+'*')
    filelist = np.sort(filelist)
    
    return filelist

def Count(filelist,keys,dir=''):
    '''
    
    Return shape of each data container if fits files and the number
    of files.

    Arguments
    filelist -- The list of files to be opened and counted
    keys -- Column names in fits file to access

    Note: Assumes pyfits table format to be: (1,DataLength,DataColumns)
    
    '''

    hdu = pyfits.open(dir+filelist[0])
    dShapes = [np.array([1.,0.,hdu[1].data[k].shape[2]],dtype='i') for k in keys]

    #Count number of samples in files    
    for i,f in enumerate(filelist):
        hdu = pyfits.open(dir+f)
        
        for ik,k in enumerate(keys):
            dShapes[ik][1] += hdu[1].data[k].shape[1]
        
        hdu.close()
        del hdu

    return dShapes,len(filelist)

def FullObs(filelist,keys,dir=''):
    '''
    Return dictionary of data from many Table Fits files

    Arguments
    filename -- Name of table fits file to access
    keys -- Array of keys to pull from table fits file

    '''

    dShapes,nFiles = Count(filelist,keys,dir=dir)

    #Setup data containers
    data = {keys[i]:np.zeros(dShapes[i]) for i in range(len(keys))}

    lastLen = np.zeros(len(keys),dtype='i')

    #Read in the data properly
    print 'TOTAL NUMBER OF FILES:',len(filelist)
    print 'OPEN: '
    for j,f in enumerate(filelist):
        print j
        hdu = pyfits.open(dir+f)

        for ik,k in enumerate(keys):
            thisLen   = hdu[1].data[k].shape[1]
            print lastLen[ik],thisLen,lastLen[ik]+thisLen
            print data[k].shape
            data[k][0,lastLen[ik]:lastLen[ik]+thisLen,:] = hdu[1].data[k][0]
            lastLen[ik]  += thisLen

        hdu.close()
        del hdu

    print ''
    print '---'
    return data


def SingleObs(filename,keys,dir=''):
    '''
    Return dictionary of data from one Table Fits file

    Arguments
    filename -- Name of table fits file to access
    keys -- Array of keys to pull from table fits file

    '''

    dShapes,nFiles = Count([filename],keys,dir='')

    #Setup data containers
    data = {keys[i]:np.zeros(dShapes[i]) for i in range(len(keys))}


    #Read in the data properly
    print 'OPEN: '
    hdu = pyfits.open(dir+filename)

    for k in keys:
        data[k][0,:,:] = hdu[1].data[k][0]

    hdu.close()
    del hdu

    print ''
    print '---'
    return data
