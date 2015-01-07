import numpy as np
import pyfits

def WriteTable(filename,data,cols,clobber=True,Header=None):
    '''
    Write FITS table data.

    Arguments:
    filename  -- Filename to write FITS table data to.
    data_list -- Contains data to write to table.
                  E.g. [dataset1,dataset2,...]
    cols -- Column names for the input data.
    
    Keyword Arguments:
    clobber -- Overwrite fits files (Default: True)
    Header -- Dictionary containing header names and values;
               E.g. Header={'CRVAL1': 88.52, 'CRVAL2': 1.44, ...}


    NOTES: Format of input data:
    - data[0].shape == (1,dataLength,nColumns)

    Even if you have only 1 data column, the format of the array must
    be as shown above, e.g. (1,dataLength,1) .

    '''

    #Setup a couple of containers for format codes
    f = np.empty(len(data),dtype='object')
    dim = np.empty(len(data),dtype='object')

    #Calculate format and dim of each data column
    for i in range(len(cols)):

        if type(data[i][0]) is str:
            #Calculate format and dim for string data *BROKEN*
            length = len(max(data[i][:],key = len))
            formatcode = 'A'
            f[i] = str(length)+formatcode
        else:
            #Calculate format and dim for numerical data
            length = data[i].shape[1]*data[i].shape[2]
            code = 'D'
            dim[i]  = '( '+str(data[i].shape[2])+', '+str(data[i].shape[1])+')'
            formatcode = str(length) + code
            f[i] = formatcode


    #Generate the table HDU
    tablelist = [pyfits.Column(name=c,format=f[i],array=data[i],dim=dim[i]) for i,c in enumerate(cols)]
    tbhdu = pyfits.new_table(tablelist) #This is deprecated, update function

    hdu = pyfits.PrimaryHDU()
        
    #Read in any header info 
    if Header:
        for key,val in Header.iteritems():
            hdu.header.set(key,val)
    
    thdulist = pyfits.HDUList([hdu,tbhdu])
    thdulist.verify('fix') #Check the header info order is correct

    print 'WRITING FILE AS: ', filename
    thdulist.writeto(filename,clobber=clobber)
    thdulist.close()


def WriteImage(filename,map,clobber=True,Header=None,Return=False):
    '''
    Write FITS formated image.

    Arguments:
    filename -- Filename to write FITS image to.
    map -- 2D array containing the fits image.

    Keyword Arguments:
    clobber -- Overwrite fits files (Default: True)
    Header -- Dictionary containing header names and values;
               E.g. Header={'CRVAL1': 88.52, 'CRVAL2': 1.44, ...}
    '''

    hdu = pyfits.PrimaryHDU()

    #Read in header info
    if Header:
        keys = Header.keys()
        for k in keys:
            hdu.header.set(k, Header[k])

    #Add image info to hdu
    hdu.data = map

    
    thdulist = pyfits.HDUList([hdu])
    thdulist.verify('fix') #Check the header info order is correct


    if Return:
        return thdulist
    else:
        print 'WRITING FILE AS: ', filename
        thdulist.writeto(filename,clobber=clobber)
