import numpy as np
import glob
import sys

def UniqueFiles(substring,dir=''):
    '''
    Return unique list of files containing substring

    Arguments
    substring -- String defining file selection criteria

    Keyword Arguments
    dir -- Directory to look for files in
    '''
    files = glob.glob(dir+substring)
    files = [f[0:-8] for f in files]

    print dir, substring
    #Find unique filenames and sort:
    ufiles = np.sort(files)    
    ufiles = np.unique(ufiles)

    ufiles = np.array([u[len(dir):] for u in ufiles])

    return ufiles
    

if __name__ == '__main__':

    #User defines directory to look in and search string:
    dir = sys.argv[1]
    substring= sys.argv[2]
    ufiles = UniqueFiles(substring,dir=dir)

    print 'Saving filenames to: FileList.txt'
    np.savetxt('FileList.txt',np.reshape(ufiles,(1,ufiles.size)),delimiter='\n',fmt="%s")
