#TextFiles.py
# Repository of text data file writing/reading

def AppendFile(filename,items,delim=' ',append=True):
    '''
    Write a line of items to a file.

    Arguments
    filename -- File to write to (append if exists, create otherwise)
    items -- List of values to write to file

    Keyword Arguments
    delim -- Delimiter separating items
    append -- If True, append to existing file, else overwrite (Default: True)
    
    '''

    #Generate the string formatter first:
    line = ''
    for i in range(len(items)):
        line += '{'+str(i)+'}'+delim
    line += '\n'

    #Open file:
    if append:
        code = 'a'
    else:
        code = 'w'
        
    f = open(filename,code)
    f.write(line.format(*items))
    f.close()
