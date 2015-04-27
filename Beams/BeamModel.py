import numpy as np
from scipy.interpolate import griddata

def BeamModelTOD(model,x,y):
    '''

    Returns TOD samples from an MFI beam model.

    Arguments
    model - hdu containing beam powers and coordinates
    x - cross-elevation coordinates (-180,180) (Radians)
    y - elevation coordinates (-90,90) (Radians)

    model should contain the rows "POWER", "AZ", "EL", power should be in dB
    
    '''

    #from matplotlib import pyplot
    #pyplot.plot(model[1].data['AZ'][0,:,0],model[1].data['EL'][0,:,0],',')
    #pyplot.figure()
    #pyplot.plot(x,y,',')
    #pyplot.show()
    
    img = griddata((-model[1].data['AZ'][0,:,0],model[1].data['EL'][0,:,0]),model[1].data['POWER'][0,:,0],(x,y),fill_value=-200)

    maximg = griddata((-model[1].data['AZ'][0,:,0],model[1].data['EL'][0,:,0]),model[1].data['POWER'][0,:,0],([0],[0]),fill_value=-200)
    img = 10**(img/10.)/ 10**(maximg/10.)
    
    return img
