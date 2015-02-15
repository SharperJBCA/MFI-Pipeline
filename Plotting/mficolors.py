import numpy as np
from matplotlib.colors import ListedColormap

def load_colourmap(filename):
    '''
    Return colour table defined in filename

    Arguments
    filename -- Contains three columns of numbers defining RGB values. From 0 -> 1.
    '''

    colourmap = ListedColormap(np.loadtxt(filename))

    colourmap.set_bad("gray")
    colourmap.set_under("white")

    return colourmap


#Colormaps:
noisemap_cmap = load_colourmap('mfi-colortable-BlWhOr.txt')
whredmap_cmap = load_colourmap('mfi-colortable.txt')
