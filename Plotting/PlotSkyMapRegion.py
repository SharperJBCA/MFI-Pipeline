import numpy as np
from matplotlib import pyplot
from matplotlib import cm
import mficolors
import sys
import healpy as hp
import astropy.visualization as vis
from astropy.visualization.mpl_normalize import ImageNormalize
from skimage import exposure
import matplotlib as mpl

import aplpy

if __name__ == "__main__":

    if len(sys.argv) == 1:
        print 'USE: python ...'
    elif len(sys.argv) == 2:
        filename = sys.argv[1]
        cmapname = 'RdBu_r'#'jet'
    else:
        filename = sys.argv[1]
        cmapname = sys.argv[2]

    if cmapname == 'RdBu_r':
        cmap = cm.get_cmap(cmapname)
        cmap.set_bad("#ADADAD")
        cmap.set_under("#EEEEEE")
    elif cmapname == 'noise':
        cmap = mficolors.noisemap_cmap
    else:
        cmap = mficolors.whredmap_cmap
        #b =  np.array(cmap.__dict__['_segmentdata']['blue'])[5:,:]
        #r =  np.array(cmap.__dict__['_segmentdata']['red'])[5:,:]
        #g =  np.array(cmap.__dict__['_segmentdata']['green'])[5:,:]

        #cmap.__dict__['_segmentdata']['blue'] = list(b)
        #cmap.__dict__['_segmentdata']['red'] = list(r)
        #cmap.__dict__['_segmentdata']['green'] = list(g)

        #help( cmap.cdict)
        #cmap = mficolors.noisemap_cmap

        
        

    gc = aplpy.FITSFigure(filename)
    gc.show_colorscale(cmap = cmap)#,stretch='arcsinh')#,vmax=150)
    gc.tick_labels.set_xformat('dd')
    gc.tick_labels.set_yformat('dd')
    gc.tick_labels.set_font(size=14)  
    gc.add_colorbar()
    gc.colorbar.set_font(size=14)
    gc.colorbar.set_axis_label_text('mK')
    gc.colorbar.set_axis_label_font(size=16)
    gc.set_nan_color('#DDDDDD')
    
    #gc.axis_labels.hide_x()
    #gc.axis_labels.hide_y()
    #gc.tick_labels.hide_x()
    #gc.tick_labels.hide_y()

    pyplot.show()

    '''
    #Plot mollweide projection:
    fig = pyplot.figure(figsize=(13,8))
    # matplotlib is doing the mollveide projection
    ax = fig.add_subplot(111,projection='mollweide')

    norm = ImageNormalize(vmin=-2.,vmax=15,stretch=vis.AsinhStretch(),clip=True)
    image = pyplot.pcolormesh(longitude[::-1], latitude, grid_map, rasterized=True, cmap=cmap,norm=norm)#, vmin=vmin, vmax=vmax)#, vmin=vmin, vmax=vmax)

    # colorbar
    cb = fig.colorbar(image, orientation='horizontal', shrink=.6, pad=0.05 )#,ticks=[vmin, vmax])
    cb.ax.xaxis.set_label_text('mK')
    cb.ax.xaxis.labelpad = -2
    # workaround for issue with viewers, see colorbar docstring
    cb.solids.set_edgecolor("face")

    # graticule
    ax.set_longitude_grid(60)
    #ax.xaxis.set_major_formatter(ThetaFormatterShiftPi(60))

    ax.set_latitude_grid(45)
    ax.set_longitude_grid_ends(90)
    ax.xaxis.set_ticklabels([])
    pyplot.grid(True,linewidth=2,color='#ADADAD',linestyle=':',alpha=0.8)
    
    pyplot.show()
    '''
