import numpy as np
from matplotlib import pyplot
from matplotlib import cm
import mficolors
import sys
import healpy as hp

if __name__ == "__main__":

    if len(sys.argv) == 1:
        print 'USE: python ...'
    elif len(sys.argv) == 2:
        filename = sys.argv[1]
        cmapname = 'jet'
    else:
        filename = sys.argv[1]
        cmapname = sys.argv[2]

    if cmapname == 'jet':
        cmap = cm.get_cmap('jet')
        cmap.set_bad("gray")
        cmap.set_under("white")
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

        
        

    #Setup mesh grid for image
    xsize = 2000
    ysize = xsize/2.

    theta = np.linspace(np.pi, 0, ysize)
    phi   = np.linspace(-np.pi, np.pi, xsize)
    longitude = np.radians(np.linspace(-180, 180, xsize))
    latitude = np.radians(np.linspace(-90, 90, ysize))

    PHI, THETA = np.meshgrid(phi, theta)

    #Import healpix map:
    m = hp.read_map(filename)
    nside = int(np.sqrt(m.size/12))

    vmin = np.min(m[m != hp.UNSEEN])
    vmax = np.max(m[m != hp.UNSEEN])

    #Get healpix pixels for cartisian grid:
    grid_pix = hp.ang2pix(nside,THETA,PHI)
    grid_map = m[grid_pix]


    #Plot mollweide projection:
    fig = pyplot.figure(figsize=(13,8))
    # matplotlib is doing the mollveide projection
    ax = fig.add_subplot(111,projection='mollweide')


    image = pyplot.pcolormesh(longitude[::-1], latitude, grid_map, rasterized=True, cmap=cmap, vmin=vmin, vmax=vmax)#, vmin=vmin, vmax=vmax)

    # colorbar
    cb = fig.colorbar(image, orientation='horizontal', shrink=.6, pad=0.05 ,ticks=[vmin, vmax])
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
    
