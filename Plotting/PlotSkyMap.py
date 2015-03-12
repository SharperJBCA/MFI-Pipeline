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

    grid_mapshape = grid_map.shape
    grid_map = np.reshape(grid_map,grid_map.shape[0]*grid_map.shape[1])
    grid_map[grid_map == 0] = np.nan
    grid_map = np.reshape(grid_map,grid_mapshape)

    #Plot mollweide projection:
    fig = pyplot.figure(figsize=(13,8))
    # matplotlib is doing the mollveide projection
    ax = fig.add_subplot(111,projection='mollweide')

    norm = ImageNormalize(vmin=-2.,vmax=20,stretch=vis.AsinhStretch(),clip=True)
    image = pyplot.pcolormesh(longitude[::-1], latitude, grid_map, rasterized=True, cmap=cmap,norm=norm)#, vmin=vmin, vmax=vmax)#, vmin=vmin, vmax=vmax)

    # colorbar
    cb = fig.colorbar(image, orientation='horizontal', shrink=.6, pad=0.05 ,ticks=np.array([-2,0,2.5,5,10,20]))#,ticks=[vmin, vmax])
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
    
