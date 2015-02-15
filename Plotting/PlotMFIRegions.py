import numpy as np
from matplotlib import pyplot
from matplotlib import cm
import mficolors
import sys
import healpy as hp
from matplotlib.patches import Rectangle

def PlotMollSquare(ax,x1,x2,y1,y2,color='#F29900',edgecolor='k' ,fill=None,hatch='/',alpha=0.5,linewidth=1):

    
    ax.plot([x1*np.pi/180.,x2*np.pi/180.,x2*np.pi/180.,x1*np.pi/180.,x1*np.pi/180.] ,[np.pi/180.*y1,np.pi/180.*y1,np.pi/180.*y2,np.pi/180.*y2,np.pi/180.*y1],color=edgecolor,linewidth=linewidth)

    w = np.abs(x2 - x1)
    h = np.abs(y2 - y1)

    Rect1=Rectangle((x1*np.pi/180.,y1*np.pi/180.),w*np.pi/180.,h*np.pi/180.,color=color,fill=fill,hatch=hatch,alpha=alpha)
    ax.add_patch(Rect1)


if __name__ == "__main__":

    if len(sys.argv) == 1:
        print 'USE: python ...'
        sys.exit()
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




    #Plot mollweide projection:
    fig = pyplot.figure(figsize=(13,8))
    # matplotlib is doing the mollveide projection
    ax = fig.add_subplot(111,projection='mollweide')


    #NOMINAL SURVEY    
    PlotMollSquare(ax,-180,180,-32,88,color='#F29900',fill=True,hatch='',linewidth=0)

    #COSMO FIELDS:
    PlotMollSquare(ax, -61.79,-10.79 , -16.57,26.57,color='#F22000',edgecolor='k',fill=True,hatch='')
    PlotMollSquare(ax, -61.79,-10.79 , -16.57,26.57,color='#4D4D4D',edgecolor='k')
    ax.text(-25.5*np.pi/180.,5.*np.pi/180.,'Cosmo 1',ha='right',weight='bold')

    PlotMollSquare(ax,-164.98,-104.98,-19.24,20.76,color='#F22000',edgecolor='k',fill=True,hatch='')    
    PlotMollSquare(ax,-164.98,-104.98,-19.24,20.76,color='#4D4D4D',edgecolor='k')
    ax.text(-134.98*np.pi/180.,0.76*np.pi/180.,'Cosmo 2',ha='center',weight='bold')

    PlotMollSquare(ax,-(277.82-360.),-(217.82-360.),42.92,72.92,color='#F22000',edgecolor='k',fill=True,hatch='')    
    PlotMollSquare(ax,-(277.82-360.),-(217.82-360.),42.92,72.92,color='#4D4D4D',edgecolor='k')
    ax.text(-(247.82-360.)*np.pi/180.,57.92*np.pi/180.,'Cosmo 3',ha='center',weight='bold')


    #Perseus Region
    PlotMollSquare(ax,-63,-43,23,43,color='#0059F2',edgecolor='k',fill=True,hatch='')    
    PlotMollSquare(ax,-63,-43,23,43,color='#4D4D4D',edgecolor='k',hatch='\ ')


    # graticule
    ax.set_longitude_grid(60)

    ax.set_latitude_grid(45)
    ax.set_longitude_grid_ends(90)
    ax.xaxis.set_ticklabels([])
    pyplot.grid(True,linewidth=2,color='#ADADAD',linestyle=':',alpha=0.4)

    
    
    pyplot.show()
    
