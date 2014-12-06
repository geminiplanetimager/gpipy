import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os


__doc__ = """ Miscellaneous Plotting Routines
by Marshall Perrin.

Utilities for formatting plots, custom color maps, and plot annotation tools. Not particularly coherently organized.


"""


# purely cosmetic: print e.g. '100' instead of '10^2' for axis labels with small exponents
class NicerLogFormatter(matplotlib.ticker.LogFormatter):
    """ A LogFormatter subclass to print nicer axis labels
        e.g. '100' instead of '10^2' 

        Parameters
        ----------
        threshhold : int
            Absolute value of log base 10 above which values are printed in exponential notation. 
            By default this is 3, which means you'll get labels like 
            10^-4, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10^4 ...

        usage:
          ax = gca()
          ax.set_yscale("log")
          ax.set_xscale("log")
          ax.xaxis.set_major_formatter(NiceLogFormatter())
    """
    def __init__(self, threshhold=3):
        self.threshhold = threshhold
    def __call__(self,val,pos=None):
        try:
            if abs(np.log10(val)) > self.threshhold:
                return "$10^{%d}$" % np.log10(val)
            elif val >= 1:
                return "%d"%val
            else: 
                return "%.1f"%val
        except:
            return str(val)



class SqrtNorm(matplotlib.colors.Normalize):
    def __call__(self, value, clip=None):
        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)

        self.real_vmin = self.vmin
        self.real_vmax = self.vmax
        self.vmin = np.sqrt(self.real_vmin)
        self.vmax = np.sqrt(self.real_vmax)
        result =  matplotlib.colors.Normalize.__call__(self,np.sqrt(value))
        self.vmin = self.real_vmin
        self.vmax = self.real_vmax
        return result


def sed_plot_base(wavelength, flux, linestyle, fluxerror=None, overplot=False,units="cgs", limit_wave=None, limit_sigma=None, limit_plotsigma=3, **kwargs):
    """
        wavelength:     in microns
        flux:           in Jansky

        overplot:       True to just overplot data
        units:          'cgs' for erg/s/cm^2, 'SI' for W/m^2

    """

    if units == 'cgs':
        units_scale=1
        units_text= 'erg s$^{-1}$ cm$^{-2}$'
    elif units =='SI':
        units_scale=0.001
        units_text= 'W m$^{-2}$'
 

    else:
        print "Unknown units!!"
        stop()

    #plt.loglog(wavelength,flux, linestyle, **kwargs)
    plt.errorbar(wavelength,flux, yerr=fluxerror, fmt=linestyle, **kwargs)

    if limit_wave is not None: 
        for lw, ls in zip(limit_wave, limit_sigma):
                dy = ls*limit_plotsigma*0.5
                plt.arrow( lw, ls*limit_plotsigma, 0, -dy, )

    ax = plt.gca()
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(NiceLogFormatter())
    #ax.yaxis.set_major_formatter(NiceLogFormatter())

    plt.xlabel('Wavelength [$\mu$m]')
    plt.ylabel('Flux $\lambda F_\lambda$ [%s]' % units_text)
 



def barOutline(xdata, ydata, *args, **kwargs):
    """ Draw an outlined-box type plot (like IDL's linestyle=10)
    Based on http://www.scipy.org/Cookbook/Matplotlib/UnfilledHistograms

    Parameters
    ----------
    xdata : array_like
        left edge of bin
    ydata : array_like
        value for bin
    """

    #(histIn, binsIn) = np.histogram(dataIn, *args, **kwargs)

    stepSize = xdata[1] - xdata[0]

    bins = np.zeros(len(xdata)*2 + 2, dtype=np.float)
    data = np.zeros(len(xdata)*2 + 2, dtype=np.float)
    for bb in range(len(xdata)):
        bins[2*bb + 1] = xdata[bb]
        bins[2*bb + 2] = xdata[min(bb+1, len(xdata)-1)]
        if bb < len(ydata):
            data[2*bb + 1] = ydata[bb]
            data[2*bb + 2] = ydata[bb]

    bins[0] = bins[1]
    bins[-1] = bins[-2]
    data[0] = data[1]
    data[-1] = data[-2]

    plt.plot(bins,data, *args, **kwargs)


    

def savepdf(filename, **kwargs):
    plt.savefig(filename, transparent=True, **kwargs)
    os.system("open "+filename)

#------- Custom Color maps ------
# A Pythoninc re-implementation of DS9's colormap 'B'.
colormap_ds9_b = matplotlib.colors.LinearSegmentedColormap('DS9_B', 
        {'red':   ((0.0, 0.0, 0.0),
                   (0.25, 0.0, 0.1),
                   (0.50, 1.0,1.0),
                   (0.75, 1.0, 1.0),
                   (1.0, 1.0, 1.0)),

        
         'green': ((0.0, 0.0, 0.0),
                   (0.25, 0.0, 0.0),
                   (0.50, 0.0, 0.0),
                   (0.75, 1.0, 1.0),
                   (1.00, 1.0, 1.0)),

         'blue':  ((0.00, 0.0, 0.0),
                   (0.25, 1.0, 1.0),
                   (0.50, 0.0, 0.0),
                   (0.75, 0.0, 0.0),
                   (1.00, 1.0, 1.0))
        })
colormap_ds9_b.set_bad('black')

# A colormap based on Jet, modified for ALICE DISks paper
colormap_jet_modified = matplotlib.colors.LinearSegmentedColormap('modified_jet', {'red':   ((0., 0.03, 0.0), (0.05, 0, 0),  (0.35, 0, 0), (0.66, 1, 1), (0.89,1, 1), (1, 0.5, 0.5)),
                        'green': ((0., 0, 0), (0.125,0, 0), (0.375,0.8, 0.8), (0.64,1, 1), (0.91,0,0), (1, 0, 0)),
                        'blue':  ((0., 0.03, 0.1), (0.05, 0, 0), (0.11, 0.3, 0.3), (0.34, 0.8, 0.8), (0.65,0, 0), (1, 0, 0))} )
colormap_jet_modified.set_bad( (0, 0, 0))

# A Pythoninc re-implementation of IDL's colormap 'blue_white'.
colormap_idl_bluewhite = matplotlib.colors.LinearSegmentedColormap('IDL_BlueWhite',
        {'red':   ((0.0, 0.0, 0.0),
                   (0.756, 0.0, 0.0),
                   (1.0, 1.0, 1.0)),

        
         'green': ((0.0, 0.0, 0.0),
                   (0.38, 0.0, 0.0),
                   (1.00, 1.0, 1.0)),

         'blue':  ((0.00, 0.0, 0.0),
                   (0.737, 1.0, 1.0),
                   (1.00, 1.0, 1.0))
        })
colormap_idl_bluewhite.set_bad('black')




#------- Things to annotate plots with -------


def compass(origin, angle=0.0, ax=None, length=5, textscale=1.4, fontsize=12, color='white', labeleast=False, **kwargs):

    if ax is None: ax = plt.gca()

    dy =  np.cos(np.deg2rad(angle))*length
    dx = -np.sin(np.deg2rad(angle))*length

    # North
    ax.arrow( origin[0], origin[1], dx, dy, color=color, **kwargs)

    ax.text( origin[0]+textscale*dx, origin[1]+textscale*dy, 'N', verticalalignment='center', horizontalalignment='center', color=color, fontsize=fontsize)

    dy =  np.cos(np.deg2rad(angle+90))*length
    dx = -np.sin(np.deg2rad(angle+90))*length

    # East
    ax.arrow( origin[0], origin[1], dx, dy, color=color, **kwargs)
    if labeleast:
        ax.text( origin[0]+textscale*dx*0.9, origin[1]+textscale*dy, 'E', verticalalignment='center', horizontalalignment='center', color=color, fontsize=fontsize)


def scalebar(origin, distance=None, pixelscale=0.0756, linewidth=3, **kwargs):
    dx = 1.0/pixelscale

    ax = pl.gca()
    auto = ax.get_autoscale_on()
    ax.set_autoscale_on(False)

    pl.plot([origin[0], origin[0]+dx], [origin[1], origin[1]], color='white', linewidth=linewidth, clip_on=False, **kwargs)
    pl.text(origin[0]+dx/2, origin[1]+3, '1 arcsec', horizontalalignment='center', color='white')

    # reset the autoscale parameter? 

    if distance is not None:
        pl.text(origin[0]+dx/2, origin[1]-3, '{0:d} AU'.format(distance), horizontalalignment='center', color='white', verticalalignment='top')


def tvcircle(radius=1, xcen=0, ycen=0, center=None,**kwargs):
    """
        draw a circle on an image.

            radius
            xcen
            ycen
            center=     tuple in (Y,X) order.
    """
    if center is not None:
        xcen=center[1]
        ycen=center[0]
    t = np.arange(0, np.pi * 2.0, 0.01)
    t = t.reshape((t.size, 1))
    x = radius * np.cos(t) + xcen
    y = radius * np.sin(t) + ycen
    plt.plot(x,y, **kwargs)

    

if __name__ == "__main__":
    pass
