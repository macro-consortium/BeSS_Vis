import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates 
import matplotlib.ticker as mticker
from astropy.time import Time
import astropy.units as u
from astropy.time import Time
from datetime import datetime, timedelta
import pytz

import file_handler as fh

def wave2vel(wave):
    return (wave - 6562.8)/6562.8 * 3e5 #km/s

def vel2wave(vel):
    return 6562.8 * (1 + vel/3e5)

def plot_bess_binned_dynamic(
    xvals,
    yvals,
    fluxes,
    xrange,
    cbarlim,
    **kwargs
):
    #ax.clear()  # Clear previous plot
    cmap=kwargs.setdefault('cmap','plasma')
    phasecheck=kwargs.setdefault('phaseplot',False)
    yformat = kwargs.setdefault('yformat','mjd')
    toplabel = kwargs.setdefault('toplabel', False)
    xvals2 = kwargs.setdefault('xvals2', None)
    yvals2 = kwargs.setdefault('yvals2', None)
    fluxes2 =kwargs.setdefault('fluxes2', None)

    """
    Plot dynamic spectrum (time vs wavelength) for BeSS data.
    """

    # --- Plot ---
    #First we will set up the plotting based on how many stars will
    #be shown
    if xvals2 is not None:
        Nfigs = 2
    else:
        Nfigs = 1
    fig, ax = plt.subplots(1,Nfigs,figsize=(8, 8),sharey=True,sharex=True)

    #There's probably a smarter/better way to do this, but it works as is
    #Basically, the algorithm should be written for plotting one or both
    #Maybe it would be better to write 2 functions. One for a single star,
    #one for two...
    if xvals2 is not None:
        ax2 = ax[1]
        ax = ax[0]

    if phasecheck:
        extentarr = [xvals[0], xvals[-1], 0, 2]
        ylabel = "Orbital Phase"
    else:
        extentarr = [xvals[0], xvals[-1], yvals[0], yvals[-1]]
        ylabel = "Observation Date (MJD)"
    
    im1=ax.imshow(
        fluxes,
        aspect='auto',
        cmap=cmap,
        extent=extentarr,
        origin='lower',
        vmin = cbarlim[0],
        vmax = cbarlim[1]
    )

    # If yformat is date: calculate where the major/minor ticks need to go, and convert them from MJD to YYYY-MM-DD when displayed 
    if (yformat == 'isot') or (yformat == 'date'):
        ax.yaxis.set_major_locator(mdates.AutoDateLocator())
        ax.yaxis.set_minor_locator(mdates.AutoDateLocator(maxticks=20)) 
        def mjd_to_iso(mjd_value, pos): 
            return Time(mjd_value, format='mjd', scale='utc').to_datetime().date()
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(mjd_to_iso)) 
        labelsize=12
        ylabel = "Observation Date"
    else:
        labelsize=12

    ax.set_xlabel("Wavelength (Å)", fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)

    ax.minorticks_on()
    
    if toplabel:
        ax.tick_params(which="both",right=True,labelsize=labelsize)
        secax = ax.secondary_xaxis('top', functions=(wave2vel, vel2wave))
        secax.minorticks_on()
        secax.set_xlabel(r'Velocity (km s$^{-1}$)', fontsize=16)
        secax.tick_params(labelsize=labelsize)
    else:
        ax.tick_params(which="both",top=True,right=True,labelsize=labelsize)

    ax.set_xlim(xrange[0],xrange[1])

    #Now we do it all again for the second plot
    if xvals2 is not None:
        if phasecheck:
            extentarr = [xvals2[0], xvals2[-1], 0, 2]
        else:
            extentarr = [xvals2[0], xvals2[-1], yvals2[0], yvals2[-1]]

        im1=ax2.imshow(
            fluxes2,
            aspect='auto',
            cmap=cmap,
            extent=extentarr,
            origin='lower',
            vmin = cbarlim[0],
            vmax = cbarlim[1]
        )
        ax2.minorticks_on()
        ax2.tick_params(labelleft=False)

        ax2.set_xlabel("Wavelength (Å)", fontsize=16)

        if toplabel:
            ax2.tick_params(which="both",right=True,labelsize=labelsize)
            secax = ax2.secondary_xaxis('top', functions=(wave2vel, vel2wave))
            secax.minorticks_on()
            secax.set_xlabel(r'Velocity (km s$^{-1}$)', fontsize=16)
            secax.tick_params(labelsize=labelsize)
        else:
            ax2.tick_params(which="both",top=True,right=True,labelsize=labelsize)
        

    cbar = fig.colorbar(im1, cmap=cmap)
    cbar.set_label("Normalized Flux", size=16)
    cbar.ax.tick_params(labelsize=12)

    
    if xvals2 is not None:
        return fig, [ax, ax2]
    else:
        return fig, ax