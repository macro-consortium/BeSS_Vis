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





def plot_bess_unbinned(
    xvals, 
    yvals, 
    fluxes, 
    **kwargs 
): 
    
    cmap=kwargs.setdefault('cmap','plasma')
    xformat = kwargs.setdefault('xformat', 'wavelength')
    vmin = kwargs.setdefault("vmin", 0.0) 
    vmax = kwargs.setdefault("vmax", 1.0) 

    if xformat == 'velocity':
        velocities = (xvals - 6562.8)/6562.8 * 3e5 #km/s
        xvals = velocities
        xlabel = r'Velocity (km s$^{-1}$)'
    else:
        xlabel = "Wavelength (Å)"

    # Create the plot using pcolormesh 
    fig, ax = plt.subplots(figsize=(8, 8))
    fig.subplots_adjust(left=0.2)
    im = ax.pcolormesh(
        xvals, 
        yvals, 
        fluxes, 
        cmap = cmap,  
        vmin = vmin, 
        vmax = vmax 
    )

    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel("Observation Date", fontsize=16)
    # ax.title(f"Dynamic Spectrum for {target_name} ({bin_days} day bins)")
    ax.minorticks_on() 
    if type(yvals[0]) == datetime: 
        ax.yaxis.set_major_locator(mdates.AutoDateLocator())
        ax.yaxis.set_minor_locator(mdates.AutoDateLocator(maxticks=20))
    labelsize = 12 
    ax.tick_params(which="both",top=True,right=True,labelsize=labelsize)
    cbar = fig.colorbar(im, cmap=cmap)
    cbar.set_label("Normalized Flux", size=16)
    cbar.ax.tick_params(labelsize=12)

    return fig, ax 





def plot_bess_binned_dynamic(
    xvals,
    yvals,
    fluxes,
    bin_days,
    target_name,
    **kwargs
):
    #ax.clear()  # Clear previous plot
    cmap=kwargs.setdefault('cmap','plasma')
    interp=kwargs.setdefault('interp','auto')
    yformat = kwargs.setdefault('yformat','mjd')
    xformat = kwargs.setdefault('xformat', 'wavelength')
    vmin = kwargs.setdefault("vmin", 0.0) 
    vmax = kwargs.setdefault("vmax", 1.0) 

    """
    Plot dynamic spectrum (time vs wavelength) for BeSS data.
    """

    # --- Plot ---
    #plt.ion()
    fig, ax = plt.subplots(figsize=(8, 8))
    fig.subplots_adjust(left=0.2)

    if xformat == 'velocity':
        velocities = (xvals - 6562.8)/6562.8 * 3e5 #km/s
        xvals = velocities
        xlabel = r'Velocity (km s$^{-1}$)'
    else:
        xlabel = "Wavelength (Å)"
        
    im=ax.imshow(
        fluxes,
        aspect='auto',
        cmap=cmap,
        extent=[xvals[0], xvals[-1], yvals[0], yvals[-1]],
        origin='lower', 
        vmin = vmin, 
        vmax = vmax, 
    )

    # If yformat is date: calculate where the major/minor ticks need to go, and convert them from MJD to YYYY-MM-DD when displayed 
    if (yformat == 'isot') or (yformat == 'date'):
        ax.yaxis.set_major_locator(mdates.AutoDateLocator())
        ax.yaxis.set_minor_locator(mdates.AutoDateLocator(maxticks=20)) 
        def mjd_to_iso(mjd_value, pos): 
            return Time(mjd_value, format='mjd', scale='utc').to_datetime().date()
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(mjd_to_iso)) 
        labelsize=12
    else:
        labelsize=12

    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel("Observation Date", fontsize=16)
    #ax.title(f"Dynamic Spectrum for {target_name} ({bin_days} day bins)")
    ax.minorticks_on()
    ax.tick_params(which="both",top=True,right=True,labelsize=labelsize)
    cbar = fig.colorbar(im, cmap=cmap)
    cbar.set_label("Normalized Flux", size=16)
    cbar.ax.tick_params(labelsize=12)

    #plt.tight_layout()
    return fig, ax