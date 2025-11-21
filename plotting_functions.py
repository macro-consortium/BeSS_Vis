import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import astropy.units as u
from astropy.time import Time
from datetime import datetime, timedelta
import pytz

import file_handler as fh

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

    """
    Plot dynamic spectrum (time vs wavelength) for BeSS data.
    """

    # --- Plot ---
    #plt.ion()
    fig, ax = plt.subplots(figsize=(8, 8))
        
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
        origin='lower'
    )


    # --- Properly aligned y-axis ticks ---
    if (yformat == 'isot') or (yformat == 'date'):
        #ax.yaxis.set_major_locator(plt.MaxNLocator(5))
        isott = Time(yvals, format='mjd', scale='utc')
        isott.format = 'isot'
        isott.out_subfmt = 'date'
        ax.set_yticks(yvals, isott.value)
        labelsize=12
    else:
        labelsize=12

    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel("Observation Date (MJD)", fontsize=16)
    #ax.title(f"Dynamic Spectrum for {target_name} ({bin_days} day bins)")
    ax.minorticks_on()
    ax.tick_params(which="both",top=True,right=True,labelsize=labelsize)
    cbar = fig.colorbar(im, cmap=cmap)
    cbar.set_label("Normalized Flux", size=16)
    cbar.ax.tick_params(labelsize=12)

    #plt.tight_layout()
    return fig, ax