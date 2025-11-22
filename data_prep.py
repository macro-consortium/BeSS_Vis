import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.time import Time
from datetime import datetime, timedelta
from scipy.interpolate import interp1d
import warnings

import file_handler as fh

def read_spectra(files, date_range, wmin, wmax, common_wave):
    hold_dates = []
    hold_mjd = []
    hold_fluxes = []
    
    for fname in files:
        try:
            wave, flux, header = fh.BeSSSpectra(fname)
            date_obs = header.get('DATE-OBS', None)
            dates_obs_T = Time(date_obs, format='isot', scale='utc')
            mjd_obs = dates_obs_T.mjd
            if not date_obs:
                continue

            try:
                date_obj = datetime.strptime(date_obs[:10], "%Y-%m-%d")
            except ValueError:
                continue

            # --- filter by date range ---
            if date_range is not None:
                if mjd_obs < date_range[0] or mjd_obs > date_range[1]:
                    continue

            # --- restrict to wavelength range ---
            mask = (wave >= wmin) & (wave <= wmax)
            if np.sum(mask) < 10:
                continue

            # --- normalize + interpolate onto common grid ---
            flux = flux / np.nanmax(flux)
            interp_flux = np.interp(
                common_wave,
                wave[mask],
                flux[mask]
            )

            hold_dates.append(date_obj)
            hold_mjd.append(mjd_obs)
            hold_fluxes.append(interp_flux)

        except Exception as e:
            print(f"⚠️ Skipping {fname}: {e}")
            continue
    return hold_dates, hold_mjd, hold_fluxes

def bin_times_and_spectra(dates, flux_arrays, bin_days):
    """
    Bin spectra into fixed time bins.

    Parameters
    ----------
    dates : list of datetime objects
    flux_arrays : list of 1D numpy arrays (already interpolated to common grid)
    bin_days : int
        Bin width in days.

    Returns
    -------
    binned_dates : list of datetime (bin centers)
    binned_fluxes : list of 1D arrays (averaged per bin)
        Missing-data bins are filled with NaNs (blank rows).
    """
    if len(dates) == 0:
        binned_dates = np.arange(50000,60000,1000)
        flux_matrix = np.empty((1000, len(binned_dates),))
        flux_matrix[:] = np.nan

        return binned_dates, flux_matrix        

    # Sort by date
    #sorted_idx = np.argsort(dates)
    #dates = np.array([dates[i] for i in sorted_idx])
    #flux_arrays = [flux_arrays[i] for i in sorted_idx]

    # Build bins
    start_date = dates[0]
    end_date = dates[-1]

    bins = []
    t = start_date
    while t <= end_date + bin_days:
        bins.append(t)
        t += bin_days

    binned_dates = []
    binned_fluxes = []

    for i in range(len(bins) - 1):
        t0 = bins[i]
        t1 = bins[i + 1]
        mid = t0 + (t1 - t0) / 2

        # indices of spectra in this bin
        idx = np.where((dates >= t0) & (dates < t1))[0]
        
        if len(idx) == 0:
            # blank bin → all NaNs
            blank = np.full_like(flux_arrays[0], np.nan)
            binned_dates.append(mid)
            binned_fluxes.append(blank)
        else:
            # average available spectra
            stacked = np.vstack([flux_arrays[k] for k in idx])
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                avg = np.nanmean(stacked, axis=0)
            binned_dates.append(mid)
            binned_fluxes.append(avg)

    flux_matrix = np.vstack(binned_fluxes)
        
    return binned_dates, flux_matrix

def data_setup(files, **kwargs):

    feature=kwargs.setdefault('feature','All')
    date_range=kwargs.setdefault('dates',None)
    interp=kwargs.setdefault('interp','auto')
    xformat = kwargs.setdefault('xformat', 'wavelength')
    
    # --- define common wavelength grid ---
    if feature == 'All':
        wmin, wmax = 6510, 6624
    elif feature == 'Hα' or feature == 'Halpha' or feature == 'H alpha':
        wmin, wmax = 6548, 6578
    elif feature == 'HeI' or feature == 'He I':
        wmin, wmax = 6648, 6708
    else:
        print("Feature not currently supported")
        return
    
    common_wave = np.linspace(wmin, wmax, 1000)

    # --- parse date range if provided ---
    if date_range is not None:
        try:
            start_date = date_range[0]
            end_date = date_range[1]
        except Exception:
            raise ValueError("❌ date_range not accepted (work on error message).")
    else:
        start_date = end_date = None

    # --- read each spectrum ---
    all_dates, all_mjd, all_fluxes = read_spectra(files, date_range, wmin, wmax, common_wave)

    # --- sort by date ---
    sorted_idx = np.argsort(all_dates)
    dates = np.array(all_dates)[sorted_idx]
    mjd = np.array(all_mjd)[sorted_idx]
    flux_arrays = np.array(all_fluxes)[sorted_idx]
    #print(mjd)
    
    return common_wave, mjd, flux_arrays