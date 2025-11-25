import marimo

__generated_with = "0.13.15"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    import os
    import numpy as np
    import glob
    import matplotlib.pyplot as plt
    import pytz
    from astropy.time import Time
    from datetime import datetime, timedelta
    return Time, datetime, glob, mo, np, os, pytz


@app.cell
def _():
    import file_handler as fh
    import data_prep as dp
    import plotting_functions as pf
    return dp, fh, pf


@app.cell
def _(Time, datetime, pytz):
    # Get the current local time and convert to UTC
    local_time = datetime.now(pytz.timezone('America/New_York'))
    utc_time = local_time.astimezone(pytz.utc)

    # Create an astropy Time object with the UTC time
    t = Time(utc_time, scale='utc')

    # Convert to MJD format
    today_mjd = t.mjd

    print(f"Today's date and time (UTC): {t.iso}")
    print(f"Modified Julian Date (MJD): {today_mjd}")

    min_time = Time("1980-01-01", scale="utc", format="isot")
    min_time_mjd = min_time.mjd
    return min_time_mjd, today_mjd


@app.cell
def _(glob, min_time_mjd, mo, os, today_mjd):
    date_slider = mo.ui.range_slider(
        start=min_time_mjd,
        stop=today_mjd,
        step=0.1,
        value=[min_time_mjd,today_mjd],
        label='Date Range',
        orientation='vertical',
        debounce=True
    )



    # List out all available data folders   
    bess_data_folder = 'BeSS_Files' 
    filenames = [os.path.basename(x) for x in glob.glob(f'{bess_data_folder}/*')] 

    filebox = mo.ui.dropdown(
        options=filenames,
        value=filenames[0],
        label='Star:'
    )

    print(f"Folders found in '{bess_data_folder}':")
    for filename in filenames: 
        print(filename) 



    featurebox = mo.ui.dropdown(
        options=['All', 'Hα', 'He I'],
        value='Hα',
        label='Feature:'
    )

    tbinbox = mo.ui.number(
        value=365,
        start=1,
        stop=1e4,
        step=1,
        label='Time bin in days'
    )

    xbutton = mo.ui.radio(
        options=["wavelength", "velocity"],
        value="wavelength",
        label="X-axis Units"
    )

    ybutton_options = ["ISO (YYYY-MM-DD)", "Modified Julian Date (MJD)"]
    ybutton = mo.ui.radio(
        options=ybutton_options,
        value=ybutton_options[0],
        label="Y-axis Units"
    )


    return date_slider, featurebox, filebox, tbinbox, xbutton, ybutton


@app.cell
def _(flux_arrays, mo, np):
    # Add slider that determines the vmin and vmax of the colorbar 

    # Use range of values in fluxes to set the default (after ignoring any points that equal exactly 0.0)
    new_fluxes = np.copy(flux_arrays) 
    new_fluxes[np.where(flux_arrays==0.0)] = np.nan    

    vrange_slider = mo.ui.range_slider(
        start = 0.0, 
        stop = 1.0, 
        step = 0.01, 
        value = [np.nanmin(new_fluxes), np.nanmax(new_fluxes)], 
        label = "Colorbar range", 
        orientation = "vertical"
    )
    return (vrange_slider,)


@app.cell
def _(Time, date_slider):
    # Create string that displays the currently selected date range in ISO (aka human-readable) format, just under the date range slider. 

    def mjd_to_iso(mjd_value): 
        return Time(mjd_value, format='mjd', scale='utc').to_datetime().date()

    date_range_str = f"{mjd_to_iso(date_slider.value[0])} --- {mjd_to_iso(date_slider.value[1])}"
    return (date_range_str,)


@app.cell
def _(filebox, os):
    file_paths = os.path.expanduser("BeSS_Files/" + filebox.value + '/*.fits')
    return (file_paths,)


@app.cell
def _(fh, file_paths, mo):
    with mo.redirect_stdout():
        bess_files_sorted, bess_mjd_sorted = fh.sort_files(file_paths, spec_feature=6562.8)
    return (bess_files_sorted,)


@app.cell
def _(bess_files_sorted, date_slider, dp, featurebox):
    common_wave, mjd, dates, flux_arrays = dp.data_setup(bess_files_sorted, feature=featurebox.value, dates=date_slider.value) 
    return common_wave, flux_arrays, mjd


@app.cell
def _(
    common_wave,
    dp,
    filebox,
    flux_arrays,
    mjd,
    pf,
    tbinbox,
    vrange_slider,
    xbutton,
):
    # --- Bin data ---
    binned_dates, flux_matrix = dp.bin_times_and_spectra(mjd, flux_arrays, tbinbox.value)

    # Original (Sean's version, requires data to be binned) 
    fig, ax = pf.plot_bess_binned_dynamic(
        common_wave, binned_dates, flux_matrix,
        target_name=filebox.value,
        yformat="mjd",
        bin_days=tbinbox.value,
        xformat = xbutton.value,
        cmap="inferno", 
        vmin = vrange_slider.value[0], 
        vmax = vrange_slider.value[1]
    ) 
    return (fig,)


@app.cell
def _():
    # # New (John's version, does not require data to be binned) 

    # # Use button to decide whether to plot MJD or ISO on y axis 
    # if ybutton.value == ybutton_options[0]: 
    #     yvals = dates 
    # if ybutton.value == ybutton_options[1]: 
    #     yvals = mjd  


    # # Create plot 
    # fig, ax = pf.plot_bess_unbinned(
    #     common_wave, 
    #     yvals, 
    #     flux_arrays, 
    #     xformat = xbutton.value, 
    #     cmap = "inferno", 
    #     vmin = vrange_slider.value[0], 
    #     vmax = vrange_slider.value[1]
    # )


    return


@app.cell
def _(
    date_range_str,
    date_slider,
    featurebox,
    fig,
    filebox,
    mo,
    tbinbox,
    vrange_slider,
    xbutton,
    ybutton,
):
    mo.hstack(
        [
            mo.mpl.interactive(fig), 
            mo.vstack([filebox, featurebox, tbinbox, "\n", xbutton, "\n", ybutton, "\n", date_slider, date_range_str, "\n", vrange_slider])
        ], 
        widths=[1,1], 
    )


    return


if __name__ == "__main__":
    app.run()
