import marimo

__generated_with = "0.15.2"
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
    return Time, datetime, glob, mo, os, pytz


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
    date_slider = mo.ui.range_slider(start=min_time_mjd,
                                    stop=today_mjd,
                                    step=0.1,
                                    value=[min_time_mjd,today_mjd],
                                    label='Date Range',
                                    orientation='vertical',
                                    debounce=True
    )
    filebox = mo.ui.dropdown(options=[os.path.basename(x) for x in glob.glob('BeSS_Files/*')],
                            value=[os.path.basename(x) for x in glob.glob('BeSS_Files/*')][0],
                            label='Star:'
    )
    featurebox = mo.ui.dropdown(options=['All', 'Hα', 'He I'],
                                value='Hα',
                                label='Feature:'
    )
    tbinbox = mo.ui.number(value=365,
                        start=1,
                        stop=1e4,
                        step=1,
                        label='Time bin in days'
    )

    xbutton = mo.ui.radio(options=["wavelength", "velocity"],
                             value="wavelength",
                             label="X-axis Units"
    )
    return date_slider, featurebox, filebox, tbinbox, xbutton


@app.cell
def _():
    hello =1 
    return


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
def _(bess_files_sorted, date_slider, dp, featurebox, tbinbox):
    x_positions, dates_arrays, flux_arrays = dp.data_setup(bess_files_sorted, feature=featurebox.value, dates=date_slider.value, bin_days=tbinbox.value)
    return dates_arrays, flux_arrays, x_positions


@app.cell
def _(dates_arrays, dp, flux_arrays, tbinbox):
    # --- Bin data ---
    binned_dates, flux_matrix = dp.bin_times_and_spectra(dates_arrays, flux_arrays, tbinbox.value)
    return binned_dates, flux_matrix


@app.cell
def _(binned_dates, filebox, flux_matrix, pf, tbinbox, x_positions, xbutton):
    fig, ax = pf.plot_bess_binned_dynamic(x_positions, binned_dates, flux_matrix,
                                        target_name=filebox.value,
                                        yformat="mjd",
                                        bin_days=tbinbox.value,
                                        xformat = xbutton.value,
                                        cmap="inferno")
    return (fig,)


@app.cell
def _(date_slider, featurebox, fig, filebox, mo, tbinbox, xbutton):
    mo.hstack([mo.mpl.interactive(fig),
               mo.vstack([filebox, featurebox, tbinbox, xbutton, date_slider])],
              widths=[1,1])
    return


@app.cell
def _(dates_arrays):
    len(dates_arrays)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
