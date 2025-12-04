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
    date_slider = mo.ui.range_slider(
        start=min_time_mjd,
        stop=today_mjd,
        step=0.1,
        value=[min_time_mjd,today_mjd],
        label='Dates',
        orientation='vertical',
        debounce=True
    )

    filenames = [os.path.basename(x) for x in glob.glob('BeSS_Files/*')] 
    filebox = mo.ui.dropdown(
        options=filenames,
        value=filenames[0],
        label='Star'
    )

    filebox2 = mo.ui.dropdown(
        options=filenames,
        value=filenames[1],
        label='Star 2'
    )

    featurebox = mo.ui.dropdown(
        options=['All', 'Hα', 'He I'],
        value='Hα',
        label='Feature:'
    )
    tbinbox = mo.ui.number(
        value=365,
        label='Time bin in days'
    )

    velaxis = mo.ui.checkbox(label="Velocity axis")

    ybutton_options = ["MJD","YYYY-MM-DD"]
    ybutton = mo.ui.radio(
        options=ybutton_options,
        value=ybutton_options[0],
        label="Y-axis Units"
    )

    comparecheck = mo.ui.checkbox(label="Compare Stars")

    phasecheck = mo.ui.checkbox(label="Plot phase")

    phasebinbox = mo.ui.number(
        value=0.1,
        start=1e-3,
        label="Phase bin"
    )

    time0box = mo.ui.number(
        value=51544.5,
        label=r'$T_0$ in MJD'
    )

    periodbox = mo.ui.number(
        value=100,
        label='Period in days'
    )

    heliocheck = mo.ui.checkbox(label="Heliocentric Correction")

    cMaps=['viridis', 'plasma', 'inferno', 'magma', 'cividis']
    cmapbox = mo.ui.dropdown(
        options=cMaps,
        value=cMaps[2],
        label="Color Map"
    )

    vlimits = mo.ui.range_slider(
        start=0,
        stop=1,
        step=0.05,
        value=[0,1],
        label='Legend',
        orientation='vertical',
        debounce=True
    )
    return (
        cmapbox,
        comparecheck,
        date_slider,
        featurebox,
        filebox,
        filebox2,
        heliocheck,
        periodbox,
        phasebinbox,
        phasecheck,
        tbinbox,
        time0box,
        velaxis,
        vlimits,
        ybutton,
    )


@app.cell
def _(filebox, os):
    file_paths = os.path.expanduser("BeSS_Files/" + filebox.value + '/*.fits')
    return (file_paths,)


@app.cell
def _(comparecheck, filebox2, os):
    if comparecheck.value:
        file_paths2 = os.path.expanduser("BeSS_Files/" + filebox2.value + '/*.fits')
    return (file_paths2,)


@app.cell
def _(fh, file_paths):
    bess_files_sorted, bess_mjd_sorted = fh.sort_files(file_paths, spec_feature=6562.8)
    return (bess_files_sorted,)


@app.cell
def _(comparecheck, fh, file_paths2):
    if comparecheck.value:
       bess_files_sorted2, bess_mjd_sorted2 = fh.sort_files(file_paths2, spec_feature=6562.8)
    return (bess_files_sorted2,)


@app.cell
def _(bess_files_sorted, date_slider, dp, featurebox, heliocheck, tbinbox):
    x_positions, dates_arrays, flux_arrays = dp.data_setup(bess_files_sorted, feature=featurebox.value, dates=date_slider.value, bin_days=tbinbox.value, heliocor=heliocheck.value)
    return dates_arrays, flux_arrays, x_positions


@app.cell
def _(mo, x_positions):
    wave_range = mo.ui.range_slider(
        debounce=True,
        start=x_positions[0],
        stop=x_positions[-1],
        step=0.1,
        value=[x_positions[0],x_positions[-1]],
        label="X-axis Limits"
    )
    return (wave_range,)


@app.cell
def _(bess_files_sorted2, comparecheck, date_slider, dp, featurebox, tbinbox):
    if comparecheck.value:
        x_positions2, dates_arrays2, flux_arrays2 = dp.data_setup(bess_files_sorted2, feature=featurebox.value, dates=date_slider.value, bin_days=tbinbox.value)
    return dates_arrays2, flux_arrays2, x_positions2


@app.cell
def _(
    dates_arrays,
    dp,
    flux_arrays,
    periodbox,
    phasebinbox,
    phasecheck,
    tbinbox,
    time0box,
):
    # --- Bin data ---
    if phasecheck.value:
        binned_y, flux_matrix = dp.phase_binning(dates_arrays, flux_arrays, time0box.value, periodbox.value, phasebinbox.value)
    else:
        binned_y, flux_matrix = dp.bin_times_and_spectra(dates_arrays, flux_arrays, tbinbox.value)
    return binned_y, flux_matrix


@app.cell
def _(
    comparecheck,
    dates_arrays2,
    dp,
    flux_arrays2,
    periodbox,
    phasebinbox,
    phasecheck,
    tbinbox,
    time0box,
):
    if comparecheck.value:
        if phasecheck.value:
            binned_y2, flux_matrix2 = dp.phase_binning(dates_arrays2, flux_arrays2, time0box.value, periodbox.value, phasebinbox.value)
        else:
            binned_y2, flux_matrix2 = dp.bin_times_and_spectra(dates_arrays2, flux_arrays2, tbinbox.value)
    return binned_y2, flux_matrix2


@app.cell
def _(
    binned_y,
    binned_y2,
    cmapbox,
    comparecheck,
    flux_matrix,
    flux_matrix2,
    pf,
    phasecheck,
    velaxis,
    vlimits,
    wave_range,
    x_positions,
    x_positions2,
    ybutton,
):
    if ybutton.value == "MJD" or phasecheck.value:
        yfor = "mjd"
    else:
        yfor = "isot"
    if comparecheck.value:
        fig1, ax1 = pf.plot_bess_binned_dynamic(x_positions, binned_y, flux_matrix, wave_range.value, vlimits.value,
                                                yformat=yfor,
                                                cmap=cmapbox.value,
                                                phaseplot=phasecheck.value,
                                                toplabel=velaxis.value,
                                                xvals2 = x_positions2,
                                                yvals2 = binned_y2,
                                                fluxes2 = flux_matrix2
                                                )
    else:
        fig1, ax1 = pf.plot_bess_binned_dynamic(x_positions, binned_y, flux_matrix, wave_range.value, vlimits.value,
                                                yformat=yfor,
                                                cmap=cmapbox.value,
                                                phaseplot=phasecheck.value,
                                                toplabel=velaxis.value
                                               )
    return (fig1,)


@app.cell
def _(Time, date_slider):
    lowerdate = Time(date_slider.value[0], format="mjd").to_value(format="isot", subfmt="date")
    upperdate = Time(date_slider.value[1], format="mjd").to_value(format="isot", subfmt="date")
    return lowerdate, upperdate


@app.cell
def _(
    cmapbox,
    comparecheck,
    date_slider,
    featurebox,
    fig1,
    filebox,
    filebox2,
    heliocheck,
    lowerdate,
    mo,
    periodbox,
    phasebinbox,
    phasecheck,
    tbinbox,
    time0box,
    upperdate,
    velaxis,
    vlimits,
    wave_range,
    ybutton,
):
    if comparecheck.value:
        if phasecheck.value:
            output = [mo.mpl.interactive(fig1),
                      mo.vstack([cmapbox,filebox, filebox2, featurebox, phasebinbox, time0box, periodbox, velaxis, phasecheck, comparecheck, heliocheck])
                     ]
        else:
            output = [mo.mpl.interactive(fig1),
                      mo.vstack([cmapbox,filebox, filebox2, featurebox, tbinbox, ybutton,wave_range,
                                 mo.hstack([date_slider,mo.vstack(["","",upperdate,"","","","",lowerdate],gap=1.15),vlimits],widths=[0,0,1]),
                                 velaxis, phasecheck, comparecheck, heliocheck
                                ])
                     ]
        widthvals = [1,1]


    else:
        if phasecheck.value:
            output = [mo.mpl.interactive(fig1),
                      mo.vstack([cmapbox,filebox, featurebox, phasebinbox, time0box, periodbox, wave_range,
                                 mo.hstack([mo.vstack([velaxis, phasecheck, comparecheck, heliocheck]),vlimits], widths=[1,3])
                                ])
                     ]
        else:
            output = [mo.mpl.interactive(fig1),
                      mo.vstack([cmapbox,filebox, featurebox, tbinbox, ybutton,wave_range,
                                 mo.hstack([date_slider,mo.vstack(["","",upperdate,"","","","",lowerdate],gap=1.15),vlimits],widths=[0,0,1]),
                                 velaxis, phasecheck, comparecheck, heliocheck
                                ])
                     ]
        widthvals = [1,1]
    return output, widthvals


@app.cell
def _(mo, output, widthvals):
    mo.hstack(output,widths=widthvals)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
