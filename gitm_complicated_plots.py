"""
Here we have some more complicated GITM plotting routines.

Currently, the data we can plot are:
    - difference between two storms
    - difference between two bandpass filtered storms

You can set:
    - column (or "TEC")
    - time range
    - oplot winds (TODO)

And the plot types are:
    - keograms
    - maps
    - the above with polar dials too (TODO)
        - option to change what is in each polar dial plot
"""

import gc
import argparse

from utility_programs.plot_help import UT_from_Storm_onset
from utility_programs import plotting_routines
from gitm_basic_plots import read_gitm_into_nparrays, make_fits
from gitm_basic_plots import remove_outliers
import datetime
import numpy as np
import pandas as pd
import os
from tqdm.auto import tqdm
import matplotlib
import scipy.integrate as integ
import matplotlib.pyplot as plt
import glob

matplotlib.use("Agg")


def main(args):

    gitm_colnames_friendly_3DALL = {
        "Rho": "Total Neutral Density",
        "[O(!U3!NP)]": "O(3P)",
        "[O!D2!N]": "O2",
        "[N!D2!N]": "N2",
        "[N(!U4!NS)]": "N(4S)",
        "[NO]": "NO",
        "[He]": "He",
        "[N(!U2!ND)]": "N(2D)",
        "[N(!U2!NP)]": "N(2P)",
        "[H]": "H",
        "[CO!D2!N]": "CO2",
        "[O(!U1!ND)]": "O(1D)",
        "Temperature": "Temperature",
        "V!Dn!N(east)": "Vn(east)",
        "V!Dn!N(north)": "Vn(north)",
        "V!Dn!N(up)": "Vn(up)",
        "V!Dn!N(up,O(!U3!NP))": "Vn(up,O(3P))",
        "V!Dn!N(up,O!D2!N)": "Vn(up,O2)",
        "V!Dn!N(up,N!D2!N)": "Vn(up,N2)",
        "V!Dn!N(up,N(!U4!NS))": "Vn(up,N(4S))",
        "V!Dn!N(up,NO)": "Vn(up,NO)",
        "V!Dn!N(up,He)": "Vn(up,He)",
        "[O_4SP_!U+!N]": "O(4Sp)+",
        "[NO!U+!N]": "NO+",
        "[O!D2!U+!N]": "O2+",
        "[N!D2!U+!N]": "N2+",
        "[N!U+!N]": "N+",
        "[O(!U2!ND)!U+!N]": "O(2D)+",
        "[O(!U2!NP)!U+!N]": "O(2P)+",
        "[H!U+!N]": "H+",
        "[He!U+!N]": "He+",
        "[e-]": "e-",
        "eTemperature": "eTemperature",
        "iTemperature": "iTemperature",
        "V!Di!N(east)": "Vi(east)",
        "V!Di!N(north)": "Vi(north)",
        "V!Di!N(up)": "Vi(up)"}

    gitm_colnames_friendly_2DANC = {
        "Local Time": "Local Time",
        "Solar Zenith Angle": "Solar Zenith Angle",
        "Vertical TEC": "Vertical TEC",
        "AltIntJouleHeating (W/m2)": "Joule Heating",
        # TODO: ADD THIS
    }

    global TWO_FILES
    if args.gitm_data_path2 is not None:
        gitm_path = [args.gitm_data_path, args.gitm_data_path2]
        TWO_FILES = True
    else:
        gitm_path = [args.gitm_data_path]
        TWO_FILES = False

    global cols, gitm_colnames_friendly
    cols = []
    gitm_colnames_friendly = {}
    twodanc = False
    threedall = False

    if type(args.variable) == str:
        args.variable = [args.variable]
    for c in args.variable:
        if c != 'tec':
            if c in gitm_colnames_friendly_3DALL.keys():
                cols.append(c)
                gitm_colnames_friendly[c] = gitm_colnames_friendly_3DALL[c]
                threedall = True
            elif c in gitm_colnames_friendly_2DANC.keys():
                cols.append(c)
                gitm_colnames_friendly[c] = gitm_colnames_friendly_2DANC[c]
                twodanc = True
            else:
                raise ValueError("Column name not recognized!")
        else:
            cols.append('Vertical TEC')
            gitm_colnames_friendly['Vertical TEC'] = 'Vertical TEC'
            twodanc = True

    CALC_TEC = False
    for idir in gitm_path:
        if threedall:
            gitm_files_3dall = np.sort(
                glob.glob(os.path.join(idir, "3DALL*")))
            if len(gitm_files_3dall) == 0:
                raise ValueError("No GITM Binary files found!",
                                 "found these files: ",
                                 os.listdir(idir))
        elif twodanc:
            gitm_files_all = np.sort(
                glob.glob(os.path.join(idir, "*.bin")))
            gitm_files_2danc = gitm_files_all[[
                '2DANC' in x for x in gitm_files_all]]

            if len(gitm_files_2danc) == 0:
                gitm_files_3dall = gitm_files_all[[
                    '3DALL' in x for x in gitm_files_all]]
                if len(gitm_files_3dall) != 0:
                    print('No 2DANC files found, but 3DALL files found.')
                    CALC_TEC = True
                    threedall = True
                    twodanc = False
                    del gitm_colnames_friendly['Vertical TEC']
                    gitm_colnames_friendly['[e-]'] = 'Vertical TEC'
                else:
                    raise ValueError("No GITM Binary files found!",
                                     "found these files: ",
                                     os.listdir(idir), "in directory: ", idir,
                                     'gitm_files_all: ', gitm_files_all,
                                     'bin' in gitm_files_all)

    if threedall and not twodanc:
        gitm_files = [gitm_files_3dall]
        if any('TEC' in x for x in cols):
            plot_cols = cols
            cols[cols.index('Vertical TEC')] = '[e-]'
    elif twodanc and not threedall:
        gitm_files = [gitm_files_2danc]
    elif threedall and twodanc:
        gitm_files = [gitm_files_3dall, gitm_files_2danc]
    else:
        raise ValueError("Something very wrong happened here!")

    if not plot_cols:
        plot_cols = cols

    # lat lims:
    global keo_lat_lim
    keo_lat_lim = args.keo_lat_lim

    # Lon to keo:
    global gitm_keo_lons
    gitm_keo_lons = args.keo_lons

    global out_path
    out_path = args.out_path+'/comp/'

    global dtime_storm_start
    dtime_storm_start = datetime.datetime.strptime(
        args.dtime_storm_start.ljust(14, '0'), '%Y%m%d%H%M%S')

    # THis presently only works if files are output at the same cadence.
    # Would be a very easy fix to add ability for more general case.
    gitm_dtimes = []
    for i in gitm_files[0]:
        yy, MM, dd = i[-17:-15], i[-15:-13], i[-13:-11]
        hr, mm, sec = i[-10:-8], i[-8:-6], i[-6:-4]
        gitm_dtimes.append(
            datetime.datetime(
                int("20" + yy), int(MM), int(dd), int(hr), int(mm), int(sec)))

    plot_start_idx = np.argmin(np.abs(np.array(gitm_dtimes)
                                      - (dtime_storm_start -
                                         datetime.timedelta(
                                             hours=args.plot_start_delta))))

    plot_end_idx = (np.argmin(np.abs(np.array(gitm_dtimes)
                                     - (dtime_storm_start +
                                        datetime.timedelta(
                                            hours=args.plot_end_delta))))
                    if args.plot_end_delta != -1
                    else -1)

    global times, gitm_grid, gitm_vars, gitm_bins
    times, gitm_grid, gitm_vars, gitm_bins = read_gitm_into_nparrays(
        gitm_files[0][plot_start_idx:plot_end_idx], cols)

    if TWO_FILES:
        global times2, gitm_grid2, gitm_vars2, gitm_bins2
        times2, gitm_grid2, gitm_vars2, gitm_bins2 = read_gitm_into_nparrays(
            gitm_files[1][plot_start_idx:plot_end_idx])

    global lats, lons, alts
    lats, lons, alts = (
        np.unique(gitm_grid["latitude"]),
        np.unique(gitm_grid["longitude"]),
        np.unique(gitm_grid["altitude"]))
    if TWO_FILES:
        lats2, lons2, alts2 = (
            np.unique(gitm_grid2["latitude"]),
            np.unique(gitm_grid2["longitude"]),
            np.unique(gitm_grid2["altitude"]))

        if lats2 != lats or lons2 != lons:
            raise ValueError("GITM grids are not the same!")
        print('alts2 has %i unique values out of %i total values'
              % (len(alts2), len(gitm_grid2["altitude"])))

    if CALC_TEC:
        global tec, tec_filt
        tec, tec_filt = get_tec_data(CALC_TEC)
        if TWO_FILES:
            global tec2, tec_filt2
            tec2, tec_filt2 = get_tec_data(CALC_TEC, file=2)

    global hrs_since_storm_onset
    hrs_since_storm_onset = np.array([(i - pd.Timestamp(dtime_storm_start))
                                      / pd.Timedelta('1 hour') for i in times])

    global fits_gitm
    print("Calculating fits. This will take a moment...")
    fits_gitm = make_fits(gitm_bins)
    if TWO_FILES:
        global fits_gitm2
        fits_gitm2 = make_fits(gitm_bins2)

    if args.gitm_alt_idxs == -1:
        gitm_alt_idxs = list(range(len(alts)))
    else:
        gitm_alt_idxs = args.gitm_alt_idxs

    if args.keogram:
        pbar = tqdm(total=len(gitm_alt_idxs) * len(gitm_keo_lons) * len(cols),
                    desc='making keograms')
        for col in plot_cols:
            for real_lon in gitm_keo_lons:
                if col == 'Vertical TEC':
                    call_keos(real_lon=real_lon,
                              namecol=col, save_or_show=args.save_or_show,
                              outliers=args.outliers,
                              figtype=args.data_to_plot,
                              TWO_FILES=TWO_FILES,)
                    pbar.update(len(gitm_alt_idxs))
                else:
                    for alt_idx in gitm_alt_idxs:
                        call_keos(alt_idx=alt_idx, real_lon=real_lon,
                                  namecol=col, save_or_show=args.save_or_show,
                                  outliers=args.outliers,
                                  OVERWRITE=args.OVERWRITE,
                                  figtype=args.data_to_plot,
                                  TWO_FILES=TWO_FILES,)
                        pbar.update(1)
        pbar.close()

    if args.map:
        pbar = tqdm(total=len(gitm_alt_idxs) * len(cols) * len(times),
                    desc='making maps')
        for alt_idx in gitm_alt_idxs:
            for col in plot_cols:
                for dtime_real in times:
                    call_maps(alt_idx=alt_idx, dtime_real=dtime_real,
                              namecol=col,
                              save_or_show=args.save_or_show,
                              outliers=args.outliers,
                              OVERWRITE=args.OVERWRITE,
                              figtype=args.data_to_plot,
                              TWO_FILES=TWO_FILES,)
                    pbar.update(1)
        pbar.close()


def call_maps(
        alt_idx,
        dtime_real=None,
        dtime_index=None,
        numcol=None,
        namecol=None,
        OVERWRITE=False,
        save_or_show="show",
        return_figs=False,
        figtype="all",
        TWO_FILES=False,
        lat_lim=[-90, 90],
        diffs=[1, 2, 3, 5, 10, 30, 50],
        outliers=False):

    # Make sure inputs are correct. either the index or actual value of the
    #    datetime and column to plot can be specified (or both).
    if numcol is None and namecol is not None:
        numcol = cols.index(namecol)
    elif namecol is None and numcol is not None:
        namecol = cols[numcol]
    elif numcol is None and namecol is None:
        raise ValueError("either namecol or numcol must be specified!")

    if dtime_real is None and dtime_index is not None:
        dtime_real = times[dtime_index]
    elif dtime_index is None and dtime_real is not None:
        dtime_index = np.argmin(np.abs(np.array(times) - dtime_real))
    elif dtime_real is None and dtime_index is None:
        raise ValueError("either dtime_index or dtime_real must be specified!")

    # get colorbar limits.
    vmin_bins = np.min(gitm_bins[:, numcol, :, :, alt_idx])
    vmax_bins = np.max(gitm_bins[:, numcol, :, :, alt_idx])

    vmin_fits = np.min(fits_gitm[:, numcol, :, :, alt_idx])
    vmax_fits = np.max(fits_gitm[:, numcol, :, :, alt_idx])

    if type(diffs) != list:

        vmin_diffs = np.min(100 * (fits_gitm[:, numcol, :, :, alt_idx]
                                   - gitm_bins[:, numcol, :, :, alt_idx])
                            / gitm_bins[:, numcol, :, :, alt_idx])
        vmax_diffs = np.max(100 * (fits_gitm[:, numcol, :, :, alt_idx]
                                   - gitm_bins[:, numcol, :, :, alt_idx])
                            / gitm_bins[:, numcol, :, :, alt_idx])

    # get data.

    # get data.
    if namecol == "Vertical TEC":
        data = tec[dtime_index, :, :] .copy()
        bandpass = tec_filt[dtime_index, :, :] .copy()
        percent = 100 * (data - bandpass) / bandpass

    else:
        data = gitm_bins[dtime_index, numcol, :, :, alt_idx].copy()
        bandpass = fits_gitm[dtime_index, numcol, :, :, alt_idx].copy()
        percent = 100 * (data - bandpass) / bandpass

    if TWO_FILES:
        if namecol == "Vertical TEC":
            data2 = tec2[dtime_index, :, :] .copy()
            bandpass2 = tec_filt2[dtime_index, :, :] .copy()
            percent2 = 100 * (data - bandpass) / bandpass
        else:
            data2 = gitm_bins2[dtime_index, numcol, :, :, alt_idx].copy()
            bandpass2 = fits_gitm2[dtime_index, numcol, :, :, alt_idx].copy()
            percent2 = 100 * (data2 - bandpass2) / bandpass2
        data = data - data2
        bandpass = bandpass - bandpass2
        percent = percent - percent2

    raw = gitm_bins[dtime_index, numcol, :, :, alt_idx].copy()
    bandpass = fits_gitm[dtime_index, numcol, :, :, alt_idx].copy()
    real_alt = alts[alt_idx]
    percent = 100 * (bandpass - raw) / raw

    if outliers:
        raw = remove_outliers(raw)
        bandpass = remove_outliers(bandpass)
        percent = remove_outliers(percent)

    made_plot = False

    # raw map
    if figtype == "all" or "raw" in figtype:
        title = (
            gitm_colnames_friendly[namecol]
            + " at "
            + str(round(float(real_alt) / 1000, 0))
            + " km at "
            + UT_from_Storm_onset(dtime_real, dtime_storm_start)
            + " from Storm Start")
        fname = os.path.join(
            out_path, 'maps'
            "raw",
            str(int(real_alt / 1000)),
            gitm_colnames_friendly[namecol],
            str(dtime_index).rjust(3, "0") + ".png",)
        cbarlims = [vmin_bins, vmax_bins]
        plotting_routines.draw_map(raw, title, cbarlims, fname=fname,
                                   save_or_show=save_or_show,
                                   plot_extent=[-180, 180, -90, 90])
        made_plot = True

    # filter map
    if figtype == "all" or "filt" in figtype:
        title = (
            gitm_colnames_friendly[namecol]
            + " at "
            + str(round(float(real_alt) / 1000, 0))
            + " km at "
            + UT_from_Storm_onset(dtime_real, dtime_storm_start)
            + " from Storm Start")
        fname = os.path.join(
            out_path, 'maps'
            "bandpass",
            str(int(real_alt / 1000)),
            gitm_colnames_friendly[namecol],
            str(dtime_index).rjust(3, "0") + ".png",)
        cbarlims = [vmin_fits, vmax_fits]
        cbar_label = "Bandpass Filtered " + gitm_colnames_friendly[namecol]
        plotting_routines.draw_map(
            bandpass,
            title,
            cbarlims,
            OVERWRITE=OVERWRITE,
            save_or_show=save_or_show,
            cbar_label=cbar_label,
            fname=fname,)
        made_plot = True

    # diffs
    if figtype == "all" or "diff" in figtype:
        title = (
            gitm_colnames_friendly[namecol]
            + " at "
            + str(round(float(real_alt) / 1000, 0))
            + " km at "
            + UT_from_Storm_onset(dtime_real, dtime_storm_start)
            + " from Storm Start")
        if type(diffs) == list:
            for v_lim in diffs:
                fname = os.path.join(
                    out_path, 'maps'
                    "diff_set_lims",
                    str(int(real_alt / 1000)),
                    gitm_colnames_friendly[namecol],
                    str(v_lim),
                    str(dtime_index).rjust(3, "0") + ".png",)
                cbarlims = [-v_lim, v_lim]
                cbar_label = "% over Background"
                plotting_routines.draw_map(
                    percent,
                    title,
                    cbarlims,
                    OVERWRITE=OVERWRITE,
                    save_or_show=save_or_show,
                    cbar_label=cbar_label,
                    fname=fname,)
        else:
            fname = os.path.join(
                out_path, 'maps'
                "diff",
                str(int(real_alt / 1000)),
                gitm_colnames_friendly[namecol],
                str(dtime_index).rjust(3, "0") + ".png",)
            cbarlims = [vmin_diffs, vmax_diffs]
            cbar_label = "% over Background"
            plotting_routines.draw_map(
                percent,
                title,
                cbarlims,
                OVERWRITE=OVERWRITE,
                save_or_show=save_or_show,
                cbar_label=cbar_label,
                fname=fname,)
            made_plot = True

    if not made_plot:
        print("No plot made. Check figtype input.")

    plt.close("all")
    gc.collect()


def call_keos(
        real_lon,
        alt_idx=None,
        numcol=None,
        namecol: str = "",
        save_or_show="show",
        return_figs=False,
        figtype="all",
        outliers=False,
        OVERWRITE=False,
        vlims=None,
        CALC_TEC=False,
        TWO_FILES=False,):

    if numcol is None and namecol != "":
        numcol = cols.index(namecol)
    elif namecol == "" and numcol is not None:
        namecol = cols[numcol]
    elif numcol is None and namecol == "":
        raise ValueError("either namecol or numcol must be specified!")

    if namecol != "Vertical TEC" and alt_idx is None:
        raise ValueError("alt_idx must be specified for non-TEC columns!")

    if vlims is None:

        vmin_bins = np.min(gitm_bins[:, numcol, :, :, alt_idx])
        vmax_bins = np.max(gitm_bins[:, numcol, :, :, alt_idx])

        vmin_fits = np.min(fits_gitm[:, numcol, :, :, alt_idx])
        vmax_fits = np.max(fits_gitm[:, numcol, :, :, alt_idx])

        vmin_diffs = np.min(100 * (fits_gitm[:, numcol, :, :, alt_idx]
                                   - gitm_bins[:, numcol, :, :, alt_idx])
                            / gitm_bins[:, numcol, :, :, alt_idx])
        vmax_diffs = np.max(100 * (fits_gitm[:, numcol, :, :, alt_idx]
                                   - gitm_bins[:, numcol, :, :, alt_idx])
                            / gitm_bins[:, numcol, :, :, alt_idx])

    else:
        vmin = -vlims
        vmax = vlims
        vmin_bins = vmin
        vmax_bins = vmax
        vmin_diffs = vmin
        vmax_diffs = vmax
        vmin_fits = vmin
        vmax_fits = vmax

    # get data.
    lon_idx = np.argmin(np.abs(lons - real_lon))
    real_lon = lons[lon_idx]
    if namecol == "Vertical TEC":
        data = tec[:, lon_idx, :].copy()
        bandpass = tec_filt[:, lon_idx, :].copy()
        percent = 100 * (data - bandpass) / bandpass
    else:
        data = gitm_bins[:, numcol, lon_idx, :, alt_idx].copy()
        bandpass = fits_gitm[:, numcol, lon_idx, :, alt_idx].copy()
        percent = 100 * (data - bandpass) / bandpass

    if TWO_FILES:
        if namecol == "Vertical TEC":
            data2 = tec2[:, lon_idx, :].copy()
            bandpass2 = tec_filt2[:, lon_idx, :].copy()
            percent2 = 100 * (data - bandpass) / bandpass
        else:
            data2 = gitm_bins2[:, numcol, lon_idx, :, alt_idx].copy()
            bandpass2 = fits_gitm2[:, numcol, lon_idx, :, alt_idx].copy()
            percent2 = 100 * (data2 - bandpass2) / bandpass2
        data = data - data2
        bandpass = bandpass - bandpass2
        percent = percent - percent2

    real_alt = alts[alt_idx]

    if outliers:
        data = remove_outliers(data)
        bandpass = remove_outliers(bandpass)
        percent = remove_outliers(percent)

    # Filter out data, bandpass, percent if keo_lat_lim:
    if keo_lat_lim:
        good_lats = []
        for n, lat in enumerate(lats):
            if np.abs(lat) <= keo_lat_lim:
                good_lats.append(n)
        good_lats = np.sort(good_lats)
        data = data[:, good_lats]
        bandpass = bandpass[:, good_lats]
        percent = percent[:, good_lats]

    hrs_start = hrs_since_storm_onset[0]
    hrs_end = hrs_since_storm_onset[-1]
    lat_start = -keo_lat_lim
    lat_end = keo_lat_lim
    plot_extent = [hrs_start, hrs_end, lat_start, lat_end]

    made_plot = False

    if figtype == "all" or "filt" in figtype:
        # plain bandpass filter
        title = "Keogram of %s along %i deg Longitude at %i km" % (
            gitm_colnames_friendly[namecol].replace(
                "(", "[").replace(")", "]"),
            real_lon,
            round(real_alt / 1000, 0),)
        color_label = "Bandpass filter"
        # print(out_path, real_alt, real_lon, namecol)
        # print(int(real_alt/1000,0), int(real_lon))
        fname = os.path.join(
            out_path, 'keo',
            "bandpass",
            str(int(real_alt / 1000)),
            "lon" + str(int(real_lon)),
            gitm_colnames_friendly[namecol] + ".png",)
        plotting_routines.make_a_keo(
            bandpass,
            plot_extent=plot_extent,
            title=title,
            OVERWRITE=OVERWRITE,
            cbarlims=(vmin_fits, vmax_fits),
            cbar_name=color_label,
            save_or_show=save_or_show,
            fname=fname,)
        made_plot = True

    if figtype == "all" or "raw" in figtype:
        # plain raw data
        title = "Keogram of %s along %i deg Longitude at %i km" % (
            gitm_colnames_friendly[namecol].replace(
                "(", "[").replace(")", "]"),
            real_lon,
            round(real_alt / 1000, 0),)
        color_label = "Raw data"
        fname = os.path.join(
            out_path, 'keo',
            "raw",
            str(int(real_alt / 1000)),
            "lon" + str(int(real_lon)),
            gitm_colnames_friendly[namecol] + ".png",)
        plotting_routines.make_a_keo(
            data,
            title,
            OVERWRITE=OVERWRITE,
            plot_extent=plot_extent,
            cbarlims=(vmin_bins, vmax_bins),
            cbar_name=color_label,
            save_or_show=save_or_show,
            fname=fname,)
        made_plot = True

    if figtype == "all" or "diff" in figtype:
        # (bandpass - raw)/ raw
        title = "Keogram of %s along %i deg Longitude at %i km" % (
            gitm_colnames_friendly[namecol].replace(
                "(", "[").replace(")", "]"),
            real_lon,
            round(real_alt / 1000, 0),)
        color_label = "% over bandpass filter"
        fname = os.path.join(
            out_path, 'keo',
            "percent-over-filter",
            str(int(real_alt / 1000)),
            "lon" + str(int(real_lon)),
            gitm_colnames_friendly[namecol] + ".png",)
        plotting_routines.make_a_keo(
            percent,
            title,
            OVERWRITE=OVERWRITE,
            plot_extent=plot_extent,
            cbarlims=(vmin_diffs, vmax_diffs),
            cbar_name=color_label,
            save_or_show=save_or_show,
            fname=fname,)
        made_plot = True

    if not made_plot:
        print("nothing made")


def get_tec_data(CALC_TEC, file=1):
    """get tec for a given lat/lon/alt

    Args:
        CALC_TEC (bool): to calculate or not, set automatically.

    Returns:
        data, bandpass, percent: np.arrays


    TODO:
    loop over lons, alts, and lats
    make this work nicely.

    """
    if file == 1:
        b = gitm_bins
    elif file == 2:
        b = gitm_bins2
    if not CALC_TEC:
        numcol = cols.index('Vertical TEC')
        data = b[:, numcol, :, :].copy()
        bandpass = b[:, numcol, :, :].copy()
    else:
        numcol = cols.index('[e-]')
        data = np.zeros_like(b[:, numcol, :, :, 0])
        for ilat in range(len(lats)):
            for ilon in range(len(lons)):
                for itime in range(len(times)):
                    vtec = integ.simps(b[
                        itime, gitm_vars.index('[e-]'), ilon, ilat, :],
                        gitm_grid['altitude'][ilon, ilat, :], "avg")
                data[itime, ilon, ilat] = vtec * 1e-16
        bandpass = make_fits(data)
    return data, bandpass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Make different types of plots of GITM data.")

    parser.add_argument(
        'dtime_storm_start',
        help='Datetime of storm start. Format YYYYMMDDHHmmss',
        action='store')

    parser.add_argument(
        '--variable', nargs='+', default="tec",
        help="which variable to plot. Default is tec.")

    parser.add_argument(
        '-gitm_data_path', type=str,
        help='Path to gitm data', default='./gitm_dir', action='store')

    parser.add_argument(
        '-gitm_data_path2', type=str, default=None, action='store',
        help='Path to another gitm data directory',)

    parser.add_argument(
        '--out_path', type=str,
        help='path to where plots are saved', default='./', action='store')

    parser.add_argument(
        '--plot_start_delta', type=int,
        action='store', default=-1, required=False)

    parser.add_argument(
        '--plot_end_delta', type=int,
        action='store', default=-1, required=False)

    parser.add_argument(
        '--keo_lons', type=float, nargs="+",
        action='store', default=[-90, 2, 90, -178], required=False,
        help='Lons to plot keograms for. Default: -90,2,90,-178')

    parser.add_argument(
        '--data_to_plot', type=str, default="diff",
        help="""What data to plot. Options are: diff, bandpass, or raw.
                Default is diff. Use raw or bandpass with two files.
                NOTE: this is the type of data plot if we choose 2 paths.""")

    parser.add_argument(
        '--gitm_alt_idxs', type=int, nargs="*",
        default=[5, 10, 15, 22, 30, 45],
        help='Which altitudes to plot. Default: 5,10,15,22,30,45')

    parser.add_argument(
        "--keo_lat_lim", type=float, default=90, action='store',
        help="limit plotted latitudes to this +/- in keos")

    parser.add_argument(
        "-m", "--map", action="store_true", help="do you want to make a map?")

    parser.add_argument(
        "-k", "--keogram", action="store_true",
        help="do you want to make a keogram?")

    parser.add_argument(
        '--save_or_show', type=str,
        action='store', default='save', required=False,
        help='Save or show plots. Default: save')

    parser.add_argument(
        "-o", "--outliers", action="store_true",
        help="do you want to remove outliers")

    parser.add_argument(
        "--OVERWRITE", action="store_true", default=False,
        help="overwrite existing files?")

    args = parser.parse_args()

    main(args)
