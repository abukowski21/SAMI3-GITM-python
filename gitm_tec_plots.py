import argparse
import datetime
# import gc
from utility_programs.read_routines import GITM
import os
import time
# from multiprocessing import Pool

import geopandas
# import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import signal
from tqdm.auto import tqdm

import scipy.integrate as integ

from utility_programs.plot_help import UT_from_Storm_onset


def main(args):

    if args.map:
        global world
        world = geopandas.read_file(
            geopandas.datasets.get_path('naturalearth_lowres'))

    # Lon to keo:
    global gitm_keo_lons
    gitm_keo_lons = args.keo_lons

    global keo_lat_lim
    keo_lat_lim = args.keo_lat_lim

    global out_path
    out_path = args.out_path

    global OVERWRITE
    OVERWRITE = args.overwrite

    global dtime_storm_start
    dtime_storm_start = datetime.datetime.strptime(
        args.dtime_storm_start.ljust(14, '0'), '%Y%m%d%H%M%S')

    # get gitm data!
    global times, gitm_grid, gitm_tec
    if args.file_type == '2DANC*.bin':
        times, gitm_grid, gitm_bins = GITM.read_gitm_into_nparrays(
            gitm_dir=args.gitm_data_path,
            gitm_file_pattern=args.file_type,
            dtime_storm_start=dtime_storm_start,
            cols=['VerticalTEC'],
            t_start_idx=args.plot_start_delta,
            t_end_idx=args.plot_end_delta)
        CALC_TEC = False

    elif args.file_type == '3DALL*.bin':
        times, gitm_grid, gitm_bins = GITM.read_gitm_into_nparrays(
            gitm_dir=args.gitm_data_path,
            gitm_file_pattern=args.file_type,
            dtime_storm_start=dtime_storm_start,
            cols=['[e-]'],
            t_start_idx=args.plot_start_delta,
            t_end_idx=args.plot_end_delta)
        CALC_TEC = True
    else:
        raise ValueError("File type not recognized",
                         "we currently only accept 2DANC*.bin or 3DALL*.bin")

    if args.gitm_data_path2:
        TWO_FILES = True
        global times2, gitm_grid2, gitm_tec2
        if not CALC_TEC:
            times2, gitm_grid2, gitm_bins2 = GITM.read_gitm_into_nparrays(
                gitm_dir=args.gitm_data_path,
                gitm_file_pattern=args.file_type,
                dtime_storm_start=dtime_storm_start,
                cols=['VerticalTEC'],
                t_start_idx=args.plot_start_delta,
                t_end_idx=args.plot_end_delta)
        if CALC_TEC:
            times2, gitm_grid2, gitm_bins2 = GITM.read_gitm_into_nparrays(
                gitm_dir=args.gitm_data_path,
                gitm_file_pattern=args.file_type,
                dtime_storm_start=dtime_storm_start,
                cols=['[e-]'],
                t_start_idx=args.plot_start_delta,
                t_end_idx=args.plot_end_delta)
    else:
        TWO_FILES = False

    global lats, lons, alts
    lats = np.unique(gitm_grid['latitude'])
    lons = np.unique(gitm_grid['longitude'])
    alts = np.unique(gitm_grid['altitude'])

    if CALC_TEC:
        # gitm_tec = np.zeros([len(times), len(lons), len(lats)])
        # print(gitm_tec.shape, len(lats), len(lons), len(alts))
        gitm_tec = integ.simps(gitm_bins[:, 0, :, :, :], alts, "avg") * 1e-16

        if TWO_FILES:
            gitm_tec2 = integ.simps(gitm_bins2[:, 0, :, :, :], alts,
                                    "avg") * 1e-16
    else:
        gitm_tec = gitm_bins.reshape(
            [len(times), len(lons), len(lats)])
        if TWO_FILES:
            gitm_tec2 = gitm_bins2.reshape(
                [len(times), len(lons), len(lats)])

    global hrs_since_storm_onset
    hrs_since_storm_onset = np.array([(i - pd.Timestamp(dtime_storm_start))
                                      / pd.Timedelta('1 hour') for i in times])

    global fits_gitm
    print("Calculating fits. This will take a moment...")
    fits_gitm = make_fits(gitm_tec)

    if TWO_FILES:
        global fits_gitm2
        fits_gitm2 = make_fits(gitm_tec2)

    if args.figtype == 'all':
        plot_types = ['raw', 'fit', 'diff']
    else:
        plot_types = [args.figtype]

    # Now for keos:
    if args.keogram:
        pbar = tqdm(total=len(gitm_keo_lons) * len(plot_types),
                    desc="Making keograms")
        for real_lon in gitm_keo_lons:
            lon_idx = np.argmin(np.abs(lons - real_lon))
            raw = gitm_tec[:, lon_idx, :].copy()
            fit = fits_gitm[:, lon_idx, :].copy()
            if TWO_FILES:
                raw2 = gitm_tec2[:, lon_idx, :].copy()
                fit2 = fits_gitm2[:, lon_idx, :].copy()
            for plot_type in plot_types:
                if plot_type == 'raw':
                    tec = raw.copy()
                    if TWO_FILES:
                        tec -= raw2.copy()
                        title = "Diff of Raw TEC at lon = {}".format(real_lon)
                    else:
                        title = "Raw TEC at lon = {}".format(real_lon)
                    cbar_lims = [np.min(tec), np.max(tec)]
                    fname = os.path.join(
                        out_path, 'keo',
                        "raw", "lon" + str(int(real_lon)),
                        'tec' + ".png")

                    make_a_keo(tec, title, cbar_lims,
                               cbar_name='Vertically Integrated TEC',
                               fname=fname, OVERWRITE=OVERWRITE)
                    pbar.update()

                if plot_type == 'fit':
                    tec = fit.copy()
                    if TWO_FILES:
                        tec -= fit2.copy()
                        title = "Diff of Fit TEC at lon = {}".format(real_lon)
                    else:
                        title = "Fit TEC at lon = {}".format(real_lon)
                    cbar_lims = [np.min(tec), np.max(tec)]
                    fname = os.path.join(
                        out_path, 'keo',
                        'fit', "lon" + str(int(real_lon)),
                        'tec' + ".png")

                    make_a_keo(tec, title, cbar_lims,
                               cbar_name='Vertically Integrated TEC',
                               fname=fname, OVERWRITE=OVERWRITE)
                    pbar.update()

                if plot_type == 'diff':
                    tec = (100*(raw.copy()
                                - fit.copy())
                           / raw.copy())
                    if TWO_FILES:
                        tec -= (100*(raw2.copy()
                                     - fit2.copy())
                                / raw2.copy())
                        title = "Diff of % over BG of TEC at lon = {}".format(
                            real_lon)
                    cbar_lims = [-5, 5]
                    fname = os.path.join(
                        out_path, 'keo',
                        'diff', "lon" + str(int(real_lon)),
                        'tec' + ".png")

                    make_a_keo(tec, title, cbar_lims,
                               cbar_name='% over BG of TEC',
                               fname=fname, OVERWRITE=OVERWRITE)
                    pbar.update()

    if args.map:
        pbar = tqdm(total=len(times) * len(plot_types),
                    desc="Making maps")
        for nt, dtime in enumerate(times):
            raw = gitm_tec[nt, :, :].copy()
            fit = fits_gitm[nt, :, :].copy()
            if TWO_FILES:
                raw2 = gitm_tec2[nt, :, :].copy()
                fit2 = fits_gitm2[nt, :, :].copy()
            for plot_type in plot_types:
                if plot_type == 'raw':
                    tec = raw.copy()
                    if TWO_FILES:
                        tec -= raw2.copy()
                        title = "Diff of Raw TEC at {} from storm onset".\
                            format(UT_from_Storm_onset(
                                dtime, dtime_storm_start))
                    else:
                        title = "Raw TEC at {} from storm onset".format(
                            UT_from_Storm_onset(dtime, dtime_storm_start))
                    cbar_lims = [np.min(gitm_tec), np.max(gitm_tec)]
                    fname = os.path.join(
                        out_path, 'map',
                        "raw", str(nt).rjust(3, '0') + ".png")

                    make_a_map(tec, title, cbar_lims,
                               cbar_label='Vertically Integrated TEC',
                               fname=fname, OVERWRITE=OVERWRITE)
                    pbar.update()

                if plot_type == 'fit':
                    tec = fit.copy()
                    if TWO_FILES:
                        tec -= fit2.copy()
                        title = "Diff of Fit TEC at {} from storm onset".\
                            format(UT_from_Storm_onset(
                                dtime, dtime_storm_start))
                    else:
                        title = "Fit TEC at {} from storm onset".format(
                            UT_from_Storm_onset(dtime, dtime_storm_start))
                    cbar_lims = [np.min(fits_gitm), np.max(fits_gitm)]
                    fname = os.path.join(
                        out_path, 'map',
                        "fit", str(nt).rjust(3, '0') + ".png")

                    make_a_map(tec, title, cbar_lims,
                               cbar_label='Vertically Integrated TEC',
                               fname=fname, OVERWRITE=OVERWRITE)
                    pbar.update()

                if plot_type == 'diff':
                    tec = (100*(raw.copy()-fit.copy())
                           / raw.copy())
                    break
                    if TWO_FILES:
                        tec -= (100*(raw2.copy()-fit2.copy())
                                / raw2.copy())
                        title = ("Diff of % over BG of TEC at " +
                                 UT_from_Storm_onset(
                                     dtime, dtime_storm_start) +
                                 " from storm onset")
                    else:
                        title = ("% over BG of TEC at {} from storm onset".
                                 format(UT_from_Storm_onset(
                                     dtime, dtime_storm_start)))
                    cbar_lims = [-5, 5]
                    fname = os.path.join(
                        out_path, 'map', "diff", str(nt).rjust(3, '0')
                        + ".png")

                    make_a_map(tec, title, cbar_lims,
                               cbar_label='% over BG of TEC',
                               fname=fname, OVERWRITE=OVERWRITE)


def make_a_keo(
        arr,
        title,
        cbarlims,
        cbar_name,
        y_label="Latitude (deg)",
        x_label="Hours since storm onset",
        save_or_show="save",
        fname=None,
        keo_lat_lim=90,
        plot_extent=None,
        OVERWRITE=False):
    """
    Inputs a data array and then generates a keogram.

    Parameters:
    -----------
    arr: np array
        The data array to be plotted. If grabbing from the gitm array,
        you do not need to transpose.
    extent: tuple/list
        The limits of the plot. [left, right, bottom, top]
    xlabel: string
        self-explanitory
    y-label: string
        self-explanitory
    title: string
        self-explanitory
    cbar limes: tuple/list
        vmin, vmax for the colorbar to be plot.
    cbar_name: string.
        Label for the colorbar.
    save_or_show: string
        Defaults to save. You can instead 'show' the plots.

    """
    if fname is not None and os.path.exists(fname) and save_or_show == "save":
        if not OVERWRITE:
            raise ValueError("We cannot overwrite the file: " + str(fname))
    fig = plt.figure(figsize=(10, 7))

    if plot_extent is None:
        hrs_start = hrs_since_storm_onset[0]
        hrs_end = hrs_since_storm_onset[-1]
        lat_start = -keo_lat_lim
        lat_end = keo_lat_lim
        plot_extent = [hrs_start, hrs_end, lat_start, lat_end]

    plt.imshow(
        arr.T,
        extent=plot_extent,
        aspect="auto",
        cmap="viridis",
        origin="lower",
        vmin=cbarlims[0],
        vmax=cbarlims[1],
    )
    plt.ylabel(y_label)
    plt.xlabel(x_label)
    plt.title(title)
    plt.colorbar(label=cbar_name)

    if save_or_show == "show":
        plt.show()
        plt.close(fig)
    elif save_or_show == "save":
        if not fname:
            raise ValueError("plot save path must be given!")
        else:
            try:
                plt.savefig(fname)
            except FileNotFoundError:
                try:
                    directory_list = os.path.join(fname).split("/")[:-1]
                    os.makedirs(os.path.join(*directory_list))
                    plt.savefig(fname)
                except PermissionError:
                    print("Permission denied. Cannot save plot.")
                    print(" tried writing to: ", fname)
            plt.close("all")
    else:
        raise ValueError(
            'save_or_show input is invalid. Accepted inputs are "save" or',
            '"show", you gave ',
            save_or_show,
        )


def make_a_map(
        data_arr,
        title,
        cbarlims,
        cbar_label=None,
        y_label="Latitude (deg)",
        x_label="Longitude (deg)",
        save_or_show="save",
        fname=None,
        plot_extent=[-180, 180, -90, 90],
        OVERWRITE=False):

    if os.path.exists(fname):
        if not OVERWRITE:
            return

    fig, ax = plt.subplots(figsize=(10, 5))
    world.plot(ax=ax, color="white", edgecolor="black", zorder=1)
    data = ax.imshow(
        data_arr.T,
        cmap="viridis",
        aspect="auto",
        extent=plot_extent,
        origin="lower",
        zorder=10,
        alpha=0.8,
        vmin=cbarlims[0],
        vmax=cbarlims[1],
        interpolation="bicubic",
        interpolation_stage="rgba",)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    if not cbar_label:
        fig.colorbar(data)
    else:
        fig.colorbar(data, label=cbar_label)

    if save_or_show == "show":
        plt.show()
        plt.close()
    elif save_or_show == "save":
        if not fname:
            raise ValueError("plot save path must be given!")
        else:
            fname = fname.replace(" ", "")
            try:
                plt.savefig(fname)
            except FileNotFoundError:
                try:
                    directory_list = os.path.join(fname).split("/")[:-1]
                    os.makedirs(os.path.join(*directory_list))
                    plt.savefig(fname)
                except FileExistsError:
                    # sometimes when we make too many plots in the same
                    # directory, it fails. this fixes that.
                    time.sleep(2)
                    try:
                        plt.savefig(fname)
                    except FileNotFoundError:
                        time.sleep(2)
                        plt.savefig(fname)

            except FileNotFoundError:
                print(fname)
                raise ValueError
            plt.close("all")
    else:
        raise ValueError(
            'save_or_show input is invalid. Accepted inputs are "save" or',
            '"show", you gave ',
            save_or_show,)


def make_filter(lowcut=100, highcut=30):
    """_summary_

    Args:
        lowcut (int, optional): lowcut of the filter. Defaults to 100.
        highcut (int, optional): highcut of the filter. Defaults to 30.

    Returns:
        scipy butterworth filter: the filter with settings defined by the user.
    """
    # Define the cutoff frequencies
    lowcut = 1 / (100 / 60)  # 100 minutes in units of sample^-1
    highcut = 1 / (30 / 60)  # 30 minutes in units of sample^-1

    # Define the Butterworth filter
    nyquist = 0.5 * 5  # 5 minutes is the sampling frequency
    low = lowcut / nyquist
    high = highcut / nyquist
    sos = signal.butter(2, [low, high], btype="bandstop", output="sos")
    return sos


def make_fits(gitm_bins):
    """
    calculate bandpass filter for all data previously read in.

    inputs: nparray of gitmdata

    returns:
    fits: np array indexed at fits[time][col][ilon][ilat][ialt]


    todo: you can thread this by splitting the alts into different threads.
    then just append the fits_full later.

    """
    sos = make_filter()

    filtered_arr = signal.sosfiltfilt(sos, gitm_bins, axis=0)
    return filtered_arr


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="This script will plot keograms of the GITM data.")

    parser.add_argument(
        'dtime_storm_start',
        help='Datetime of storm start. Format YYYYMMDDHHmmss',
        action='store')

    parser.add_argument(
        '-gitm_data_path', type=str,
        help='Path to gitm data', default='./gitm_dir', action='store')

    parser.add_argument(
        '-gitm_data_path2', type=str,
        help='Path to gitm data', default=None, action='store')

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
        '--save_or_show', type=str,
        action='store', default='save', required=False,
        help='Save or show plots. Default: save')

    parser.add_argument(
        '--keo_lons', type=float, nargs="+",
        action='store', default=[-90, 2, 90, -178], required=False,
        help='Lons to plot keograms for. Default: -90,2,90,-178')

    parser.add_argument(
        "--keo_lat_lim", type=float, default=90, action='store',
        help="limit plotted latitudes to this +/- in keos")

    parser.add_argument(
        '--figtype', type=str, action='store', default='all',
        help='Which type of plot to make.' +
        'Options: raw, filt, diffs. Default: all')

    parser.add_argument(
        "-f", "--file-type", type=str,
        default="3DALL*.bin",
        help="which filetype to plot, e.g. 3DALL* or 2DANC* (3DALL default.)",)

    parser.add_argument(
        "-k", "--keogram", action="store_true",
        help="do you want to make a keogram?")

    parser.add_argument(
        "-m", "--map", action="store_true", help="do you want to make a map?")

    parser.add_argument(
        "-o", "--overwrite", action="store_true",
        help="overwrite existing files?")

    args = parser.parse_args()

    main(args)
