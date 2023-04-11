""" 
Make plots of any column from GITM - veritcally integrated!

- Currently only supports a single column


"""


import argparse
import datetime
from utility_programs.read_routines import GITM, aether_read_routines
import os
import time
import glob
import geopandas
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import signal
from tqdm.auto import tqdm
from utility_programs import filters
from utility_programs.plotting_routines import make_a_keo, draw_map
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

    # Test input column:
    col = args.data_var
    flist = np.sort(glob.glob(os.path.join(
        args.gitm_data_path, args.file_type)))
    f = aether_read_routines.read_gitm_file(flist[0])
    if col not in f['vars']:
        raise ValueError('The column you requested is not available',
                         'in the file you chose.\n',
                         'you gave me %s and the available columns are:\n\n' % col, f.vars)

    col_nice = col.replace('[', '').replace(']', '').replace('!', '').replace('(', '').replace(
        ')', '').replace('!', '')  # for writing files.

    # get gitm data!
    global times, gitm_grid, gitm_tec
    if args.file_type == '2DANC*.bin':
        times, gitm_grid, gitm_bins = GITM.read_gitm_into_nparrays(
            gitm_dir=args.gitm_data_path,
            gitm_file_pattern=args.file_type,
            dtime_storm_start=dtime_storm_start,
            cols=[col],
            t_start_idx=args.plot_start_delta,
            t_end_idx=args.plot_end_delta)
        CALC_TEC = False

    elif args.file_type == '3DALL*.bin':
        times, gitm_grid, gitm_bins = GITM.read_gitm_into_nparrays(
            gitm_dir=args.gitm_data_path,
            gitm_file_pattern=args.file_type,
            dtime_storm_start=dtime_storm_start,
            cols=[col],
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
                gitm_dir=args.gitm_data_path2,
                gitm_file_pattern=args.file_type,
                dtime_storm_start=dtime_storm_start,
                cols=[col],
                t_start_idx=args.plot_start_delta,
                t_end_idx=args.plot_end_delta)
        if CALC_TEC:
            times2, gitm_grid2, gitm_bins2 = GITM.read_gitm_into_nparrays(
                gitm_dir=args.gitm_data_path2,
                gitm_file_pattern=args.file_type,
                dtime_storm_start=dtime_storm_start,
                cols=[col],
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
        if col == '[e-]':
            gitm_tec = integ.simps(
                gitm_bins[:, 0, :, :, :], alts, "avg") * 1e-16
        else:
            gitm_tec = integ.simps(gitm_bins[:, 0, :, :, :], alts, "avg")

        if TWO_FILES:
            if col == '[e-]':
                gitm_tec2 = integ.simps(gitm_bins2[:, 0, :, :, :], alts,
                                        "avg") * 1e-16
            else:
                gitm_tec2 = integ.simps(gitm_bins2[:, 0, :, :, :], alts,
                                        "avg")
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
    fits_gitm = filters.make_fits(gitm_tec)

    if TWO_FILES:
        global fits_gitm2
        fits_gitm2 = filters.make_fits(gitm_tec2)

    if args.figtype == 'all':
        plot_types = ['raw', 'fit', 'diff']
    else:
        plot_types = [args.figtype]

    ylims = [-args.lat_lim, args.lat_lim]

    cbar_lims_dict = {
        'TWO_FILES': {
            'raw': [np.min(gitm_tec-gitm_tec2), np.max(gitm_tec-gitm_tec2)],
            'fit': [np.min(gitm_tec-gitm_tec2), np.max(gitm_tec-gitm_tec2)],
            'diff': [-2, 2]},
        'ONE_FILE': {
            'raw': [np.min(gitm_tec), np.max(gitm_tec)],
            'fit': [np.min(gitm_tec), np.max(gitm_tec)],
            'diff': [-2, 2]}}

    # Now for keos:
    if args.keogram:
        k_extent = [-args.plot_start_delta, args.plot_end_delta,
                    -90, 90]
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
                        tec -= raw2
                        title = "Diff of Raw "+col + \
                            " at lon = " + str(int(real_lon))
                        cbar_lims = cbar_lims_dict['TWO_FILES']['raw']
                    else:
                        title = "Raw  "+col+" at lon = " + str(int(real_lon))
                        cbar_lims = cbar_lims_dict['ONE_FILE']['raw']
                    fname = os.path.join(
                        out_path, 'keo',
                        "raw", "lon" + str(int(real_lon)),
                        col_nice + ".png")

                    make_a_keo(tec, title=title, cbarlims=cbar_lims,
                               cbar_name='Vertically Integrated ' + col,
                               extent=k_extent,
                               fname=fname, OVERWRITE=OVERWRITE)
                    pbar.update()

                if plot_type == 'fit':
                    tec = fit.copy()
                    if TWO_FILES:
                        tec -= fit2
                        title = "Diff of Fit "+col + \
                            " at lon = " + str(int(real_lon))
                        cbar_lims = cbar_lims_dict['TWO_FILES']['fit']
                    else:
                        title = "Fit "+col+" at lon = " + str(int(real_lon))
                        cbar_lims = cbar_lims_dict['ONE_FILE']['fit']
                    cbar_lims = [np.min(tec), np.max(tec)]
                    fname = os.path.join(
                        out_path, 'keo',
                        'fit', "lon" + str(int(real_lon)),
                        col_nice + ".png")

                    make_a_keo(tec, title=title, cbarlims=cbar_lims,
                               cbar_name='Vertically Integrated ' + col,
                               extent=k_extent
                               fname=fname, OVERWRITE=OVERWRITE)
                    pbar.update()

                if plot_type == 'diff':
                    tec = (100*(raw.copy()
                                - fit.copy())
                           / raw.copy())
                    if TWO_FILES:
                        tec -= (100*(raw2 - fit2) / raw2.copy())
                        title = "Diff of % over BG of "+col+" at lon = " \
                                + str(int(real_lon))
                        cbar_lims = cbar_lims_dict['TWO_FILES']['diff']
                    else:
                        title = "% over BG of "+col + \
                            " at lon = " + str(int(real_lon))
                        cbar_lims = cbar_lims_dict['ONE_FILE']['diff']
                    fname = os.path.join(
                        out_path, 'keo',
                        'diff', "lon" + str(int(real_lon)),
                        col_nice + ".png")

                    make_a_keo(tec, title=title, cbarlims=cbar_lims,
                               cbar_name='Vertically Integrated ' + col,
                               extent=k_extent,
                               fname=fname, OVERWRITE=OVERWRITE)

                    pbar.update()
        pbar.close()

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
                        title = "Diff of Raw  %s at %s from storm onset = " % (
                            col, UT_from_Storm_onset(
                                dtime, dtime_storm_start))
                        cbar_lims = cbar_lims_dict['TWO_FILES']['raw']
                    else:
                        title = "Raw %s at %s from storm onset = " % (
                            col, UT_from_Storm_onset(
                                dtime, dtime_storm_start))
                        cbar_lims = cbar_lims_dict['ONE_FILE']['raw']
                    fname = os.path.join(
                        out_path, 'map', "raw", col_nice, str(nt).rjust(3, '0')
                        + ".png")

                    draw_map(tec, title=title, cbarlims=cbar_lims,
                             cbar_label='Vertically Integrated '+col,
                             fname=fname, OVERWRITE=OVERWRITE)
                    pbar.update()

                if plot_type == 'fit':
                    tec = fit.copy()
                    if TWO_FILES:
                        tec -= fit2.copy()
                        title = "Diff of Fit %s at %s from storm onset = " % (
                            col, UT_from_Storm_onset(
                                dtime, dtime_storm_start))
                        cbar_lims = cbar_lims_dict['TWO_FILES']['fit']
                    else:
                        title = "Fit %s at %s from storm onset = " % (
                            col, UT_from_Storm_onset(
                                dtime, dtime_storm_start))
                        cbar_lims = cbar_lims_dict['ONE_FILE']['fit']
                    fname = os.path.join(
                        out_path, 'map', "fit", col_nice, str(nt).rjust(3, '0')
                        + ".png")

                    draw_map(tec, title=title, cbarlims=cbar_lims,
                             cbar_label='Vertically Integrated '+col,
                             fname=fname, OVERWRITE=OVERWRITE)
                    pbar.update()

                if plot_type == 'diff':
                    tec = (100*(raw.copy()-fit.copy())
                           / raw.copy())
                    if TWO_FILES:
                        tec -= (100*(raw2-fit2) / raw2)
                        title = "Diff of percent over BG of %s at %s from storm onset = " % (
                            col, UT_from_Storm_onset(
                                dtime, dtime_storm_start))
                        cbar_lims = cbar_lims_dict['TWO_FILES']['diff']
                    else:
                        title = "% over BG of %s at %s from storm onset = " % (
                            col, UT_from_Storm_onset(
                                dtime, dtime_storm_start))
                        cbar_lims = cbar_lims_dict['ONE_FILE']['diff']
                    fname = os.path.join(
                        out_path, 'map', "diff", col_nice, str(
                            nt).rjust(3, '0')
                        + ".png")

                    draw_map(tec, title=title, cbarlims=cbar_lims,
                             cbar_label='Vertically Integrated '+col,
                             fname=fname, OVERWRITE=OVERWRITE)
                    pbar.update()
        pbar.close()


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
        '--data_var', type=str, default='[e-]', action='store',
        help='Which variable to plot (In GITM output format)')

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
        "--lat_lim", type=float, default=90, action='store',
        help="limit plotted latitudes to this (type=float)")

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