
import pandas as pd
import numpy as np
from multiprocessing import Pool
import math
import os
import shutil
from tqdm.auto import tqdm
import matplotlib.pyplot as plt
import time
import datetime
from scipy import spatial, signal

import sys
import geopandas
import gc

from scipy.interpolate import LinearNDInterpolator, interp1d, griddata
from utility_programs.plot_help import UT_from_Storm_onset
from utility_programs.read_routines import SAMI

import argparse


def main(args):
    data_path = args.sami_data_path

    global dtime_storm_start
    dtime_storm_start = datetime.datetime.strptime(
        args.dtime_storm_start.ljust(14, '0'), '%Y%m%d%H%M%S')

    global dtime_sim_start
    dtime_sim_start = datetime.datetime.strptime(
        args.dtime_sim_start.ljust(14, '0'), '%Y%m%d%H%M%S')

    global times, sami_data
    sami_data, times = SAMI.read_sami_data(sami_data_path=data_path,
                                           dtime_sim_start=dtime_sim_start,
                                           dtime_storm_start=dtime_storm_start,
                                           t_start_idx=args.plot_start_delta,
                                           t_end_idx=args.plot_end_delta,
                                           pbar=True)

    if args.cols == 'all':
        cols_to_plot = sami_data['data'].keys()
    else:
        cols_to_plot = args.cols

    if args.plot_type == 'all':
        plot_type = ['raw', 'bandpass', 'diff']
    else:
        plot_type = [args.plot_type]

    if args.map:
        global world
        world = geopandas.read_file(
            geopandas.datasets.get_path('naturalearth_lowres'))


def call_maps(col, real_time, real_alt, figtype='all', save_or_show='show'):

    itime = np.where(times == real_time)[0][0]
    ialt = np.where(alts == real_alt)[0][0]

    raw = data_dict[col][itime, :, :, ialt]
    fit = fits[col][itime, :, :, ialt]

    diff = 100*(raw - fit)/fit

    fname = os.path.join(sami_map_save_path, col, 'plot_type', str(
        int(real_alt)), str(itime).rjust(3, '0') + '.png')

    ut_diff = UT_from_Storm_onset(real_time)

    plotted = False

    if figtype == 'all' or 'raw' in figtype:
        title = '%s at %s from storm onset \n %i km altitude' % (
            col, ut_diff, int(real_alt))
        make_a_plot(raw, title=title, save_or_show=save_or_show,
                    fname=fname.replace('plot_type', 'raw'))
        plotted = True

    if figtype == 'all' or 'bandpass' in figtype:
        title = '%s at %s from storm onset \n %i km altitude' % (
            col, ut_diff, int(real_alt))
        make_a_plot(fit, title=title, save_or_show=save_or_show,
                    fname=fname.replace('plot_type', 'bandpass'))
        plotted = True

    if figtype == 'all' or 'diff' in figtype:
        title = '%s at %s from storm onset \n %i km altitude' % (
            col, ut_diff, int(real_alt))
        for v in diff_vs:
            make_a_plot(diff, title=title, save_or_show=save_or_show,
                        fname=fname.replace('plot_type', 'filt-' + str(v)),
                        cbar_label='% over background', cbar_lims=[-v, v])

        plotted = True
    plt.close('all')
    gc.collect()
    if not plotted:
        raise ValueError(
            "no plots were made. options are 'diff', 'bandpass', 'raw',",
            "or 'all', you gave %s" % (figtype))


def thread_call_maps(arg_arr):
    call_maps(arg_arr[0], arg_arr[1], arg_arr[2],
              figtype=arg_arr[3], save_or_show=arg_arr[4])


def call_keos(col, real_lon, real_alt, figtype='all', save_or_show='show',
              cbar_lims=None):

    ilon = np.argmin(np.abs(lons - real_lon))
    ialt = np.argmin(np.abs(alts - real_alt))

    real_lon = lons[ilon]
    real_alt = alts[ialt]

    raw = data_dict[col][:, :, ilon, ialt]
    fit = fits[col][:, :, ilon, ialt]

    diff = 100*(raw - fit)/fit

    fname = os.path.join(sami_keo_save_path, col, 'plot_type', str(
        int(real_alt)), 'lon' + str(int(real_lon)) + '.png')

    plotted = False

    if figtype == 'all' or 'raw' in figtype:
        title = '%s at %i (deg) glon \n %i km altitude' % (
            col, int(real_lon), int(real_alt))
        make_a_plot(raw, title=title, save_or_show=save_or_show,
                    fname=fname.replace('plot_type', 'raw'))
        plotted = True

    if figtype == 'all' or 'bandpass' in figtype:
        title = '%s at %i (deg) glon \n %i km altitude' % (
            col, int(real_lon), int(real_alt))
        make_a_plot(fit, title=title, save_or_show=save_or_show,
                    fname=fname.replace('plot_type', 'bandpass'))
        plotted = True

    if figtype == 'all' or 'diff' in figtype:
        title = '%s at %i (deg) glon \n %i km altitude' % (
            col, int(real_lon), int(real_alt))
        for v in diff_vs:
            make_a_plot(diff, title=title, save_or_show=save_or_show,
                        fname=fname.replace('plot_type', 'filt-' + str(v)),
                        cbar_label='% over background', cbar_lims=[-v, v])
        plotted = True

    if not plotted:
        raise ValueError(
            "no plots were made. options are 'diff', 'bandpass', 'raw',",
            "or 'all', you gave %s" % (figtype))


def thread_call_keos(arg_arr):
    call_keos(arg_arr[0], arg_arr[1], arg_arr[2],
              figtype=arg_arr[3], save_or_show=arg_arr[4])


def loop_maps(cols=cols, times=times, alts=alts, thread=True, plottype='all',
              save_or_show='save'):

    if thread:
        arg_arr = []
        for col in cols:
            for itime in times:
                for alt in alts:
                    arg_arr.append([col, itime, alt, plottype, 'save'])

        with Pool(num_pool_workers) as pool:
            with tqdm(desc='making maps!', total=len(arg_arr)) as pbar:
                for _ in pool.imap_unordered(thread_call_maps, arg_arr):
                    pbar.update(1)

    else:
        for col in cols:
            for itime in times:
                for alt in alts:
                    call_maps(col, itime, alt, figtype=plottype,
                              save_or_show=save_or_show)


def loop_keos(cols=cols, lon_keos=lon_keos, alts=alts, thread=True,
              plottype='all'):

    if thread:
        arg_arr = []
        for col in cols:
            for lon in lon_keos:
                for alt in alts:
                    arg_arr.append([col, lon, alt, plottype, 'save'])

        with Pool(num_pool_workers) as pool:
            with tqdm(desc='making keos!', total=len(arg_arr)) as pbar:
                for _ in pool.imap_unordered(thread_call_keos, arg_arr):
                    pbar.update(1)

    else:
        for col in cols:
            for lon in lon_keos:
                for alt in alts:
                    call_keos(col, lon, alt, figtype=plottype,
                              save_or_show=save_or_show)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'dtime_storm_start',
        help='Datetime of storm start. Format YYYYMMDDHHmmss',
        action='store')

    parser.add_argument(
        'dtime_sim_start',
        help='Datetime of the start of the silulations',
        action='store')

    parser.add_argument(
        '-sami_data_path', type=str,
        help='Path to sami data', default='./sami_dir', action='store')

    parser.add_argument(
        '--out_path', type=str,
        help='path to where plots are saved', default='./', action='store')

    parser.add_argument(
        '--cols', nargs="+", type=str,
        help='Which columns to plot. Default: all', default='all')

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
        '--figtype', type=str, action='store', default='all',
        help='Which type of plot to make.' +
        'Options: raw, filt, diffs. Default: all')

    parser.add_argument(
        "--keo_lat_lim", type=float, default=90, action='store',
        help="limit plotted latitudes to this +/- in keos")

    parser.add_argument(
        '--threading',
        help='Use threading. Default: False',
        action='store_true', default=False, required=False)

    parser.add_argument(
        '--num_workers', type=int,
        help='Number of workers to use. Default: 48',
        action='store', default=os.cpu_count(), required=False)

    parser.add_argument(
        "-f", "--file-type", type=str, nargs="+",
        default="3DALL*",
        help="which filetype to plot, e.g. 3DALL* or 2DANC*",)

    parser.add_argument(
        "-o", "--outliers", action="store_true",
        help="do you want to remove outliers")

    parser.add_argument(
        "-d", "--diff", type=int, nargs="+", default=None,
        help="Specify the colorbar limits for 'diff' plots.")

    parser.add_argument(
        "-k", "--keogram", action="store_true",
        help="do you want to make a keogram?")

    parser.add_argument(
        "-m", "--map", action="store_true", help="do you want to make a map?")
