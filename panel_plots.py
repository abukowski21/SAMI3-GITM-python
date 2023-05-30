import argparse
import glob
import os
from datetime import timedelta

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import utility_programs.filters as filters
from utility_programs.plotting_routines import panel_plot
from utility_programs.read_routines.GITM import auto_read as auto_read_gitm
from utility_programs.read_routines.GITM import gitm_times_from_filelist
from utility_programs.read_routines.SAMI import auto_read as auto_read_sami


def get_fit(array, lowcut=80, highcut=40):

    return xr.apply_ufunc(filters.make_fits, array.load(), lowcut, highcut)


def get_diffs(array):
    return 100 * (array - get_fit(array)) / array


def main(
        directories,
        prefixes,
        cuts,
        samicols,
        gitmcols,
        start_plotting,
        end_plotting,
        outdir=None,
        progress=True,
        use_dask=False,):

    # read in data
    data = {}
    dirs_read = 0

    # make times, get start & end idx.
    if isinstance(directories, str):
        directories = [directories]
    all_files = glob.glob(os.path.join(directories[0], 'GITM*.nc'))
    if len(all_files) == 0:
        all_files = glob.glob(directories[0] + 'SAMI_REGRID*.nc')
    if len(all_files) == 0:
        raise ValueError('No netCDF files found in %s' % directories[0])
    all_times = np.sort(gitm_times_from_filelist(all_files))

    t_start = all_times[0] + timedelta(hours=start_plotting)
    t_end = all_times[0] + timedelta(hours=end_plotting)
    start_idx = np.argmin(np.abs(all_times - t_start))
    end_idx = np.argmin(np.abs(all_times - t_end))

    # check if GITM or SAMI, then read all the data in to a dict

    if prefixes is not None:
        if not isinstance(prefixes, list):
            prefixes = [prefixes]
        for fpath in directories:
            for prefix in prefixes:

                if 'sami' in prefix.lower():
                    data[str(dirs_read) + '-sami'] = auto_read_sami(
                        fpath + prefix, cols=samicols,
                        start_idx=start_idx, end_idx=end_idx,
                        use_dask=use_dask)
                elif 'gitm' in prefix.lower():
                    data[str(dirs_read) + '-gitm'] = auto_read_gitm(
                        fpath + prefix, cols=gitmcols,
                        start_idx=start_idx, end_idx=end_idx,
                        use_dask=use_dask)
                else:
                    raise ValueError('Prefix %s is not valid' % prefix)
            dirs_read += 1

    else:
        for fpath in directories:
            data[str(dirs_read) + '-sami'] = auto_read_sami(
                fpath, cols=samicols,
                start_idx=start_idx, end_idx=end_idx, use_dask=use_dask)
            data[str(dirs_read) + '-gitm'] = auto_read_gitm(
                fpath, cols=gitmcols,
                start_idx=start_idx, end_idx=end_idx, use_dask=use_dask)
            dirs_read += 1

    # loop thru all plots
    # Set up locations to plot:
    alts = [225, 400, 660, 850]
    lons = np.linspace(0, 360, 7)[:-1]
    lats = np.linspace(-45, 45, 6)
    afewtimes = [10, 20, 30, 40, 50, 60, 70, 80]
    lonlatalt = [
        [0, 50, 550],
        [180 + 50, 40, 150],
        [180, 35, 250],
        [60, 40, 650],
        [300, -20, 450],
        [300, 20, 450],
        [300, 40, 500],
        [300, -40, 500],
        [320, 45, 200],
        [320, 45, 400],
        [320, 45, 600],
        [320, -45, 200],
        [320, -45, 400],
        [320, -45, 600]]

    # Lon keos:
    if cuts == ['all'] or 'lon' in cuts:
        for a in alts:
            out_name = os.path.join(outdir, 'lon-keos-alt_%i' % int(a))
            # check if out_name directory exists:
            if not os.path.isdir(out_name):
                os.makedirs(name=out_name)

            # Now do the plots
            for model_key in data.keys():
                if 'sami' in model_key:
                    for s_col in samicols:
                        panel_plot(get_diffs(
                            data[model_key][s_col].sel(
                                alt=a, method='nearest')),
                            x='time', y='lat', wrap_col='lon',
                            col_wrap=2, plot_vals=lons, suptitle=model_key,
                            out_fname=out_name + '/%s-%s.png' %
                            (model_key, s_col,))
                if 'gitm' in model_key and a < 700:
                    for g_col in gitmcols:
                        panel_plot(data[model_key][g_col].sel(
                            alt=a, method='nearest'),
                            x='time', y='lat', wrap_col='lon',
                            col_wrap=2, plot_vals=lons, suptitle=model_key,
                            out_fname=out_name + '/%s-%s.png' %
                            (model_key, g_col,))

    # Lat keos:
    if cuts == ['all'] or 'lat' in cuts:
        for a in alts:
            out_name = os.path.join(outdir, 'lat-keos-alt_%i' % int(a))
            # check if out_name directory exists:
            if not os.path.isdir(out_name):
                os.makedirs(name=out_name)

            # Now do the plots
            for model_key in data.keys():
                if 'sami' in model_key:
                    for s_col in samicols:
                        panel_plot(get_diffs(
                            data[model_key][s_col].sel(
                                alt=a, method='nearest')),
                            x='time', y='lon', wrap_col='lat',
                            col_wrap=2, plot_vals=lats, suptitle=model_key,
                            out_fname=out_name + '/%s-%s.png' %
                            (model_key, s_col,))
                if 'gitm' in model_key and a < 700:
                    for g_col in gitmcols:
                        panel_plot(data[model_key][g_col].sel(
                            alt=a, method='nearest'),
                            x='time', y='lon', wrap_col='lat',
                            col_wrap=2, plot_vals=lats, suptitle=model_key,
                            out_fname=out_name + '/%s-%s.png' %
                            (model_key, g_col,))

    # Maps:
    if cuts == ['all'] or 'alt' in cuts:
        for a in alts:
            out_name = os.path.join(outdir, 'maps-alt_%i' % int(a))
            # check if out_name directory exists:
            if not os.path.isdir(out_name):
                os.makedirs(name=out_name)

            # Now do the plots
            for model_key in data.keys():
                if 'sami' in model_key:
                    for s_col in samicols:
                        panel_plot(get_diffs(
                            data[model_key][s_col].sel(
                                alt=a, method='nearest')),
                            x='lon', y='lat', wrap_col='time',
                            do_map=True,
                            plot_vals=afewtimes,
                            suptitle=model_key,
                            out_fname=out_name + '/%s-%s.png' %
                            (model_key, s_col,))
                if 'gitm' in model_key and a < 700:
                    for g_col in gitmcols:
                        panel_plot(data[model_key][g_col].sel(
                            alt=a, method='nearest'),
                            x='lon', y='lat', wrap_col='time',
                            do_map=True,
                            plot_vals=afewtimes,
                            suptitle=model_key,
                            out_fname=out_name + '/%s-%s.png' %
                            (model_key, g_col,))

    # Single Point plots:

    nrows = len(directories)

    for plotvar in (gitmcols + samicols):
        if plotvar in samicols:
            model = 'sami'
        elif plotvar in gitmcols:
            model = 'gitm'
        else:
            raise ValueError(plotvar, 'not found in model outputs')

        for pt in lonlatalt:
            ilon, ilat, ialt = pt
            fig, axs = plt.subplots(nrows=nrows, ncols=2,
                                    figsize=(11, 1 + nrows * 4),)
            out_name = os.path.join(outdir, 'single-point-%s-%s-%i-%i-%i' %
                                    (model, plotvar, ilon, ilat, ialt))

            for row in range(nrows):
                get_diffs(data[str(row) + '-' + model][plotvar].sel(
                    lon=ilon, alt=ialt, lat=ilat, method='nearest')).plot(
                    ax=axs[row, 0], label=str(row) + ' (diff)')

                get_fit(data[str(row) + '-' + model][plotvar].sel(
                    lon=ilon, alt=ialt, lat=ilat, method='nearest')).plot(
                    ax=axs[row, 1], label=str(row) + ' (fit)')
                data[str(row) + '-' + model][plotvar].sel(
                    lon=ilon, alt=ialt, lat=ilat, method='nearest').plot(
                    ax=axs[row, 1], label=str(row) + ' (raw)')

            for ax in axs.flatten():
                ax.legend()
            fig.tight_layout()

            if out_name is None:
                plt.show()
                plt.close('all')
            else:
                plt.savefig(out_name)
                plt.close('all')

    # exit?

    return


if __name__ == '__main__':

    args = argparse.ArgumentParser()

    args.add_argument('--dirs', type=str, nargs='*',
                      help='Directories of postprocessed data to read from')

    args.add_argument(
        '--out_dir',
        type=str,
        default=None,
        help='Directory to save plots to.'
        'Default: None (MUST BE SET.) (TODO: MAKE THIS MANDATORY.)')

    args.add_argument(
        '--prefixes',
        default=None,
        nargs='*',
        help='Prefixes of postprocessed data to read from'
        'default=None, which will read all files in the directory')

    args.add_argument(
        '--cuts',
        type=str,
        default='all', nargs='*',
        help="Which cuts to use."
        " Options are 'all', 'alt', 'lon', 'lat', 'lonlatalt'"
        " lon/lat cut will make keos, alt makes maps, lonlatalt makes"
        " plots of a single pt throughout time. Default: 'all'"
        " Change cut values in source code.")

    args.add_argument('--samicols', type=str, nargs='*',
                      default=['edens'],
                      help="Which columns to plot from SAMI data."
                      "Default: ['edens',] (only reads SAMI_REGRID)")
    args.add_argument('--gitmcols', type=str, nargs='*',
                      default=['rho'],
                      help="Which columns to plot from GITM."
                      " Default: ['rho']  NOTE: cannot handle altitude"
                      " integrated quantities correctly.")

    args.add_argument(
        '-s', '--start_plotting',
        type=int,
        default=24,
        help='hours after simulation start to begin plotting. Default: 24'
        ' NOTE: these options, the times (sim runtime & relative storm start'
        ' time) must be the same across all directories')
    args.add_argument(
        '-e', '--end_plotting',
        type=int,
        default=48,
        help='hours after simulation start to stop plotting. Default: 48')

    args.add_argument(
        '--dask', action='store_true',
        help='Use dask to read in data. Default: False')

    args = args.parse_args()

    main(directories=args.dirs,
         outdir=args.out_dir,
         prefixes=args.prefixes,
         cuts=args.cuts,
         samicols=args.samicols,
         gitmcols=args.gitmcols,
         start_plotting=args.start_plotting,
         end_plotting=args.end_plotting,
         use_dask=args.dask,)
