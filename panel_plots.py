import argparse
import glob
import os
from datetime import timedelta

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd

import utility_programs.filters as filters
from utility_programs.utils import str_to_ut
from utility_programs.plotting_routines import panel_plot
from utility_programs.read_routines.GITM import auto_read as auto_read_gitm
from utility_programs.read_routines.GITM import gitm_times_from_filelist
from utility_programs.read_routines.SAMI import auto_read as auto_read_sami


def get_fit(array, lowcut=80, highcut=40):
    """Wrapper for filters.make_fits to accept xarrays

    Parameters
    ----------
    array : xarray dataset/dataarray
        Array to perform the filter on
    lowcut : int , optional, default 80
        Lowcut of the bandpass filter (in minutes)
    highcut : int, optional, default 40
        Highcut of the bandpass (in minutes)

    Returns
    -------
    xarray.dataset/dataarray
        Xarray object with filter applied.

        """

    return xr.apply_ufunc(filters.make_fits, array.load(), lowcut, highcut)


def get_diffs(array, lowcut=80, highcut=40):
    """Returns the % difference of the input array over background 
    calculated with a bandpass filter.

    Parameters
    ----------
    array : xarray dataset/dataarray
        array to calculate on.
    lowcut : int , optional, default 80
        Lowcut of the bandpass filter (in minutes)
    highcut : int, optional, default 40
        Highcut of the bandpass (in minutes)

    Returns
    -------
    xarray dataset/dataarray
        Dataset/array after performing the bandpass fit.
        """

    return 100 * (array - get_fit(array, lowcut, highcut)) / array


def main(
        directories,
        cuts,
        start_plotting,
        end_plotting,
        samicols=['edens'],
        gitmcols=['Rho'],
        prefixes=None,
        outdir=None,
        progress=True,
        use_dask=False,
        map_start_times=None,
        cache_data=False, # Used for debugging and to speed up reads to later.
        ):

    # read in data
    data = {}
    dirs_read = 0

    #check if cache is requested (and exists)
    data_read_from_cache = False
    if cache_data and len(directories) == 1:
        if not directories[0].endswith('/'):
            directory = directories[0] + '/'
        for file in os.listdir(directory):
            if use_dask:
                data[file] = xr.open_mfdataset(directory+file)
            else:
                data[file] = xr.open_dataset(directory+file)
        data_read_from_cache = True
    # print(directories)
        
    if cache_data and not data_read_from_cache:
        save_data_to = directories[-1]
        directories = directories[:-1]

    if not data_read_from_cache:
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
    if prefixes is not None and not data_read_from_cache:
        if not isinstance(prefixes, list):
            prefixes = [prefixes]
        for fpath in directories:
            for prefix in prefixes:

                if 'sami' in prefix.lower():
                    data[str(dirs_read) + '-sami'] = auto_read_sami(
                        fpath + prefix, cols=samicols,
                        start_idx=start_idx, end_idx=end_idx,
                        use_dask=use_dask, progress_bar=progress)
                elif 'gitm' in prefix.lower():
                    data[str(dirs_read) + '-gitm'] = auto_read_gitm(
                        fpath + prefix, cols=gitmcols,
                        start_idx=start_idx, end_idx=end_idx,
                        use_dask=use_dask, progress_bar=progress)
                else:
                    raise ValueError('Prefix %s is not valid' % prefix)
            dirs_read += 1

    elif not data_read_from_cache:
        for fpath in directories:
            data[str(dirs_read) + '-sami'] = auto_read_sami(
                fpath, cols=samicols, progress_bar=progress,
                start_idx=start_idx, end_idx=end_idx, use_dask=use_dask)
            data[str(dirs_read) + '-gitm'] = auto_read_gitm(
                fpath, cols=gitmcols, progress_bar=progress,
                start_idx=start_idx, end_idx=end_idx, use_dask=use_dask)
            dirs_read += 1


    if cache_data and not data_read_from_cache:
        for file in data.keys():
            data[file].to_netcdf(os.path.join(save_data_to, file),
                encoding={'time': {"dtype": "int16"}})

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

        if map_start_times is None:
            map_times = [afewtimes for i in range(len(data.keys()))]
        else:
            map_times = []
            if not data_read_from_cache:
                all_files_started_read = start_idx
                selected_times = all_times[start_idx:end_idx]
            else:
                start_idx = 0
                # Get times from first value of data.keys()
                selected_times = [pd.Timestamp(xr.decode_cf(
                    data[list(data.keys())[0]]).time.values[i]) for i in range(
                    len(data[list(data.keys())[0]].time))]
                selected_times = np.array(selected_times)

            
            for plot_group in range(len(map_start_times)):
                time_here = str_to_ut(map_start_times[plot_group])
                idx_start = np.argmin(np.abs(time_here - selected_times))
                map_times.append(idx_start + np.array(afewtimes))

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
                            plot_vals=map_times[int(model_key[0])],
                            suptitle=model_key,
                            out_fname=out_name + '/%s-%s.png' %
                            (model_key, s_col,))
                if 'gitm' in model_key and a < 700:
                    for g_col in gitmcols:
                        panel_plot(data[model_key][g_col].sel(
                            alt=a, method='nearest'),
                            x='lon', y='lat', wrap_col='time',
                            do_map=True,
                            plot_vals=map_times[int(model_key[0])],
                            suptitle=model_key,
                            out_fname=out_name + '/%s-%s.png' %
                            (model_key, g_col,))

    # Single Point plots:
    ## Number of individual int's in data.keys()
    nrows = np.unique([int(s) for s in str(data.keys()) if s.isdigit()])
    # nrows = len(directories)
    if cuts == ['all'] or 'lonlatalt' in cuts:
        if not os.path.isdir(os.path.join(outdir, 'single-point')):
            os.makedirs(os.path.join(outdir, 'single-point'))

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
                out_name = os.path.join(
                    outdir, 'single-point', '%s-%s-%i-%i-%i' %
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
        default=['all'], nargs='*',
        help="Which cuts to use."
        " Options are 'all', 'alt', 'lon', 'lat', 'lonlatalt'"
        " lon/lat cut will make keos, alt makes maps, lonlatalt makes"
        " plots of a single pt throughout time. Default: 'all'"
        " Change cut values in source code.")

    args.add_argument(
        '--map_start_times',
        type=str,
        nargs='*',
        help='When making maps out of data in several directories, you can'
        ' specify a datetime (in YYYMMDDHHMMSS format [will be padded with 0s'
        ' so 20110521 is valid]) to be included,'
        ' otherwise, plots will be made according to parameters set in here.')

    args.add_argument('--samicols', type=str, nargs='*',
                      default=['edens'],
                      help="Which columns to plot from SAMI data."
                      "Default: ['edens',] (only reads SAMI_REGRID)")
    args.add_argument('--gitmcols', type=str, nargs='*',
                      default=['Rho'],
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

    args.add_argument(
        '--cache_data',action='store_true',
        help='Cache data to make reads faster. Set directories to a single'
        ' string to read from cached data. If multiple dirs are given,'
        ' data will be written to the LAST listed directory!')

    args = args.parse_args()

    main(directories=args.dirs,
         outdir=args.out_dir,
         prefixes=args.prefixes,
         cuts=args.cuts,
         samicols=args.samicols,
         gitmcols=args.gitmcols,
         start_plotting=args.start_plotting,
         end_plotting=args.end_plotting,
         use_dask=args.dask,
         map_start_times=args.map_start_times,
         cache_data=args.cache_data)
