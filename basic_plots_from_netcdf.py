import argparse
import glob
import os

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from tqdm.auto import tqdm
import pandas as pd

from utility_programs.filters import make_fits
from utility_programs.utils import get_var_names, str_to_ut


def run_processing_options(ds,
                           process_options,):
    """Process Dataset for Plots:

    Args:
        ds (xarray.Dataset):
            Dataset to be processed
        process_options (str or list):
            Options to be applied to the dataset.
            Currently supported options:
                'alt_int': integrate over altitude
                'bandpass': apply bandpass filter
                'transpose': transpose the dataset

    Returns:
        xarray.Dataset: Dataset with processing options applied

    """
    if 'alt_int' in process_options:
        ds = ds.mean(dim='alt')

    if 'bandpass' in process_options:
        ds = make_fits(ds)

    if 'transpose' in process_options:
        ds = ds.transpose()

    return ds


def autoplot(
        file_list,
        columns_to_plot,
        output_dir=None,
        show_map=False,
        time_lims=[0, -1],
        cut_dict={},
        lim_dict={},
        loop_var='time',
        process_options=None,
        plot_arg_dict=None,
        concat_dim='time'):

    # Check validity of cuts & args... :
    if 'alt' in cut_dict:
        if 'alt_int' in process_options:
            raise ValueError('Cannot select altitude when using alt_int.')
    if show_map:
        if 'lon' in cut_dict.keys() or 'lat' in cut_dict.keys():
            raise ValueError('Cannot make maps with lon/lat cuts.'
                             ' Try running again without map flag on.')

    file_list = np.sort(file_list)

    # Only grab data for the requested column(s)
    if isinstance(columns_to_plot, str):
        columns_to_plot = [columns_to_plot]
    drops = []

    if len(file_list) >= 0:  # only need to do this if we're not
        # using single_file outputs...
        ds0 = xr.open_dataset(file_list[0])
        for v in ds0.data_vars:
            if v not in columns_to_plot:
                drops.append(v)
        del ds0  # save memory

    # open & read the files, drop variables we don't want
    print('Reading in {} files...'.format(len(file_list)))
    ds = [xr.open_dataset(f, drop_variables=drops) for f in file_list]
    ds = xr.concat(ds, dim=concat_dim)
    print('Done reading files.')

    # process the plotlims now:
    if lim_dict is not None:
        if 'alt' in lim_dict:
            alt_lim = lim_dict.pop('alt')
            if alt_lim > 90:
                ds = ds.sel(alt=alt_lim, method='nearest')
            else:
                ds = ds.isel(alt=alt_lim)
        # now make sure it's not empty:
        if len(lim_dict) > 0:
            ds = ds.sel(lim_dict, method='nearest')

    if time_lims[0] > 100000 or time_lims[1] > 100000:
        # it's probably a datetime string

        if time_lims[0] != 0:
            dtime_lim_0 = str_to_ut(str(int(time_lims[0])))
            ds = ds.where(ds.time > pd.Timestamp(dtime_lim_0), drop=True)
            time_lims[0] = 0  # for isel-ing later
        if time_lims[1] != -1:
            dtime_lim_1 = str_to_ut(str(int(time_lims[1])))
            ds = ds.where(ds.time < pd.Timestamp(dtime_lim_1), drop=True)
            time_lims[1] = -1  # for isel-ing later

    ds = ds.isel(time=slice(time_lims[0], time_lims[1]))

    # check if output dir exists:
    a = ''
    if len(cut_dict) > 0:
        for i in cut_dict.keys():
            a += str(i) + '-' + str(int(cut_dict[i])) + '_'
        a = a[:-1]
    out_dir = os.path.join(output_dir, a)
    out_dir = out_dir.replace('//', '/')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        print('created directory: ', out_dir)
    print('Saving plots as: ', os.path.join(out_dir, 'var_####.png'))
    # Now plot the data with the cuts specified.

    for var in columns_to_plot:

        # look thru the process options:
        if process_options is not None:
            ds[var] = run_processing_options(ds[var], process_options)
            # With Dask it shouldn't matter if we trim the ds before or after

        for nplot in tqdm(range(ds[loop_var].shape[0]),
                          desc='%s loop for %s plots: '
                          % (loop_var, var)):
            out_fname = os.path.join(out_dir,
                                     var + '_' + str(nplot))

            if show_map:
                p = ds[var].isel({loop_var: nplot}).sel(
                    cut_dict, method='nearest').plot(
                        transform=ccrs.PlateCarree(),
                        subplot_kws={"projection": ccrs.PlateCarree()},
                        x='lon', y='lat',
                        **plot_arg_dict,)
                p.axes.coastlines()
                p.axes.gridlines(alpha=0.6)
                plt.savefig(out_fname)
                plt.close()

            else:
                p = ds[var].isel({loop_var: nplot}).sel(
                    cut_dict, method='nearest').plot(
                        **plot_arg_dict,)
                plt.savefig(out_fname)
                plt.close()

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""
        Make plots from netCDF files.
        --------------------------
        Requires files to be post-processed with "PostProcessModelResults.py"
        Will allow creation of keograms and maps of various variables.

        Lat/lon cuts are performed with method='nearest' in xarray.
            Note: models probably use lon limits from 0-360, not -180-180.
            Errors will not be raised if you select things wrong.
        Alt cut can either be 'nearest' or by index (auto-selected)
            ('nearest' is used when alt_cut is below 90 km)
        """,
        usage='',
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('data_dir', type=str,
                        help='Directory where (postprocessed) model outputs'
                        ' are located.\n')

    parser.add_argument('-col', type=str, nargs='+', required=False,
                        help='List of columns to plot. \n '
                        "Program will automatically grab data from the "
                        "correct model, depending on your input. \n"
                        " If output exists in both models, it will plot both."
                        "Specify --model to select only one model. \n"
                        " NOTE: to plot TEC or something like that, "
                        " use 'SAMI_REGRID' as the model, 'edens' as the"
                        " column, and set --process_options to 'alt_int'\n")

    parser.add_argument('--model', type=str, nargs='*',
                        default=['SAMI_REGRID', 'GITM'],
                        help='Model to plot. Currently only "SAMI" or "GITM".'
                        '\n  Default is to look inside of both model outputs'
                        ' and find the variable requested. \n'
                        'NOTE: When plotting SAMI outputs, this program only '
                        'handles SAMI_REGRID. Use #TODO to plot magnetic'
                        'grid\n')

    parser.add_argument('-out_dir', type=str, default='./',
                        help='Directory where plots will be saved.'
                        '\n  Default is current directory.\n')

    parser.add_argument('-m', '--show_map', action='store_true',
                        help='Plot data on a map.\n')

    parser.add_argument('-t', '--time_lims', type=float, nargs='*',
                        default=[0, -1],
                        help="Time limits to make plots for.\n"
                        "Usage can either be index of time array or datetime"
                        " (as str: YYYYMMDDHHMMSS)\n"
                        "Default is to plot all times [0,-1].\n "
                        "Set to 0 or -1 to plot the first (or last) time.\n "
                        "You can set to any two values to plot a range of"
                        " times (e.g. '-t 25 30' will make plots "
                        "from timestep 25 to 30)\n "
                        "Can also use datetime strings (YYYYMMDDHHMMSS).\n"
                        )

    parser.add_argument('-lon', '--lon_cut', type=float, nargs='*',
                        help="Longitude cut to make (for keograms).\n"
                        "If two values are given, they are taken to be\n "
                        "the limits of the plots made (maps, etc.).\n "
                        "example: -lon 50 will plot along longitude=50.\n "
                        "         -lon 50 100 will set plot limits between 50"
                        " and 100.\n "
                        )

    parser.add_argument('-lat', '--lat_cut', type=float, nargs='*',
                        help="Latitude cut to make (for keograms).\n"
                        "If two values are given, they are taken to be\n "
                        "the limits of the plots made (maps, etc.).\n "
                        "example: -lat 50 will plot along latitude=50.\n "
                        "         -lat 50 100 will set plot limits between 50"
                        " and 100.\n "
                        )

    parser.add_argument('-a', '--alt_cut', type=float, nargs='*',
                        help="Altitude cut to make (for keograms).\n"
                        "If two values are given, they are taken to be\n "
                        "the limits of the plots made (maps, etc.).\n "
                        "example: -alt 350 will plot along altitude=350km.\n "
                        "         -alt 150 700 will set plot limits between"
                        "         150 and 750.\n "
                        )

    parser.add_argument('--loop_var', type=str, nargs='*', default='time',
                        help='Dimension to loop over when plotting.\n'
                        ' Example: --loop_var time (default) will make a'
                        ' single plot for each time step.\n'
                        '          --loop_var alt will make a single plot'
                        ' for each altitude.\n'
                        'Not required (will just make one plot).\n')

    parser.add_argument('--col_help', action='store_true',
                        help='Prints a list of available columns to plot'
                        ' from all available data files in data_dir.\n')

    parser.add_argument('--process_option', type=str, nargs='*',
                        default=None,
                        help='Apply any post-processing options to the data.\n'
                        'Currently, only "alt_int", "transpose", or "bandpass"'
                        ' are supported. \nMore can be added (easily) by '
                        'modifying run_processing_options() in this script.\n')

    parser.add_argument('--plot_args', type=str, nargs='*', default=None,
                        help='Extra arguments to pass the plotting function.'
                        '\n Example: --plot_args cmap=jet vmin=-3 vmax=3')

    args = parser.parse_args()

    # check if user needs help:
    if args.col_help:
        get_var_names(args.data_dir, args.model)

    # format plot_arguments. Cannot pass NoneType to plotting functions.
    if args.plot_args is not None:
        args.plot_args = dict(x.split('=') for x in args.plot_args)
        for k in args.plot_args.keys():
            # change str to int or float if possible
            try:
                args.plot_args[k] = int(args.plot_args[k])
            except ValueError:
                try:
                    args.plot_args[k] = float(args.plot_args[k])
                except ValueError:
                    pass

    else:
        args.plot_args = {}

    # set up cut & plot lim dicts.
    plot_lims = {}
    plot_cuts = {}
    if args.lon_cut is not None:
        if len(args.lon_cut) == 1:
            plot_cuts['lon'] = args.lon_cut[0]
        elif len(args.lon_cut) == 2:
            plot_lims['lon'] = args.lon_cut
        else:
            raise ValueError('lon_cut must be either 1 or 2 values.'
                             ' If you mean to run multiple lons, set '
                             'loop_var=lon')

    if args.lat_cut is not None:
        if len(args.lat_cut) == 1:
            plot_cuts['lat'] = args.lat_cut[0]
        elif len(args.lat_cut) == 2:
            plot_lims['lat'] = args.lat_cut
        else:
            raise ValueError('lat_cut must be either 1 or 2 values.'
                             ' If you mean to run multiple lats, set '
                             'loop_var=lat')

    if args.alt_cut is not None:
        if len(args.alt_cut) == 1:
            plot_cuts['alt'] = args.alt_cut[0]
        elif len(args.alt_cut) == 2:
            plot_lims['alt'] = args.alt_cut
        else:
            raise ValueError('alt_cut must be either 1 or 2 values.'
                             ' If you mean to run multiple alts, set '
                             'loop_var=alt')

    if not isinstance(args.loop_var, str):
        if len(args.loop_var) > 1:
            raise ValueError('Can only loop over one variable at a time.')
        else:
            loop_var = args.loop_var[0]
    else:
        loop_var = args.loop_var

    """
    This loop will look for the column requested in the data files,
    and then make the requested plots.

    - Files are listed out twice (here & in autoplot), this is so that
        autoplot can be run interactively and leaves the time stuff
        transparent to the user.
    """

    if not isinstance(args.col, list):
        cols = [args.col]
    else:
        cols = args.col

    if not isinstance(args.model, list):
        models = [args.model]
    else:
        models = args.model

    made_plots = {}
    for col in cols:
        if ('AltInt' not in col) and (len(plot_cuts) == 0)\
                and ('alt_int' not in args.process_option):
            raise ValueError('Must specify at least one cut (lon, lat, alt)'
                             ' when not using AltInt variable.')
        made_plots[col] = 0
        for model in models:
            files = glob.glob(
                os.path.join(args.data_dir, '*' + model + '.nc'))
            ds0 = xr.open_dataset(files[0])
            if col in ds0:
                autoplot(files,
                         columns_to_plot=col,
                         output_dir=args.out_dir,
                         show_map=args.show_map,
                         time_lims=args.time_lims,
                         cut_dict=plot_cuts,
                         lim_dict=plot_lims,
                         loop_var=loop_var,
                         process_options=args.process_option,
                         plot_arg_dict=args.plot_args)

                made_plots[col] += 1

    for col in made_plots.keys():
        if made_plots[col] == 0:
            print('Column {} not found in any data files.'.format(col))
            print('You may want to run --col_help to see available columns.')
        else:
            print('Made plots from {} models for column {}.'.format(
                made_plots[col], col))
