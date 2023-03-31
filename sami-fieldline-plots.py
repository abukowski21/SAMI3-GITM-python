"""SAMI3 fieldline plots.

Created March 13 2023 by Aaron Bukowski.

"""
import datetime
from utility_programs.read_routines import SAMI
from utility_programs.plot_help import UT_from_Storm_onset
from scipy.interpolate import LinearNDInterpolator, interp1d
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
import time
import aacgmv2
from tqdm.auto import tqdm


def make_filter(params=None):
    # Define the cutoff frequencies
    lowcut = 1/(100/60)  # 100 minutes in units of sample^-1
    highcut = 1/(30/60)  # 30 minutes in units of sample^-1

    # Define the Butterworth filter
    nyquist = 0.5 * 5  # 5 minutes is the sampling frequency
    low = lowcut / nyquist
    high = highcut / nyquist
    sos = signal.butter(2, [low, high], btype='bandstop', output='sos')
    return sos


def remove_background(time_series, sos, axis = 0):

    # Apply the filter to the time series
    filtered_data = signal.sosfiltfilt(sos, time_series, axis)

    return filtered_data


def draw_field_line_plot(x, y, z, title, interpolate=False,
                         x_label='Mlat (Deg)', y_label='Altitude (km)',
                         x_lims=[-65, 65], y_lims=[0, 1200],
                         cbar_label=None, cbar_lims=None,
                         fname=None, save_or_show='save',
                         fpeak_col=None):

    if cbar_lims is None:
        cbar_lims = [min(z), max(z)]

    plot_mask = (x >= x_lims[0]) & (x <= x_lims[1]) & \
        (y >= y_lims[0]) & (y <= y_lims[1])

    x = x[plot_mask]
    y = y[plot_mask]
    z = z[plot_mask]
    if fpeak_col is not None:
        fpeak_col = fpeak_col[plot_mask]

    if interpolate:

        grid_x, grid_y = np.meshgrid(
            np.linspace(x_lims[0], x_lims[1], 100),
            np.linspace(y_lims[0], y_lims[1], 150))
        in_x, in_y = np.meshgrid(
            x,y)

        loc_grid = list(zip(x, y))
        interp = LinearNDInterpolator(
            np.array([in_x.flatten(), in_y.flatten()]).T, z,
            rescale=True)

        znew = interp(list(zip(grid_x.flatten(), grid_y.flatten())))

        plt.imshow(znew.reshape(grid_x.shape), origin='lower',
                   vmin=cbar_lims[0], vmax=cbar_lims[1],
                   # interpolation='bilinear', interpolation_stage='rgba',
                   extent=[x_lims[0], x_lims[1], y_lims[0], y_lims[1]],
                   aspect='auto')

    else:
        plt.scatter(x, y, c=z, vmin=cbar_lims[0], vmax=cbar_lims[1])

    plt.ylim(y_lims)
    plt.xlim(x_lims)
    plt.colorbar(label=cbar_label)

    plt.figtext(.5, .9, title, fontsize=12, ha='center')
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    if fpeak_col is not None:
        xs = []
        ys = []
        bins = np.arange(np.min(x), np.max(x), 3)

        for i in range(len(bins)-1):
            mlat_mask = (x >= bins[i]) & (x <= bins[i+1])
            max_here = (np.argmax(fpeak_col[mlat_mask]))
            xs.append(x[mlat_mask][max_here])
            ys.append(y[mlat_mask][max_here])

        cubic_interpolation_model = interp1d(xs, ys, kind="quadratic")
        newx = np.linspace(min(xs), max(xs), 100)
        newy = cubic_interpolation_model(newx)

        plt.plot(newx, newy, c='k')

    if save_or_show == 'show':
        plt.show()
        plt.close()
    elif save_or_show == 'save':
        if fname is None:
            raise ValueError('plot save path must be given!')
        else:
            # fname = fname.replace(' ','')
            try:
                plt.savefig(fname)
            except FileNotFoundError:
                try:
                    directory_list = os.path.join(fname).split('/')[:-1]
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
                raise ValueError('File not found error. Check your path.')
            plt.close()

    else:
        raise ValueError(
            'save_or_show input is invalid. Accepted inputs are "save" or',
            '"show", you gave ', save_or_show)


def main(args):
    data_path = args.sami_data_path

    global dtime_storm_start
    dtime_storm_start = datetime.datetime.strptime(
        args.dtime_storm_start.ljust(14, '0'), '%Y%m%d%H%M%S')

    global dtime_sim_start
    dtime_sim_start = datetime.datetime.strptime(
        args.dtime_sim_start.ljust(14, '0'), '%Y%m%d%H%M%S')

    sami_data, times = SAMI.read_sami_data(sami_data_path=data_path,
                                           dtime_sim_start=dtime_sim_start,
                                           dtime_storm_start=dtime_storm_start,
                                           t_start_idx=args.plot_start_delta,
                                           t_end_idx=args.plot_end_delta,
                                           cols=args.cols,
                                           pbar=True)

    mlons = np.unique(sami_data['grid']['mlon'].round(2))

    mlons_to_plot = mlons[::10]

    if args.cols == 'all':
        cols_to_plot = sami_data['data'].keys()
    else:
        cols_to_plot = args.cols

    if args.plot_type == 'all':
        plot_type = ['raw', 'bandpass', 'diff']
    else:
        plot_type = [args.plot_type]

    pbar = tqdm(total=len(mlons_to_plot)*len(cols_to_plot) *
                len(plot_type)*len(times))

    for mlon in mlons_to_plot:
        mask = (sami_data['grid']['mlon'].round(2) == mlon)
        x = sami_data['grid']['mlat'][mask]
        y = sami_data['grid']['alt'][mask]

        for col in cols_to_plot:

            for type in plot_type:
                if type == 'raw':
                    data_src = sami_data['data'][col][mask, :].copy()
                    cbar_lims = [np.min(data_src), np.max(data_src)]
                elif type == 'bandpass':
                    data_raw = sami_data['data'][col][mask, :].copy()
                    sos = make_filter()
                    data_src = remove_background(data_raw, sos, axis = 1)
                    cbar_lims = [np.min(data_src), np.max(data_src)]
                elif type == 'diff':
                    data_raw = sami_data['data'][col][mask, :].copy()
                    sos = make_filter()
                    data_src = remove_background(data_raw, sos, axis = 1)
                    data_src = 100*(data_raw - data_src)/data_raw
                    cbar_lims = [-3, 3]
                else:
                    raise ValueError('Invalid plot type')

                for nt, dtime in enumerate(times):

                    if args.fpeak:
                        fpeak_col = sami_data['data']['edens'][mask, nt]
                    else:
                        fpeak_col = None

                    mlt = aacgmv2.wrapper.convert_mlt(
                        mlon, dtime, m2a=False).round(2)[0]

                    title = '%s, MLON = %s, MLT = %s' % (
                        col, str(mlon), str(mlt) + '\n UT from storm onset:' +
                        UT_from_Storm_onset(dtime, dtime_storm_start))

                    fname = os.path.join(args.out_path, col, type,
                                         str(int(mlon)),
                                         str(nt).rjust(3, '0'))

                    draw_field_line_plot(
                        x, y, data_src[:, nt], interpolate=args.interpolate,
                        cbar_label=col, title=title, fname=fname,
                        cbar_lims=cbar_lims,
                        save_or_show='save', fpeak_col=fpeak_col)

                    pbar.update()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='SAMI3 fieldline plots.')

    parser.add_argument(
        'dtime_storm_start',
        help='Datetime of storm start. Format YYYYMMDDHHmmss',
        action='store')

    parser.add_argument(
        'dtime_sim_start',
        help='Datetime of simulation start. Format YYYYMMDDHHmmss',
        action='store')

    parser.add_argument(
        'sami_data_path', type=str,
        help='Path to sami data', action='store')

    parser.add_argument(
        'out_path', type=str,
        help='path to where plots will be saved', action='store')

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
        '--cols', nargs="+", type=str,
        help='Which columns to plot. Default: all', default='all')

    parser.add_argument(
        '--interpolate', action='store_true',
        help='whether to make the plots interpolated or not. Default: False')

    parser.add_argument(
        '--fpeak', action='store_true',
        help='whether to plot the fpeak or not. Default: False')

    parser.add_argument(
        '--plot_type', type=str,
        default='diff',
        help='which type of plot? raw, bandpass, diff, or all. Default: diff')

    args = parser.parse_args()

    main(args)
