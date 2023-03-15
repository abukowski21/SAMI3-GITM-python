"""This will make GITM polar dial plots and a map of the data.

User selects which GITM file(s) to plot, and the program will make a
nice and pretty graph.  The user can also select which variable to
put in which plot.

created mar 15 2023 by aaron
# """
# # %%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from utility_programs.read_routines import GITM

from utility_programs import filters

import argparse
from tqdm.auto import tqdm
from utility_programs.plot_help import UT_from_Storm_onset
import time
import os

import geopandas
world = geopandas.read_file(geopandas.datasets.get_path("naturalearth_lowres"))


# # %%
# dials_var = 'AltIntJouleHeating(W/m2)'
# map_var = 'VerticalTEC'

# # %%


def main(args):

    # let's make sure the user doesn't need help with vars:
    if args.var_help:
        from aetherpy.io.read_routines import read_gitm_file
        import glob
        var3d = read_gitm_file(glob.glob(
            os.path.join(args.gitm_data_path, '3DALL*.bin'))[0])['vars']
        var1d = read_gitm_file(glob.glob(
            os.path.join(args.gitm_data_path, '2DANC*.bin'))[0])['vars']
        found = False
        if args.polar_var is not None:
            if args.polar_var in var3d:
                print(args.polar_var, ' was found as a 3D variable')
                found = True
            if args.polar_var in var1d:
                print(args.polar_var, ' was found as a 1D variable')
                found = True
        if args.map_var is not None:
            if args.map_var in var3d:
                print(args.map_var, ' was found as a 3D variable')
                found = True
            if args.map_var in var1d:
                print(args.map_var, ' was found as a 1D variable')
                found = True
        if not found:
            print('Variable not found in any GITM files')
        print('avalibale variables are: \n 3D: \n', var3d,
              '\n 1D: \n', var1d)
        return

    # Read in files, if cut alt, use 3DALL. make sure all vars r there
    # If not, read in 2DANC.
    need_2d = True
    need_3d = False
    if args.alt_cut is not None:
        need_2d = False
        times3, gitm_grid3, gitm_f3, gitm_vars3 = GITM.read_gitm_into_nparrays(
            args.gitm_data_path, args.dtime_start,
            gitm_file_pattern='3DALL*.bin',
            t_start_idx=args.plot_start_delta,
            t_end_idx=args.plot_end_delta, return_vars=True)

        if args.polar_var not in gitm_vars3:
            need_2d = True
        else:
            need_3d = True
            polar_in_3d = True
        if args.map_var not in gitm_vars3:
            need_2d = True
        else:
            need_3d = True
            map_in_3d = True
        if need_3d:
            gitm_fits3 = filters.make_fits(gitm_f3)
            percent3 = 100*(gitm_f3 - gitm_fits3)/gitm_f3

    if need_2d:
        times2, gitm_grid2, gitm_f2, gitm_vars2 = GITM.read_gitm_into_nparrays(
            args.gitm_data_path, args.dtime_start,
            gitm_file_pattern='2DANC*.bin',
            t_start_idx=args.plot_start_delta,
            t_end_idx=args.plot_end_delta, return_vars=True)

        if args.polar_var not in gitm_vars2:
            print('Polar variable %s not found in 3D or 2D files'
                  % args.polar_var, gitm_vars3, gitm_vars2)
            return
        if args.map_var not in gitm_vars2:
            print('Map variable %s not found in 3D or 2D files'
                  % args.map_var, gitm_vars3, gitm_vars2)
            return
        gitm_fits2 = filters.make_fits(gitm_f2)
        percent2 = 100*(gitm_f2 - gitm_fits2)/gitm_f2

    # alt cuts?
    if args.alt_cut is not None and need_3d is False:
        raise ValueError('Alt cut requested, but no 3D files needed')
    elif args.alt_cut is not None and need_3d is True:
        if args.alt_cut < np.min(gitm_grid3['altitude']):
            alt_cuts = args.alt_cut
        else:
            alt_cuts = []
            for a in args.alt_cuts:
                alt_cuts.append(np.abs(gitm_grid3['altitude']/1e3 - a))
    else:
        alt_cuts = [None]

    # Set up plot_types:
    if args.figtype_map == 'all':
        figtype_map = ['raw', 'fit', 'diff']
    else:
        figtype_map = [args.figtype_map]
    if args.figtype_polar == 'all':
        figtype_polar = ['raw', 'fit', 'diff']
    else:
        figtype_polar = [args.figtype_polar]

    # make sure lats and lons match between 2d and 3d files
    if need_2d:
        lats2 = np.unique(gitm_grid2['latitude'])
        lons2 = np.unique(gitm_grid2['longitude'])
    if need_3d:
        lats3 = np.unique(gitm_grid3['latitude'])
        lons3 = np.unique(gitm_grid3['longitude'])
    if need_2d and need_3d:
        if lats2 != lats3 or lons2 != lons3:
            raise ValueError(
                'Lats and Lons do not match between 2D and 3D files')
        elif times2 != times3:
            raise ValueError(
                'Times do not match between 2D and 3D files')
        else:
            lats = lats2, lons = lons2, times = times2

    maskNorth = lats > 45
    maskSouth = lats < -45

    pbar = tqdm(total=len(times)*len(alt_cuts)*len(figtype_polar)
                * len(figtype_map), desc='making plots')

    for a in alt_cuts:
        for nt, dtime in times:
            for p_fig in figtype_polar:
                for m_fig in figtype_map:
                    # get data.
                    if polar_in_3d:
                        if p_fig == 'raw':
                            poldata = gitm_f3[
                                nt, gitm_vars3.index(args.polar_var), :, :, a
                            ].copy()
                            pol_str = ''
                            minP = np.min(gitm_f3[
                                :, gitm_vars3.index(args.polar_var), :, :, a])
                            maxP = np.max(gitm_f3[
                                :, gitm_vars3.index(args.polar_var), :, :, a])
                        elif p_fig == 'fit':
                            poldata = gitm_fits3[
                                nt, gitm_vars3.index(args.polar_var), :, :, a
                            ].copy()
                            pol_str = ' (fit)'
                            minP = np.min(gitm_fits3[
                                :, gitm_vars3.index(args.polar_var), :, :, a])
                            maxP = np.max(gitm_fits3[
                                :, gitm_vars3.index(args.polar_var), :, :, a])
                        elif p_fig == 'diff':
                            poldata = percent3[
                                nt, gitm_vars3.index(args.polar_var), :, :, a
                            ].copy()
                            pol_str = ' (% over BG)'
                            minP = np.min(percent3[
                                :, gitm_vars3.index(args.polar_var), :, :, a])
                            maxP = np.max(percent3[
                                :, gitm_vars3.index(args.polar_var), :, :, a])

                    else:
                        if p_fig == 'raw':
                            poldata = gitm_f2[
                                nt, gitm_vars2.index(args.polar_var)].copy()
                            pol_str = ''
                            minP = np.min(gitm_f2[
                                :, gitm_vars2.index(args.polar_var), :, :, a])
                            maxP = np.max(gitm_f2[
                                :, gitm_vars2.index(args.polar_var), :, :, a])
                        elif p_fig == 'fit':
                            poldata = gitm_fits2[
                                nt, gitm_vars2.index(args.polar_var)].copy()
                            pol_str = ' (fit)'
                            minP = np.min(gitm_fits2[
                                :, gitm_vars2.index(args.polar_var), :, :, a])
                            maxP = np.max(gitm_fits2[
                                :, gitm_vars2.index(args.polar_var), :, :, a])
                        elif p_fig == 'diff':
                            poldata = percent2[
                                nt, gitm_vars2.index(args.polar_var)].copy()
                            pol_str = ' (% over BG)'
                            minP = np.min(percent2[
                                :, gitm_vars2.index(args.polar_var), :, :, a])
                            maxP = np.max(percent2[
                                :, gitm_vars2.index(args.polar_var), :, :, a])

                    if map_in_3d:
                        if m_fig == 'raw':
                            mapdata = gitm_f3[
                                nt, gitm_vars3.index(args.map_var), :, :, a
                            ].copy()
                            map_str = ''
                            minM = np.min(gitm_f3[
                                :, gitm_vars3.index(args.map_var), :, :, a])
                            maxM = np.max(gitm_f3[
                                :, gitm_vars3.index(args.map_var), :, :, a])
                        elif m_fig == 'fit':
                            mapdata = gitm_fits3[
                                nt, gitm_vars3.index(args.map_var), :, :, a
                            ].copy()
                            map_str = ' (fit)'
                            minM = np.min(gitm_fits3[
                                :, gitm_vars3.index(args.map_var), :, :, a])
                            maxM = np.max(gitm_fits3[
                                :, gitm_vars3.index(args.map_var), :, :, a])
                        elif m_fig == 'diff':
                            mapdata = percent3[
                                nt, gitm_vars3.index(args.map_var), :, :, a
                            ].copy()
                            map_str = ' (% over BG)'
                            minM = np.min(percent3[
                                :, gitm_vars3.index(args.map_var), :, :, a])
                            maxM = np.max(percent3[
                                :, gitm_vars3.index(args.map_var), :, :, a])

                    else:
                        if m_fig == 'raw':
                            mapdata = gitm_f2[
                                nt, gitm_vars2.index(args.map_var)].copy()
                            map_str = ''
                            minM = np.min(gitm_f2[
                                :, gitm_vars2.index(args.map_var), :, :, a])
                            maxM = np.max(gitm_f2[
                                :, gitm_vars2.index(args.map_var), :, :, a])
                        elif m_fig == 'fit':
                            mapdata = gitm_fits2[
                                nt, gitm_vars2.index(args.map_var)].copy()
                            map_str = ' (fit)'
                            minM = np.min(gitm_fits2[
                                :, gitm_vars2.index(args.map_var), :, :, a])
                            maxM = np.max(gitm_fits2[
                                :, gitm_vars2.index(args.map_var), :, :, a])
                        elif m_fig == 'diff':
                            mapdata = percent2[
                                nt, gitm_vars2.index(args.map_var)].copy()
                            map_str = ' (% over BG)'
                            minM = np.min(percent2[
                                :, gitm_vars2.index(args.map_var), :, :, a])
                            maxM = np.max(percent2[
                                :, gitm_vars2.index(args.map_var), :, :, a])

                    # make fig.
                    fig = plt.figure(figsize=(10, 8.5), layout='tight')

                    fig.suptitle('%s UT from storm onset \n (%s) UT'
                                 % (UT_from_Storm_onset(
                                     dtime, args.dtime_storm_start),
                                    str(dtime)), fontsize=17)

                    gs1 = GridSpec(nrows=2, ncols=2, wspace=.1, hspace=.1)
                    ax0 = fig.add_subplot(gs1[0, 0], projection='polar')
                    ax1 = fig.add_subplot(gs1[0, 1], projection='polar')
                    ax2 = fig.add_subplot(gs1[1, :2])

                    # add in plots. polar left, polar right, map
                    r, theta = np.meshgrid(90-lats[maskNorth], lons)
                    ax0.pcolor(np.deg2rad(theta), r, poldata[:, maskNorth],
                               vmin=minP, vmax=maxP)
                    ylabels = ['80', '70', '60', '50']
                    ax0.set_yticklabels(ylabels)
                    ax0.set_xticks(np.arange(0, 2*np.pi, np.pi/2))
                    ax0.set_yticks(np.arange(10, 50, 10))
                    ax0.set_title('North ' + args.polar_var + pol_str)

                    r, theta = np.meshgrid(lats[maskSouth], lons)
                    cb = ax1.pcolor(np.deg2rad(theta), r,
                                    poldata[:, maskSouth],
                                    vmin=minP, vmax=maxP)
                    ylabels = ['-80', '-70', '-60', '-50']
                    ax1.set_yticklabels(ylabels)
                    ax1.set_xticks(np.arange(0, 2*np.pi, np.pi/2))
                    ax1.set_yticks(np.arange(10, 50, 10))
                    ax1.set_title('South ' + args.polar_var + pol_str)
                    fig.colorbar(cb)

                    mapping(data_arr=mapdata, lats=lats, lons=lons,
                            map_var=args.mapvar + map_str,
                            ax=ax2, zorder=1, origin='lower', alpha=0.8,
                            vmin=minM, vmax=maxM,
                            title=m_fig + ' ' + args.map_var)

                    if args.save_or_show == "show":
                        plt.show()
                        plt.close()
                    elif args.save_or_show == "save":
                        fname = os.path.join(args.out_path, str(a),
                                             'map_' + args.map_var + m_fig,
                                             'polar_' + args.polar_var + p_fig,
                                             str(nt).rjust(3, '0') + '.png')
                        fname = fname.replace(" ", "")
                        try:
                            plt.savefig(fname)
                        except FileNotFoundError:
                            try:
                                directory_list = os.path.join(
                                    fname).split("/")[:-1]
                                os.makedirs(os.path.join(*directory_list))
                                plt.savefig(fname)
                            except FileExistsError:
                                # sometimes when we make too many plots in
                                # the same directory, it fails.
                                # this fixes that.
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
                            'save_or_show input is invalid. Accepted inputs',
                            'are "save" or "show", you gave ',
                            args.save_or_show)

                    pbar.update()


def mapping(data_arr, lats, lons, map_var, title=None,
            ax=None, vmin=None, vmax=None,
            **kwargs):
    """This will make a map of the data.

    Parameters
    ----------
    data : np.array
        The data to be plotted.
    lats : np.array
        The latitudes of the data.
    lons : np.array
        The longitudes of the data.
    map_var : str
        The name of the variable being plotted.
    ax : matplotlib.axes.Axes
        The axes to plot on.
    cbar_label : str, optional
        The label for the colorbar, by default None

    Returns
    -------
    matplotlib.axes.Axes
    """

    ax = ax or plt.gca()
    world.plot(ax=ax, color="white", edgecolor="black", zorder=1)
    data = ax.imshow(
        data_arr.T,
        cmap="viridis",
        aspect="auto",
        extent=[min(lons), max(lons), min(lats), max(lats)],
        **kwargs)
    plt.colorbar(data, ax=ax, label=map_var)

    return ax


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Make plots of gitm vars with polar dials!")

    parser.add_argument(
        '-dtime_storm_start',
        help='Datetime of storm start. Format YYYYMMDDHHmmss',
        action='store')

    parser.add_argument(
        '-polar_var', type=str,
        help='Which variable to plot on the polar dials')

    parser.add_argument(
        '-map_var', type=str,
        help='which variable to plot on the map')

    parser.add_argument(
        '--alt_cut', type=int, default=None, nargs='+',
        help='Whether or not to make an alt cut. Default: None' +
        '(no cut) You can specify an alt-idx or an alt in km')

    parser.add_argument(
        '-gitm_data_path', type=str,
        help='Path to gitm data', default='./gitm_dir', action='store')

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
        '--figtype_polar', type=str, action='store', default='raw',
        help='Which type of plot to make (in polar dials).' +
        'Options: raw, filt, diffs. Default: raw')

    parser.add_argument(
        '--figtype_map', type=str, action='store', default='all',
        help='Which type of plot to make (in map view).' +
        'Options: all, raw, filt, diffs. Default: all')

    parser.add_argument(
        '-o', '--overwrite', action='store_true',
        help='overwrite existing files?')

    parser.add_argument(
        '--var_help', action='store_true',
        help='print out available variables and exit' +
        'OR you can specify some vars for plots and check if they work')

    args = parser.parse_args()

    main(args)
