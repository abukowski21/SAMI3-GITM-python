"""This will make GITM polar dial plots and a map of the data.

User selects which GITM file(s) to plot, and the program will make a
nice and pretty graph.  The user can also select which variable to
put in which plot.

created mar 15 2023 by aaron
 """

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from utility_programs.read_routines import GITM, SAMI

from cartopy import crs as ccrs
from cartopy.feature.nightshade import Nightshade
from utility_programs import filters

import argparse
from tqdm.auto import tqdm
from utility_programs.plot_help import UT_from_Storm_onset
import time
import os
import datetime


def var_help(args):
    var3d = GITM.read_bin_to_nparrays(args.gitm_data_path, '3D*.bin',
                                 end_idx=1, return_vars=True)['gitmvars']
    var2d = GITM.read_bin_to_nparrays(args.gitm_data_path, '2D*.bin', 
                                 end_idx=1, return_vars=True)['gitmvars']

    if args.polar_var is not None:
        found = False
        if args.polar_var in var3d:
            print(args.polar_var, ' was found as a 3D variable')
            found = True
        if args.polar_var in var2d:
            print(args.polar_var, ' was found as a 2D variable')
            found = True
        if not found:
            print('Polar variable: %s not found in any GITM files'
                 %args.polar_var)
    if args.map_var is not None:
        found = False
        if args.map_var in var3d:
            print(args.map_var, ' was found as a 3D variable')
            found = True
        if args.map_var in var2d:
            print(args.map_var, ' was found as a 2D variable')
            found = True
        if not found:
            print('Map variable: %s not found in any GITM files'
                 %args.map_var)

    print('available variables are: \n 3D: \n', var3d,
          '\n 2D: \n', var2d)
    return


        
def main(args):

    # let's make sure the user doesn't need help with vars:
    if args.var_help:
        var_help(args)
        return


    # Read in files, if cut alt, use 3DALL. make sure all vars are there
    # If not, read in 2DANC.
    found_map = False
    found_polar = False
    
    times3, gitm_grid3, gitm_f3, gitm_vars3 = GITM.read_bin_to_nparrays(
            gitm_dir=args.gitm_data_path,
            start_idx=args.time_idx,
            end_idx=args.time_idx+1,
            return_vars=True)
    if args.polar_var in gitm_vars3:
        found_polar = True
    if args.map_var in gitm_vars3:
        found_map = True
    
    if not found_map or not found_polar:
        times2, gitm_grid2, gitm_f2, gitm_vars2 = GITM.read_bin_to_nparrays(
            gitm_dir=args.gitm_data_path,
            start_idx=args.time_idx,
            end_idx=args.time_idx+1,
            return_vars=True)

        if args.polar_var in gitm_vars2:
            found_polar = True
        if args.map_var in gitm_vars2:
            found_map = True
            
    if not found_map or not found_polar:
        if args.sami_path is None:
            print(found_map, found_polar, args.sami_path)
            var_help(args)
            return

    # alt cuts?
    if args.alt_cut is not None and need_3d is False and \
            args.sami_path is not None:
        raise ValueError('Alt cut requested, but no 3D files needed')
    elif args.alt_cut is not None and need_3d is True:
        if args.alt_cut < np.min(gitm_grid3['altitude']):
            alt_cut = args.alt_cut
        else:
            alt_cut = np.abs(gitm_grid3['altitude'] /
                             1e3 - args.alt_cut).argmin()
    else:
        alt_cut = None

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
    if need_2d and not need_3d:
        lats = np.unique(gitm_grid2['latitude'])
        lons = np.unique(gitm_grid2['longitude'])
        times = times2
        new_shape = [len(times), len(lons), len(lats)]
    elif need_2d:
        lats2 = np.unique(gitm_grid2['latitude'])
        lons2 = np.unique(gitm_grid2['longitude'])
        new_shape = [len(times), len(lons), len(lats)]

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

    # pull out the data. makes everything much easier.
    if polar_in_3d:
        polar_data = gitm_f3[:, gitm_vars2.index(
            args.polar_var), :, :, alt_cut]
    else:
        polar_data = gitm_f2[:, gitm_vars2.index(
            args.polar_var), :, :].copy().reshape(new_shape)
    polar_fits = filters.make_fits(polar_data)
    polar_percents = 100*(polar_data - polar_fits)/polar_data

    if args.sami_path is None:
        if map_in_3d:
            map_data = gitm_f3[:, gitm_vars2.index(
                args.map_var), :, :, alt_cut].copy()

        else:
            map_data = gitm_f2[:, gitm_vars2.index(
                args.map_var), :, :].copy().reshape(new_shape)
        map_fits = filters.make_fits(map_data)
        map_percents = 100*(map_data - map_fits)/map_data

    else:
        if args.dtime_sim_start != '':
            dtime_sim_start = datetime.datetime.strptime(
                args.dtime_sim_start.ljust(14, '0'),
                '%Y%m%d%H%M%S')
        else:
            raise ValueError("""Simulation start time must
                             be given to read SAMI data""")
        f, tectimes = SAMI.read_sami_dene_tec(
            args.sami_path,dtime_sim_start=dtime_sim_start,
            reshape=True)
        tectimes = list(tectimes)

        if args.plot_start_delta != 0:
            start_idx = tectimes.index(
                dtime_storm_start - datetime.timedelta(
                    hours=args.plot_start_delta))
            f['data'][args.map_var] = f['data'][args.map_var][start_idx:]
            tectimes = tectimes[start_idx:]

        if args.plot_end_delta != -1:
            end_idx = tectimes.index(
                dtime_storm_start + datetime.timedelta(
                    hours=args.plot_end_delta))
            f['data'][args.map_var] = f['data'][args.map_var][:end_idx]
            tectimes = tectimes[:end_idx]

        if alt_cut is not None:
            sami_alt_cut = np.argmin(np.abs(f['data']['alt'] - args.alt_cut))
            map_data = f['data'][args.map_var][:, :, sami_alt_cut, :]
            map_fits = filters.make_fits(map_data)
            map_percents = 100*(map_data - map_fits)/map_data
        else:
            map_data = f['data'][args.map_var]
            map_fits = filters.make_fits(map_data)
            map_percents = (map_data - map_fits)

    # Set up masks
    maskNorth = ((lats > 45))
    maskSouth = ((lats < -45))

    # Set colorbar limits:
    minP = dict()
    maxP = dict()
    minM = dict()
    maxM = dict()
    poldata = {}
    mapdata = {}
    for p_fig in figtype_polar:
        if p_fig == 'raw':
            minP[p_fig] = np.min(polar_data)
            maxP[p_fig] = np.max(polar_data)
            poldata[p_fig] = polar_data
        elif p_fig == 'fit':
            minP[p_fig] = np.min(polar_fits)
            maxP[p_fig] = np.max(polar_fits)
            poldata[p_fig] = polar_fits
        elif p_fig == 'diff':
            minP[p_fig] = np.min(polar_percents)
            maxP[p_fig] = np.max(polar_percents)
            poldata[p_fig] = polar_percents
        else:
            raise ValueError('Unknown polar figure type %s' % p_fig)

        for m_fig in figtype_map:
            if m_fig == 'raw':
                minM[m_fig] = np.min(map_data)
                maxM[m_fig] = np.max(map_data)
                mapdata[m_fig] = map_data
            elif m_fig == 'fit':
                minM[m_fig] = np.min(map_fits)
                maxM[m_fig] = np.max(map_fits)
                mapdata[m_fig] = map_fits
            elif m_fig == 'diff':
                if args.diff_lim is not None:
                    minM[m_fig] = -args.diff_lim
                    maxM[m_fig] = args.diff_lim
                else:
                    minM[m_fig] = np.min(map_percents)
                    maxM[m_fig] = np.max(map_percents)
                mapdata[m_fig] = map_percents
            else:
                raise ValueError('Unknown map figure type %s' % m_fig)

    pbar = tqdm(total=len(times)*len(figtype_polar)
                * len(figtype_map), desc='making plots')

    # make plots
    for nt, dtime in enumerate(times):
        for p_fig in figtype_polar:
            for m_fig in figtype_map:

                # make fig.
                fig = plt.figure(figsize=(10, 8.5))

                fig.suptitle('%s UT from storm onset \n (%s) UT'
                             % (UT_from_Storm_onset(
                                 dtime, dtime_storm_start),
                                str(dtime)), fontsize=15)

                gs1 = GridSpec(nrows=2, ncols=2, wspace=.1, hspace=.1)
                ax0 = fig.add_subplot(gs1[0, 0], projection='polar')
                ax1 = fig.add_subplot(gs1[0, 1], projection='polar')
                ax2 = fig.add_subplot(gs1[1, :2])

                # add in plots. polar left, polar right, map
                r, theta = np.meshgrid(90-lats[maskNorth], lons)
                ax0.pcolor(np.deg2rad(theta), r,
                           poldata[p_fig][nt, :, maskNorth].T.copy(),
                           vmin=np.min(poldata[p_fig]),
                           vmax=np.max(poldata[p_fig]))
                ylabels = ['80', '70', '60', '50']
                ax0.set_xticks(np.arange(0, 2*np.pi, np.pi/2))
                ax0.set_yticks(np.arange(10, 50, 10))
                ax0.set_yticklabels(ylabels)
                ax0.set_title('North')

                r, theta = np.meshgrid(lats[maskSouth], lons)
                cb = ax1.pcolor(np.deg2rad(theta), r,
                                poldata[p_fig][nt, :, maskSouth].T.copy(),
                                vmin=minP[p_fig], vmax=maxP[p_fig])
                ylabels = ['-80', '-70', '-60', '-50']
                ax1.set_yticklabels(ylabels)
                ax1.set_title('South')
                fig.colorbar(cb, ax=ax1,
                             label=p_fig + ' ' + args.polar_var)

                mapping(data_arr=mapdata[m_fig][nt, :, :], lats=lats,
                        lons=lons,
                        map_var=m_fig + ' ' + args.map_var,
                        ax=ax2, zorder=1, origin='lower', alpha=0.8,
                        vmin=minM[m_fig], vmax=maxM[m_fig],
                        title=m_fig + ' ' + args.map_var)

                if args.save_or_show == "show":
                    plt.show()
                    plt.close()
                elif args.save_or_show == "save":
                    pvar = args.polar_var.replace(
                        '/', '-').replace('(', '-').replace(')', '-')
                    mvar = args.map_var.replace(
                        '/', '-').replace('(', '-').replace(')', '-')

                    fname = os.path.join(
                        args.out_path, str(alt_cut),
                        'map_' + mvar + m_fig + str(args.diff_lim),
                        'polar_' + pvar + p_fig,
                        str(nt).rjust(3, '0') + '.png')

                    fname = fname.replace(" ", "")
                    try:
                        plt.savefig(fname)
                    except FileNotFoundError:
                        try:
                            last_slash = fname.rfind('/')
                            os.makedirs(fname[:last_slash])
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
    print(fname)


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

    ax = ax or plt.gca(projection=ccrs.PlateCarree())
    
    data = ax.imshow(
        data_arr.T,
        cmap="viridis",
        vmin=vmin,
        vmax=vmax,
        extent=[min(lons), max(lons), min(lats), max(lats)],
        transform=ccrs.PlateCarree(),
        **kwargs)
    
    ax.coastlines(zorder=3, color='black', alpha=1)
    ax.gridlines(color='black', linestyle='--',
                                  alpha=0.6, )
    plt.colorbar(data, ax=ax, label=map_var)

    if title:
        ax.set_title(title)
    return ax


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Make plots of gitm vars with polar dials!")

    parser.add_argument(
        '-time_idx', type=int,
        help='Index of timestamp to plot',
        action='store')

    parser.add_argument(
        '-polar_var', type=str,
        help='Which variable to plot on the polar dials')

    parser.add_argument(
        '-map_var', type=str,
        help='which variable to plot on the map')
    
    parser.add_argument(
        '-dtime_sim_start',type=str,default='',
        help='Start datetime (as YYMMDDhhmmss). Reqd if using SAMI')

    parser.add_argument(
        '--alt_cut', type=int, default=None,
        help='Whether or not to make an alt cut. Default: None' +
        '(no cut) You can specify an alt-idx or an alt in km')

    parser.add_argument(
        '-gitm_data_path', type=str,
        help='Path to gitm data', default='./gitm_dir', action='store')

    parser.add_argument(
        '-sami_path', type=str,
        help='Path to SAMI data (if we want it)', default=None, action='store')

    parser.add_argument(
        '--out_path', type=str,
        help='path to where plots are saved', default='./', action='store')


    parser.add_argument(
        '--save_or_show', type=str,
        action='store', default='save', required=False,
        help='Save or show plots. Default: save'
        ' (To run woth "show", you need to run locally, or enable X11 forwarding.)')

    parser.add_argument(
        '--figtype_polar', type=str, action='store', default='raw',
        help='Which type of plot to make (in polar dials).' +
        'Options: raw, filt, diff. Default: raw')

    parser.add_argument(
        '--figtype_map', type=str, action='store', default='all',
        help='Which type of plot to make (in map view).' +
        'Options: all, raw, filt, diff. Default: all')

    parser.add_argument(
        '-o', '--overwrite', action='store_true',
        help='overwrite existing files?')

    parser.add_argument(
        '--var_help', action='store_true',
        help='print out available variables and exit' +
        'OR you can specify some vars for plots and check if they work')

    parser.add_argument(
        '--diff_lim', type=float, default=None,
        help='Set the limits of the colobar on diff maps.')

    args = parser.parse_args()

    main(args)
