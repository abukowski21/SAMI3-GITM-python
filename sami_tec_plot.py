import argparse
import datetime
# import gc
from utility_programs.read_routines import SAMI
import os
import time
# from multiprocessing import Pool

import geopandas
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import signal
from tqdm.auto import tqdm

import scipy.integrate as integ

from utility_programs.plot_help import UT_from_Storm_onset
from utility_programs import plotting_routines
from utility_programs import filters


def main(args):

    if args.map:
        global world
        world = geopandas.read_file(
            geopandas.datasets.get_path('naturalearth_lowres'))

    # Lon to keo:
    global sami_keo_lons
    sami_keo_lons = args.keo_lons

    global lat_lim
    lat_lim = args.lat_lim

    global out_path
    out_path = args.out_path

    global OVERWRITE
    OVERWRITE = args.overwrite

    global dtime_storm_start
    dtime_storm_start = datetime.datetime.strptime(
        args.dtime_storm_start.ljust(14, '0'), '%Y%m%d%H%M%S')

    global dtime_sim_start
    dtime_sim_start = datetime.datetime.strptime(
        args.dtime_sim_start.ljust(14, '0'), '%Y%m%d%H%M%S')

    print('reading data...')

    sami_data, times = SAMI.read_sami_dene_tec(args.sami_data_path,
                                               dtime_sim_start, reshape=True)

    if args.sami_data_path2:
        TWO_FILES = True
        global times2, sami_data2
        sami_data2, times2 = SAMI.read_dene_tec(
            args.sami_data_path2, dtime_sim_start,
            reshape=True)

    else:
        TWO_FILES = False

    global glats, glons
    glats = []
    glons = []
    nlons = sami_data['grid']['glon'].shape[0]
    nlats = sami_data['grid']['glat'].shape[2]
    for lat in range(nlats):
        for lon in range(nlons):
            glon_i = sami_data['grid']['glon'][lon,0,lat]
            if glon_i > 180:
                glons.append(glon_i - 360)
            else:
                glons.append(glon_i)
            glats.append(sami_data['grid']['glat'][lon,0,lat])

    glons = np.array(glons)
    glats = np.array(glats)

    global hrs_since_storm_onset
    hrs_full = [(i - pd.Timestamp(dtime_storm_start))
                / pd.Timedelta('1 hour') for i in times]

    ins = []
    hrs_since_storm_onset = []
    t_start = dtime_storm_start - datetime.timedelta(hours= args.plot_start_delta)
    t_end = dtime_storm_start + datetime.timedelta(hours=args.plot_end_delta)
    new_times = []
    for index, t in enumerate(times):
        if t > t_start and t < t_end:
            ins.append(index)
            new_times.append(t)
            hrs_since_storm_onset.append((pd.Timestamp(t) -
                                           dtime_storm_start)/
                                         pd.Timedelta(1, 'hour'))
    times = new_times

    print("Calculating fits....")

    if args.altitudes is not None:
        given_alts = np.asarray(args.altitudes)
        nalts = sami_data['grid']['glon'].shape[1]
        alts = sami_data['grid']['alt']
        alt_idxs = []
        #find closest alts...
        for a in given_alts:
            alt_idxs.append(np.argmin(np.abs(a - alts)))
        print('\n=======\n you gave %s alt inputs. \n' 
            %(str(given_alts)),
            'available alts are: \n',
            alts,'\nThe closest I found are:\n',
            alts[alt_idxs])
        sami_edens = []
        actual_alts = []
        for a in alt_idxs:
            sami_edens.append(sami_data['data']["edens"][ins,:,a,:].reshape(
                [len(times), nlons, nlats]))
            actual_alts.append(alts[a])
        # for plotting:
        DO_ALT_PLOTS = True
        alt_num = 0
    
    else:
        sami_tec = sami_data['data']["tec"][ins].reshape(
            [len(times), nlons, nlats])
        
        fits_sami = filters.make_fits(sami_tec)
        DO_ALT_PLOTS = False

    if TWO_FILES:
        if DO_ALT_PLOTS:
            raise ValueError("Not yet compatible with alt plots.")

        sami_tec2 = sami_data2['data'][tec][ins].reshape(
            [len(times), len(lons), len(lats)])
        global fits_sami2
        fits_sami2 = filters.make_fits(sami_tec2)

    if DO_ALT_PLOTS:
        print('Done! Shape of edens data: ', len(sami_edens),
            sami_edens[0].shape,
            '\nLaunching plotting routine now.')

    else:
        print('Done! Shape of TEC data: ', sami_tec.shape,
            '\nLaunching plotting routine now.')

    if args.figtype == 'all':
        plot_types = ['raw', 'fit', 'diff']
    else:
        plot_types = args.figtype

    try:
        cbar_lims_dict = {
        'ONE_FILE': {
            'raw': [0, .7*np.max(sami_edens)],
            'fit': [0, .7*np.max(sami_edens)],
            'diff': [-4, 4]}}
    except UnboundLocalError:
        cbar_lims_dict = {
            'TWO_FILES': {
                'raw': [-5, 5], 'fit': [-5, 5], 'diff': [-5, 5]},
            'ONE_FILE': {
                'raw': [0, 70], 'fit': [0, 70], 'diff': [-.1, .1]}}

        


    DOING_PLOTS = True

    while DOING_PLOTS:

        if DO_ALT_PLOTS:
            sami_tec = sami_edens[alt_num]
            fits_sami = filters.make_fits(sami_tec)
            alt_here = round(actual_alts[alt_num], -1)

        # Now for keos:
        if args.keogram:
            pbar = tqdm(total=len(sami_keo_lons) * len(plot_types),
                        desc="Making keograms")
            for real_lon in sami_keo_lons:

                sel_pts = np.where(
                    (np.abs(glons - real_lon) < 3) &
                    (np.abs(glats) < lat_lim) )[0]

                raw = []
                fit = []
                for t in range(len(times)):
                    raw.append(
                        sami_tec[t,:,:].copy().T.flatten()[sel_pts])
                    fit.append(
                        fits_sami[t,:,:].copy().T.flatten()[sel_pts])
                raw = np.array(raw)
                fit = np.array(fit)

                if TWO_FILES:
                    raw2 = []
                    fit2 = []
                    for t in range(len(times)):
                        raw2.append(
                            sami_tec2[t,:,:].copy().T.flatten()[sel_pts])
                        fit2.append(
                            fits_sami2[t,:,:].copy().T.flatten()[sel_pts])

                for plot_type in plot_types:
                    if plot_type == 'raw':
                        tec = raw.copy()
                        if TWO_FILES:
                            tec -= raw2
                            title = "Diff of Raw TEC at lon = {}".format(real_lon)
                            cbar_lims = cbar_lims_dict['TWO_FILES']['raw']
                        else:
                            title = "Raw TEC at lon = {}".format(real_lon)
                            cbar_lims = cbar_lims_dict['ONE_FILE']['raw']

                        fname = os.path.join(
                            out_path, 'keo',
                            "raw", "lon" + str(int(real_lon)),
                            'tec' + ".png")

                        cbar_name = 'Vertically Integrated TEC'

                        if DO_ALT_PLOTS:
                            title = title.replace('TEC','edens at %i km' 
                                %(alt_here))
                            fname = fname.replace('keo','keo/edens-at-%i-km' 
                                %(alt_here))
                            cbar_name = 'Electron Density'


                        data, extent = plotting_routines.interpolate_2d_plot(
                            hrs_since_storm_onset, glats[sel_pts], tec, 
                            len(hrs_since_storm_onset), 80)

                        plotting_routines.make_a_keo(data.T, title, cbar_lims,
                                   cbar_name=cbar_name,
                                   fname=fname, OVERWRITE=OVERWRITE,
                                   extent = extent)
                        pbar.update()

                    if plot_type == 'fit':
                        tec = fit.copy()
                        if TWO_FILES:
                            tec -= fit2
                            title = "Diff of Fit TEC at lon = {}".format(real_lon)
                            cbar_lims = cbar_lims_dict['TWO_FILES']['fit']
                        else:
                            title = "Fit TEC at lon = {}".format(real_lon)
                            cbar_lims = cbar_lims_dict['ONE_FILE']['fit']
                        
                        cbar_name = 'Vertically Integrated TEC'

                        if DO_ALT_PLOTS:
                            title = title.replace('TEC','edens at %i km' 
                                %(alt_here))
                            fname = fname.replace('keo','keo/edens-at-%i-km' 
                                %(alt_here))
                            cbar_name = 'Electron Density'


                        data, extent = plotting_routines.interpolate_2d_plot(
                            hrs_since_storm_onset, glats[sel_pts], tec, 
                            len(hrs_since_storm_onset), 80)

                        plotting_routines.make_a_keo(data.T, title, cbar_lims,
                                   cbar_name=cbar_name,
                                   fname=fname, OVERWRITE=OVERWRITE,
                                   extent = extent)
                        pbar.update()

                    if plot_type == 'diff':
                        tec = raw.copy() - fit.copy()
                        if TWO_FILES:
                            tec -= (100*(raw2.copy() - fit2.copy()) 
                                        / raw2.copy())
                            title = "Diff of % over BG of TEC at lon = {}".format(
                                real_lon)
                            cbar_lims = cbar_lims_dict['TWO_FILES']['diff']
                        else:
                            title = "TEC over Background at lon = {}".format(
                                real_lon)
                            cbar_lims = cbar_lims_dict['ONE_FILE']['diff']
                        fname = os.path.join(
                            out_path, 'keo',
                            'diff', "lon" + str(int(real_lon)),
                            'tec' + ".png")

                        cbar_name = 'Vertically Integrated TEC'

                        if DO_ALT_PLOTS:
                            tec = 100 * tec / (raw.copy())
                            title = title.replace('TEC','edens at %i km' 
                                %(alt_here))
                            fname = fname.replace('keo','keo/edens-at-%i-km' 
                                %(alt_here))
                            cbar_name = '% over Background Electron Density'


                        data, extent = plotting_routines.interpolate_2d_plot(
                            hrs_since_storm_onset, glats[sel_pts], tec, 
                            len(hrs_since_storm_onset), 80)

                        plotting_routines.make_a_keo(data.T, title, cbar_lims,
                                   cbar_name=cbar_name,
                                   fname=fname, OVERWRITE=OVERWRITE,
                                   extent = extent)
                        pbar.update()

            pbar.close()

        if args.map:
            pbar = tqdm(total=len(times) * len(plot_types),
                        desc="Making maps")

            for nt, dtime in enumerate(times):
                raw = sami_tec[nt, :, :].copy()
                fit = fits_sami[nt, :, :].copy()
                if TWO_FILES:
                    raw2 = sami_tec2[nt, :, :].copy()
                    fit2 = fits_sami2[nt, :, :].copy()
                for plot_type in plot_types:
                    if plot_type == 'raw':
                        tec = raw.copy()
                        if TWO_FILES:
                            tec -= raw2
                            title = "Diff of Raw TEC at {} from storm onset".\
                                format(UT_from_Storm_onset(
                                    dtime, dtime_storm_start))
                            cbar_lims = cbar_lims_dict['TWO_FILES']['raw']
                        else:
                            title = "Raw TEC at {} from storm onset".format(
                                UT_from_Storm_onset(dtime, dtime_storm_start))
                            cbar_lims = cbar_lims_dict['ONE_FILE']['raw']
                        fname = os.path.join(
                            out_path, 'map',
                            "raw", str(nt).rjust(3, '0') + ".png")

                        data, extent = plotting_routines.interpolate_2d_plot(
                            glons, glats, tec, 80, 100, map = True)

                        cbar_name = 'Vertically Integrated TEC'

                        if DO_ALT_PLOTS:
                            title = title.replace('TEC','edens at %i km' 
                                %(alt_here))
                            fname = fname.replace('map','map/edens-at-%i-km' 
                                %(alt_here))
                            cbar_name = 'Electron Density'

                        plotting_routines.draw_map(data.T, title, cbar_lims,
                                   cbar_label=cbar_name,
                                   fname=fname, OVERWRITE=OVERWRITE,
                                   ylims = (-lat_lim, lat_lim), plot_extent = extent, )
                        pbar.update()

                    if plot_type == 'fit':
                        tec = fit.copy()
                        if TWO_FILES:
                            tec -= fit2
                            title = "Diff of Fit TEC at {} from storm onset".\
                                format(UT_from_Storm_onset(
                                    dtime, dtime_storm_start))
                            cbar_lims = cbar_lims_dict['TWO_FILES']['fit']
                        else:
                            title = "Fit TEC at {} from storm onset".format(
                                UT_from_Storm_onset(dtime, dtime_storm_start))
                            cbar_lims = cbar_lims_dict['ONE_FILE']['fit']
                        fname = os.path.join(
                            out_path, 'map',
                            "fit", str(nt).rjust(3, '0') + ".png")

                        data, extent = plotting_routines.interpolate_2d_plot(
                            glons, glats, tec, 80, 100, map = True)

                        cbar_name = 'Vertically Integrated TEC'

                        if DO_ALT_PLOTS:
                            title = title.replace('TEC','edens at %i km' 
                                %(alt_here))
                            fname = fname.replace('map','map/edens-at-%i-km' 
                                %(alt_here))
                            cbar_name = 'Electron Density'

                        plotting_routines.draw_map(data.T, title, cbar_lims,
                                   cbar_label=cbar_name,
                                   fname=fname, OVERWRITE=OVERWRITE,
                                   ylims = (-lat_lim, lat_lim), plot_extent = extent, )
                        pbar.update()

                    if plot_type == 'diff':
                        tec = raw.copy() - fit.copy()
                        if TWO_FILES:
                            tec -= (raw2-fit2)
                            title = ("Diff of % over BG of TEC at " +
                                     UT_from_Storm_onset(
                                         dtime, dtime_storm_start) +
                                     " from storm onset")
                            cbar_lims = cbar_lims_dict['TWO_FILES']['diff']
                        else:
                            title = ("TEC over background at {} from storm onset".
                                     format(UT_from_Storm_onset(
                                         dtime, dtime_storm_start)))
                            cbar_lims = cbar_lims_dict['ONE_FILE']['diff']

                        fname = os.path.join(
                            out_path, 'map', "diff", str(nt).rjust(3, '0')
                            + ".png")

                        cbar_name = 'Vertically Integrated TEC'

                        if DO_ALT_PLOTS:
                            tec = 100 * tec / raw.copy()
                            title = title.replace('TEC','Electron Density % over Background at %i km' 
                                %(alt_here))
                            fname = fname.replace('map','map/edens-at-%i-km' 
                                %(alt_here))
                            cbar_name = '% over Background Electron Density'

                        data, extent = plotting_routines.interpolate_2d_plot(
                            glons, glats, tec, 80, 100, map = True)

                        plotting_routines.draw_map(data.T, title, cbar_lims,
                                   cbar_label=cbar_name,
                                   fname=fname, OVERWRITE=OVERWRITE,
                                   ylims = (-lat_lim, lat_lim), plot_extent = extent, )

                        pbar.update()

            pbar.close()

        if DO_ALT_PLOTS:
            alt_num +=1
            print('done with %i/%i loops' %(alt_num, len(alt_idxs)))
            if alt_num == len(alt_idxs):
                DOING_PLOTS = False
        else:
            DOING_PLOTS = False




if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="This script will plot keograms of the post-processed sami tec data.")

    parser.add_argument(
        'dtime_storm_start',
        help='Datetime of storm start. Format YYYYMMDDHHmmss',
        action='store')

    parser.add_argument(
        'dtime_sim_start',
        help='Datetime of simulation start. Format YYYYMMDDHHmmss',
        action='store')

    parser.add_argument(
        '-sami_data_path', type=str,
        help='Path to sami data', default='./sami_dir', action='store')

    parser.add_argument(
        '-sami_data_path2', type=str,
        help='Path to sami data', default=None, action='store')

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
        "--lat_lim", type=float, default=90, action='store',
        help="limit plotted latitudes to this +/- in all made plots")

    parser.add_argument(
        '--figtype', type=str, action='store', nargs = '+', default='all',
        help='Which type of plot to make.' +
        'Options: raw, filt, diff. Default: all')

    parser.add_argument(
        "-k", "--keogram", action="store_true",
        help="do you want to make a keogram?")

    parser.add_argument(
        "-m", "--map", action="store_true", help="do you want to make a map?")

    parser.add_argument(
        "-o", "--overwrite", action="store_true",
        help="overwrite existing files?")

    parser.add_argument(
        '--altitudes', type=float, nargs="+",
        action='store', default=None, required=False,
        help='Altitudes to plot electron density. Note: setting this will'+
        'not plot TEC anymore, just electron density at the closest '+
        'given altitudes.')

    args = parser.parse_args()

    main(args)



