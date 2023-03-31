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
    global fits_sami
    sami_tec = sami_data['data']["tec"][ins].reshape(
        [len(times), nlons, nlats])
    
    fits_sami = filters.make_fits(sami_tec)

    if TWO_FILES:
        sami_tec2 = sami_data2['data'][tec][ins].reshape(
            [len(times), len(lons), len(lats)])
        global fits_sami2
        fits_sami2 = filters.make_fits(sami_tec2)

    print('Done! Shape of TEC data: ', sami_tec.shape,
        '\nLaunching plotting routine now.')

    if args.figtype == 'all':
        plot_types = ['raw', 'fit', 'diff']
    else:
        plot_types = [args.figtype]

    cbar_lims_dict = {
        'TWO_FILES': {
            'raw': [-5, 5], 'fit': [-5, 5], 'diff': [-5, 5]},
        'ONE_FILE': {
            'raw': [0, 70], 'fit': [0, 70], 'diff': [-.2, .2]}}

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

                    data, extent = plotting_routines.interpolate_2d_plot(
                        hrs_since_storm_onset, glats[sel_pts], tec, 
                        len(hrs_since_storm_onset), 80)

                    plotting_routines.make_a_keo(data.T, title, cbar_lims,
                               cbar_name='Vertically Integrated TEC',
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
                    
                    # cbar_lims = [np.min(tec), np.max(tec)]
                    fname = os.path.join(
                        out_path, 'keo',
                        'fit', "lon" + str(int(real_lon)),
                        'tec' + ".png")

                    data, extent = plotting_routines.interpolate_2d_plot(
                        hrs_since_storm_onset, glats[sel_pts], tec, 
                        len(hrs_since_storm_onset), 80)

                    plotting_routines.make_a_keo(data.T, title, cbar_lims,
                               cbar_name='Vertically Integrated TEC',
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
                        title = "% over BG of TEC at lon = {}".format(
                            real_lon)
                        cbar_lims = cbar_lims_dict['ONE_FILE']['diff']
                    fname = os.path.join(
                        out_path, 'keo',
                        'diff', "lon" + str(int(real_lon)),
                        'tec' + ".png")

                    data, extent = plotting_routines.interpolate_2d_plot(
                        hrs_since_storm_onset, glats[sel_pts], tec, 
                        len(hrs_since_storm_onset), 80)

                    plotting_routines.make_a_keo(data.T, title, cbar_lims,
                               cbar_name='Vertically Integrated TEC',
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
                        glons, glats, tec, 100, 80, map = True)

                    plotting_routines.draw_map(data.T, title, cbar_lims,
                               cbar_label='Vertically Integrated TEC',
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
                        glons, glats, tec, 100, 80, map = True)

                    plotting_routines.draw_map(data.T, title, cbar_lims,
                               cbar_label='Vertically Integrated TEC',
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
                        title = ("% over BG of TEC at {} from storm onset".
                                 format(UT_from_Storm_onset(
                                     dtime, dtime_storm_start)))
                        cbar_lims = cbar_lims_dict['ONE_FILE']['diff']
                    fname = os.path.join(
                        out_path, 'map', "diff", str(nt).rjust(3, '0')
                        + ".png")

                    data, extent = plotting_routines.interpolate_2d_plot(
                        glons, glats, tec, 80, 100, map = True)

                    plotting_routines.draw_map(data.T, title, cbar_lims,
                               cbar_label='Vertically Integrated TEC',
                               fname=fname, OVERWRITE=OVERWRITE,
                               ylims = (-lat_lim, lat_lim), plot_extent = extent, )

                    pbar.update()
        pbar.close()


# def make_a_keo(
#         arr,
#         title,
#         cbarlims,
#         cbar_name,
#         y_label="Latitude (deg)",
#         x_label="Hours since storm onset",
#         save_or_show="save",
#         fname=None,
#         keo_lat_lim=90,
#         plot_extent=None,
#         OVERWRITE=False):
#     """
#     Inputs a data array and then generates a keogram.

#     Parameters:
#     -----------
#     arr: np array
#         The data array to be plotted. If grabbing from the sami array,
#         you do not need to transpose.
#     extent: tuple/list
#         The limits of the plot. [left, right, bottom, top]
#     xlabel: string
#         self-explanitory
#     y-label: string
#         self-explanitory
#     title: string
#         self-explanitory
#     cbar limes: tuple/list
#         vmin, vmax for the colorbar to be plot.
#     cbar_name: string.
#         Label for the colorbar.
#     save_or_show: string
#         Defaults to save. You can instead 'show' the plots.

#     """
#     if fname is not None and os.path.exists(fname) and save_or_show == "save":
#         if not OVERWRITE:
#             raise ValueError("We cannot overwrite the file: " + str(fname))
#     fig = plt.figure(figsize=(10, 7))

#     if plot_extent is None:
#         hrs_start = hrs_since_storm_onset[0]
#         hrs_end = hrs_since_storm_onset[-1]
#         lat_start = -keo_lat_lim
#         lat_end = keo_lat_lim
#         plot_extent = [hrs_start, hrs_end, lat_start, lat_end]

#     plt.imshow(
#         arr.T,
#         extent=plot_extent,
#         aspect="auto",
#         cmap="viridis",
#         origin="lower",
#         vmin=cbarlims[0],
#         vmax=cbarlims[1],
#     )
#     plt.ylabel(y_label)
#     plt.xlabel(x_label)
#     plt.title(title)
#     plt.colorbar(label=cbar_name)

#     if save_or_show == "show":
#         plt.show()
#         plt.close(fig)
#     elif save_or_show == "save":
#         if not fname:
#             raise ValueError("plot save path must be given!")
#         else:
#             try:
#                 plt.savefig(fname)
#             except FileNotFoundError:
#                 try:
#                     last_slash = fname.rfind('/')
#                     os.makedirs(fname[:last_slash])
#                     plt.savefig(fname)
#                 except PermissionError:
#                     print("Permission denied. Cannot save plot.")
#                     print(" tried writing to: ", fname)
#             plt.close("all")
#     else:
#         raise ValueError(
#             'save_or_show input is invalid. Accepted inputs are "save" or',
#             '"show", you gave ',
#             save_or_show,
#         )


# def make_a_map(
#         data_arr,
#         title,
#         cbarlims,
#         cbar_label=None,
#         y_label="Latitude (deg)",
#         x_label="Longitude (deg)",
#         save_or_show="save",
#         fname=None,
#         plot_extent=[-180, 180, -90, 90],
#         OVERWRITE=False):

#     if os.path.exists(fname):
#         if not OVERWRITE:
#             return

#     fig, ax = plt.subplots(figsize=(10, 5))
#     world.plot(ax=ax, color="white", edgecolor="black", zorder=1)
#     data = ax.imshow(
#         data_arr.T,
#         cmap="viridis",
#         aspect="auto",
#         extent=plot_extent,
#         origin="lower",
#         zorder=10,
#         alpha=0.8,
#         vmin=cbarlims[0],
#         vmax=cbarlims[1],
#         interpolation="bicubic",
#         interpolation_stage="rgba",)
#     plt.title(title)
#     plt.xlabel(x_label)
#     plt.ylabel(y_label)

#     if not cbar_label:
#         fig.colorbar(data)
#     else:
#         fig.colorbar(data, label=cbar_label)

#     if save_or_show == "show":
#         plt.show()
#         plt.close()
#     elif save_or_show == "save":
#         if not fname:
#             raise ValueError("plot save path must be given!")
#         else:
#             fname = fname.replace(" ", "")
#             try:
#                 plt.savefig(fname)
#             except FileNotFoundError:
#                 try:
#                     last_slash = fname.rfind('/')
#                     os.makedirs(fname[:last_slash])
#                     plt.savefig(fname)
#                 except FileExistsError:
#                     # sometimes when we make too many plots in the same
#                     # directory, it fails. this fixes that.
#                     time.sleep(2)
#                     try:
#                         plt.savefig(fname)
#                     except FileNotFoundError:
#                         time.sleep(2)
#                         plt.savefig(fname)

#             except FileNotFoundError:
#                 print(fname)
#                 raise ValueError
#             plt.close("all")
#     else:
#         raise ValueError(
#             'save_or_show input is invalid. Accepted inputs are "save" or',
#             '"show", you gave ',
#             save_or_show,)


# def make_filter(lowcut=100, highcut=30):
#     """_summary_

#     Args:
#         lowcut (int, optional): lowcut of the filter. Defaults to 100.
#         highcut (int, optional): highcut of the filter. Defaults to 30.

#     Returns:
#         scipy butterworth filter: the filter with settings defined by the user.
#     """
#     # Define the cutoff frequencies
#     lowcut = 1 / (100 / 60)  # 100 minutes in units of sample^-1
#     highcut = 1 / (30 / 60)  # 30 minutes in units of sample^-1

#     # Define the Butterworth filter
#     nyquist = 0.5 * 5  # 5 minutes is the sampling frequency
#     low = lowcut / nyquist
#     high = highcut / nyquist
#     sos = signal.butter(2, [low, high], btype="bandstop", output="sos")
#     return sos


# def make_fits(gitm_bins):
#     """
#     calculate bandpass filter for all data previously read in.

#     inputs: nparray of gitmdata

#     returns:
#     fits: np array indexed at fits[time][col][ilon][ilat][ialt]


#     todo: you can thread this by splitting the alts into different threads.
#     then just append the fits_full later.

#     """
#     sos = make_filter()

#     filtered_arr = signal.sosfiltfilt(sos, gitm_bins, axis=0)
#     return filtered_arr


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="This script will plot keograms of the post-calc sami tec data.")

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
        '--figtype', type=str, action='store', default='all',
        help='Which type of plot to make.' +
        'Options: raw, filt, diffs. Default: all')

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