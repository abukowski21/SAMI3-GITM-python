import argparse
import datetime
import gc
import glob
import os
import time
from multiprocessing import Pool

import geopandas
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from aetherpy.io import read_routines
from scipy import signal
from tqdm.auto import tqdm

from utility_programs.plot_help import UT_from_Storm_onset

matplotlib.use("Agg")

np.seterr(divide='ignore')


def main(args):
    """Main plotting routine`.

    Args:
        args (namespace): alll the args

    Raises:
        ValueError: _description_
        ValueError: _description_
    """
    # Set variables
    global world
    world = geopandas.read_file(
        geopandas.datasets.get_path("naturalearth_lowres"))

    global gitm_colnames_friendly
    gitm_colnames_friendly = {
        "Rho": "Total Neutral Density",
        "[O(!U3!NP)]": "O(3P)",
        "[O!D2!N]": "O2",
        "[N!D2!N]": "N2",
        "[N(!U4!NS)]": "N(4S)",
        "[NO]": "NO",
        "[He]": "He",
        "[N(!U2!ND)]": "N(2D)",
        "[N(!U2!NP)]": "N(2P)",
        "[H]": "H",
        "[CO!D2!N]": "CO2",
        "[O(!U1!ND)]": "O(1D)",
        "Temperature": "Temperature",
        "V!Dn!N(east)": "Vn(east)",
        "V!Dn!N(north)": "Vn(north)",
        "V!Dn!N(up)": "Vn(up)",
        "V!Dn!N(up,O(!U3!NP))": "Vn(up,O(3P))",
        "V!Dn!N(up,O!D2!N)": "Vn(up,O2)",
        "V!Dn!N(up,N!D2!N)": "Vn(up,N2)",
        "V!Dn!N(up,N(!U4!NS))": "Vn(up,N(4S))",
        "V!Dn!N(up,NO)": "Vn(up,NO)",
        "V!Dn!N(up,He)": "Vn(up,He)",
        "[O_4SP_!U+!N]": "O(4Sp)+",
        "[NO!U+!N]": "NO+",
        "[O!D2!U+!N]": "O2+",
        "[N!D2!U+!N]": "N2+",
        "[N!U+!N]": "N+",
        "[O(!U2!ND)!U+!N]": "O(2D)+",
        "[O(!U2!NP)!U+!N]": "O(2P)+",
        "[H!U+!N]": "H+",
        "[He!U+!N]": "He+",
        "[e-]": "e-",
        "eTemperature": "eTemperature",
        "iTemperature": "iTemperature",
        "V!Di!N(east)": "Vi(east)",
        "V!Di!N(north)": "Vi(north)",
        "V!Di!N(up)": "Vi(up)", }

    # files
    gitm_files = np.sort(
        glob.glob(os.path.join(args.gitm_data_path, args.file_type)))
    if len(gitm_files) == 0:
        raise ValueError("No GITM Binary files found!",
                         "found these files: ",
                         os.listdir(args.gitm_data_path))

    # Columns to plot

    global cols
    if args.cols == 'all':
        cols = gitm_colnames_friendly.keys()
    else:
        cols = []
        for c in args.cols:
            if c in gitm_colnames_friendly.keys():
                cols.append(c)
            else:
                raise ValueError('col %s not found in: \n' %
                                 c, gitm_colnames_friendly.keys())

    # Lon to keo:
    global gitm_keo_lons
    gitm_keo_lons = args.keo_lons

    global keo_lat_lim
    keo_lat_lim = args.keo_lat_lim

    global out_path
    out_path = args.out_path

    global dtime_storm_start
    dtime_storm_start = datetime.datetime.strptime(
        args.dtime_storm_start.ljust(14, '0'), '%Y%m%d%H%M%S')

    gitm_dtimes = []
    for i in gitm_files:
        yy, MM, dd = i[-17:-15], i[-15:-13], i[-13:-11]
        hr, mm, sec = i[-10:-8], i[-8:-6], i[-6:-4]
        gitm_dtimes.append(
            datetime.datetime(
                int("20" + yy), int(MM), int(dd), int(hr), int(mm), int(sec)))

    plot_start_idx = np.argmin(np.abs(np.array(gitm_dtimes)
                                      - (dtime_storm_start -
                                         datetime.timedelta(
                                             hours=args.plot_start_delta))))

    plot_end_idx = (np.argmin(np.abs(np.array(gitm_dtimes)
                                     - (dtime_storm_start +
                                        datetime.timedelta(
                                            hours=args.plot_end_delta))))
                    if args.plot_end_delta != -1
                    else -1)

    # Read in grid & data files.
    global times, gitm_grid, gitm_vars, gitm_bins
    times, gitm_grid, gitm_vars, gitm_bins = read_gitm_into_nparrays(
        gitm_files[plot_start_idx:plot_end_idx], cols)

    print(gitm_vars, cols, gitm_bins.shape)

    global hrs_since_storm_onset
    hrs_since_storm_onset = np.array([(i - pd.Timestamp(dtime_storm_start))
                                      / pd.Timedelta('1 hour') for i in times])

    global lats, lons, alts
    lats, lons, alts = (
        np.unique(gitm_grid["latitude"]),
        np.unique(gitm_grid["longitude"]),
        np.unique(gitm_grid["altitude"]))

    if args.gitm_alt_idxs == -1:
        gitm_alt_idxs = list(range(len(alts)))
    else:
        gitm_alt_idxs = args.gitm_alt_idxs

    global fits_gitm
    print("Calculating fits. This will take a moment...")
    fits_gitm = make_fits(gitm_bins)

    # Start plotting.
    if args.keogram:
        print("Making keogram")
        if args.threading:
            arg_arr = []
            for namecol in cols:
                numcol = cols.index(namecol)
                for nalt in gitm_alt_idxs:
                    for nlon in gitm_keo_lons:
                        arg_arr.append([numcol, namecol, nalt, nlon,
                                        args.outliers])
            with Pool(args.num_workers) as pool:
                with tqdm(desc='threading keo making',
                          total=len(arg_arr)) as pbar:
                    for _ in pool.imap_unordered(thread_call_keos, arg_arr):
                        pbar.update(1)

        else:
            pbar = tqdm(
                total=len(cols) * len(gitm_alt_idxs) * len(gitm_keo_lons),
                desc='keogram making')
            for alt_idx in gitm_alt_idxs:
                for real_lon in gitm_keo_lons:
                    for col in cols:
                        call_keos(alt_idx=alt_idx, real_lon=real_lon,
                                  namecol=col, save_or_show=args.save_or_show,
                                  outliers=args.outliers,
                                  keo_lat_lim=args.keo_lat_lim,
                                  figtype=args.figtype)
                        pbar.update()
            pbar.close()

    if args.map:
        print("Making map")
        if args.threading:
            arg_arr = []
            for namecol in cols:
                numcol = cols.index(namecol)
                for nalt in gitm_alt_idxs:
                    for dtime_real in times:
                        arg_arr.append([nalt, namecol, dtime_real,
                                        args.outliers])
            with Pool(args.num_workers) as pool:
                with tqdm(desc='threading map making',
                          total=len(arg_arr)) as pbar:
                    for _ in pool.imap_unordered(thread_call_maps, arg_arr):
                        pbar.update(1)
        else:
            pbar = tqdm(total=len(gitm_alt_idxs) * len(times) * len(cols),
                        desc='map making')

            for col in cols:
                numcol = cols.index(col)
                for nalt in gitm_alt_idxs:
                    for dtime_real in times:
                        call_maps(nalt, dtime_real=dtime_real,
                                  numcol=numcol,
                                  save_or_show=args.save_or_show,
                                  figtype=args.figtype,
                                  outliers=args.outliers)
                        pbar.update()
            pbar.close()


def read_gitm_into_nparrays(flist, cols):
    """reads a list of gitm filenames and returns a few numpy arrays.

    Parameters
    ----------
    flist: list
        List of gitm filenames to read in.

    Returns
    -------
    gitmtimes = list
        Datetimes corresponding to the times of the gitm files.
    gitmgrid = dict.
        Holds the gitm grid for reference. lons/lats have been
        changed to degrees from rads. Ghost cells removed.
        Keys are ['longitude', 'latitude', 'altitude']
        Index with gitmgrid[time][key][lon][lat][alt]
    gitmvars = list
        The gitm variables.
    gitmbins = numpy array.
        All of the gitm data (except grid)
        index with gitmbins[time,varnumber,lat,lon,alt]

    """

    f = read_routines.read_gitm_file(flist[0])
    gitmgrid = {f["vars"][k].lower(): f[k][2:-2, 2:-2, 2:-2]
                for k in [0, 1, 2]}
    nlons, nlats, nalts = np.array(f[0].shape) - 4  # ghost cells
    gitmvars = [i for i in f["vars"][3:] if i in cols]

    gitmtimes = []
    gitmbins = np.zeros([len(flist), len(cols), nlons, nlats, nalts])

    for ifile, file_name in enumerate(tqdm(flist)):
        f = read_routines.read_gitm_file(file_name)

        gitmtimes.append(f["time"])  # times

        for num_var, real_var in enumerate(gitmvars):
            num_v_src = f["vars"].index(real_var)
            gitmbins[ifile, num_var] = f[num_v_src][2:-2, 2:-2, 2:-2]

    gitmgrid["latitude"] = np.rad2deg(gitmgrid["latitude"])
    gitmgrid["longitude"] = np.rad2deg(gitmgrid["longitude"])

    # Fix the ordering of the longitudes and go from -180-180 not 0->360
    newlons_for_order = []
    for ilon in range(len(gitmgrid["longitude"])):
        oldlon = gitmgrid["longitude"][ilon, 0, 0]
        if oldlon <= 180:
            newlons_for_order.append(int(oldlon))

        else:
            newlons_for_order.append(int(oldlon) - 360)
            gitmgrid["longitude"][ilon] = gitmgrid["longitude"][ilon] - 360

    new_lons_sorted = np.sort(newlons_for_order)
    new_order = np.array(
        [newlons_for_order.index(new_lons_sorted[i])
         for i in range(len(new_lons_sorted))])

    gitmbins = gitmbins[:, :, new_order, :, :]
    gitmgrid["longitude"] = np.sort(gitmgrid["longitude"], axis=0)

    # # enforce lat_lim:
    # if global_lat_lim != 90:
    #     lat_mask = (gitmgrid["latitude"] >= -global_lat_lim) & (
    #         gitmgrid["latitude"] <= global_lat_lim)
    #     gitmgrid["latitude"] = gitmgrid["latitude"][lat_mask]
    #     num_lats_in = len(np.unique(gitmgrid["latitude"]))
    #     gitmbins = gitmbins[:, :, lat_mask].reshape(
    #         len(flist),
    #         len(gitm_colnames_friendly.keys()),
    #         nlons,
    #         num_lats_in,
    #         nalts)

    #     gitmgrid["longitude"] = gitmgrid["longitude"][lat_mask].reshape(
    #         nlons, num_lats_in, nalts)
    #     gitmgrid["altitude"] = gitmgrid["altitude"][lat_mask].reshape(
    #         nlons, num_lats_in, nalts)
    #     gitmgrid['latitude'] = gitmgrid['latitude'].reshape(
    #         nlons, num_lats_in, nalts)

    return gitmtimes, gitmgrid, gitmvars, gitmbins


def make_filter(lowcut=100, highcut=30):
    """_summary_

    Args:
        lowcut (int, optional): lowcut of the filter. Defaults to 100.
        highcut (int, optional): highcut of the filter. Defaults to 30.

    Returns:
        scipy butterworth filter: the filter with settings defined by the user.
    """
    # Define the cutoff frequencies
    lowcut = 1 / (100 / 60)  # 100 minutes in units of sample^-1
    highcut = 1 / (30 / 60)  # 30 minutes in units of sample^-1

    # Define the Butterworth filter
    nyquist = 0.5 * 5  # 5 minutes is the sampling frequency
    low = lowcut / nyquist
    high = highcut / nyquist
    sos = signal.butter(2, [low, high], btype="bandstop", output="sos")
    return sos


def make_fits(gitm_bins):
    """
    calculate bandpass filter for all data previously read in.

    inputs: nparray of gitmdata

    returns:
    fits: np array indexed at fits[time][col][ilon][ilat][ialt]


    todo: you can thread this by splitting the alts into different threads.
    then just append the fits_full later.

    """
    sos = make_filter()

    filtered_arr = signal.sosfiltfilt(sos, gitm_bins, axis=0)
    return filtered_arr


def remove_outliers(array):
    arr2 = array.copy()
    # calculate mean, standard deviation, and median over all elements
    mean, std, median = np.mean(arr2), np.std(arr2), np.median(arr2)
    # set outlier threshold (in terms of number of standard deviations)
    outlier_threshold = 5
    outliers = np.logical_or(
        arr2 < mean - outlier_threshold * std, arr2 > mean +
        outlier_threshold * std)  # find outliers
    arr2[outliers] = median  # set outliers to median
    return arr2

# KEO MAKING FUNCTIONS:


def make_a_keo(
        arr,
        title,
        cbarlims,
        cbar_name,
        y_label="Latitude (deg)",
        x_label="Hours since storm onset",
        save_or_show="save",
        fname=None,
        keo_lat_lim=90,
        plot_extent=None,
        OVERWRITE=False):
    """
    Inputs a data array and then generates a keogram.

    Parameters:
    -----------
    arr: np array
        The data array to be plotted. If grabbing from the gitm array,
        you do not need to transpose.
    extent: tuple/list
        The limits of the plot. [left, right, bottom, top]
    xlabel: string
        self-explanitory
    y-label: string
        self-explanitory
    title: string
        self-explanitory
    cbar limes: tuple/list
        vmin, vmax for the colorbar to be plot.
    cbar_name: string.
        Label for the colorbar.
    save_or_show: string
        Defaults to save. You can instead 'show' the plots.

    """
    if fname is not None and os.path.exists(fname) and save_or_show == "save":
        if not OVERWRITE:
            raise ValueError("We cannot overwrite the file: " + str(fname))
    fig = plt.figure(figsize=(10, 7))

    if plot_extent is None:
        hrs_start = hrs_since_storm_onset[0]
        hrs_end = hrs_since_storm_onset[-1]
        lat_start = -keo_lat_lim
        lat_end = keo_lat_lim
        plot_extent = [hrs_start, hrs_end, lat_start, lat_end]

    plt.imshow(
        arr.T,
        extent=plot_extent,
        aspect="auto",
        cmap="viridis",
        origin="lower",
        vmin=cbarlims[0],
        vmax=cbarlims[1],
    )
    plt.ylabel(y_label)
    plt.xlabel(x_label)
    plt.title(title)
    plt.colorbar(label=cbar_name)

    if save_or_show == "show":
        plt.show()
        plt.close(fig)
    elif save_or_show == "save":
        if not fname:
            raise ValueError("plot save path must be given!")
        else:
            try:
                plt.savefig(fname)
            except FileNotFoundError:
                try:
                    directory_list = os.path.join(fname).split("/")[:-1]
                    os.makedirs(os.path.join(*directory_list))
                    plt.savefig(fname)
                except PermissionError:
                    print("Permission denied. Cannot save plot.")
                    print(" tried writing to: ", fname)
            plt.close("all")
    else:
        raise ValueError(
            'save_or_show input is invalid. Accepted inputs are "save" or',
            '"show", you gave ',
            save_or_show,
        )


def call_keos(
        alt_idx,
        real_lon,
        numcol=None,
        namecol: str = "",
        save_or_show="show",
        return_figs=False,
        figtype="all",
        keo_lat_lim=90,
        outliers=False,
        vlims=None):

    if numcol is None and namecol != "":
        numcol = cols.index(namecol)
    elif namecol == "" and numcol is not None:
        namecol = cols[numcol]
    elif numcol is None and namecol == "":
        raise ValueError("either namecol or numcol must be specified!")

    if vlims is None:

        vmin_bins = np.min(gitm_bins[:, numcol, :, :, alt_idx])
        vmax_bins = np.max(gitm_bins[:, numcol, :, :, alt_idx])

        vmin_fits = np.min(fits_gitm[:, numcol, :, :, alt_idx])
        vmax_fits = np.max(fits_gitm[:, numcol, :, :, alt_idx])

        vmin_diffs = np.min(100 * (fits_gitm[:, numcol, :, :, alt_idx]
                                   - gitm_bins[:, numcol, :, :, alt_idx])
                            / gitm_bins[:, numcol, :, :, alt_idx])
        vmax_diffs = np.max(100 * (fits_gitm[:, numcol, :, :, alt_idx]
                                   - gitm_bins[:, numcol, :, :, alt_idx])
                            / gitm_bins[:, numcol, :, :, alt_idx])

    else:
        vmin = -vlims
        vmax = vlims
        vmin_bins = vmin
        vmax_bins = vmax
        vmin_diffs = vmin
        vmax_diffs = vmax
        vmin_fits = vmin
        vmax_fits = vmax

    # get data.
    lon_idx = np.argmin(np.abs(lons - real_lon))
    real_lon = lons[lon_idx]
    data = gitm_bins[:, numcol, lon_idx, :, alt_idx].copy()
    bandpass = fits_gitm[:, numcol, lon_idx, :, alt_idx].copy()
    if np.sum(data) == 0:
        raise ValueError("No data at this altitude and longitude!",
                         data.shape, real_lon, namecol, numcol)
    real_alt = alts[alt_idx]
    percent = 100 * (data - bandpass) / bandpass

    if outliers:
        data = remove_outliers(data)
        bandpass = remove_outliers(bandpass)
        percent = remove_outliers(percent)

    made_plot = False

    if figtype == "all" or "filt" in figtype:
        # plain bandpass filter
        title = "Keogram of %s along %i deg Longitude at %i km" % (
            gitm_colnames_friendly[namecol].replace(
                "(", "[").replace(")", "]"),
            real_lon,
            round(real_alt / 1000, 0),)
        color_label = "Bandpass filter"
        # print(out_path, real_alt, real_lon, namecol)
        # print(int(real_alt/1000,0), int(real_lon))
        fname = os.path.join(
            out_path, 'keo',
            "bandpass",
            str(int(real_alt / 1000)),
            "lon" + str(int(real_lon)),
            gitm_colnames_friendly[namecol] + ".png",)
        make_a_keo(
            bandpass,
            title,
            cbarlims=(vmin_fits, vmax_fits),
            cbar_name=color_label,
            keo_lat_lim=keo_lat_lim,
            save_or_show=save_or_show,
            fname=fname,)
        made_plot = True

    if figtype == "all" or "raw" in figtype:
        # plain raw data
        title = "Keogram of %s along %i deg Longitude at %i km" % (
            gitm_colnames_friendly[namecol].replace(
                "(", "[").replace(")", "]"),
            real_lon,
            round(real_alt / 1000, 0),)
        color_label = "Raw data"
        fname = os.path.join(
            out_path, 'keo',
            "raw",
            str(int(real_alt / 1000)),
            "lon" + str(int(real_lon)),
            gitm_colnames_friendly[namecol] + ".png",)
        make_a_keo(
            data,
            title,
            keo_lat_lim=keo_lat_lim,
            cbarlims=(vmin_bins, vmax_bins),
            cbar_name=color_label,
            save_or_show=save_or_show,
            fname=fname,)
        made_plot = True

    if figtype == "all" or "diff" in figtype:
        # (bandpass - raw)/ raw
        title = "Keogram of %s along %i deg Longitude at %i km" % (
            gitm_colnames_friendly[namecol].replace(
                "(", "[").replace(")", "]"),
            real_lon,
            round(real_alt / 1000, 0),)
        color_label = "% over bandpass filter"
        fname = os.path.join(
            out_path, 'keo',
            "percent-over-filter",
            str(int(real_alt / 1000)),
            "lon" + str(int(real_lon)),
            gitm_colnames_friendly[namecol] + ".png",)
        make_a_keo(
            percent,
            title,
            keo_lat_lim=keo_lat_lim,
            cbarlims=(vmin_diffs, vmax_diffs),
            cbar_name=color_label,
            save_or_show=save_or_show,
            fname=fname,)
        made_plot = True

    if not made_plot:
        print("nothing made")


def thread_call_keos(args):
    call_keos(
        args[2],
        args[3],
        args[0],
        args[1],
        outliers=args[4],
        save_or_show="save",
        return_figs=False,
        figtype="all")


def draw_map(
        data_arr,
        title,
        cbarlims,
        cbar_label=None,
        y_label="Latitude (deg)",
        x_label="Longitude (deg)",
        save_or_show="save",
        fname=None,
        plot_extent=[-180, 180, -90, 90],
        OVERWRITE=False):

    if os.path.exists(fname):
        if not OVERWRITE:
            return

    fig, ax = plt.subplots(figsize=(10, 5))
    world.plot(ax=ax, color="white", edgecolor="black", zorder=1)
    data = ax.imshow(
        data_arr.T,
        cmap="viridis",
        aspect="auto",
        extent=plot_extent,
        origin="lower",
        zorder=10,
        alpha=0.8,
        vmin=cbarlims[0],
        vmax=cbarlims[1],
        interpolation="bicubic",
        interpolation_stage="rgba",)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    if not cbar_label:
        fig.colorbar(data)
    else:
        fig.colorbar(data, label=cbar_label)

    if save_or_show == "show":
        plt.show()
        plt.close()
    elif save_or_show == "save":
        if not fname:
            raise ValueError("plot save path must be given!")
        else:
            fname = fname.replace(" ", "")
            try:
                plt.savefig(fname)
            except FileNotFoundError:
                try:
                    directory_list = os.path.join(fname).split("/")[:-1]
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
                raise ValueError
            plt.close("all")
    else:
        raise ValueError(
            'save_or_show input is invalid. Accepted inputs are "save" or',
            '"show", you gave ',
            save_or_show,)


def call_maps(
        alt_idx,
        dtime_real=None,
        dtime_index=None,
        numcol=None,
        namecol=None,
        save_or_show="show",
        return_figs=False,
        figtype="all",
        lat_lim=[-90, 90],
        diffs=[1, 2, 3, 5, 10, 30, 50],
        outliers=False):

    # Make sure inputs are correct. either the index or actual value of the
    #    datetime and column to plot can be specified (or both).
    if numcol is None and namecol is not None:
        numcol = cols.index(namecol)
    elif namecol is None and numcol is not None:
        namecol = cols[numcol]
    elif numcol is None and namecol is None:
        raise ValueError("either namecol or numcol must be specified!")

    if dtime_real is None and dtime_index is not None:
        dtime_real = times[dtime_index]
    elif dtime_index is None and dtime_real is not None:
        dtime_index = np.argmin(np.abs(np.array(times) - dtime_real))
    elif dtime_real is None and dtime_index is None:
        raise ValueError("either dtime_index or dtime_real must be specified!")

    # get colorbar limits.
    vmin_bins = np.min(gitm_bins[:, numcol, :, :, alt_idx])
    vmax_bins = np.max(gitm_bins[:, numcol, :, :, alt_idx])

    vmin_fits = np.min(fits_gitm[:, numcol, :, :, alt_idx])
    vmax_fits = np.max(fits_gitm[:, numcol, :, :, alt_idx])

    if type(diffs) != list:

        vmin_diffs = np.min(100 * (fits_gitm[:, numcol, :, :, alt_idx]
                                   - gitm_bins[:, numcol, :, :, alt_idx])
                            / gitm_bins[:, numcol, :, :, alt_idx])
        vmax_diffs = np.max(100 * (fits_gitm[:, numcol, :, :, alt_idx]
                                   - gitm_bins[:, numcol, :, :, alt_idx])
                            / gitm_bins[:, numcol, :, :, alt_idx])

    # get data.
    raw = gitm_bins[dtime_index, numcol, :, :, alt_idx].copy()
    bandpass = fits_gitm[dtime_index, numcol, :, :, alt_idx].copy()
    real_alt = alts[alt_idx]
    percent = 100 * (bandpass - raw) / raw

    if outliers:
        raw = remove_outliers(raw)
        bandpass = remove_outliers(bandpass)
        percent = remove_outliers(percent)

    made_plot = False

    # raw map
    if figtype == "all" or "raw" in figtype:
        title = (
            gitm_colnames_friendly[namecol]
            + " at "
            + str(round(float(real_alt) / 1000, 0))
            + " km at "
            + UT_from_Storm_onset(dtime_real, dtime_storm_start)
            + " from Storm Start")
        fname = os.path.join(
            out_path, 'maps'
            "raw",
            str(int(real_alt / 1000)),
            gitm_colnames_friendly[namecol],
            str(dtime_index).rjust(3, "0") + ".png",)
        cbarlims = [vmin_bins, vmax_bins]
        draw_map(raw, title, cbarlims, fname=fname, save_or_show=save_or_show)
        made_plot = True

    # filter map
    if figtype == "all" or "filt" in figtype:
        title = (
            gitm_colnames_friendly[namecol]
            + " at "
            + str(round(float(real_alt) / 1000, 0))
            + " km at "
            + UT_from_Storm_onset(dtime_real, dtime_storm_start)
            + " from Storm Start")
        fname = os.path.join(
            out_path, 'maps'
            "bandpass",
            str(int(real_alt / 1000)),
            gitm_colnames_friendly[namecol],
            str(dtime_index).rjust(3, "0") + ".png",)
        cbarlims = [vmin_fits, vmax_fits]
        cbar_label = "Bandpass Filtered " + gitm_colnames_friendly[namecol]
        draw_map(
            bandpass,
            title,
            cbarlims,
            save_or_show=save_or_show,
            cbar_label=cbar_label,
            fname=fname,)
        made_plot = True

    # diffs
    if figtype == "all" or "diff" in figtype:
        title = (
            gitm_colnames_friendly[namecol]
            + " at "
            + str(round(float(real_alt) / 1000, 0))
            + " km at "
            + UT_from_Storm_onset(dtime_real, dtime_storm_start)
            + " from Storm Start")
        if type(diffs) == list:
            for v_lim in diffs:
                fname = os.path.join(
                    out_path, 'maps'
                    "diff_set_lims",
                    str(int(real_alt / 1000)),
                    gitm_colnames_friendly[namecol],
                    str(v_lim),
                    str(dtime_index).rjust(3, "0") + ".png",)
                cbarlims = [-v_lim, v_lim]
                cbar_label = "% over Background"
                draw_map(
                    percent,
                    title,
                    cbarlims,
                    save_or_show=save_or_show,
                    cbar_label=cbar_label,
                    fname=fname,)
        else:
            fname = os.path.join(
                out_path, 'maps'
                "diff",
                str(int(real_alt / 1000)),
                gitm_colnames_friendly[namecol],
                str(dtime_index).rjust(3, "0") + ".png",)
            cbarlims = [vmin_diffs, vmax_diffs]
            cbar_label = "% over Background"
            draw_map(
                percent,
                title,
                cbarlims,
                save_or_show=save_or_show,
                cbar_label=cbar_label,
                fname=fname,)
            made_plot = True

    if not made_plot:
        print("No plot made. Check figtype input.")

    plt.close("all")
    gc.collect()


def thread_call_maps(args):
    call_maps(
        alt_idx=args[0],
        namecol=args[1],
        dtime_real=args[2],
        outliers=args[3],
        save_or_show="save",
        return_figs=False,
        figtype="all",)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Make plots of GITM data at a given altitude and time.")

    parser.add_argument(
        'dtime_storm_start',
        help='Datetime of storm start. Format YYYYMMDDHHmmss',
        action='store')

    parser.add_argument(
        '-gitm_data_path', type=str,
        help='Path to gitm data', default='./gitm_dir', action='store')

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
        action='store', default=48, required=False)

    parser.add_argument(
        "-f", "--file-type", type=str, nargs="+",
        default="3DALL*",
        help="which filetype to plot, e.g. 3DALL* or 2DANC*",)

    parser.add_argument(
        "-o", "--outliers", action="store_true",
        help="do you want to remove outliers")

    parser.add_argument(
        "-d", "--diff", type=int, nargs="+",
        help="do you want to plot the difference between the reg and raw ",)

    parser.add_argument(
        "-k", "--keogram", action="store_true",
        help="do you want to make a keogram?")

    parser.add_argument(
        '--gitm_alt_idxs', type=int, nargs="*",
        default=[5, 10, 15, 22, 30, 45],
        help='Which altitudes to plot. Default: 5,10,15,22,30,45')

    parser.add_argument(
        "-m", "--map", action="store_true", help="do you want to make a map?")

    args = parser.parse_args()

    main(args)
