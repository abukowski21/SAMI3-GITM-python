import os
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import time
from cartopy import crs as ccrs
from cartopy.feature.nightshade import Nightshade
from utility_programs import utils
import numpy as np
import pandas as pd
from utility_programs.utils import ut_to_lt, add_lt_to_dataset
from scipy.interpolate import LinearNDInterpolator, interp1d
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable


def make_a_keo(arr,
               title=None,
               cbarlims=[None, None],
               cbar_name=None,
               y_label="Latitude (deg)",
               x_label="Index of Time-step (relative to plot start)",
               p_extent=None,
               save_or_show="save", cmap='viridis',
               fname=None,
               ax=None,
               ylims=None,
               OVERWRITE=False,
               **kwargs):
    """
    Generate a keogram plot from a data array.

    Parameters
    ----------
    arr : xarray DataArray
        Data array to be plotted
    title : str, optional
        Title of plot, defaults to None
    cbarlims : list, optional
        Colorbar limits, Defaults to [None, None] (automatic)
    cbar_name : str, optional
        Label for colorbar, defaults to None
    y_label : str, optional
        Y-axis label, defaults to "Latitude (deg)"
    x_label : str, optional
        Label for x-axis, defaults to "Index of Time-step (relative to plot start)"
    p_extent : list-like, optional
        Plot extent. Set to 0 for auto.
    save_or_show : str, optional
        Save, show, or return plots? Defaults to "save"
    cmap : str, optional
        Colormap to use to plot, defaults to 'viridis'
    fname : str, optional
        Filename to save plot to, defaults to None
    ax : Matplotlib axes, optional
        If plotting on an existing axis, specify it here, defaults to None
    ylims : int, optional
        Limits of y axis, defaults to None
    OVERWRITE : bool, optional
        Overwrite existing plot, defaults to False

    Raises
    ------
    ValueError
        If fname exists and OVERWRITE is False
    ValueError
        If save_or_show is not "save", "show", or "return"
    ValueError
        If fname is not given when save_or_show is "save"

    Returns
    -------
    Matplotlib imshow object
        Data from imshow
    """
    if fname is not None and os.path.exists(fname) and save_or_show == "save":
        if not OVERWRITE:
            raise ValueError("We cannot overwrite the file: " + str(fname))

    if save_or_show != 'return':
        fig, ax = plt.subplots(figsize=(10, 7))

    if p_extent is None:
        p_extent = [0, len(arr.T) if type(arr) is np.ndarray
                    else len(arr.time.values),
                    -90, 90]

    data = ax.imshow(
        arr.T,
        cmap,
        aspect="auto",
        origin="lower",
        vmin=cbarlims[0],
        vmax=cbarlims[1],
        interpolation="bicubic",
        interpolation_stage="rgba",
        extent=p_extent,
        **kwargs)

    if save_or_show != "return":
        plt.title(title)
        plt.xlabel(x_label)
        plt.ylabel(y_label)

    if ylims is not None:
        plt.ylim(ylims)

    if cbar_name is not None and save_or_show != "return":
        fig.colorbar(data)
    elif save_or_show != "return":
        fig.colorbar(data, label=cbar_name)

    if save_or_show == "return":
        return data

    elif save_or_show == "show":
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
                    last_slash = fname.rfind('/')
                    os.makedirs(fname[:last_slash])
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


def draw_map(data_arr,
             cbarlims=[None, None],
             title=None, ylims=None,
             cbar_label=None,
             y_label="Latitude (deg)",
             x_label="Longitude (deg)",
             save_or_show="save",
             cmap='viridis',
             ax=None,
             fname=None,
             OVERWRITE=False,
             **kwargs):
    """
    Draw a map of the data array.

    Parameters
    ----------
    data_arr : xarray DataArray
        Data array to be plotted
    cbarlims : list, optional
        Colorbar limits as [vmin, vmax], defaults to None (automatic limits)
    title : str, optional
        Title to draw, defaults to None
    ylims : list, optional
        Limits of y-axis [ymin, ymax], defaults to None
    cbar_label : str, optional
        Label of colorbar, defaults to None
    y_label : str, optional
        Label of y-axis, defaults to "Latitude (deg)"
    x_label : str, optional
        Label of x-axis, defaults to "Longitude (deg)"
    save_or_show : str, optional
        Save or show plots? Defaults to "save"
    cmap : str, optional
        Which matplotlib colormap to use, defaults to 'viridis'
    ax : Matplotlib axes, optional
        If plotting on an existing axis, specify it here, defaults to None
    fname : str, optional
        Filename when saving, defaults to None
    OVERWRITE : bool, optional
        Overwrite existing files when saving, defaults to False

    Raises
    ------
    ValueError
        If fname exists and OVERWRITE is False
    ValueError
        If save_or_show is not "save", "show", or "return"
    ValueError
        If fname is not given when save_or_show is "save"

    Returns
    -------
    Matplotlib imshow object
        Data from imshow
    """

    if save_or_show == "save":
        if os.path.exists(fname):
            if not OVERWRITE:
                raise ValueError("We cannot overwrite the file: " + str(fname))

    if ax is None and save_or_show != "return":
        fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()},
                               figsize=(10, 5))

    elif ax is None and save_or_show == "return":
        raise ValueError("Cannot return figure if ax is not given.")

    elif ax is not None:
        fig = ax.figure

    # Create a world map background
    ax.coastlines(zorder=3, color='black', alpha=1)
    ax.gridlines(color='black', linestyle='--', alpha=0.6)

    data = ax.imshow(data_arr.T, cmap, zorder=10, alpha=0.8,
                     extent=[-180, 180, -90, 90], vmin=cbarlims[0],
                     vmax=cbarlims[1], interpolation="bicubic",
                     interpolation_stage="rgba", **kwargs)

    if ylims is not None:
        plt.ylim(ylims)

    if not cbar_label and save_or_show != "return":
        plt.colorbar(data)
    elif save_or_show != "return":
        plt.colorbar(data, label=cbar_label)

    if save_or_show == "return":
        return data

    elif save_or_show == "show":
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
                    last_slash = fname.rfind('/')
                    os.makedirs(fname[:last_slash])
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

            except PermissionError:
                print(fname, 'Has permission errors')
            plt.close("all")

    else:
        if save_or_show != "return":
            raise ValueError('save_or_show input is invalid. '
                             'Accepted inputs are "save", "show"'
                             'or "return", you gave ', save_or_show)


def draw_field_line_plot(x, y, z, title=None, interpolate=False,
                         x_label='Mlat (Deg)', y_label='Altitude (km)',
                         x_lims=[-65, 65], y_lims=[0, 1200],
                         cbar_label=None, cbar_lims=None,
                         ax=None, cmap='viridis',
                         fname=None, save_or_show='save',
                         fpeak_col=None):
    """
    Draw a plot of data along a single field line (longitude).

    Parameters
    ----------
    x : numpy array
        X-axis data, or latitude
    y : numpy array
        Y-axis data, or altitude
    z : numpy array
        Z-axis (color) data, or data to be plotted
    title : str, optional
        Tiitle of generated plot, by default None
    interpolate : bool, optional
        Interpolate the data (to the x-y grid), by default False
    x_label : str, optional
        Label of x-axis, by default 'Mlat (Deg)'
    y_label : str, optional
        Label of y-axis, by default 'Altitude (km)'
    x_lims : list, optional
        Limits of x-axis, by default [-65, 65]
    y_lims : list, optional
        Limits of y-axis, by default [0, 1200]
    cbar_label : str, optional
        Label of colorbar, by default None
    cbar_lims : list, optional
        Limits of colorbar, by default None
    ax : Matplotlib axes, optional
        If plotting on an existing axis, specify it here, by default None
    cmap : str, optional
        Which matplotlib colormap to use, by default 'viridis'
    fname : str, optional
        Filename when saving, by default None
    save_or_show : str, optional
        Save or show plots? Defaults to "save", by default 'save'
    fpeak_col : numpy array, optional
        To plot the F-peak location, give the latitude-altitude coordinates here, by default None

    Raises
    ------
    ValueError
        If fname exists and OVERWRITE is False
    ValueError
        If save_or_show is not "save", "show", or "return"
    ValueError
        If fname is not given when save_or_show is "save"

    Returns
    -------
    Matplotlib imshow object
        Data from imshow
    """

    if cbar_lims is None:
        cbar_lims = [min(z), max(z)]

    if ax is None and save_or_show != "return":
        fig, ax = plt.subplots(figsize=(10, 5))

    elif ax is None and save_or_show == "return":
        raise ValueError("Cannot return figure if ax is not given.")

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
            x, y)

        loc_grid = list(zip(x, y))
        interp = LinearNDInterpolator(
            loc_grid, z,
            rescale=True)

        znew = interp(list(zip(grid_x.flatten(), grid_y.flatten())))

        data = ax.imshow(znew.reshape(grid_x.shape), origin='lower',
                         vmin=cbar_lims[0], vmax=cbar_lims[1], cmap=cmap,
                         extent=[x_lims[0], x_lims[1], y_lims[0], y_lims[1]],
                         aspect='auto')

    else:
        data = ax.scatter(x, y, c=z, vmin=cbar_lims[0], vmax=cbar_lims[1],
                          cmap=cmap)

    plt.ylim(y_lims)
    plt.xlim(x_lims)

    if save_or_show != 'return':
        plt.colorbar(label=cbar_label)
        plt.figtext(.5, .9, title, fontsize=12, ha='center')
        plt.xlabel(x_label)
        plt.ylabel(y_label)

    if fpeak_col is not None:
        xs = []
        ys = []
        bins = np.arange(np.min(x), np.max(x), 3)

        for i in range(len(bins) - 1):
            mlat_mask = (x >= bins[i]) & (x <= bins[i + 1])
            max_here = (np.argmax(fpeak_col[mlat_mask]))
            xs.append(x[mlat_mask][max_here])
            ys.append(y[mlat_mask][max_here])

        cubic_interpolation_model = interp1d(xs, ys, kind="quadratic")
        newx = np.linspace(min(xs), max(xs), 100)
        newy = cubic_interpolation_model(newx)

        ax.plot(newx, newy, c='k')

    if save_or_show == 'return':
        return data

    elif save_or_show == 'show':
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
            plt.close()

    else:
        raise ValueError(
            'save_or_show input is invalid. Accepted inputs are "save" or',
            '"show", you gave ',
            save_or_show,)


def panel_plot(da,
               x='time',
               y='lat',
               wrap_col='lon',
               plot_vals=[0, 45, 90, 135, 180, 225, 270, 315],
               do_map=False,
               col_wrap=4,
               suptitle=None,
               vlims=None,
               cmap='bwr',
               out_fname=None,
               isel_plotvals=False,
               tight_layout=False,
               cbar_label=None
               ):
    """
    Plot a panel plot of a data array.

    Parameters
    ----------
    da : xr.DataArray
        Data array to be plotted
    x : str, optional
        X-axis data, by default 'time'
    y : str, optional
        Y-axis data, by default 'lat'
    wrap_col : str, optional
        Column to wrap around, by default 'lon'
    plot_vals : list, optional
        Values (of wrap_col) to plot, by default
        [0, 45, 90, 135, 180, 225, 270, 315]
    do_map : bool, optional
        Plot on a map?, by default False
    col_wrap : int, optional
        How many columns to make, by default 4
    suptitle : str, optional
        Title of plot, by default None
    vlims : int, optional
        Absolute value of colorbar limits, by default None (auto).
    cmap : str, optional
        Which matplotlib colormap to use, by default 'bwr'
    out_fname : str, optional
        Filename when saving, by default None
    isel_plotvals : bool, optional
        Set to true if the plot_vals are indices
        rather than values, by default False
    tight_layout : bool, optional
        Use tight layout?, by default False
    cbar_label : str, optional
        Label of colorbar, by default None

    Returns
    -------
    None
    """

    if suptitle is None:
        add_cbar = False
    elif cbar_label is not None:
        add_cbar = False
    else:
        add_cbar = True

    if do_map:
        if isel_plotvals:
            p = da.isel({wrap_col: plot_vals}).plot(
                x=x, y=y, col=wrap_col,
                transform=ccrs.PlateCarree(),
                subplot_kws={"projection": ccrs.PlateCarree(),
                             },
                col_wrap=col_wrap, vmin=-vlims, vmax=vlims,
                add_colorbar=add_cbar,
                cmap=cmap, aa=True)
        else:
            if vlims is not None:
                p = da.sel({wrap_col: plot_vals}, method='nearest').plot(
                    x=x, y=y, col=wrap_col,
                    transform=ccrs.PlateCarree(),
                    subplot_kws={"projection": ccrs.PlateCarree(),
                                 },
                    add_colorbar=add_cbar,
                    col_wrap=col_wrap, vmin=-vlims, vmax=vlims,
                    cmap=cmap, aa=True)
            else:
                p = da.sel({wrap_col: plot_vals}, method='nearest').plot(
                    x=x, y=y, col=wrap_col,
                    transform=ccrs.PlateCarree(),
                    add_colorbar=add_cbar,
                    subplot_kws={"projection": ccrs.PlateCarree(),
                                 },
                    col_wrap=col_wrap,
                    cmap=cmap, aa=True)
        for ax in p.axs.flatten():
            ax.coastlines(alpha=0.6)
            ax.gridlines(color='black', alpha=0.5, linestyle='--')

    elif not isel_plotvals:
        if vlims is not None:
            p = da.sel({wrap_col: plot_vals}, method='nearest').plot(
                x=x, y=y, col=wrap_col,
                add_colorbar=add_cbar,
                col_wrap=col_wrap, vmin=-vlims, vmax=vlims,
                cmap=cmap, aa=True)
        else:
            p = da.sel({wrap_col: plot_vals}, method='nearest').plot(
                x=x, y=y, col=wrap_col,
                add_colorbar=add_cbar,
                col_wrap=col_wrap,
                cmap=cmap, aa=True)

    else:
        p = da.isel({wrap_col: plot_vals}).plot(
            x=x, y=y, col=wrap_col,
            add_colorbar=add_cbar,
            col_wrap=col_wrap, vmin=-vlims, vmax=vlims,
            cmap=cmap, aa=True)

    p.fig.suptitle(suptitle)

    if tight_layout:
        p.fig.tight_layout()

    if not add_cbar:
        p.add_colorbar()
    if cbar_label is not None:
        p.cbar.set_label(cbar_label)

    if out_fname is None:
        plt.show()
        plt.close('all')
    else:
        plt.savefig(out_fname)
        plt.close('all')


def panel_with_lt(da,
                  x='time',
                  y='lat',
                  wrap_col='lon',
                  lons=[0, 45, 90, 135, 180, 225, 270, 315],
                  col_wrap=4,
                  suptitle=None,
                  figsize=None,
                  vlims=None,
                  cmap='bwr',
                  out_fname=None,
                  tight_layout=False,):
    """
    Plot a panel plot of a data array, with local time on the x-axis.

    Parameters
    ----------
    da : DataArray
        Data array to be plotted
    x : str, optional
        X-axis data, defaults to 'time'
    y : str, optional
        Y-axis data, defaults to 'lat'
    wrap_col : str, optional
        Column to wrap around, defaults to 'lon'
    lons : list, optional
        Values (of wrap_col) to plot, defaults to
        [0, 45, 90, 135, 180, 225, 270, 315]
    col_wrap : int, optional
        How many columns to make, defaults to 4
    suptitle : str, optional
        Title of plot, defaults to None
    figsize : list, optional
        Size of figure, defaults to None (automatic)
    vlims : int, optional
        Absolute value of colorbar limits, defaults to None (automatic)
    cmap : str, optional
        Which matplotlib colormap to use, defaults to 'bwr'
    out_fname : str, optional
        Filename when saving, defaults to None
    tight_layout : bool, optional
        Use tight layout? Defaults to False
    """

    lons = np.array(lons).reshape(int(len(lons) / col_wrap), col_wrap)

    fig, axs = plt.subplots(
        nrows=lons.shape[0], ncols=lons.shape[1], sharey=True)
    if figsize is None:
        fig.set_size_inches([3 * col_wrap + 1, axs.shape[0] * 3])
    else:
        fig.set_size_inches(figsize)

    for i in np.ndindex(lons.shape):
        j = da.sel({wrap_col: lons[i]}, method='nearest').copy()
        if vlims is None:
            j.plot(x=x, y=y, ax=axs[i], cmap=cmap)
        else:
            j.plot(x=x, y=y, ax=axs[i],
                   vmin=-vlims, vmax=vlims,
                   cmap=cmap)

        axs[i].set_xlabel('Local Time (Hours)')
        tick_locs = axs[i].get_xticks()
        axs[i].set_xticks(tick_locs)
        try:
            lts = utils.ut_to_lt([pd.Timestamp(t) for t in da.time.values[
                ::int(np.floor(len(da.time.values) / len(tick_locs)))]],
                lons[i])
            axs[i].set_xticklabels(lts.astype('int'))
        except ValueError:
            lts = utils.ut_to_lt([pd.Timestamp(t) for t in da.time.values[
                ::int(np.ceil(len(da.time.values) / len(tick_locs)))]],
                lons[i])
            axs[i].set_xticklabels(lts.astype('int'))
        axs[i].set_title('Glon=%i°' % int(lons[i]))

        if i[1] != 0:
            axs[i].set_ylabel('')
        if i[0] != lons.shape[0] - 1:
            axs[i].set_xlabel('')

    fig.suptitle(suptitle)
    # 'Alt=%i km, Date=%s\nRun=%s' %(int(400), str(j.time.dt.date.values[0]),
    #       'Full storm'))
    if tight_layout:
        fig.tight_layout()

    if out_fname is None:
        plt.show()
        plt.close('all')
    else:
        plt.savefig(out_fname)
        plt.close('all')


def panel_of_dials(da, hemi_titles, times, time_titles=None, title=None,
                   mask_dials=False, lon_labels=False, **plotargs):
    """
    Plot a panel of polar dials.

    Parameters
    ----------
    da : xarray DataArray
        Data array to be plotted.
    hemi_titles : list, optional
        Titles of each hemisphere.
    times : list, optional
        Times to plot.
    time_titles : list, optional
        Titles of each time, defaults to None.
    title : str, optional
        Title of plot, defaults to None.
    mask_dials : int, optional
        Mask values below this value, defaults to False.
    lon_labels : bool, optional
        Show longitude labels, defaults to False.

    Returns
    -------
    Matplotlib imshow object
        Data from imshow (a figure).

    Raises
    ------
    ValueError
        If hemi_titles and times do not match in shape.
    """

    if time_titles is not None:
        time_titles = np.array(time_titles)
    hemi_titles = np.array(hemi_titles)
    times = np.array(times)

    if hemi_titles.shape != times.shape:
        raise ValueError('Title specs must match!')

    fig = plt.figure(
        figsize=[
            4 * hemi_titles.shape[1],
            4 * hemi_titles.shape[0]])

    axs = []

    for i, ax in enumerate(
            np.ndindex(
            hemi_titles.shape[0], hemi_titles.shape[1])):

        time_here = pd.Timestamp(str(times[ax])).to_pydatetime()

        # Find central longitude (midnight Local time)
        lons = np.arange(0, 360)
        lts = ut_to_lt([time_here], lons)
        central_lon = lons[np.argmin(np.abs(24 - lts))]

        if hemi_titles[ax] == 'North':
            axs.append(
                fig.add_subplot(
                    hemi_titles.shape[0],
                    hemi_titles.shape[1],
                    i + 1,
                    projection=ccrs.Orthographic(
                        central_lon,
                        90)))

        elif hemi_titles[ax] == 'South':
            axs.append(
                fig.add_subplot(
                    hemi_titles.shape[0],
                    hemi_titles.shape[1],
                    i + 1,
                    projection=ccrs.Orthographic(
                        central_lon - 180,
                        -90)))

        if not mask_dials:
            da.sel(time=times[ax],
                   method='nearest').load().plot(x='lon',
                                                 y='lat',
                                                 ax=axs[-1],
                                                 transform=ccrs.PlateCarree(),
                                                 cbar_kwargs={'label': ""},
                                                 **plotargs)
        else:
            da.where(np.abs(da) >= mask_dials).sel(time=times[ax],
                                                   method='nearest').plot(
                x='lon', y='lat', ax=axs[-1], transform=ccrs.PlateCarree(),
                cbar_kwargs={'label': ""}, **plotargs)

        if time_titles is None:
            axs[-1].set_title('%s (%s UT)' %
                              (hemi_titles[ax], str(time_here.time())[:5]))
        else:
            axs[-1].set_title('%s, %s Storm (%s UT)' % (hemi_titles[ax],
                              time_titles[ax], str(time_here.time())[:5]))

        axs[-1].coastlines(zorder=3, color='black', alpha=1)
        if lon_labels:
            axs[-1].gridlines(color='black', linestyle='--',
                              alpha=0.6, draw_labels=True)
        else:
            axs[-1].gridlines(color='black', linestyle='--', alpha=0.6)

        axs[-1].add_feature(Nightshade(time_here), alpha=0.3)

    if title is not None:
        size = plt.Text.get_size(fig.suptitle('foo'))
        fig.suptitle(title, fontsize=size + 4)

    fig.tight_layout()

    return fig


def loop_panels(da,
                ncols,
                start_time,
                time_delta='1 hour',
                sel_criteria=None,
                suptitle=None,
                title=None,
                col_names=None,
                save=None,
                lon_labels=False,
                mask_dials=False,
                **plotargs):
    """
    When making panel plots, this script is easier to interface with.

    Parameters
    ----------
    da : xarray DataArray
        Data array to be plotted.
    ncols : int
        Number of columns to plot.
    start_time : str
        Start time of plot to be converted to pd.Timestamp. Format as 'YYYY-MM-DD HH:MM:SS'; '2011-05-20' is acceptable (00:00 UT).
    time_delta : str, optional
        Time between each plot, defaults to '1 hour'.
    sel_criteria : dict, optional
        Criteria to select data (i.e. plotting multiple altitudes), defaults to None.
    suptitle : str, optional
        Title of plot, defaults to None.
    title : str, optional
        Title of each subplot, defaults to None.
    col_names : list of str, optional
        Name(s) of column(s) to plot, defaults to None (use this when plotting a DataSet instead of a DataArray).
    save : str, optional
        Filename to save plot to, defaults to None.
    lon_labels : bool, optional
        Show longitude labels, defaults to False.
    mask_dials : int, optional
        Mask values below this value, defaults to False.

    Returns
    -------
    Matplotlib imshow object
        Data from imshow (a figure).
    """

    hemititles = [['North', 'South'] for i in range(ncols)]
    hemititles = np.array(hemititles).T
    times = []
    for i in range(ncols):
        times.append(pd.Timestamp(start_time) + pd.Timedelta(time_delta) * i)

    times = np.array([times, times])
    timetitles = None
    if col_names is not None:
        timetitles = np.array([col_names, col_names])

    if sel_criteria is None:
        fig = panel_of_dials(da,
                             hemi_titles=hemititles,
                             times=times,
                             time_titles=timetitles,
                             title=title,
                             lon_labels=lon_labels,
                             mask_dials=mask_dials,
                             **plotargs)
        if suptitle is not None:
            fig.suptitle(suptitle)
        if save is not None:
            plt.savefig(save)
        else:
            return fig
    else:
        if isinstance(sel_criteria, dict):
            fig = panel_of_dials(
                da.sel(
                    sel_criteria,
                    method='nearest'),
                hemi_titles=hemititles,
                times=times,
                mask_dials=mask_dials,
                time_titles=timetitles,
                lon_labels=lon_labels,
                title=title +
                ' at %s = %i' %
                (list(
                    sel_criteria.keys())[0],
                    list(
                    sel_criteria.values())[0]) ** plotargs)
            if suptitle is not None:
                fig.suptitle(suptitle)
            if save is not None:
                plt.savefig(save)
            else:
                return fig
        else:
            for entry in sel_criteria:
                fig = panel_of_dials(
                    da.sel(
                        entry,
                        method='nearest'),
                    hemi_titles=hemititles,
                    times=times,
                    mask_dials=mask_dials,
                    lon_labels=lon_labels,
                    time_titles=timetitles,
                    title=title +
                    ' at %s = %i' %
                    (list(
                        entry.keys())[0],
                        int(
                        list(
                            entry.values())[0])),
                    ** plotargs)
                if suptitle is not None:
                    fig.suptitle(suptitle)
                if save is not None:
                    plt.savefig(save)
                else:
                    return fig


def custom_panels_keos(da, numplots=8, sel_col='localtime', max_per_row=4,
                       suptitle=None, vmin=None, vmax=None, sharex=True,
                       sharey=True, x='time', cmap='rainbow',
                       one_colorbars=True, colorbar_label=''):
    """
    A script to make a panel of keogram-like plots.

    Parameters
    ----------
    da : xarray DataArray
        Data array to be plotted.
    numplots : int, optional
        Number of plots to make, by default 8.
    sel_col : str, optional
        Column to select data from, by default 'localtime'.
    max_per_row : int, optional
        Maximum number of plots per row, by default 4.
    suptitle : str, optional
        Title of plot, by default None.
    vmin : int, optional
        Minimum value of colorbar, by default None.
    vmax : int, optional
        Maximum value of colorbar, by default None.
    sharex : bool, optional
        Share x-axis?, by default True.
    sharey : bool, optional
        Share y-axis?, by default True.
    x : str, optional
        X-axis data, by default 'time'.
    cmap : str, optional
        Which matplotlib colormap to use, by default 'rainbow'.
    one_colorbars : bool, optional
        Use one colorbar for all plots?, by default True.
    colorbar_label : str, optional
        Label of colorbar, by default ''.

    Returns
    -------
    Matplotlib imshow object
        Data from imshow (a figure).
    """

    if sel_col == 'localtime' and sel_col not in da.coords:
        da = add_lt_to_dataset(da, localtimes=90)

    nrows = int(np.ceil(numplots / max_per_row))
    ncols = max_per_row

    f, axs = plt.subplots(nrows, ncols,
                          figsize=(5 * nrows, 1.3 * ncols
                                   if suptitle is not None
                                   else 1 * ncols),
                          sharey=sharey, sharex=sharex)

    sel_list = np.linspace(da[sel_col].min().values,
                           da[sel_col].max().values,
                           numplots + 1)[:-1]

    if one_colorbars:
        if vmin is None:
            vmin = da.min().compute()
        if vmax is None:
            vmax = da.max().compute()

    for a, ax in enumerate(axs.flatten()):
        ims = da.sel({sel_col: sel_list[a]}, method='nearest').plot(
            x=x, ax=ax, cmap=cmap, vmin=vmin, vmax=vmax,
            add_colorbar=not one_colorbars)

    if one_colorbars:
        divider = make_axes_locatable(axs[nrows - 1, ncols - 1])
        cax = divider.append_axes('right', size='5%', pad=0.05, in_layout=True)
        f.colorbar(ims, cax=cax, orientation='vertical', label=colorbar_label)

        divider = make_axes_locatable(axs[0, ncols - 1])
        cax = divider.append_axes('right', size='5%', pad=0.05, in_layout=True)
        f.colorbar(ims, cax=cax, orientation='vertical', label=colorbar_label)

    if suptitle is not None:
        f.suptitle(suptitle)
    f.tight_layout()

    return f


def map_and_dials(dial_da,
                  total,
                  map_da=None,
                  max_per_row=3,
                  isel_dials=None,
                  sel_dials=None,
                  isel_map=None,
                  sel_map=None,
                  quiver_map_cols=None,
                  suptitle=None,
                  suptitlesize='large',
                  time_start=None,
                  time_delta='1 hour',
                  save=None,
                  mask_dials=0.001,
                  mask_maps=False,
                  dial_cmap='rainbow',
                  map_cmap='rainbow',
                  vmin_dial=None,
                  vmax_dial=None,
                  dial_kwargs={},
                  map_kwargs={},
                  vmin_map=None,
                  vmax_map=None,
                  several_datasets=False,
                  times_datasets=None,
                  latlon_labeled=False):
    """
    This script will make a plot of a map and polar dials for several timesteps.

    Parameters
    ----------
    dial_da : xarray DataArray
        Data array to be plotted on dials
    total : int
        Number of plots to make
    map_da : xarray DataArray, optional
        Data array to be plotted on map, defaults to None
    max_per_row : int, optional
        Maximum number of plots per row, defaults to 3
    isel_dials : dict, optional
        Indices to select dial data from, defaults to None
    sel_dials : dict, optional
        Criteria to select dial data from, defaults to None
    isel_map : dict, optional
        Indices to select map data from, defaults to None
    sel_map : dict, optional
        Criteria to select map data from, defaults to None
    quiver_map_cols : str, optional
        Columns to plot on map as quiver, defaults to None (no quivers)
    suptitle : str, optional
        Title of plot, defaults to None
    suptitlesize : str, optional
        Size of suptitle, defaults to 'large'
    time_start : str, optional
        Start time of plot, defaults to None
    time_delta : str, optional
        Time between each plot, defaults to '1 hour'
    save : str, optional
        Filename to save plot to, defaults to None (don't save)
    mask_dials : float, optional
        Mask values below this value, defaults to 0.001
    mask_maps : bool, optional
        Mask values below this value, defaults to False
    dial_cmap : str, optional
        Which matplotlib colormap to use for dials, defaults to 'rainbow'
    map_cmap : str, optional
        Which matplotlib colormap to use for maps, defaults to 'rainbow'
    vmin_dial : int, optional
        Minimum value of dial colorbar, defaults to None
    vmax_dial : int, optional
        Maximum value of dial colorbar, defaults to None
    dial_kwargs : dict, optional
        Keyword arguments for dial plot, defaults to {}
    map_kwargs : dict, optional
        Keyword arguments for map plot, defaults to {}
    vmin_map : int, optional
        Minimum value of map colorbar, defaults to None
    vmax_map : int, optional
        Maximum value of map colorbar, defaults to None
    several_datasets : bool, optional
        Plot several datasets on the same plot, defaults to False
    times_datasets : list of str, optional
        Times to plot for each dataset, defaults to None
    latlon_labeled : bool, optional
        Show latitude and longitude labels, defaults to False

    Returns
    -------
    Matplotlib imshow object
        Data from imshow (a figure)
    """

    # setup data:
    if max_per_row is None:
        max_per_row = 3

    if isel_map is None and sel_map is None:
        times = [pd.Timestamp(time_start) +
                 pd.Timedelta(time_delta) * i for i in range(total)]
        sel_map = {'time': times}

    if isel_dials is None and sel_dials is None:
        times = [pd.Timestamp(time_start) +
                 pd.Timedelta(time_delta) * i for i in range(total)]
        sel_dials = {'time': times}

    # If no colorbar limits are defined, make them for the user.

    # setup figure
    nrows = int(np.ceil(total / max_per_row))
    ncols = max_per_row

    # if not latlon_labeled else (15*ncols, 15*nrows))
    fig = plt.figure(figsize=(6 * ncols, 6 * nrows))

    gs0 = gridspec.GridSpec(nrows, ncols,
                            figure=fig,
                            wspace=.2)

    axs = []

    subgrids = []

    # keep track of times
    times_plotted = []

    for i in range(total):
        # make figure
        subgrids.append(gs0[i].subgridspec(2, 2))  # , wspace=0.5, hspace=0.3))

        if several_datasets:
            dial_data = dial_da[list(dial_da.data_vars)[i]]
            map_data = map_da[list(map_da.data_vars)[i]]
            # for local time:
            time_here = pd.Timestamp(dial_data.time.values)

        else:
            # get data (first dial then map):
            if isel_dials is None:
                new = {}
                new[list(sel_dials.keys())[0]] = list(sel_dials.values())[0][i]
                dial_data = dial_da.sel(new, method='nearest')
            else:
                new = {}
                new[list(isel_dials.keys())[0]] = list(
                    isel_dials.values())[0][i]
                dial_data = dial_da.sel(new)

            if isel_map is None:
                new = {}
                new[list(sel_map.keys())[0]] = list(sel_map.values())[0][i]
                map_data = map_da.sel(new, method='nearest')
            else:
                new = {}
                new[list(isel_map.keys())[0]] = list(isel_map.values())[0][i]
                map_data = map_da.sel(new)

        # for local time:
        if times_datasets is None:
            time_here = pd.Timestamp(dial_data.time.values)
        else:
            time_here = pd.Timestamp(times_datasets[i])
        # Find central longitude (midnight Local time)
        lons = np.arange(0, 360)
        lts = ut_to_lt([time_here], lons)
        central_lon = lons[np.argmin(np.abs(24 - lts))]

        if mask_dials:
            dial_data = dial_data.where(np.abs(dial_data) > mask_dials)
        if mask_maps:
            map_data = map_data.where(np.abs(map_data) > mask_maps)

        axs.append(fig.add_subplot(
            subgrids[-1][0, 0], projection=ccrs.Orthographic(central_lon, 90)))
        dial_data.plot(ax=axs[-1], x='lon', transform=ccrs.PlateCarree(),
                       cmap=dial_cmap, cbar_kwargs={'label': "",
                                                    'pad': 0.12},
                       vmin=vmin_dial, vmax=vmax_dial, **dial_kwargs)
        axs[-1].set_title('')
        axs[-1].add_feature(Nightshade(time_here), alpha=0.3)

        axs.append(fig.add_subplot(
            subgrids[-1][0, 1],
            projection=ccrs.Orthographic(central_lon - 180, -90)))
        dial_data.plot(ax=axs[-1], x='lon', transform=ccrs.PlateCarree(),
                       cmap=dial_cmap, cbar_kwargs={'label': "",
                                                    'pad': 0.12},
                       vmin=vmin_dial, vmax=vmax_dial, **dial_kwargs)
        axs[-1].set_title('')
        axs[-1].add_feature(Nightshade(time_here), alpha=0.3)

        axs.append(fig.add_subplot(
            subgrids[-1][1, :], projection=ccrs.PlateCarree()))

        if quiver_map_cols is None:
            map_data.plot(ax=axs[-1], x='lon', transform=ccrs.PlateCarree(),
                          cmap=map_cmap, cbar_kwargs={'label': "", },
                          vmin=vmin_map, vmax=vmax_map, **map_kwargs)
        else:
            if quiver_map_cols[2] == 'amp':
                map_data[quiver_map_cols[2]] = np.sqrt(
                    map_data[quiver_map_cols[0]]**2
                    + map_data[quiver_map_cols[1]]**2)
            map_data[quiver_map_cols[2]].plot(ax=axs[-1],
                                              x='lon',
                                              transform=ccrs.PlateCarree(),
                                              cmap=map_cmap,
                                              cbar_kwargs={'label': ""},
                                              vmin=vmin_map,
                                              vmax=vmax_map,
                                              *map_kwargs)
            map_data.coarsen(lat=3,
                             lon=2).mean().where(
                                 np.abs(map_data.lat) < 84).plot.quiver(
                                     x='lon', y='lat',
                                     u=quiver_map_cols[0],
                                     v=quiver_map_cols[1],
                                     ax=axs[-1],
                                     transform=ccrs.PlateCarree(),
                                     add_guide=False)

        if not several_datasets:
            axs[-1].set_title(str(time_here))
        else:
            axs[-1].set_title(list(map_da.data_vars)[i])
        axs[-1].add_feature(Nightshade(time_here), alpha=0.3)

        # if subtitle is None:
        #     axs[-1].set_title(str(time_here))
        # else:
        #     axs[-1].set_title(subtitle[i])

        times_plotted.append(time_here)

    times_plotted = np.array([times_plotted for i in range(total)]).flatten()
    for t, ax in enumerate(axs):
        ax.coastlines(zorder=3, color='black', alpha=1)
        ax.gridlines(
            color='black',
            linestyle='--',
            alpha=0.6,
            draw_labels='x' if latlon_labeled else False,
            xlabel_style={
                'fontsize': 'xx-small'} if latlon_labeled else None)

    plt.suptitle(suptitle, fontsize=suptitlesize)

    fig.tight_layout()

    if save is not None:
        plt.savefig(save)
        plt.close()
    else:
        return fig
