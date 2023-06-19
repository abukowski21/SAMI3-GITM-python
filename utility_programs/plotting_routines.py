import os
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import geopandas
import time
from cartopy import crs as ccrs
from cartopy.feature.nightshade import Nightshade
from utility_programs import utils
import numpy as np
import pandas as pd
from utility_programs.utils import ut_to_lt

world = geopandas.read_file(geopandas.datasets.get_path("naturalearth_lowres"))


def make_a_keo(
        arr,
        title,
        cbarlims,
        cbar_name,
        y_label="Latitude (deg)",
        x_label="Hours since storm onset",
        save_or_show="save",
        fname=None,
        plot_extent=None,
        OVERWRITE=False):
    """
    Inputs a data array and then generates a keogram.

    Parameters:
    -----------
    arr: np array
        The data array to be plotted. If grabbing from the gitm array,
        you do not need to transpose.
    plot_extent: tuple/list
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

    plt.imshow(
        arr.T,
        extent=plot_extent if plot_extent else None,
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


def panel_plot(da,
               x='time',
               y='lat',
               wrap_col='lon',
               plot_vals=[0, 45, 90, 135, 180, 225, 270, 315],
               do_map=False,
               col_wrap=4,
               suptitle=None,
               vlims=2,
               cmap='bwr',
               out_fname=None,
               isel_plotvals=False,
               tight_layout=False,
               ):

    if do_map:
        if isel_plotvals:
            p = da.isel({wrap_col: plot_vals}).plot(
                x=x, y=y, col=wrap_col,
                transform=ccrs.PlateCarree(),
                subplot_kws={"projection": ccrs.PlateCarree(),
                             },
                col_wrap=col_wrap, vmin=-vlims, vmax=vlims,
                cmap=cmap, aa=True)
        else:
            if vlims is not None:
                p = da.sel({wrap_col: plot_vals}, method='nearest').plot(
                    x=x, y=y, col=wrap_col,
                    transform=ccrs.PlateCarree(),
                    subplot_kws={"projection": ccrs.PlateCarree(),
                                 },
                    col_wrap=col_wrap, vmin=-vlims, vmax=vlims,
                    cmap=cmap, aa=True)
            else:
                p = da.sel({wrap_col: plot_vals}, method='nearest').plot(
                    x=x, y=y, col=wrap_col,
                    transform=ccrs.PlateCarree(),
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
                col_wrap=col_wrap, vmin=-vlims, vmax=vlims,
                cmap=cmap, aa=True)
        else:
            p = da.sel({wrap_col: plot_vals}, method='nearest').plot(
                x=x, y=y, col=wrap_col,
                col_wrap=col_wrap,
                cmap=cmap, aa=True)

    else:
        p = da.isel({wrap_col: plot_vals}).plot(
            x=x, y=y, col=wrap_col,
            col_wrap=col_wrap, vmin=-vlims, vmax=vlims,
            cmap=cmap, aa=True)

    p.fig.suptitle(suptitle)

    if tight_layout:
        p.fig.tight_layout()

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
                  do_map=False,
                  col_wrap=4,
                  suptitle=None,
                  figsize=None,
                  vlims=None,
                  cmap='bwr',
                  out_fname=None,
                  isel_plotvals=False,
                  tight_layout=False,):

    lons = np.array(lons).reshape(int(len(lons)/col_wrap), col_wrap)

    fig, axs = plt.subplots(
        nrows=lons.shape[0], ncols=lons.shape[1], sharey=True)
    if figsize is None:
        fig.set_size_inches([3*col_wrap + 1, axs.shape[0]*3])
    else:
        fig.set_size_inches(figsize)

    for i in np.ndindex(lons.shape):
        j = da.sel({wrap_col: lons[i]}, method='nearest').copy()
        if vlims is None:
            im = j.plot(x=x, y=y, ax=axs[i], cmap=cmap)
        else:
            im = j.plot(x=x, y=y, ax=axs[i], vmin=-
                        vlims, vmax=vlims, cmap=cmap)

        axs[i].set_xlabel('Local Time (Hours)')
        tick_locs = axs[i].get_xticks()
        axs[i].set_xticks(tick_locs)
        try:
            lts = utils.ut_to_lt([pd.Timestamp(t) for t in da.time.values[
                ::int(np.floor(len(da.time.values)/len(tick_locs)))]], lons[i])
            axs[i].set_xticklabels(lts.astype('int'))
        except ValueError:
            lts = utils.ut_to_lt([pd.Timestamp(t) for t in da.time.values[
                ::int(np.ceil(len(da.time.values)/len(tick_locs)))]], lons[i])
            axs[i].set_xticklabels(lts.astype('int'))
        axs[i].set_title('Glon=%i°' % int(lons[i]))

        if i[1] != 0:
            axs[i].set_ylabel('')
        if i[0] != lons.shape[0]-1:
            axs[i].set_xlabel('')

    fig.suptitle(suptitle)
    # 'Alt=%i km, Date=%s\nRun=%s' %(int(400), str(j.time.dt.date.values[0]), 'Full storm'))
    if tight_layout:
        fig.tight_layout()

    if out_fname is None:
        plt.show()
        plt.close('all')
    else:
        plt.savefig(out_fname)
        plt.close('all')


def panel_of_dials(da, hemi_titles, times,
                   time_titles=None,
                   title=None,
                   mask_dials=False,
                   lon_labels=False,
                   **plotargs):

    if time_titles is not None:
        time_titles = np.array(time_titles)
    hemi_titles = np.array(hemi_titles)
    times = np.array(times)

    if hemi_titles.shape != times.shape:
        raise ValueError('Title specs must match!')

    fig = plt.figure(figsize=[4*hemi_titles.shape[1], 4*hemi_titles.shape[0]])

    axs = []

    for i, ax in enumerate(np.ndindex(hemi_titles.shape[0], hemi_titles.shape[1])):

        time_here = pd.Timestamp(str(times[ax])).to_pydatetime()

        # Find central longitude (midnight Local time)
        lons = np.arange(0, 360)
        lts = ut_to_lt([time_here], lons)
        central_lon = lons[np.argmin(np.abs(24-lts))]

        if hemi_titles[ax] == 'North':
            axs.append(fig.add_subplot(
                hemi_titles.shape[0], hemi_titles.shape[1], i+1, projection=ccrs.Orthographic(central_lon, 90)))

        elif hemi_titles[ax] == 'South':
            axs.append(fig.add_subplot(
                hemi_titles.shape[0], hemi_titles.shape[1], i+1, projection=ccrs.Orthographic(central_lon-180, -90)))

        if mask_dials == False:
            da.sel(time=times[ax], method='nearest').load().plot(x='lon', y='lat', ax=axs[-1],
                                                                 transform=ccrs.PlateCarree(),
                                                                 cbar_kwargs={
                'label': ""},
                **plotargs)
        else:
            da.where(np.abs(da) >= mask_dials).sel(time=times[ax], method='nearest').plot(x='lon', y='lat', ax=axs[-1],
                                                                                          transform=ccrs.PlateCarree(),
                                                                                          cbar_kwargs={
                'label': ""},
                **plotargs)

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
        fig.suptitle(title, fontsize=size+4)

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

    hemititles = [['North', 'South'] for i in range(ncols)]
    hemititles = np.array(hemititles).T
    times = []
    for i in range(ncols):
        times.append(pd.Timestamp(start_time) + pd.Timedelta(time_delta)*i)

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
        if type(sel_criteria) == dict:
            fig = panel_of_dials(da.sel(sel_criteria, method='nearest'), hemi_titles=hemititles,
                                 times=times, mask_dials=mask_dials,
                                 time_titles=timetitles,
                                 lon_labels=lon_labels,
                                 title=title + ' at %s = %i' % (list(sel_criteria.keys())[0],
                                                                list(sel_criteria.values())[0])
                                 ** plotargs)
            if suptitle is not None:
                fig.suptitle(suptitle)
            if save is not None:
                plt.savefig(save)
            else:
                return fig
        else:
            for entry in sel_criteria:
                fig = panel_of_dials(da.sel(entry, method='nearest'), hemi_titles=hemititles,
                                     times=times, mask_dials=mask_dials,
                                     lon_labels=lon_labels,
                                     time_titles=timetitles,
                                     title=title + ' at %s = %i' % (list(entry.keys())[0],
                                                                    int(list(entry.values())[0])),
                                     ** plotargs)
                if suptitle is not None:
                    fig.suptitle(suptitle)
                if save is not None:
                    plt.savefig(save)
                else:
                    return fig

                
def custom_panels_keos(da,
                       numplots=8,
                       sel_col='localtime',
                       max_per_row=4,
                       suptitle=None,
                       vmin=None,
                       vmax=None,
                       sharex=True,
                       sharey=True,
                       x='time',
                       cmap='rainbow',
                       one_colorbars=True,
                       colorbar_label=''):

    if sel_col == 'localtime' and sel_col not in da.coords:
        # print('adding')
        da = add_lt_to_dataset(da, localtimes=90)
    # print(da.coords)
    # return da

    nrows = int(np.ceil(numplots/max_per_row))
    ncols = max_per_row

    f, axs = plt.subplots(nrows,
                          ncols,
                          figsize=(5*nrows, 1.3*ncols if suptitle is not None else 1*ncols),
                          sharey=sharey,
                          sharex=sharex)

    sel_list = np.linspace(da[sel_col].min().values,
                           da[sel_col].max().values,
                           numplots+1)[:-1]
    
    if one_colorbars:
        if vmin is None:
            vmin = da.min().compute()
        if vmax is None:
            vmax = da.max().compute()
            
    for a, ax in enumerate(axs.flatten()):
        ims = da.sel({sel_col: sel_list[a]}, method='nearest').plot(
            x=x, ax=ax, cmap=cmap,vmin=vmin, vmax=vmax, add_colorbar=not no_colorbar)
        
    if one_colorbars:
        divider = make_axes_locatable(axs[nrows-1, ncols-1])
        cax = divider.append_axes('right', size='5%', pad=0.05, in_layout=True)
        f.colorbar(ims, cax=cax, orientation='vertical', label=colorbar_label)

        divider = make_axes_locatable(axs[0, ncols-1])
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
                  time_start=None,
                  time_delta='1 hour',
                  save=None,
                  mask_dials=0.001,
                  dial_cmap='rainbow',
                  dial_JH_defaults=True,
                  map_cmap='rainbow',
                  vmin_dial=None,
                  vmax_dial=None,
                  vmin_map=None,
                  vmax_map=None,
                  several_datasets=False,
                  times_datasets=None):

    # setup data:
    if max_per_row is None:
        max_per_row = 3

    if isel_map is None and sel_map is None:
        times = [pd.Timestamp(time_start) +
                 pd.Timedelta(time_delta)*i for i in range(total)]
        sel_map = {'time': times}

    if isel_dials is None and sel_dials is None:
        times = [pd.Timestamp(time_start) +
                 pd.Timedelta(time_delta)*i for i in range(total)]
        sel_dials = {'time': times}

    # If no colorbar limits are defined, make them for the user.

    # setup figure
    nrows = int(np.ceil(total/max_per_row))
    ncols = max_per_row

    fig = plt.figure(figsize=(5*ncols, 5*nrows))

    gs0 = gridspec.GridSpec(nrows, ncols,
                            figure=fig)

    axs = []

    subgrids = []

    # keep track of times
    times_plotted = []

    for i in range(total):
        # make figure
        subgrids.append(gs0[i].subgridspec(2, 2))

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
        central_lon = lons[np.argmin(np.abs(24-lts))]

        if mask_dials != False:
            dial_data = dial_data.where(np.abs(dial_data) > mask_dials)

        axs.append(fig.add_subplot(
            subgrids[-1][0, 0], projection=ccrs.Orthographic(central_lon, 90)))
        dial_data.plot(ax=axs[-1], x='lon', transform=ccrs.PlateCarree(),
                       cmap=dial_cmap, cbar_kwargs={'label': "", },
                       vmin=vmin_dial, vmax=vmax_dial)
        axs[-1].set_title('')
        axs[-1].add_feature(Nightshade(time_here), alpha=0.3)

        axs.append(fig.add_subplot(
            subgrids[-1][0, 1], projection=ccrs.Orthographic(central_lon-180, -90)))
        dial_data.plot(ax=axs[-1], x='lon', transform=ccrs.PlateCarree(),
                       cmap=dial_cmap, cbar_kwargs={'label': "", },
                       vmin=vmin_dial, vmax=vmax_dial)
        axs[-1].set_title('')
        axs[-1].add_feature(Nightshade(time_here), alpha=0.3)

        axs.append(fig.add_subplot(
            subgrids[-1][1, :], projection=ccrs.PlateCarree()))

        if quiver_map_cols is None:
            map_data.plot(ax=axs[-1], x='lon', transform=ccrs.PlateCarree(),
                          cmap=map_cmap, cbar_kwargs={'label': "", },
                          vmin=vmin_map, vmax=vmax_map)
        else:
            if quiver_map_cols[2] == 'amp':
                map_data[quiver_map_cols[2]] = np.sqrt(
                    map_data[quiver_map_cols[0]]**2
                    + map_data[quiver_map_cols[1]]**2)
            map_data[quiver_map_cols[2]].plot(ax=axs[-1], x='lon', transform=ccrs.PlateCarree(),
                                              cmap=map_cmap,
                                              cbar_kwargs={'label': ""},
                                              vmin=vmin_map, vmax=vmax_map)
            map_data.coarsen(lat=3, lon=2).mean().where(np.abs(map_data.lat) < 84)\
                .plot.quiver(x='lon', y='lat', u=quiver_map_cols[0], v=quiver_map_cols[1],
                             ax=axs[-1], transform=ccrs.PlateCarree(),
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
        ax.gridlines(color='black', linestyle='--', alpha=0.6)

    plt.suptitle(suptitle)

    fig.tight_layout()

    if save is not None:
        plt.savefig(save)
        plt.close()
    else:
        return fig