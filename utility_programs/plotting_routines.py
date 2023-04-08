import os
from matplotlib import pyplot as plt
import geopandas
import time
from scipy.interpolate import LinearNDInterpolator, interp1d
import numpy as np

world = geopandas.read_file(geopandas.datasets.get_path("naturalearth_lowres"))


def make_a_keo(
        arr,
        title = None,
        cbarlims = None,
        cbar_name = None,
        y_label="Latitude (deg)",
        x_label="Hours since storm onset",
        save_or_show="save",
        cmap='viridis',
        fname=None,
        ax=None,
        ylims = None,
        OVERWRITE=False,
        **kwargs):
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
        Defaults to save. You can instead 'show' or 'return' the plots.

    """
    if fname is not None and os.path.exists(fname) and save_or_show == "save":
        if not OVERWRITE:
            raise ValueError("We cannot overwrite the file: " + str(fname))
            
    if save_or_show != 'return':
        fig, ax = plt.subplots(figsize=(10, 7))

    data = ax.imshow(
        arr.T,
        cmap,
        aspect="auto",
        origin="lower",
        vmin=cbarlims[0],
        vmax=cbarlims[1],
        interpolation="bicubic",
        interpolation_stage="rgba",
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
        fig.colorbar(data, label=cbar_lname)

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


def draw_map(
        data_arr,
        cbarlims,
        title = None,
        ylims = None,
        cbar_label=None,
        y_label="Latitude (deg)",
        x_label="Longitude (deg)",
        save_or_show="save",
        cmap='viridis',
        ax=None,
        fname=None,
        plot_extent=[-180, 180, -90, 90],
        OVERWRITE=False,
        **kwargs):

    if save_or_show == "save":
        if os.path.exists(fname):
            if not OVERWRITE:
                return

    if ax is None and save_or_show != "return":
        fig, ax = plt.subplots(figsize=(10, 5))
        
    elif ax is None and save_or_show == "return":
        raise ValueError("Cannot return figure if ax is not given.")
    
    world.plot(ax=ax, color="white", edgecolor="black", zorder=1)
    data = ax.imshow(
        data_arr.T,
        cmap,
        aspect="auto",
        extent=plot_extent,
        origin="lower",
        zorder=10,
        ax=ax,
        alpha=0.8,
        vmin=cbarlims[0],
        vmax=cbarlims[1],
        interpolation="bicubic",
        interpolation_stage="rgba",
        **kwargs)
        
    if save_or_show != "return":
        plt.title(title)
        plt.xlabel(x_label)
        plt.ylabel(y_label)

    if ylims is not None:
        plt.ylim(ylims)

    if not cbar_label and save_or_show != "return":
        fig.colorbar(data)
    elif save_or_show != "return":
        fig.colorbar(data, label=cbar_label)

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
            raise ValueError(
                'save_or_show input is invalid. Accepted inputs are "save",',
                '"show" or "return", you gave ',
                save_or_show,)
    


def interpolate_2d_plot(x,y,c,nx_out, ny_out, map = False):

    x = np.array(x)
    y = np.array(y)
    c = np.array(c)

    if map:
        in_x = x
        in_y = y
    else:
        in_x, in_y = np.meshgrid(x,y)
    
    out_x, out_y = np.meshgrid(
        np.linspace(min(x), max(x), nx_out),
        np.linspace(min(y), max(y), ny_out))


    interp = LinearNDInterpolator(
        np.array([in_x.flatten(), in_y.flatten()]).T, c.T.flatten(),
        rescale=True)
    znew = interp(list(zip(out_x.flatten(), out_y.flatten())))
    znew = znew.reshape(out_x.shape)

    p_extent = [min(x), max(x), min(y), max(y)]

    return znew, p_extent



def draw_field_line_plot(x, y, z, title=None, interpolate=False,
                         x_label='Mlat (Deg)', y_label='Altitude (km)',
                         x_lims=[-65, 65], y_lims=[0, 1200],
                         cbar_label=None, cbar_lims=None,
                         ax = None, cmap = 'viridis',
                         fname=None, save_or_show='save',
                         fpeak_col=None):
    
    import aacgmv2

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
            x,y)

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

        for i in range(len(bins)-1):
            mlat_mask = (x >= bins[i]) & (x <= bins[i+1])
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

            except FileNotFoundError:
                print(fname)
                raise ValueError('File not found error. Check your path.')
            plt.close()

    else:
        raise ValueError(
            'save_or_show input is invalid. Accepted inputs are "save" or',
            '"show", you gave ', save_or_show)
