import os
from matplotlib import pyplot as plt
import geopandas
import time

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
        raise ValueError(
            'save_or_show input is invalid. Accepted inputs are "save" or',
            '"show", you gave ',
            save_or_show,)
