# import aacgmv2, time
import pandas as pd
import numpy as np
from multiprocessing import Pool
import math
import os
import shutil
from tqdm.auto import tqdm
import matplotlib.pyplot as plt
import time
# from scipy.io import readsav
# import pymap3d as pm
# import glob
import datetime
# from aetherpy.io import read_routines
# from math import cos, radians, sin, sqrt
from scipy import spatial, signal

# from spacepy.coordinates import Coords
# from spacepy.time import Ticktock
# import fnmatch
#
import gc

import sys
from mpl_toolkits.basemap import Basemap
import geopandas

from scipy.interpolate import LinearNDInterpolator, interp1d, griddata


def make_a_plot(data, x_label=None, y_label=None, title=None,
                cbar_lims=None, cbar_label=None, save_or_show='show',
                fname=None, plot_extent=None):
    # Clean inputs:
    if cbar_lims == None:
        cbar_lims = [np.min(data), np.max(data)]

    # Plots have to be done in weird ways.
    # this will make calling things easier...
    ismap = False
    iskeo = False

    if plot_extent is None:
        if data.shape == (65, 75):
            plot_extent = [-180, 180, np.min(lats), np.max(lats)]
            x_label = 'Longitude(deg)'
            y_label = 'Latitude (deg)'
            ismap = True
        if data.shape == (144, 65):
            plot_extent = [np.min(hrs), np.max(
                hrs), np.min(lats), np.max(lats)]
            x_label = 'Hours from storm onset'
            y_label = 'Latitude (deg)'
            iskeo = True

    if ismap:

        # Need to get the data from -180-180 not 0-360...

        # Fix the ordering of the longitudes and go from -180-180 not 0->360
        newlons_for_order = []
        newlons = np.zeros_like(lons)
        for ilon in lons:
            oldlon = ilon
            if oldlon <= 180:
                newlons_for_order.append(int(oldlon))

            else:
                newlons_for_order.append(int(oldlon)-360)

        new_lons_sorted = np.sort(newlons_for_order)
        new_order = np.array([newlons_for_order.index(
            new_lons_sorted[i]) for i in range(len(new_lons_sorted))])

        fig, ax = plt.subplots(figsize=(10, 5))
        world.plot(ax=ax, color='white', edgecolor='black', zorder=1)

        p = ax.imshow(data[:, new_order], extent=plot_extent,
                      vmin=cbar_lims[0], vmax=cbar_lims[1], origin='lower',
                      aspect='auto', interpolation='bicubic',
                      interpolation_stage='rgba', label=cbar_label, zorder=10,
                      alpha=0.9)
        fig.colorbar(p, label=cbar_label)

    elif iskeo:
        plt.figure(figsize=(10, 5))
        plt.imshow(data.T, extent=plot_extent, vmin=cbar_lims[0],
                   vmax=cbar_lims[1], origin='lower', aspect='auto',
                   interpolation='bicubic', interpolation_stage='rgba',
                   label=cbar_label)
        plt.colorbar()

    plt.ylabel(y_label)
    plt.xlabel(x_label)

    plt.title(title)

    if save_or_show == 'show':
        plt.show()
        plt.close()
        if fname:
            print(fname)
    elif save_or_show == 'save':
        if not fname:
            raise ValueError('plot save path must be given!')
        else:
            # fname = fname.replace(' ','')
            try:
                plt.savefig(fname)
            except FileNotFoundError:
                try:
                    directory_list = os.path.join(fname).split('/')[:-1]
                    os.makedirs('/'+os.path.join(*directory_list))
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

            except:
                print(fname)
                raise ValueError

    else:
        raise ValueError(
            'save_or_show input is invalid. Accepted inputs are "save"',
            'or "show", you gave ', save_or_show)
    plt.close()
