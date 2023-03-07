#!/usr/bin/env python
# coding: utf-8


import aacgmv2
import time
import pandas as pd
import numpy as np
from multiprocessing import Pool
import math
import os
import shutil
from tqdm.auto import tqdm
import matplotlib.pyplot as plt
from scipy.io import readsav
import pymap3d as pm
import glob
import datetime
import statistics
from aetherpy.io import read_routines
from math import cos, radians, sin, sqrt
from scipy import spatial, signal

from spacepy.coordinates import Coords
from spacepy.time import Ticktock
import fnmatch

import gc

import sys
from mpl_toolkits.basemap import Basemap
import geopandas

from scipy.interpolate import LinearNDInterpolator, interp1d, griddata


# ## Settings


dtime_storm_start = datetime.datetime(2011, 5, 21, 13, 40)

dtime_sim_start = datetime.datetime(2011, 5, 20)

t_step_minutes = 5  # minutes


# hours before storm onset to start making plots. set to -1 to run the whole time
plot_start_delta = 3
# hours after storm onset to end plots. Set to -1 to run for the whole time
plot_end_delta = 9

sami_data_path = "/home/axb170054/scratch/GITM-testing/test_folders/step_function_driving/SAMI3-stretch/"


global_lat_lim = 65  # will limit SAMI output data to this GEO latitude


OVERWRITE = True  # be careful!


# SHOULD be the same as t_step_minutes above. If not, something will be wrong. Won't throw errors. Gotta pay attention.
sample_rate_min = 5
low_cut = 80  # min, lowest freq wave the filter will allow thru
high_cut = 40  # min, highest freq the filter will allow thru

# Do we want to run the filtering? if not, this will need to be done later.
to_filter = False  # I'll adjust this later..
# make plots at the end? I would leave this as False. At the bottom there are some simple plotting routines.
debug = False


cols = ['edens', 'hplusdens', 'oplusdens', 'noplusdens', 'o2plusdens', 'heplusdens',
        'n2plusdens', 'nplusdens', 'hdens', 'odens', 'nodens', 'o2dens', 'hedens', 'n2dens', 'ndens']

# above is all cols (that I care about), below is just edens.

# cols = ['edens']


# Available columns (for now) are:
#
# ['edens', 'hplusdens', 'oplusdens', 'noplusdens', 'o2plusdens', 'heplusdens', 'n2plusdens', 'nplusdens', 'hdens', 'odens', 'nodens', 'o2dens', 'hedens', 'n2dens', 'ndens']

# In[ ]:


data_out_path = '/home/axb170054/scratch/pickles/SimStormPaper/simstorm_sami_files/'

save_raw = False  # Save original SAMI data in data_out_path?

# In[ ]:


thread = True
num_workers = len(cols)  # adjust this to fit your system.


# In[ ]:


out_grid_lats = 65
out_grid_lons = 75


#   ==========================   #


# Now it's processing time!


# out_alts = np.linspace(out_alt_min, out_alt_max, num_out_grid_alts)
out_lats = np.linspace(-global_lat_lim, global_lat_lim, out_grid_lats)
# we don't need both 0 & 360. the 1's deal with that.
out_lons = np.linspace(0, 360, out_grid_lons + 1)[1:]
out_alts = np.array([250, 300, 350, 400, 450, 500, 550,
                    600, 700, 800, 840, 880, 900, 1000])


print('\n out lats are ', *out_lats.round(2), sep=' ')
print('\n out lons are ', *out_lons.round(2), sep=' ')
print('\n out alts are ', *out_alts, sep=' ')
print('\n lengths are ', len(out_lats), len(out_lons), len(out_alts))


# In[ ]:


# ## Constants:

# In[ ]:


geo_grid_files = {'glat': 'glatu.dat', 'glon': 'glonu.dat', 'alt': 'zaltu.dat',
                  'mlat': 'blatu.dat', 'mlon': 'blonu.dat', 'malt': 'baltu.dat'}


data_files = {'edens': 'deneu.dat', 'hplusdens': 'deni1u.dat', 'oplusdens': 'deni2u.dat',
              'noplusdens': 'deni3u.dat', 'o2plusdens': 'deni4u.dat',
              'heplusdens': 'deni5u.dat', 'n2plusdens': 'deni6u.dat',
              'nplusdens': 'deni7u.dat', 'hdens': 'denn1u.dat', 'odens': 'denn2u.dat',
              'nodens': 'denn3u.dat', 'o2dens': 'denn4u.dat', 'hedens': 'denn5u.dat',
              'n2dens': 'denn6u.dat', 'ndens': 'denn7u.dat'}

time_file = 'time.dat'


# In[ ]:


# ## Define Functions

# In[ ]:


def get_grid_elems_from_parammod(data_dir):
    """
    Will look for: words = ['nz0','nf','nl'] in SAMI files.

    inputs:
    ------
    sami path

    outputs:
    -------
    nz,nf,nlt,nt :
    - nz  = num points along field line
    - nf  = num field lines along each mag lon
    - nlt = num mag lons
    - nt  = num times

    """

    # Make sure that we only grab the first instance of each var in the file.
    # SOmetimes they repeat and we don't want them
    returns = [False, False, [False, False], False]

    with open(data_dir + 'parameter_mod.f90', 'r') as fp:
        # read all lines in a list
        lines = fp.readlines()
        for line in lines:
            # check if string present on a current line

            if not returns[0]:
                if line.find('nz0') != -1:
                    nz0 = []
                    for l in line:
                        if l.isdigit():
                            nz0.append(l)
                    if len(nz0[1:4]) == 3:
                        nz = int(''.join(nz0[1:4]))
                        returns[0] = True

            if not returns[1]:
                if line.find('nf') != -1:
                    nf = []
                    for l in line:
                        if l.isdigit():
                            nf.append(l)
                    nf = int(''.join(nf))
                    returns[1] = True

            if not returns[2][0]:
                if line.find('nl ') != -1:
                    nl = []
                    for l in line:
                        if l.isdigit():
                            nl.append(l)
                    nl = int(''.join(nl))
                    returns[2][0] = True

            if not returns[2][1]:
                if line.find('numwork ') != -1:
                    numwork = []
                    for l in line:
                        if l.isdigit():
                            numwork.append(l)
                    numwork = int(''.join(numwork))
                    returns[2][1] = True

    # time
    with open(data_dir + 'time.dat', 'r') as fp:
        lines = fp.readlines()
        nt = len(lines) - 1

    return nz, nf, numwork*(nl - 2), nt


# In[ ]:


def make_times(t0, nt, plot_start_delta=None, plot_end_delta=None):
    times = []
    hrs_since_storm_start = []

    for t in range(nt):
        time_here = pd.Timestamp(dtime_sim_start) + \
            t * pd.Timedelta(5, 'minutes')
        times.append(time_here.to_pydatetime())
        hrs = (time_here - dtime_storm_start)/pd.Timedelta(1, 'hour')
        hrs_since_storm_start.append(hrs)

    times_df = pd.read_fwf(os.path.join(sami_data_path, 'time.dat'),
                           names=['istep', 'hour', 'minute', 'second', 'hrdelta'], infer_nrows=115)
    times_df.pop('istep')

    times_list = []
    for hr in times_df['hrdelta']:
        times_list.append(dtime_sim_start + datetime.timedelta(hours=hr))

    truths = np.array([pd.Timestamp(times_list[t]).round(
        'T') == times[t] for t in range(len(times))])
    if truths.sum() != len(truths):
        raise ValueError(
            'The times are wrong! Somehow this needs to be fixed. probably outputting fake files again. Take a look and debug before proceeding.')

    # maybe chop the time lists, depending on if the plot start/end are given.
    # adjusted to allow for -1 in plot start/end deltas (plot all times)

    if plot_start_delta and plot_end_delta:
        if plot_start_delta != -1:
            start_idx = np.argmin(np.abs(np.array(times)
                                         - (dtime_storm_start - pd.Timedelta(plot_start_delta, 'hour'))))
        else:
            start_idx = 0

        if plot_end_delta != -1:
            end_idx = np.argmin(np.abs(np.array(times)
                                       - (dtime_storm_start + pd.Timedelta(plot_end_delta, 'hour'))))
        elif plot_end_delta == -1:
            end_idx = len(times)
        else:
            end_idx = len(times)

        times = times[start_idx:end_idx]
        hrs_since_storm_start = hrs_since_storm_start[start_idx:end_idx]
        times_list = times_list[start_idx:end_idx]

        return times, hrs_since_storm_start, (start_idx, end_idx)

    elif plot_start_delta != plot_end_delta:
        raise ValueError('You cannot specify one and not the other!')

    return times, hrs_since_storm_start


# In[ ]:


def UT_from_Storm_onset(itime):
    """input a datetime

    returns the UT as HH:MM from storm onset, as a string"""
    l = (pd.Timestamp(itime) - dtime_storm_start) / pd.Timedelta(
        '1 minute')  # get pd datetime of this iter, find minute diff from storm start
    if l > 0:
        hrs = np.floor(l/60)
        hrs = str(int(hrs)).rjust(2, '0')
    else:
        hrs = np.ceil(l/60)
        hrs = '-' + str(np.abs(int(hrs))).rjust(2, '0')
    mins = str(int(np.abs(l) % 60)).rjust(2, '0')
    ut = hrs + ':' + mins
    return ut


# In[ ]:


def get_sami_grid(sami_data_path=sami_data_path, geo_grid_files=geo_grid_files):

    grid = {}

    for f in geo_grid_files:
        file = open(os.path.join(sami_data_path, geo_grid_files[f]), 'rb')
        raw = np.fromfile(file, dtype='float32')[1:-1].copy()
        file.close()

        grid[f] = raw.reshape(nlt, nf, nz).copy()
    return grid


# In[ ]:


def read_sami_data(cols, nts):
    """
    Read in sami data for the specified columns and return sama data dict

    inputs:
    -------
    cols: list-like
        - Columns you want data for. Does not have to be everything

    nts: int OR tuple/list
        - either nt (number of times) if you want all sami data from simulation or:
        - nts (start_time, end_time) if you want plots from a select time period.

    """
    sami_data = {}

    # handle cut time list and full time list
    if type(nts) != int:
        t_start = nts[0]
        t_end = nts[1]
        ntimes = t_end - t_start
    else:
        t_start = 0
        t_end = nt
        ntimes = nt

    pbar = tqdm(total=len(cols) * ntimes, desc='reading SAMI data')

    for f in cols:

        sami_data[f] = np.zeros((nlt, nf, nz, ntimes))

        file = open(os.path.join(sami_data_path, data_files[f]), 'rb')
        for t in range(t_end):
            raw = np.fromfile(file, dtype='float32', count=(nz*nf*nlt)+2)[1:-1]
            if t >= t_start:
                sami_data[f][:, :, :, t -
                             t_start] = raw.reshape(nlt, nf, nz).copy()
                pbar.update(1)
        file.close()
    pbar.close()

    return sami_data


# In[ ]:


def make_filter(params=None):
    # Define the cutoff frequencies
    lowcut = 1/(100/60)  # 100 minutes in units of sample^-1
    highcut = 1/(30/60)  # 30 minutes in units of sample^-1

    # Define the Butterworth filter
    nyquist = 0.5 * 5  # 5 minutes is the sampling frequency
    low = lowcut / nyquist
    high = highcut / nyquist
    sos = signal.butter(2, [low, high], btype='bandstop', output='sos')
    return sos


def remove_background(time_series, sos):

    # Apply the filter to the time series
    filtered_data = signal.sosfiltfilt(sos, time_series)

    return filtered_data


# In[ ]:


def geo_to_cart(lats_, lons_, alts_, ntime):
    """
    get cartesian from a grid slice
    """

    coord_arr = []

    if type(alts_) != list:
        if type(alts_) != np.array:
            if type(alts_) != np.ndarray:
                alts_ = [alts_]

    for lat in lats_:
        for lon in lons_:
            for alt in alts_:
                coord_arr.append([(alt + 6371)/6371, lat, lon])

    dtime = times[ntime]

    coords = Coords(coord_arr, 'GEO', 'sph')
    coords.ticks = Ticktock([dtime for k in range(len(coord_arr))])

    newcoords = coords.convert('GEO', 'car')

    return newcoords


# In[ ]:


def grid_to_cartesian(grid, ntime):

    coords0 = list(zip(grid['malt'].flatten() / 6371,
                   grid['mlat'].flatten(), grid['mlon'].flatten()))
    dtime = times[ntime]

    coords = Coords(coords0, 'CDMAG', 'sph')
    coords.ticks = Ticktock([dtime for k in range(len(coords0))])

    newcoords = coords.convert('GEO', 'car')

    return newcoords


# ## Read in data


nz, nf, nlt, nt = get_grid_elems_from_parammod(sami_data_path)


# In[ ]:


times, hrs, times_list, nts = make_times(
    dtime_sim_start, nt, plot_start_delta, plot_end_delta)
new_nt = np.diff(nts)[0]


# In[ ]:


# OR.... to run for all times:
# times, hrs, times_list = make_times(dtime_sim_start, nt)


# In[ ]:


print('nlt, nf, nz, nts, new_nt = \n', nlt, nf, nz, nts, new_nt)
nt = new_nt


# In[ ]:


grid = get_sami_grid()


# In[ ]:


for g in grid.keys():
    print(g, grid[g].shape)


# In[ ]:


# sami grid
sami_data = read_sami_data(cols, nts)
# or
# sami_data = read_sami_data(cols, nt)

print('sami data shape: ', sami_data[cols[0]].shape)


# In[ ]:


# if to_filter:
#     print('Calculating fits. This will take a moment...')
#     fits_sami = make_fits(sami_data)


# ## Some setup

# In[ ]:


mlons = np.unique(grid['mlon'].round(2))


# ## I think this is the best way to do it... Need to make sure we can resolve LSTIDS though.

# In[ ]:


norm_alts = (grid['alt'].flatten() < (max(out_alts) + 300)
             ) & (grid['alt'].flatten() > (min(out_alts)-75))


# In[ ]:


pbar = tqdm(total=len(times),
            desc='making preds... pbar is a very rough estimate! ')

# this will be very messy. clean up after the processing is done.


def interp_grid(col):
    preds = {}
    preds[col] = np.zeros(
        [len(times), len(out_lats), len(out_lons),  len(out_alts)])

    for ntime, dt in enumerate(times):

        # print('cart grid')
        grid_cart = grid_to_cartesian(grid, ntime)

        xs = grid_cart.data[:, 0][norm_alts]
        ys = grid_cart.data[:, 1][norm_alts]
        zs = grid_cart.data[:, 2][norm_alts]

        # print('zipping')

        loc_grid = list(zip(xs, ys, zs))

        # print('making output grid')

        out_grid = geo_to_cart(out_lats, out_lons, out_alts, ntime)

        out_xs = out_grid.data[:, 0]
        out_ys = out_grid.data[:, 1]
        out_zs = out_grid.data[:, 2]

        datas = sami_data[col][:, :, :, ntime].flatten()[norm_alts]

        # print('building interpolator')

        interp = LinearNDInterpolator(loc_grid, datas, rescale=True)
        # print('interpolating')
        pred = interp(list(zip(out_xs, out_ys, out_zs)))

        # break

        preds[col][ntime] = pred.reshape(
            [len(out_lats), len(out_lons), len(out_alts)])

        pbar.update(1)
    return preds


# In[ ]:


# Thread the above function, if thread option is set.
if thread:
    with Pool(num_workers) as pool:

        pred_inter = pool.map(interp_grid, cols)


else:
    pred_inter = []
    for col in cols:
        pred_inter.append(interp_grid[col])


preds_cleaned = {}
for p in pred_inter:
    for k in p.keys():
        preds_cleaned[k] = p[k]
preds = preds_cleaned
del preds_cleaned


# In[ ]:


if debug:
    a = 5

    # COmpare maps
    plt.imshow(preds[cols[0]][0, :, :, a], origin='lower', extent=[
               min(out_lons), max(out_lons), min(out_lats), max(out_lats)], aspect='auto')

    plt.show()
    plt.close()

    alt_mask_plot = (np.abs(grid['alt'].flatten() - out_alts[a]) < 10)
    plt.scatter(grid['glon'].flatten()[alt_mask_plot], grid['glat'].flatten()[
                alt_mask_plot], c=sami_data['edens'][:, :, :, ntime].flatten()[alt_mask_plot])
    plt.ylim(min(out_lats), max(out_lats))

    plt.show()
    plt.close()


# In[ ]:


if to_filter:
    sos = make_filter()
    filtered = {}
    for col in cols:
        filtered[col] = signal.sosfiltfilt(sos, preds[col], axis=0)

    if debug:
        plt.imshow(preds[cols[0]][0, :, :, a], origin='lower', extent=[
                   min(out_lons), max(out_lons), min(out_lats), max(out_lats)], aspect='auto')

        plt.show()
        plt.close()

        plt.imshow(100*(preds[cols[0]][25, :, :, a] - filtered[cols[0]][25, :, :, a])/preds[cols[0]][1, :, :, a], origin='lower',
                   extent=[min(out_lons), max(out_lons), min(out_lats), max(out_lats)], aspect='auto', vmin=-5, vmax=5)
        plt.colorbar()
        plt.show()
        plt.close()

        for l in range(0, len(out_lons), 4):
            plt.figure(figsize=(8, 4))
            plt.imshow((100*(preds[cols[0]][:, :, l, a] - filtered[cols[0]][:, :, l, a])/preds[cols[0]][:, :, l, a]).T,
                       extent=[min(hrs), max(hrs), min(out_lats), max(out_lats)], aspect='auto', vmin=-4, vmax=4)
            plt.colorbar(label='edends % over background at ' +
                         str(out_alts[a]) + ' km')
            plt.title('glon = ' + str(out_lons[l].round(2)))
            plt.show()
            plt.close()


# In[ ]:


# ## Save everything to files.

# In[ ]:


# In[ ]:
print('writing files...')

np.array(times).tofile(data_out_path + 'times', format='%s', sep=',')
out_lats.tofile(data_out_path + 'out-lats', sep=',')
out_lons.tofile(data_out_path + 'out-lons', sep=',')
out_alts.tofile(data_out_path + 'out-alts', sep=',')

np.array(preds[cols[0]].shape).tofile(sami_data_path + 'out-shape', sep=',')

for col in cols:

    preds[col].tofile(data_out_path + 'preds-' + col, sep=',')

    if to_filter:
        filtered[col].file(data_out_path + 'bp_filtered-' + col, sep=',')

if save_raw:
    for col in cols:

        sami_data[col].tofile(data_out_path + 'raw-' + col, sep=',')

    np.array(sami_data[cols[0]].shape).tofile(
        sami_data_path + 'raw-shape', sep=',')

# In[ ]:


print('written! Exiting. ')

# In[ ]:
