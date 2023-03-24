"""This will make GITM polar dial plots and a map of the data.

User selects which GITM file(s) to plot, and the program will make a
nice and pretty graph.  The user can also select which variable to
put in which plot.

created mar 15 2023 by aaron
"""
# %%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from utility_programs.read_routines import GITM
from datetime import datetime

from utility_programs import plotting_routines as pr

# %%
import geopandas
world = geopandas.read_file(geopandas.datasets.get_path("naturalearth_lowres"))


# %%
dials_var = 'AltIntJouleHeating(W/m2)'
map_var = 'VerticalTEC'

# %%

times, gitm_grid, gitm_f0, gitm_vars = GITM.read_gitm_into_nparrays(
    '/petastore/phil/GITM/cheyenne_runs/FullAmp/data',
    datetime(2011, 5, 21, 12), gitm_file_pattern='2DANC*.bin',
    t_start_idx=1, t_end_idx=1, return_vars=True)


# %
f_now = gitm_f0[-1].copy()
lats = np.unique(gitm_grid['latitude'])
lons = np.unique(gitm_grid['longitude'])
alts = np.unique(gitm_grid['altitude'])

# %%

data = f_now[gitm_vars.index(map_var)]
newshape = len(lons), len(lats)
data = data.reshape(newshape)
gitm_grid['latitude'] = gitm_grid['latitude'].reshape(newshape)
gitm_grid['longitude'] = gitm_grid['longitude'].reshape(newshape)

minI = np.min(data)
max_I = np.max(data)

# %%


def mapping(data_arr, lats, lons, cbar_label=None, title=None,
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

    ax = ax or plt.gca()
    world.plot(ax=ax, color="white", edgecolor="black", zorder=1)
    data = ax.imshow(
        data_arr.T,
        cmap="viridis",
        aspect="auto",
        extent=[min(lons), max(lons), min(lats), max(lats)],
        **kwargs)
    plt.colorbar(data, ax=ax, label=map_var)

    return ax


# %%
fig = plt.figure(figsize=(10, 8.5), layout='tight')

gs1 = GridSpec(nrows=2, ncols=2, wspace=.1, hspace=.1)
# gs = GridSpec(nrows=2, ncols=2, wspace=0.0, left=0.0, right=0.9)

ax0 = fig.add_subplot(gs1[0, 0], projection='polar', title='North')
ax1 = fig.add_subplot(gs1[0, 1], projection='polar', title='South')
ax2 = fig.add_subplot(gs1[1, :2])


mapping(data_arr=data, lats=lats, lons=lons, cbar_label='TEC',
        ax=ax2, zorder=1, origin='lower', alpha=0.8, vmin=minI, vmax=max_I,
        title='Vertically Integrated TEC')

maskNorth = ((lats > 45) & (lats < 90.0))
maskSouth = ((lats > -90.0) & (lats < -45))


# north
r, theta = np.meshgrid(90-lats[maskNorth], lons)
ax0.pcolor(np.deg2rad(theta), r, data[:, maskNorth], vmin=minI, vmax=max_I)
xlabels = ['', '12', '18', '00']
ylabels = ['80', '70', '60', '50']
ax0.set_xticklabels(xlabels)
ax0.set_yticklabels(ylabels)
pi = 3.14159
ax0.set_xticks(np.arange(0, 2*pi, pi/2))
ax0.set_yticks(np.arange(10, 50, 10))
ax0.set_title('North TEC')

# South
r, theta = np.meshgrid(lats[maskSouth], lons)
ax1.pcolor(np.deg2rad(theta), r, data[:, maskSouth], vmin=minI, vmax=max_I)
# xlabels = ['', '12', '18', '00']
ylabels = ['-80', '-70', '-60', '-50']
# ax1.set_xticklabels(xlabels)
ax1.set_yticklabels(ylabels)
ax1.set_title('South TEC')

plt.savefig('test')

# %


# # %%
# ntime = 4

# f_now = gitm_f0[ntime].copy()
# lats = np.unique(gitm_grid['latitude'])
# lons = np.unique(gitm_grid['longitude'])
# alts = np.unique(gitm_grid['altitude'])

# data = f_now[gitm_vars.index(map_var)]
# newshape = len(lons), len(lats)
# data = data.reshape(newshape)
# gitm_grid['latitude'] = gitm_grid['latitude'].reshape(newshape)
# gitm_grid['longitude'] = gitm_grid['longitude'].reshape(newshape)

# minI = np.min(data)
# max_I = np.max(data)

# fig = plt.figure(figsize=(10, 8.5), layout='tight')

# gs1 = GridSpec(nrows=2, ncols=2, wspace=.1, hspace=.2)
# # gs = GridSpec(nrows=2, ncols=2, wspace=0.0, left=0.0, right=0.9)

# ax0 = fig.add_subplot(gs1[0, 0], projection='polar', title='North')
# ax1 = fig.add_subplot(gs1[0, 1], projection='polar', title='South')
# ax2 = fig.add_subplot(gs1[1, :2])


# mapping(data_arr=data, lats=lats, lons=lons, cbar_label='TEC',
#         ax=ax2, zorder=1, origin='lower', alpha=0.8, vmin=minI, vmax=max_I,
#         title='Vertically Integrated TEC')

# maskNorth = ((lats > 45) & (lats < 90.0))
# maskSouth = ((lats > -90.0) & (lats < -45))

# minP = np.min(data[:, np.ma.mask_or(maskSouth, maskNorth)])
# maxP = np.max(data[:, np.ma.mask_or(maskSouth, maskNorth)])

# # north
# r, theta = np.meshgrid(90-lats[maskNorth], lons)
# ax0.pcolor(np.deg2rad(theta), r, data[:, maskNorth], vmin=minP, vmax=maxP)
# # xlabels = ['', '12', '18', '00']
# ylabels = ['80', '70', '60', '50']
# # ax0.set_xticklabels(xlabels)
# ax0.set_yticklabels(ylabels)
# pi = 3.14159
# ax0.set_xticks(np.arange(0, 2*pi, pi/2))
# ax0.set_yticks(np.arange(10, 50, 10))
# ax0.set_title('North TEC')

# # South
# r, theta = np.meshgrid(lats[maskSouth], lons)
# cb = ax1.pcolor(np.deg2rad(theta), r, data[:, maskSouth], vmin=minP,
# vmax=maxP)
# # xlabels = ['', '12', '18', '00']
# ylabels = ['-80', '-70', '-60', '-50']
# # ax1.set_xticklabels(xlabels)
# ax1.set_yticklabels(ylabels)
# ax1.set_title('South TEC')
# fig.colorbar(cb)


# fig.suptitle('From GITM run on %s' % (str(times[ntime])),
#              fontsize=17)

# plt.savefig('test')


# %
