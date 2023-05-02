# Postprocess SAMI data to be read into xarrays.

# Data will be interpolated to a standard grographic grid when possible.


import xarray as xr
import pandas as pd

from glob import glob
from tqdm import tqdm
import numpy as np

from datetime import datetime
from scipy.spatial import KDTree

from utility_programs.read_routines import SAMI
import argparse

data2, times = SAMI.read_to_nparray(sami_data_path, dtime_sim_start,
                                    dtime_event_start, 1, 2, cols='edens')


def latlonalt_to_cart(lat, lon, radius):
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)
    x = radius * np.cos(lat) * np.cos(lon)
    y = radius * np.cos(lat) * np.sin(lon)
    z = radius * np.sin(lat)
    return np.array([x, y, z])


def make_weights(latout,
                 lonout,
                 altout,
                 sami_data_path,
                 dtime_sim_start):

    nz, nf, nlt, nt = SAMI.get_grid_elems_from_parammod(sami_data_path)

    old_shape = [nlt, nf, nz]

    grid2 = SAMI.get_sami_grid(sami_data_path, nlt, nf, nz)
    grid = {}

    for k in grid2.keys():
        grid[k] = grid2[k].flatten()
        grid2[k] = grid2[k]
        print(k, grid[k].shape, grid2[k].shape)
    in_cart = latlonalt_to_cart(grid['glat'], grid['glon'], grid['malt'])

    out_lats = []
    out_lons = []
    out_alts = []
    for a in latout:
        for o in lonout:
            for l in altout:
                out_lats.append(a)
                out_lons.append(o)
                out_alts.append(l)
    out_cart = latlonalt_to_cart(out_lats, out_lons, np.array(out_alts)+6371)

    tree = KDTree(in_cart.T)
    dists, nearest = tree.query(out_cart.T)

    zero1dsrc = np.empty([len(out_cart[0]), 8])
    weights = {'weight': zero1dsrc.copy(),
               'srcidxs1d': zero1dsrc.copy(), }

    for n, i in enumerate(tqdm(nearest)):
        l, f, z = np.unravel_index(i, old_shape)

        if l == old_shape[0]-1:
            l2 = 0
        else:
            l2 = l+1

        if z == 0:
            z2 = -1
        else:
            z2 = z-1

        f2 = f+1
        if f == old_shape[1]-1:
            continue

        cs = [[l,  f,  z],
              [l,  f,  z2],
              [l2, f,  z2],
              [l2, f,  z],
              [l,  f2, z],
              [l,  f2, z],
              [l2, f2, z2],
              [l2, f2, z2], ]

        id_pt = []
        xs = []
        ys = []
        zs = []
        for c in cs:
            try:
                index = np.ravel_multi_index(c, old_shape)
            except:
                break
            id_pt.append(index)

            xs.append(in_cart[0, index])
            ys.append(in_cart[1, index])
            zs.append(in_cart[2, index])

        d = np.sqrt((xs-out_cart[0, n])**2 +
                    (ys-out_cart[1, n])**2 + (zs-out_cart[2, n])**2)
        dtot = np.sum(d)

        try:
            weights['weight'][n] = 1/(d)
        except ValueError:
            weights['weight'][n] = 0
        weights['srcidxs1d'][n] = id_pt
    return weights


def apply_weights(weights,
                  data_dict,
                  times,
                  altout,
                  latout,
                  lonout,
                  **coarsen_args=None):

    ds = xr.Dataset(coords={
        'time': (['time'], times),
        'alt': (['alt'], altout),
        'lat': (['lat'], latout),
        'lon': (['lon'], lonout)},)

    outv = []
    for var in data_dict['data'].keys():
        t0 = data2['data'][var]
        for t in range(t0.shape[-1]):
            outv.append(
                np.sum(weights['weight'] * np.take(t0[:, :, :, t],
                                                   weights['srcidxs1d'].astype(int)), axis=1)
                / np.sum(weights['weight'], axis=1))

        ds[var] = np.array(outv).reshape(
            len(times), len(latout), len(lonout), len(altout))

    if **coarsen_args is not None:
        ds = ds.coarsen(**coarsen_args).mean()

    return ds


def main(args):

    weights = make_weights(latout,
                           lonout,
                           altout,
                           sami_data_path,
                           dtime_sim_start)

    if save_weights:
        # TODO: SAVE WEIGHTS

    ds = apply_weights(weights,
                       data_dict,
                       times,
                       altout,
                       latout,
                       lonout,)
    ds.to_cdf(os.path.join(out_path, 'sami-regridded'))
    
    return