# Postprocess SAMI data to be read into xarrays.

# Data will be interpolated to a standard grographic grid when possible.


import xarray as xr

from tqdm import tqdm
import numpy as np

from scipy.spatial import KDTree
import os
from utility_programs.read_routines import SAMI
import argparse


def latlonalt_to_cart(lat, lon, radius):
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)
    x = radius * np.cos(lat) * np.cos(lon)
    y = radius * np.cos(lat) * np.sin(lon)
    z = radius * np.sin(lat)
    return np.array([x, y, z])


def is_inside_cube(x, y, z, xs, ys, zs):
    if (x > min(xs)) and (x < max(xs)) and (y > min(ys)) and (y < max(ys)) \
            and (z > min(zs)) and (z < max(zs)):
        return True
    else:
        return False


def generate_interior_points(in_cart, old_shape):

    nlt, nf, nz = old_shape
    centers = []
    coords = []
    pbar = tqdm(total=np.product(old_shape), desc='Generating interior points')
    badbadbad = []

    for lt in range(nlt):
        for f in range(nf):
            for z in range(nz):

                if lt == old_shape[0]-1:
                    l2 = 0
                else:
                    l2 = lt+1

                if z == 0:
                    z2 = -1
                else:
                    z2 = z-1

                f2 = f+1
                if f == old_shape[1]-1:
                    badbadbad.append([lt, f, z])
                    pbar.update()
                    continue

                cs = [[lt,  f,  z],
                      [lt,  f,  z2],
                      [l2, f,  z2],
                      [l2, f,  z],
                      [lt,  f2, z],
                      [lt,  f2, z],
                      [l2, f2, z2],
                      [l2, f2, z2]]

                for c in cs:
                    id_pt = []
                    xs = []
                    ys = []
                    zs = []
                    for c in cs:
                        try:
                            index = np.ravel_multi_index(c, old_shape)
                        except ValueError:
                            pbar.update()
                            break
                        id_pt.append(index)

                        xs.append(in_cart[0, index])
                        ys.append(in_cart[1, index])
                        zs.append(in_cart[2, index])

                center = np.sum(xs)/8, np.sum(ys)/8, np.sum(zs)/8

                centers.append(center)
                coords.append(cs)
                pbar.update()
    print('From %i grid points we generated %i cubes'
          % (len(in_cart[0]), len(centers)))

    return centers, coords


def make_weights(in_cart, out_cart, nearest, old_shape, coords):

    weights = np.zeros([len(out_cart[0]), 8])
    src_idxs = np.zeros([len(out_cart[0]), 8])
    for n, pt in enumerate(tqdm(nearest)):

        idxs = [np.ravel_multi_index(c, old_shape) for c in coords[pt]]

        xs = in_cart[0, idxs]
        ys = in_cart[0, idxs]
        zs = in_cart[2, idxs]

        val = is_inside_cube(out_cart[0, n], out_cart[1, n], out_cart[2, n],
                             xs, ys, zs)
        if val:
            d = np.sqrt((xs-out_cart[0, n])**2 + (ys -
                        out_cart[1, n])**2 + (zs-out_cart[2, n])**2)

            weights[n] = 1/(d)

            src_idxs[n] = idxs


def find_pairs(centers, out_cart):
    tree = KDTree(centers)
    _, nearest = tree.query(out_cart.T)
    return nearest


def apply_weights(weights,
                  data_dict,
                  times,
                  altout,
                  latout,
                  lonout):

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

    return ds


def main(
        sami_data_path,
        out_path=None,
        save_weights=False,
        use_saved_weights=True,
        cols='all',
        lat_step=1,
        lon_step=4,
        alt_step=50,
        minmax_alt=[100, 2200],
        lat_finerinterps=1,
        lon_finerinterps=1,
        alt_finerinterps=1):

    if out_path is None:
        out_path = sami_data_path

    # Read in the data
    nz, nf, nlt, nt = SAMI.get_grid_elems_from_parammod(sami_data_path)
    old_shape = [nlt, nf, nz]

    grid2 = SAMI.get_sami_grid(sami_data_path, nlt, nf, nz)

    grid = {}
    for k in grid2.keys():
        grid[k] = grid2[k].flatten()
        grid2[k] = grid2[k]

    in_cart = latlonalt_to_cart(grid['glat'], grid['glon'], grid['malt'])

    latout = np.arange(-90, 90, lat_step/lat_finerinterps)
    lonout = np.arange(0, 360, lon_step/lon_finerinterps)
    altout = np.arange(200, 2200, alt_step/alt_finerinterps)

    out_lats = []
    out_lons = []
    out_alts = []

    for a in latout:
        for o in lonout:
            for l1 in altout:
                out_lats.append(a)
                out_lons.append(o)
                out_alts.append(l1)

    out_cart = latlonalt_to_cart(out_lats, out_lons, np.array(out_alts)+6371)

    centers, coords = generate_interior_points(in_cart, old_shape)

    nearest = find_pairs(centers, out_cart)

    weights, idxs = make_weights(in_cart, out_cart, nearest, old_shape, coords)

    if save_weights:
        # TODO: SAVE WEIGHTS
        weights.tofile(os.path.join(out_path, 'weights'))
        idxs.tofile(os.path.join(out_path, 'indexes'))

    ds = apply_weights(weights,
                       data_dict,
                       times,
                       altout,
                       latout,
                       lonout,)
    ds.to_cdf(os.path.join(out_path, 'sami-regridded'))

    return
