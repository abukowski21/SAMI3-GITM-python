# Postprocess SAMI data to be read into xarrays.

# Data will be interpolated to a standard grographic grid when possible.


import xarray as xr

from tqdm.auto import tqdm
import numpy as np

from scipy.spatial import KDTree
import os
from utility_programs.read_routines import SAMI
from utility_programs.utils import str_to_ut
import argparse


def latlonalt_to_cart(lat, lon, radius):
    """Convert lat, lon, alt to cartesian coordinates.

    Args:
        lat (numpy array): latitude (in degrees)
        lon (numpy array): longitude (in degrees)
        radius (numpy array): radius in Re from center of Earth

    Returns:
        numpy array: 3xN array of cartesian coordinates
    """
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)
    x = radius * np.cos(lat) * np.cos(lon)
    y = radius * np.cos(lat) * np.sin(lon)
    z = radius * np.sin(lat)
    return np.array([x, y, z])


def is_inside_cube(x, y, z, xs, ys, zs):
    """Check if a point is inside a cube.

    Args:
        x (float): x coordinate of point
        y (float): y coordinate of point
        z (float): z coordinate of point
        xs (list-like): cartesian x coordinates of cube vertices
        ys (list-like): cartesian y coordinates of cube vertices
        zs (list-like): cartesian z coordinates of cube vertices

    Returns:
        bool: True if point is inside cube, False otherwise
    """
    if (x > min(xs)) and (x < max(xs)) and (y > min(ys)) and (y < max(ys)) \
            and (z > min(zs)) and (z < max(zs)):
        return True
    else:
        return False


def generate_interior_points(in_cart, old_shape):
    """Generates interior points for each cube in the old grid.

    Args:
        in_cart (numpy array): Cartesian coordinates of the old grid
        old_shape (list): (nlt, nf, nz) of original sami run.

    Returns:
        centers (numpy array): cartesian coordinates of the center of each
            cube
        coords (list): list of lists of the coordinates of each corner
            of each cube.
    """

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
                            break
                        id_pt.append(index)

                        xs.append(in_cart[0, index])
                        ys.append(in_cart[1, index])
                        zs.append(in_cart[2, index])

                center = np.sum(xs)/8, np.sum(ys)/8, np.sum(zs)/8

                centers.append(center)
                coords.append(cs)
                pbar.update()
    pbar.close()

    print('From %i grid points we generated %i cubes'
          % (len(in_cart[0]), len(centers)))

    return centers, coords


def make_weights(in_cart, out_cart, nearest, old_shape, coords):
    """Generates weights for each point in the new grid.

    Args:
        in_cart (list-like): Cartesian coordinates of the old grid
        out_cart (list-like): Cartesian coordinates of the new grid
        nearest (list-like): list of nearest cube center in the old grid to
            each point in the new grid
        old_shape (list-like): (nlt, nf, nz) of original sami run.
        coords (list-like): coordinates of each corner of each cube in the old
            grid.

    Returns:
        weights (numpy array): Weights to multiply the old grid by to get
            the new grid
        src_idxs (numpy array): Indices of the old grid that contribute to
            each point in the new grid.
    """
    weights = np.zeros([len(out_cart[0]), 8])
    src_idxs = np.zeros([len(out_cart[0]), 8])
    print('Calculating weights...')
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
    print('Done, found %i valid points' % np.sum(weights > 0))

    return weights, src_idxs


def find_pairs(centers, out_cart):
    """Find nearest point between two sets of coordinates.
            (only works well for cartesian coordinates)

    Args:
        centers (numpy array): Centers of cubes in the old grid
        out_cart (numpy array): Grid to search from.

    Returns:
        numpy array: indices of the nearest point in the old grid to each
            point in the new grid (out_cart).
    """
    tree = KDTree(centers)
    _, nearest = tree.query(out_cart.T)
    return nearest


def do_apply_weights(weights,
                     src_idxs,
                     data_dict,
                     times,
                     altout,
                     latout,
                     lonout):
    """Apply weights generated by make_weights to the old grid to get the new
        grid.

    Args:
        weights (numpy array): Weights to multiply the old grid by to get
            Values of the new grid.
        src_idxs (numpy array): Indices of the old grid that contribute to
            each point in the new grid.
        data_dict (dict): Dictionary of data from the old grid.
        times (numpy array): List of times.
        altout (numpy array): Altitudes out.
        latout (numpy array): Latitudes out.
        lonout (numpy array): Longitudes out.

    Returns:
        xarray Dataset: Dataset of the new grid.
    """

    ds = xr.Dataset(coords={
        'time': (['time'], times),
        'alt': (['alt'], altout),
        'lat': (['lat'], latout),
        'lon': (['lon'], lonout)},)

    outv = []
    for var in data_dict['data'].keys():
        t0 = data_dict['data'][var]
        for t in range(t0.shape[-1]):
            outv.append(
                np.sum(weights * np.take(
                    t0[:, :, :, t], src_idxs.astype(int)), axis=1)
                / np.sum(weights, axis=1))

        ds[var] = np.array(outv).reshape(
            len(times), len(latout), len(lonout), len(altout))

    return ds


def main(
        sami_data_path,
        out_path=None,
        save_weights=True,
        use_saved_weights=True,
        apply_weights=True,
        cols='all',
        dtime_sim_start=None,
        lat_step=1,
        lon_step=4,
        alt_step=50,
        minmax_alt=[100, 2200],
        lat_finerinterps=1,
        lon_finerinterps=1,
        alt_finerinterps=1,
        split_by_var=True,
        single_file=False):

    if out_path is None:
        out_path = sami_data_path

    if not os.path.exists(out_path):
        print('making outpath %s' % out_path)
        os.makedirs(out_path)

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

    if apply_weights:

        print('reading SAMI data...')
        data, times = SAMI.read_to_nparray(
            sami_data_path, dtime_sim_start, cols=cols, pbar=True)

        ds = do_apply_weights(weights, idxs,
                              data, times,
                              altout, latout, lonout)

        if 1 not in [lat_finerinterps, lon_finerinterps, alt_finerinterps]:

            if lat_finerinterps > 1:
                ds = ds.coarsen(lat=lat_finerinterps, boundary='trim').mean()
            if lon_finerinterps > 1:
                ds = ds.coarsen(lon=lon_finerinterps, boundary='trim').mean()
            if alt_finerinterps > 1:
                ds = ds.coarsen(alt=alt_finerinterps, boundary='trim').mean()

        if split_by_var:
            for varname in ds.data_vars:
                ds[varname].to_cdf(
                    os.path.join(out_path, '%s' % varname))
        elif single_file:
            ds.to_cdf(os.path.join(out_path, 'sami-regridded'))

    return


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('sami_data_path', type=str, help='path to sami data')

    parser.add_argument('-o', '--out_path', type=str, default=None,
                        help='path to save regridded data,'
                        'deafults to sami_data_path')
    parser.add_argument('-d', '--dtime_sim_start', type=str, default=None,
                        help='datetime string of simulation start time.'
                        'Required to read raw SAMI files')
    parser.add_argument('--save_weights', action='store_true',
                        help='save weights and indexes')
    parser.add_argument('--reuse_weights', action='store_true',
                        help='use saved weights and indexes'
                        '(helpful if you have already run on other columns)')
    parser.add_argument('-a', '--apply_weights', action='store_true',
                        help='apply weights to data. Generates new data file')
    parser.add_argument('--cols', type=str, default='all', nargs='*',
                        help='columns to read from sami data. Defaults to all'
                        'input as a list with spaces between.')
    parser.add_argument('--split_by_var', action='store_false', default=True,
                        help='Split output files by variable? Default: True')
    parser.add_argument('--single_file', action='store_true', default=False,
                        help='Write all output datasets to single file?'
                        'Default: False')
    parser.add_argument('--custom_grid', action='store_true',
                        help='Launches interactive script to set custom grid.'
                        'Default grid is 4 x 1 deg lonxlat, 50 alts.'
                        'minimum alt is 100, max is 2200, global lon x lat')
    parser.add_argument('--lat_step', type=float, default=1,
                        help='latitude step size in degrees')
    parser.add_argument('--lon_step', type=float, default=4,
                        help='longitude step size in degrees')
    parser.add_argument('--alt_step', type=float, default=50,
                        help='altitude step size in km')
    parser.add_argument('--lat_finerinterps', type=int, default=1,
                        help='Interpolate at a finer resolution in latitude?'
                        'defaults to 1 (no finer), set to any value > 1 to'
                        'interpolate at high resolution and then'
                        'coarsen when writing the files out.')
    parser.add_argument('--lon_finerinterps', type=int, default=1,
                        help='Interpolate at a finer resolution in longitude?'
                        'defaults to 1 (no finer), set to any value > 1 to'
                        'interpolate at high resolution and then'
                        'coarsen when writing the files out.')
    parser.add_argument('--alt_finerinterps', type=int, default=1,
                        help='Interpolate at a finer resolution in altitude?'
                        'defaults to 1 (no finer), set to any value > 1 to'
                        'interpolate at high resolution and then'
                        'coarsen when writing the files out.')
    parser.add_argument('--minalt', type=float, default=100,
                        help='Minimum altitude in km')
    parser.add_argument('--maxalt', type=float, default=2200,
                        help='Maximum altitude in km')

    args = parser.parse_args()

    if args.custom_grid:
        print('We are going to set a custom grid for you. ')
        print('Press enter to accept the default values in parentheses')
        latstep = input('latitude step size in degrees (1):', default=1)
        lonstep = input('longitude step size in degrees: (4):', default=4)
        altstep = input('altitude step size in km (50):', default=50)
        minalt = input('minimum altitude in km (100):', default=100)
        maxalt = input('maximum altitude in km (2200):', default=2200)
        print('Now for the options to interpolate at a finer resolution'
              ' and then coarsen afterwards. If you dont know what this'
              ' means you can run with 1s and it will be faster. if you'
              ' see weird artifacts in your outputs you can try '
              ' adjusting this. Number given multiplies the step size')
        latfiner = input('interpolate a finer resolution in latitude? (1):',
                         default=1)
        lonfiner = input('interpolate a finer resolution in longitude? (1):',
                         default=1)
        altfiner = input('interpolate a finer resolution in altitude? (1):',
                         default=1)
    else:
        latstep = args.lat_step
        lonstep = args.lon_step
        altstep = args.alt_step
        latfiner = args.lat_finerinterps
        lonfiner = args.lon_finerinterps
        altfiner = args.alt_finerinterps
        minalt = args.minalt
        maxalt = args.maxalt

    if args.dtime_sim_start is None:
        if args.apply_weights:
            raise ValueError('You must specify a simulation start time'
                             'when using apply_weights to read data')
    else:
        dtime_sim_start = str_to_ut(args.dtime_sim_start)

    if args.reuse_weights:
        # check if out_weight file exists
        if not os.path.isfile(os.path.join(args.out_path, 'weights')):
            raise ValueError('You specified to reuse weights, but no weights'
                             'file was found in the out_path')

    if args.split_by_var and args.single_file:
        raise ValueError(
            'You cannot specify both split_by_var and single_file')

    main(args.sami_data_path,
         out_path=args.out_path,
         save_weights=args.save_weights,
         use_saved_weights=args.reuse_weights,
         apply_weights=args.apply_weights,
         cols=args.cols,
         dtime_sim_start=dtime_sim_start,
         lat_step=latstep,
         lon_step=lonstep,
         alt_step=altstep,
         minmax_alt=[minalt, maxalt],
         lat_finerinterps=latfiner,
         lon_finerinterps=lonfiner,
         alt_finerinterps=altfiner,
         split_by_var=args.split_by_var,
         single_file=args.single_file,)
