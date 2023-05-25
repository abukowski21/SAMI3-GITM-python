"""
Postprocess SAMI data to be read into xarrays.

- Data will be interpolated to a standard grographic grid when possible.
- SAMI model results can be regridded to geographic grid or written
	with the same indexing as raw files.
- Can specify to write one file for the whole run or to split_by_(var/time)

"""


import xarray as xr

from tqdm.auto import tqdm
import numpy as np

from scipy.spatial import KDTree
import os
from utility_programs.read_routines import SAMI
from utility_programs.utils import str_to_ut, make_ccmc_name
import argparse


def latlonalt_to_cart(lat, lon, radius):
    """Convert lat, lon, alt to cartesian coordinates.

    Args:
        lat (numpy.ndarray): latitude (in degrees)
        lon (numpy.ndarray): longitude (in degrees)
        radius (numpy.ndarray): radius in Re from center of Earth

    Returns:
        numpy.ndarray: 3xN array of cartesian coordinates
    """
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)
    x = radius * np.cos(lat) * np.cos(lon)
    y = radius * np.cos(lat) * np.sin(lon)
    z = radius * np.sin(lat)
    return np.array([x, y, z])


# def is_inside_cube(x, y, z, xs, ys, zs):
#     """Check if a point is inside a cube.

#     Args:
#         x (float): x coordinate of point
#         y (float): y coordinate of point
#         z (float): z coordinate of point
#         xs (list-like): cartesian x coordinates of cube vertices
#         ys (list-like): cartesian y coordinates of cube vertices
#         zs (list-like): cartesian z coordinates of cube vertices

#     Returns:
#         bool: True if point is inside cube, False otherwise
#     """
#     if (x > min(xs)) and (x < max(xs)) and (y > min(ys)) and (y < max(ys)) \
#             and (z > min(zs)) and (z < max(zs)):
#         return True
#     else:
#         return False


def generate_interior_points(in_cart, old_shape):
    """Generates interior points for each cube in the old grid.

    Args:
        in_cart (numpy.ndarray): Cartesian coordinates of the old grid
        old_shape (list): (nlt, nf, nz) of original sami run.

    Returns:
        centers (numpy.ndarray): cartesian coordinates of the center of each
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

                if lt == old_shape[0] - 1:
                    l2 = 0
                else:
                    l2 = lt + 1

                if z == 0:
                    z2 = -1
                else:
                    z2 = z - 1

                f2 = f + 1
                if f == old_shape[1] - 1:
                    badbadbad.append([lt, f, z])
                    continue

                cs = [[lt, f, z],
                      [lt, f, z2],
                      [l2, f, z2],
                      [l2, f, z],
                      [lt, f2, z],
                      [lt, f2, z],
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

                center = np.sum(xs) / 8, np.sum(ys) / 8, np.sum(zs) / 8

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
        weights (numpy.ndarray): Weights to multiply the old grid by to get
            the new grid
        src_idxs (numpy.ndarray): Indices of the old grid that contribute to
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

        d = np.sqrt((xs - out_cart[0, n])**2 + (ys -
                    out_cart[1, n])**2 + (zs - out_cart[2, n])**2)

        weights[n] = 1 / (d)
        src_idxs[n] = idxs

    print('Done, found %i valid points' % np.sum(weights > 0))

    return weights, src_idxs


def find_pairs(centers, out_cart):
    """Find nearest point between two sets of coordinates.
            (only works well for cartesian coordinates)

    Args:
        centers (numpy.ndarray): Centers of cubes in the old grid
        out_cart (numpy.ndarray): Grid to search from.

    Returns:
        numpy.ndarray: indices of the nearest point in the old grid to each
            point in the new grid (out_cart).
    """
    tree = KDTree(centers)
    _, nearest = tree.query(out_cart.T)
    return nearest


try:
    from numba import jit
    from numba import prange
    print('Numba available! This will drastically speed up applying weights,')

    @jit(nopython=True)
    def numba_do_apply_weights(t0, src_idxs, weights, outv):
        """Speed up applying weight function.

        Args:
            t0 (numpy.ndarray): data (all times)
            src_idxs (numpy.ndarray): indexes from dst to src grid
            weights (numpy.ndarray): generated weights.
        """

        for t in prange(t0.shape[-1]):
            outv[t] += np.sum(weights * np.take(
                t0[:, :, :, t], src_idxs), axis=1) \
                / np.sum(weights, axis=1)

        return outv


except ImportError:
    print('No Numba available')

    def numba_do_apply_weights(t0, src_idxs, weights, outv):
        """Speed up applying weight function.

        Args:
            t0 (numpy.ndarray): data (all times)
            src_idxs (numpy.ndarray): indexes from dst to src grid
            weights (numpy.ndarray): generated weights.
        """

        for t in range(t0.shape[-1]):
            outv[t] += np.sum(weights * np.take(
                t0[:, :, :, t], src_idxs), axis=1) \
                / np.sum(weights, axis=1)

        return outv


def main(
        sami_data_path,
        out_path=None,
        save_weights=True,
        use_saved_weights=True,
        apply_weights=True,
        cols='all',
        dtime_sim_start=None,
        lat_step=2,
        lon_step=5,
        alt_step=50,
        minmax_alt=[100, 2200],
        lat_finerinterps=3,
        lon_finerinterps=3,
        alt_finerinterps=2,
        out_coord_file=None,
        split_by_var=False,
        single_file=False,
        split_by_time=True,
        use_ccmc=False,
        numba=False,
        skip_time_check=False):
    """Interpolate SAMI3 outputs to a new grid.

    Args:
        sami_data_path (str, path-like):
            Path to SAMI3 output files.
        out_path (str; path-like, optional):
            Location to save output files. Defaults to "./"
        save_weights (bool, optional):
            Whether or not to save weight and index files. Defaults to True.
        use_saved_weights (bool, optional):
            Read in existing weight & index files (if they exist). Will not
                error if files don't exist, just create them.
                Defaults to True.
        apply_weights (bool, optional):
            Whether or not to apply weights. Otherwise weight files are made
                and the program exits. Defaults to True.
        cols (str or list-like, optional):
            Columns to read data from. Can be string or list-like.
                Defaults to 'all'.
        dtime_sim_start (datetime.datetime, optional):
            Datetime of simulation start. Required to read raw SAMI outputs.
                Defaults to None.
        lat_step (int, optional):
            Integer step to use in output grid (deg). Defaults to 2.
        lon_step (int, optional):
            Integer step to use in output grid (deg). Defaults to 5.
        alt_step (int, optional):
            Integer step to use in output grid (in km). Defaults to 50.
        minmax_alt (list, optional):
            Min & Max altitude to output in km. Defaults to [100, 2200].
        lat_finerinterps (int, optional):
            Finer interpolation factor to use in latitude.
            Will by coarsened to the desired output resolution before writing.
            Defaults to 3.
        lon_finerinterps (int, optional):
            Finer interpolation factor to use in longitude.
            Defaults to 3.
        alt_finerinterps (int, optional):
            Finer interpolation factor to use in altitude. Defaults to 2.
        out_coord_file (_type_, optional):
            Output coordinates from a file instead of using the default grid.
            Cannot be used with finerinterps.
            Defaults to None (no coordinate file).
        split_by_var (bool, optional):
            Write each variable to a separate file. Defaults to False.
        single_file (bool, optional):
            Write all variables to a single file. Defaults to False.
        split_by_time (bool, optional):
            Write each time to a separate file. Defaults to True.
        use_ccmc (bool, optional):
            Whether or not to use CCMC naming conventions for outputs. Must be
                used with split_by_time. Defaults to False.
        numba (bool, optional):
            Use numba to speed up applying weights. Defaults to False.

    Raises:
        ValueError: If out_coord_file does not have the reuired variables.
        ValueError: If split_by_var and single_file are both True.
    """

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

    if out_coord_file is None:
        latout = np.arange(-90, 90, lat_step / lat_finerinterps)
        lonout = np.arange(0, 360, lon_step / lon_finerinterps)
        altout = np.arange(200, 2200, alt_step / alt_finerinterps)

        out_lats = []
        out_lons = []
        out_alts = []

        for a in latout:
            for o in lonout:
                for l1 in altout:
                    out_lats.append(a)
                    out_lons.append(o)
                    out_alts.append(l1)

        out_cart = latlonalt_to_cart(
            out_lats, out_lons, np.array(out_alts) + 6371)

    else:
        ds_out = xr.open_dataset(out_coord_file)
        print('Read coordinate file.')
        try:
            latout = ds_out['lat'].values
            lonout = ds_out['lon'].values
            altout = ds_out['alt'].values
        except ValueError:
            raise ValueError('Coordinate file must have lat, lon, alt vars.')

        out_cart = latlonalt_to_cart(
            altout, lonout, np.array(altout) + 6371)
        lon_finerinterps = 1
        lat_finerinterps = 1
        alt_finerinterps = 1

    if use_saved_weights:
        weights = np.fromfile(os.path.join(out_path, 'weights'))
        idxs = np.fromfile(os.path.join(out_path, 'indexes'))
        weights = weights.reshape([int(len(weights) / 8), 8])
        idxs = idxs.reshape([int(len(idxs) / 8), 8])
        idxs = idxs.astype(int)
        print('using weights from %s' % out_path)

    else:
        centers, coords = generate_interior_points(in_cart, old_shape)
        nearest = find_pairs(centers, out_cart)

        weights, idxs = make_weights(in_cart, out_cart,
                                     nearest, old_shape, coords)
        idxs = idxs.astype(int)

        if save_weights:
            weights.tofile(os.path.join(out_path, 'weights'))
            idxs.tofile(os.path.join(out_path, 'indexes'))

        del centers, coords

    if apply_weights:

        sami_og_vars = SAMI.sami_og_vars
        if cols == 'all':
            cols = sami_og_vars.values()
        elif isinstance(cols, str):
            cols = [cols]

        # First make sure all cols are valid
        data_files = {}
        for ftype in sami_og_vars:
            if sami_og_vars[ftype] in cols:
                data_files[ftype] = sami_og_vars[ftype]
        if len(data_files.keys()) == 0:
            print('the available data files are: \n', sami_og_vars.keys(),
                  'you gave: ', cols)
            raise ValueError('no data files to read in')

        # then read thru & apply weights.
        for ftype in data_files:

            varname = data_files[ftype]
            print('reading in %s' % ftype)

            data, times = SAMI.read_to_nparray(
                sami_data_path, dtime_sim_start, cols=varname, pbar=False,
                skip_time_check=skip_time_check)

            ds = xr.Dataset(coords={
                'time': (['time'], times),
                'alt': (['alt'], altout),
                'lat': (['lat'], latout),
                'lon': (['lon'], lonout)},)
            varname = list(data['data'].keys())[0]

            outv = np.zeros([data['data'][varname].shape[-1], len(weights)])
            outv = numba_do_apply_weights(data['data'][varname].copy(),
                                          idxs, weights, outv)

            print(
                'received weights from numba function, writing & continuing')

            ds[varname] = (('time', 'lat', 'lon', 'alt'),
                           np.array(outv).reshape(
                len(times), len(latout), len(lonout), len(altout)))

            if not all([i == 1 for i in [lat_finerinterps,
                       lon_finerinterps, alt_finerinterps]]):
                print('coarsening dataset')
                ds = ds.coarsen(lat=lat_finerinterps,
                                lon=lon_finerinterps,
                                alt=alt_finerinterps,
                                boundary='trim').mean()

            if split_by_var:
                for varname in ds.data_vars:
                    ds.to_netcdf(
                        os.path.join(out_path, '%s' % varname))
            elif single_file or out_coord_file is not None:
                ds.to_netcdf(os.path.join(
                    out_path, 'sami-regridded'), mode='a')
            elif split_by_time:
                for t in ds.time.values:
                    fname = make_ccmc_name('SAMI', t, data_type='REGRID')
                    ds.sel(time=t).to_netcdf(
                        os.path.join(out_path, fname), mode='a')

            del data

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
    parse.add_argument('--skip_time_check', action='store_true',
                        help='Skip verifying accuracy of times. Useful when'
                        ' SAMI has been configured to skip some outputs '
                        '(hrpr != 0)')
    parser.add_argument('--split_by_var', action='store', default=True,
                        type=bool,
                        help='Split output files by variable? Default: True')
    parser.add_argument('--single_file', action='store_true', default=False,
                        help='Write all output datasets to single file?'
                        'Default: False')
    parser.add_argument('--custom_grid', action='store_true', default=False,
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
         single_file=args.single_file,
         numba=args.numba,
         skip_time_check=args.skip_time_check)
