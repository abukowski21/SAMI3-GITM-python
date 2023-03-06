
import argparse
# import re
import datetime

import numpy as np
from spacepy.coordinates import Coords
from spacepy.time import Ticktock

from utility_programs.read_routines import SAMI

try:
    from tqdm import tqdm
    pbar = True
except ImportError:
    pbar = False

from scipy.interpolate import LinearNDInterpolator

# def get_command_line_args(argv):

#     args = {'sami_path': '', dtime_storm_start: '', dtime_sim_start: '0',
#             t_start_idx: 0, t_end_idx: -1, cols: 'all', help: False,
#      thread: 48}

#     for arg in argv:


def convert_to_cart(lats_, lons_, alts_, dtime, coord_source):
    """Converts from GEO to cartesian coordinates at a single datetime.

    Args:
        lats_ (_type_):
            Latitudes to perform the conversion at.
        lons_ (_type_):
            Longitudes to perform the conversion at.
        alts_ (int or list-like):
            Altitude(s) to perform the conversion at. When using magnetic
            sources, this needs to be 'magnetic altitude'. In geographic
            coordinates, this is the altitude above the surface of the Earth.
        dtime (datetime.datetime):
            Datetime object at which to perform the conversion.
        coord_source (str):
            Which type of coordinates to convert from. Options are currently:
            'GEO' or 'MAG'.

    Returns:
        array: 3D array of cartesian coordinates
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
    if coord_source == 'MAG':
        coords = Coords(coord_arr, 'MAG', 'sph')
    elif coord_source == 'GEO':
        coords = Coords(coord_arr, 'GEO', 'sph')
    else:
        raise ValueError('coord_source must be either "MAG" or "GEO"')
    coords.ticks = Ticktock([dtime for k in range(len(coord_arr))])
    newcoords = coords.convert('GEO', 'car')

    return newcoords


def do_interpolating(out_lats, out_lons, out_alts, times,
                     columns='all'):
    """_summary_

    Args:
        out_lats (_type_): _description_
        out_lons (_type_): _description_
        out_alts (_type_): _description_
        times (_type_): _description_
        columns (str, optional): _description_. Defaults to 'all'.
        threading (bool, optional): _description_. Defaults to False.

    Raises:
        ValueError: _description_

    Returns:
        _type_: _description_
    """

    # check columns
    if columns == 'all':
        columns = sami_data.keys()
    elif type(columns) == str:
        columns = [columns]
    for col in columns:
        if col not in sami_data.keys():
            raise ValueError('column {} not in sami_data'.format(col))

    preds = {}
    if pbar:
        progress = tqdm(total=len(times),
                        desc='making preds... pbar is a very rough estimate! ')

    norm_alts = (sami_data['grid']['alt'].flatten() < (max(out_alts) + 300)
                 ) & (sami_data['grid']['alt'].flatten() > (min(out_alts)-75))

    for col in columns:
        preds[col] = np.zeros(
            [len(times), len(out_lats), len(out_lons), len(out_alts)])

    for ntime, dt in enumerate(times):
        grid_cart = convert_to_cart(sami_data['grid']['mlat'],
                                    sami_data['grid']['mlon'],
                                    sami_data['grid']['malt'], ntime, 'MAG')

        xs = grid_cart.data[:, 0][norm_alts]
        ys = grid_cart.data[:, 1][norm_alts]
        zs = grid_cart.data[:, 2][norm_alts]

        loc_grid = list(zip(xs, ys, zs))

        out_grid = convert_to_cart(out_lats, out_lons, out_alts, ntime, 'GEO')

        out_xs = out_grid.data[:, 0]
        out_ys = out_grid.data[:, 1]
        out_zs = out_grid.data[:, 2]

        datas = sami_data[col][:, :, :, ntime].flatten()[norm_alts]

        interp = LinearNDInterpolator(loc_grid, datas, rescale=True)
        # print('interpolating')
        pred = interp(list(zip(out_xs, out_ys, out_zs)))

        # break

        preds[col][ntime] = pred.reshape(
            [len(out_lats), len(out_lons), len(out_alts)])

        if pbar:
            progress.update(1)
    return preds


def set_up_interpolations(out_lats=None, num_out_lats=90, specific_lat=None,
                          lat_lim=60, out_lons=None, num_out_lons=90,
                          lon_lim=180, specific_lon=None, out_alts=None,
                          num_out_alts=10, min_alt=200, max_alt=1000,
                          specific_alts=None, args=None):

    if args is not None:
        out_lats = args.out_lats
        num_out_lats = args.num_out_lats
        specific_lat = args.lat
        lat_lim = args.lat_lim
        out_lons = args.out_lons
        num_out_lons = args.num_out_lons
        lon_lim = args.lon_lim
        specific_lon = args.lon
        out_alts = args.out_alts
        num_out_alts = args.num_out_alts
        min_alt = args.min_alt
        max_alt = args.max_alt
        specific_alts = args.alts

    print(
        out_lats, num_out_lats, specific_lat, lat_lim, out_lons, num_out_lons,
        lon_lim, specific_lon, out_alts, num_out_alts, min_alt, max_alt,
        specific_alts)

    if out_lats is None:
        out_lats = np.linspace(-lat_lim, lat_lim, num_out_lats)

    if specific_lat:
        out_lats = np.append(out_lats, specific_lat)

    if out_lons is None:
        out_lons = np.linspace(-lon_lim, lon_lim, num_out_lons)

    if specific_lon:
        out_lons = np.append(out_lons, specific_lon)

    if out_alts is None:
        out_alts = np.linspace(min_alt, max_alt, num_out_alts)

    if specific_alts:
        out_alts = np.append(out_alts, specific_alts)

    return out_lats, out_lons, out_alts


def main(args):
    """_summary_

    Args:
        sami_data_path (_type_): _description_
        dtime_sim_start (_type_): _description_
        dtime_storm_start (_type_): _description_

    Returns:
        _type_: _description_
    """

    # FORMAT TIMES
    dtime_sim_start = datetime.datetime.strptime(
        args.dtime_sim_start.ljust(14, '0'), '%Y%m%d%H%M%S')
    dtime_storm_start = datetime.datetime.strptime(
        args.dtime_storm_start.ljust(14, '0'), '%Y%m%d%H%M%S')

    # Read in SAMI data
    global sami_data
    sami_data, times_array = SAMI.read_sami_data(
        sami_data_path=args.sami_data_path,
        dtime_sim_start=dtime_sim_start,
        dtime_storm_start=dtime_storm_start,
        t_start_idx=args.t_start_idx, t_end_idx=args.t_end_idx,
        pbar=pbar, cols=args.columns,
        help=args.info)

    cols = sami_data['data'].keys()

    out_lats, out_lons, out_alts = set_up_interpolations(args=args)

    preds = do_interpolating(out_lats, out_lons, out_alts, times_array,
                             columns='all')

    print('writing files...')

    np.array(times_array).tofile(args.data_out_path + 'times',
                                 format='%s', sep=',')
    out_lats.tofile(args.data_out_path + 'out-lats', sep=',')
    out_lons.tofile(args.data_out_path + 'out-lons', sep=',')
    out_alts.tofile(args.data_out_path + 'out-alts', sep=',')

    np.array(preds[cols[0]].shape).tofile(
        args.data_out_path + 'out-shape', sep=',')

    for col in cols:

        preds[col].tofile(args.data_out_path + 'preds-' + col, sep=',')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Transforms SAMI data into a normal geographic grid')

    parser.add_argument(
        'dtime_sim_start',
        help='Datetime of simulation start. Format YYYYMMDDHHmmss',
        action='store')
    parser.add_argument(
        'dtime_storm_start',
        help='Datetime of storm start. Format YYYYMMDDHHmmss',
        action='store')
    parser.add_argument(
        '-sami_data_path', type=str,
        help='Path to SAMI data', default='./data_dir', action='store')
    parser.add_argument(
        '--out_lats', type=str,
        help='List of latitudes to interpolate to. Format: [lat1, lat2, lat3]',
        action='store', default=None, required=False)
    parser.add_argument(
        '--num_out_lats', type=int,
        help='Number of latitudes to interpolate to. Default: 90',
        action='store', default=90, required=False)
    parser.add_argument(
        '--lat', type=float,
        help='Specific latitude to interpolate to. Default: None',
        action='append', default=None, required=False)
    parser.add_argument(
        '--lat_lim', type=float,
        help='Latitude limit for interpolation. Default: 60',
        action='store', default=65, required=False)
    parser.add_argument(
        '--out_lons', type=str,
        help='List of longitudes to interpolate. Format: [lon1, lon2, lon3]',
        action='store', default=None, required=False)
    parser.add_argument(
        '--lon_lim', type=float,
        help='Longitude limit for interpolation. Default: 180',
        action='store', default=180, required=False)
    parser.add_argument(
        '--num_out_lons', type=int,
        help='Number of longitudes to interpolate to. Default: 90',
        action='store', default=90, required=False)
    parser.add_argument(
        '--lon',
        help='Specific longitude to interpolate to. Default: None',
        action='append', default=None, required=False)
    parser.add_argument(
        '--alts', type=str,
        help='List of altitudes to interpolate to. Format: [alt1, alt2, alt3]',
        action='store', default=None, required=False)
    parser.add_argument(
        '--num_out_alts', type=int,
        help='Number of altitudes to interpolate to. Default: 10',
        action='store', default=10, required=False)
    parser.add_argument(
        '--out_alts', type=str,
        help='List of altitudes to interpolate to. Format: [alt1, alt2, alt3]',
        action='store', default=None, required=False)
    parser.add_argument(
        '--min_alt', type=int,
        help='Minimum altitude to interpolate to. Default: 200',
        action='store', default=200, required=False)
    parser.add_argument(
        '--max_alt', type=int,
        help='Maximum altitude to interpolate to. Default: 1000',
        action='store', default=1000, required=False)
    parser.add_argument(
        '--alt',
        help='Specific altitudes to interpolate. Format: [alt1, alt2, alt3]',
        action='append', default=None, required=False)
    parser.add_argument(
        '--columns', type=str, help='Columns to interpolate. Default: all',
        action='store', default='all', required=False)
    parser.add_argument(
        '--threading',
        help='Use threading. Default: False',
        action='store_true', default=False, required=False)
    parser.add_argument(
        '--num_workers', type=int,
        help='Number of workers to use. Default: 48',
        action='store', default=48, required=False)
    parser.add_argument(
        '--data_out_path', type=str,
        help='Path to save data to. Default: ./data_dir',
        action='store', default='./data_dir', required=False)
    parser.add_argument(
        '--info', help='Print run information. Will be a lot',
        action='store_true', required=False, default=False)
    parser.add_argument(
        '--t_start_idx', type=int,
        action='store', default=None, required=False)
    parser.add_argument(
        '--t_end_idx', type=int,
        action='store', default=None, required=False)

    args = parser.parse_args()

    main(args)
