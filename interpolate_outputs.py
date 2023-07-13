"""
Code to interpolate any model data to either a:
- user-defined grid
- collection of satellite points


It's easiest to interface with the  `interpolate` function,
    or use the `main` function from the command line...

Can read output coords from cdf/csv or by user-defined grid.
> Default output grid is 5deg lon x 2deg lat x 50km alt

"""

import xarray as xr

from tqdm.auto import tqdm
import numpy as np

import os
from utility_programs.read_routines import SAMI, GITM
from utility_programs.utils import str_to_ut, make_ccmc_name
# import argparse
import glob
import pickle
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
import math


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

def gps_to_ecef_custom(lon, lat, alt, degrees=True, alt_in_m=False):
    a = 6378137.0
    finv = 298.257223563
    f = 1 / finv
    e2 = 1 - (1 - f) * (1 - f)
    
    if not alt_in_m:
        alt = np.asarray(alt) * 1000
    
    if degrees:
        lat = np.deg2rad(lat)
        lon = np.deg2rad(lon)
        
    v = a / np.sqrt(1 - e2 * np.sin(lat) * np.sin(lat))

    x = (v + alt) * np.cos(lat) * np.cos(lon)
    y = (v + alt) * np.cos(lat) * np.sin(lon)
    z = (v * (1 - e2) + alt) * np.sin(lat)

    return np.array([x, y, z])


def do_interpolations(
    sami_data_path=None,
    dtime_sim_start=None,
    gitm_data_path=None,
    gitm_output_each_var=True,
    gitm_output_each_time=False,
    out_lat_lon_alt=None,
    out_path=None,
    out_runname='',
    save_delauney=False,
    max_alt=None,
    cols='all',
    show_progress=False,
    engine='h5netcdf',
    return_ds_too=False,
):
    """Interpolate SAMI (GITM functionality not done yet) to either a
        standard geographic grid or to user-defined points.

    Args:
        sami_data_math (string): path to sami data.
        dtime_sim_start (string/datetime.datetime): Start time of simulation.
            Required to read SAMI data. Can be str (YYYYMMDD) or a pre-computed
            datetime object.
        gitm_data_path (string): path to gitm data.
        gitm_output_each_var (bool): If True, output each variable to a
            separate file. Requires looping through the GITM output files
            multiple times. If False, gitm_output_each_time must be True.
        gitm_output_each_time (bool): If True, output each time to a separate
            file. Will run faster for all variables than gitm_output_each_var,
            but will include variables the user does not care for.
        out_lat_lon_alt (numpy.array): Coordinates to interpolate to.
            Must have dimenstions 3xN, where N is number of points.
            Will be converted to cartesian coordinates.
            Lon and Lat in degrees, Alt in km above earth surface.
        out_path (str): Path to save regridded to.
            Default is same as MODEL_data_path.
        out_runname (str): Descriptive name for output file. Appended to
            out_path + SAMI_REGRID.
        save_delauney (bool): Option to save/read delauney weights from file.
            It takes a while to compute them. Weight file is saved to the
            sami_data_path path with a specification of the max_alt.
            Setting to True allows the program to read weights as well as save
            them.
        max_alt (int): specify maximum altitude of data grid to feed in to
            delauney calculations. Useful if you don't want to recalculate
            weights and the interpoaltion is different from one already done.
        cols (str/list-like): Which variables to interpolate. Default is 'all'.
            Can be any str from
            utility_programs.read_routines.SAMI.sami_og_vars.
        show_progress (bool): Show progressbars? Default: False
            Requires tqdm.
        engine (str): Which engine to use when writing netcdf files.
            Default is 'h5netcdf' but can cause some issues on some systems and
            some python environments. Set to None to use default xarray engine.
        return_ds_too (bool): Set to True to also return the interpolated
            dataset.
            !! Does not support multiple variables.
            !! ONLY works for SAMI (currently).


    Returns:
        Nothing. Unless return_ds_too == True
            - The interpolated data is written to a file.
    """
    if sami_data_path is not None and gitm_data_path is not None:
        raise ValueError('Only one of sami_data_path or gitm_data_path can be'
                         ' specified at a time.')

    if out_path is None:
        if sami_data_path is not None:
            out_path = sami_data_path
        else:
            out_path = gitm_data_path

    # outputs...
    if out_lat_lon_alt is None:
        latout = np.arange(-90, 90, 2)
        lonout = np.arange(0, 360, 5)
        if sami_data_path is not None:
            altout = np.arange(200, 2200, 50)
        else:
            altout = np.arange(120, 670, 50)
        out_lats = []
        out_lons = []
        out_alts = []

        for a in latout:
            for o in lonout:
                for l1 in altout:
                    out_lats.append(a)
                    out_lons.append(o)
                    out_alts.append(l1)

        out_lat_lon_alt = latlonalt_to_cart(
            out_lats, out_lons, np.array(out_alts) + 6371)
    else:
<<<<<<< HEAD
        out_lat_lon_alt = latlonalt_to_cart(
            out_lat_lon_alt[0], out_lat_lon_alt[1],
            np.array(out_lat_lon_alt[2]) + 6371)
=======
        latout = np.array(out_lat_lon_alt[0])
        lonout = np.array(out_lat_lon_alt[1])
        altout = np.array(out_lat_lon_alt[2])

    out_lats = []
    out_lons = []
    out_alts = []

    for a in latout:
        for o in lonout:
            for l1 in altout:
                out_lats.append(a)
                out_lons.append(o)
                out_alts.append(l1)

    out_lat_lon_alt = gps_to_ecef_custom(
        out_lons, out_lats, out_alts)
>>>>>>> 548a59d664ddd861776d5dfd38fafee4f0fddf7f

    # deal with sami first
    if sami_data_path is not None:

        if isinstance(dtime_sim_start, str):
            dtime_sim_start = str_to_ut(dtime_sim_start)

        nz, nf, nlt, nt = SAMI.get_grid_elems_from_parammod(sami_data_path)
        # old_shape = [nlt, nf, nz]
        grid = SAMI.get_sami_grid(sami_data_path, nlt, nf, nz)

        # specify max alt to build delauney at
        if max_alt is None:
            if out_lat_lon_alt is not None:
                # alt will be the biggest coord:
                max_alt = np.max(out_lat_lon_alt) + 300
                # add 300 to make sure we have enough points above top
            else:
                max_alt = 2500

        mask = np.where(grid['alt'] < max_alt)
        grid2 = {}
        for k in grid.keys():
            grid2[k] = grid[k][mask].flatten()
        del grid

        in_cart = gps_to_ecef_custom(grid2['glon'],
                                    grid2['glat'],
                                    grid2['alt']).T

        if os.path.exists(os.path.join(sami_data_path,
                                       'delauney_max-%i.pkl' % max_alt)):
            if save_delauney:
                print('attempting to reuse existing triangulation file')
                with open(os.path.join(
                        sami_data_path, 'delauney_max-%i.pkl' % max_alt),
                        'rb') as f:
                    tri = pickle.load(f)
            else:
                print('Found existing triangulation file. Recalculating...',
                      '\n(Specify save_delauney=True to reuse)')
                tri = Delaunay(in_cart)
        else:
            print('Calculating Delauney Triangulation..')
            tri = Delaunay(in_cart)
            if save_delauney:
                print('Saving')
                with open(os.path.join(sami_data_path,
                                       'delauney_max-%i.pkl' % max_alt),
                          'wb') as f:
                    pickle.dump(tri, f)

        # format 'cols' variable
        if cols == 'all':
            cols = SAMI.sami_og_vars.items()
        else:
            if isinstance(cols, str):
                cols = [cols]
            else:
                cols = np.asarray(cols)

        if show_progress:
            pbar = tqdm(total=len(cols) * nt,
                        desc='Reading in SAMI data')

        first = True  # for choosing which mode to write
        for data_var in cols:
            interpd = []
            data, times = SAMI.read_to_nparray(
                sami_data_path, dtime_sim_start,
                cols=data_var,
                skip_time_check=True)

            if show_progress:
                pbar.set_description('interpolating')
            for t in tqdm(range(len(times))):
                interp = LinearNDInterpolator(
                    tri,
                    data['data'][data_var][:, :, :, t][mask].flatten())
                interpd.append(interp(out_lat_lon_alt.T))
                if show_progress:
                    pbar.update()
            if show_progress:
                pbar.set_description('writing Dataset...')
            ds = xr.Dataset(coords={
                'time': (['time'], times),
                'alt': (['alt'], altout),
                'lat': (['lat'], latout),
                'lon': (['lon'], lonout)},)
            ds[data_var] = (('time', 'lat', 'lon', 'alt'),
                            np.array(interpd).reshape(
                len(times),
                len(latout),
                len(lonout),
                len(altout)))
            if out_runname != '':
                out_runname = '_' + out_runname + '_'
            ds.to_netcdf(os.path.join(
                out_path, 'SAMI_REGRID' + out_runname + '.nc'),
                engine=engine,
                mode='w' if first else 'a',
                encoding={'time': {'dtype': float}})
            first = False
            if return_ds_too:
                return ds
            del ds, interpd, data  # clean up memory

    if gitm_data_path is not None:
        doraw = False
        if len(glob.glob(os.path.join(gitm_data_path, '*.bin'))) == 0:
            if len(glob.glob(os.path.join(gitm_data_path, '*.nc'))) != 0:
                raise NotImplementedError('NetCDF files not yet supported')

            raise ValueError(
                'No GITM files found in {}'.format(gitm_data_path),
                'Go run pGITM and rerun this with the .bin'
                ' files.')
        else:
            doraw = True

        if doraw:
            if cols == 'all':
                cols = ['all']
            elif isinstance(cols, str):
                cols = [cols]
            else:
                cols = np.asarray(cols)

            # Double check varnames are correct.
            f0 = GITM.read_bin_to_nparrays(
                gitm_dir=gitm_data_path,
                start_idx=0,
                end_idx=1,
                return_vars=True)
            if cols != ['all']:
                for varname in cols:
                    if varname not in f0['gitmvars']:
                        raise ValueError(
                            '{} not found in GITM files'.format(varname),
                            '\nData variables are: \n{}'.format(f0.keys()))
            else:
                cols = f0['gitmvars']

            # make/load Delauney triangulation
            if os.path.exists(os.path.join(gitm_data_path,
                                           'delauney.pkl')):
                if save_delauney:
                    print('attempting to reuse existing triangulation file')
                    with open(os.path.join(
                            gitm_data_path, 'delauney.pkl'),
                            'rb') as f:
                        tri = pickle.load(f)
                else:
                    print(
                        'Found existing triangulation file. Recalculating...',
                        '\n(Specify save_delauney=True to reuse)')
                    tri = Delaunay(in_cart)
            else:
                print('Calculating Delauney Triangulation..')

                in_lat = f0['gitmgrid']['latitude'].flatten()
                in_lon = f0['gitmgrid']['longitude'].flatten()
                in_alt = f0['gitmgrid']['altitude'].flatten()

                in_cart = latlonalt_to_cart(in_lon, in_lat, in_alt).T

                tri = Delaunay(in_cart)
                if save_delauney:
                    print('Saving')
                    with open(os.path.join(gitm_data_path,
                                           'delauney.pkl'), 'wb') as f:
                        pickle.dump(tri, f)

            # Now start grabbing model outputs...
            files = glob.glob(os.path.join(gitm_data_path, '*.bin'))
            times = GITM.gitm_times_from_filelist(files)
            numfiles = len(files)

            if gitm_output_each_var:
                if show_progress:
                    pbar = tqdm(total=len(cols) * numfiles)
                for varname in cols:
                    if show_progress:
                        pbar.set_description('Reading in GITM data')
                    darr = GITM.read_bin_to_nparrays(
                        gitm_dir=gitm_data_path,
                        cols=[varname],
                        progress_bar=False)
                    interpd = []
                    for t in range(numfiles):
                        interp = LinearNDInterpolator(
                            tri,
                            darr['gitmbins'][t, 0, :, :].T.flatten())
                        interpd.append(interp(out_lat_lon_alt.T))
                        if show_progress:
                            pbar.update()
                    if show_progress:
                        pbar.set_description('writing Dataset...')
                    ds = xr.Dataset(coords={
                        'time': (['time'], times),
                        'alt': (['alt'], altout),
                        'lat': (['lat'], latout),
                        'lon': (['lon'], lonout)},)
                    ds[varname] = (('time', 'lat', 'lon', 'alt'),
                                   np.array(interpd).reshape(
                        len(times),
                        len(latout),
                        len(lonout),
                        len(altout)))
                    ds.to_netcdf(os.path.join(
                        out_path,
                        'GITM_INTERP%s_%s.nc' % ('_' + out_runname if
                                                 out_runname != '' else '',
                                                 varname)),
                                 engine=engine,
                                 mode='w',
                                 encoding={'time': {'dtype': float}})
                    del ds, interpd, darr  # clean up memory

            elif gitm_output_each_time:
                for t in range(numfiles):
                    if show_progress:
                        pbar.set_description('Reading in GITM data')
                    darr = GITM.read_bin_to_nparrays(
                        gitm_dir=gitm_data_path,
                        cols=cols,
                        start_idx=t,
                        end_idx=t + 1,
                        return_vars=True,
                        progress_bar=False)
                    interpd = []
                    ds = xr.Dataset(coords={
                        'time': (['time'], [times[t]]),
                        'alt': (['alt'], altout),
                        'lat': (['lat'], latout),
                        'lon': (['lon'], lonout)},)
                    for varnum, varname in enumerate(cols):
                        # print(darr['gitmbins'][t, varnum, :, :].shape,
                        #       darr['gitmbins'][
                        #           t, varnum, :, :].T.flatten().shape)
                        interp = LinearNDInterpolator(
                            tri,
                            darr['gitmbins'][0, varnum, :, :].T.flatten())

                        ds[varname] = (
                            ('time', 'lat', 'lon', 'alt'), np.array(
                                interp(out_lat_lon_alt.T)).reshape(
                                    1,  # single time value
                                    len(latout),
                                    len(lonout),
                                    len(altout)))

                        if show_progress:
                            pbar.update()

                    if show_progress:
                        pbar.set_description('writing Dataset...')
                    fname = make_ccmc_name('GITM_REGRID',
                                           times[t],
                                           out_runname if out_runname != ''
                                           else '')
                    ds.to_netcdf(os.path.join(out_path, fname),
                                 engine=engine,
                                 mode='w',
                                 encoding={'time': {'dtype': float}})
                    del ds, interpd, darr  # clean up memory
