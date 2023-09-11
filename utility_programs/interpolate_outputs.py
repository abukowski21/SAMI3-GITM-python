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
import pandas as pd
import glob
import pickle
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator


def gps_to_ecef_custom(lon, lat, alt, degrees=True, alt_in_m=False):
    """Convert GPS coordinates to ECEF coordinates.
    Units are degrees and km unless specified otherwise.

    :param lon: Longitude
    :type lon: float or numpy.array
    :param lat: Latitude
    :type lat: float or numpy.array
    :param alt: Altitude
    :type alt: float or numpy.array
    :param degrees: Is input Lat/Lon in degrees? Defaults to True
    :type degrees: bool, optional
    :param alt_in_m: Is input Alt in km? Defaults to False
    :type alt_in_m: bool, optional
    :return: ECEF coordinates
    :rtype: numpy.array
    """
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
    skip_time_check=False,
    out_lat_lon_alt=None,
    out_path=None,
    out_runname=None,
    sat_times=None,
    cols='all',
    show_progress=False,
    gitm_data_path=None,
    gitm_output_each_var=True,
    gitm_output_each_time=False,
    is_grid=False,
    aarons_mods=False,
    sami_mintime=0,
    save_delauney=False,
    max_alt=None,
    engine='h5netcdf',
    return_ds_too=False,
):
    """Interpolate SAMI (GITM functionality not fully tested) to either a
    standard geographic grid or to user-defined points.

    :param sami_data_path: path to sami data.
    :type sami_data_path: str, optional
    :param dtime_sim_start: Start time of simulation. Required to read SAMI
        data. Can be str (YYYYMMDD) or a pre-computed datetime object.
    :type dtime_sim_start: str or datetime, optional
    :param skip_time_check: (bool) If True, skip the check to make sure the
        SAMI times are self-consistent (not always true...).
    :type skip_time_check: bool, optional
    :param out_lat_lon_alt: (numpy.array) Coordinates to interpolate to. Must
        have dimenstions 3xN, where N is number of points. Will be converted to
        cartesian coordinates. Lon and Lat in degrees, Alt in km above earth
        surface.
    :type out_lat_lon_alt: numpy.array, optional
    :param out_path: (str) Path to save regridded to. Default is same as
        MODEL_data_path.
    :type out_path: str or os.PathLike, optional
    :param out_runname: (str) Descriptive name for output file. Saved as:
        out_path/{out_runname + "SAMI_REGRID.nc"}.
    :type out_runname: str, optional
    :param sat_times: (list) List of times to interpolate to. Must be a list
        of (python or pandas) datetime objects.
    :type sat_times: list, optional
    :param cols: (str/list-like) Which variables to interpolate. Default is
        'all'. Can be any str from
        utility_programs.read_routines.SAMI.sami_og_vars.
    :type cols: str or list-like, optional
    :param show_progress: (bool) Show progressbars? Default: False.
        Requires tqdm.
    :type show_progress: bool, optional
    :param gitm_data_path: (string) path to gitm data.
    :type gitm_data_path: str, optional
    :param gitm_output_each_var: (bool) If True, output each variable to a
        separate file. Requires looping through the GITM output files multiple
        times. If False, gitm_output_each_time must be True.
    :type gitm_output_each_var: bool, optional
    :param gitm_output_each_time: (bool) If True, output each time to a
        separate file. Will run faster for all variables than
        gitm_output_each_var, but will include variables the user does not
        care for.
    :type gitm_output_each_time: bool, optional
    :param save_delauney: (bool) Option to save/read delauney weights from
        file. It takes a while to compute them. Weight file is saved to the
        sami_data_path path with a specification of the max_alt. Setting to
        True allows the program to read weights as well as save them.
    :type save_delauney: bool, optional
    :param max_alt: (int) specify maximum altitude of data grid to feed in to
        delauney calculations. Useful if you don't want to recalculate weights
        and the interpoaltion is different from one already done.
    :type max_alt: int, optional
    :param engine: (str) Which engine to use when writing netcdf files. Default
        is 'h5netcdf' but can cause some issues on some systems and some python
        environments. Set to None to use default xarray engine.
    :type engine: str, optional
    :param return_ds_too: (bool) Set to True to also return the interpolated
        dataset. !! Does not support multiple variables. !! ONLY works for SAMI
        (currently).
    :type return_ds_too: bool, optional

    :raises ValueError: Only one of sami_data_path or gitm_data_path can be
        specified at a time.
    :raises ValueError: Must specify sat_times if not using a grid
    :raises NotImplementedError: Interpolating from NetCDF files not yet
        supported
    :raises ValueError: No GITM files found in gitm_data_path. Go run pGITM and
        rerun this with the .bin files.
    :raises ValueError: Invalid column requested.
    :return: Interpolated data. Optional, only returned if return_ds_too=True.
    :rtype: xarray.Dataset
    """

    if sami_data_path is not None and gitm_data_path is not None:
        raise ValueError(
            'Only one of sami_data_path or gitm_data_path can be specified'
            ' at a time.')

    if out_path is None:
        if sami_data_path is not None:
            out_path = sami_data_path
        else:
            out_path = gitm_data_path

    # outputs...
    if out_lat_lon_alt is None:
        if aarons_mods:
            latout = np.arange(-90, 90, 0.5)
            lonout = np.arange(0, 360, 1)
            if sami_data_path is not None:
                altout = np.arange(150, 2200, 25)
            else:
                altout = np.arange(120, 670, 25)
        else:
            latout = np.arange(-90, 90, 2)
            lonout = np.arange(0, 360, 4)
            if sami_data_path is not None:
                altout = np.arange(200, 2200, 25)
            else:
                altout = np.arange(120, 670, 25)

        out_lats, out_lons, out_alts = np.meshgrid(latout, lonout, altout)

    else:
        out_lats = out_lat_lon_alt[0]
        out_lons = out_lat_lon_alt[1]
        out_alts = out_lat_lon_alt[2]

        if sat_times is None and not is_grid:
            raise ValueError('Must specify sat_times if not using a grid')

        if is_grid:
            latout = np.unique(out_lats)
            lonout = np.unique(out_lons)
            altout = np.unique(out_alts)

    out_lon_lat_alt = gps_to_ecef_custom(
        out_lons.flatten(), out_lats.flatten(), out_alts.flatten()).T

    # deal with sami first
    if sami_data_path is not None:
        if isinstance(dtime_sim_start, str):
            dtime_sim_start = str_to_ut(dtime_sim_start)

        nz, nf, nlt, nt = SAMI.get_grid_elems_from_parammod(sami_data_path)
        # old_shape = [nlt, nf, nz]
        grid = SAMI.get_sami_grid(sami_data_path, nlt, nf, nz)

        # specify max alt to build delauney at
        if max_alt is None:
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
            cols = SAMI.sami_og_vars.values()
        else:
            if isinstance(cols, str):
                cols = [cols]
            else:
                cols = np.asarray(cols)

        if show_progress:
            pbar = tqdm(total=len(cols) * (nt - sami_mintime),
                        desc='Reading in SAMI data')
        else:
            print('Interpolating SAMI data...')

        first = True  # for choosing which mode to write
        numcol_for_pbar = 1
        for data_var in cols:

            if out_runname is None:
                out_runname = data_var

            data, times = SAMI.read_to_nparray(
                sami_data_path, dtime_sim_start,
                cols=data_var,
                skip_time_check=skip_time_check)

            if show_progress:
                pbar.set_description('interpolating %s (%i/%i)'
                                     % (data_var, numcol_for_pbar, len(cols)))

            if is_grid:  # holds interpolated grid
                ds = xr.Dataset(coords={
                    'time': (['time'], times[sami_mintime:]),
                    'alt': (['alt'], altout),
                    'lat': (['lat'], latout),
                    'lon': (['lon'], lonout)})
                ds[data_var] = (('time', 'lon', 'lat', 'alt'),
                                np.zeros([len(times) - sami_mintime,
                                          len(lonout),
                                          len(latout),
                                          len(altout)]))
            # If we have a list of points that aren't a grid, do it this way
            else:
                ds = xr.Dataset(
                    coords={
                        'sami_time': (['sami_time'], times),
                        'glat': (['sat_step'], out_lats),
                        'glon': (['sat_step'], out_lons),
                        'alt': (['sat_step'], out_alts),
                        'sat_time': (['sat_step'], [
                            pd.Timestamp(i) for i in sat_times])})
                ds[data_var] = (('sami_time', 'sat_step'),
                                np.zeros([len(times), len(sat_times)]))

            for t in range(len(times) - sami_mintime):
                interp = LinearNDInterpolator(
                    tri,
                    data['data'][data_var][:, :, :, t + sami_mintime][mask]
                    .flatten())
                if is_grid:
                    ds[data_var][t] = interp(out_lon_lat_alt).reshape(
                        len(lonout),
                        len(latout),
                        len(altout))
                else:
                    ds[data_var][t] = interp(out_lon_lat_alt)

                if show_progress:
                    pbar.update()
            if show_progress:
                pbar.set_description('writing Dataset...')

            if aarons_mods:
                ds2 = ds.coarsen(lat=4, alt=2, lon=8,
                                 boundary='pad').mean()
                del ds
                ds = ds2
                del ds2

                ds = ds.interp(
                    lon=np.linspace(0, 360, 90),
                    lat=np.linspace(-90, 90, 90),
                    alt=np.arange(altout[1], altout[-2], 50))

            ds.to_netcdf(
                os.path.join(
                    out_path,
                    out_runname +
                    '_SAMI-REGRID.nc' if is_grid else
                    out_runname + '_SAMI-INTERP.nc'),
                engine=engine,
                mode='w' if first else 'a')
            first = False
            numcol_for_pbar += 1
            if return_ds_too:
                return ds
            del ds, data  # clean up memory

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

                in_cart = gps_to_ecef_custom(in_lon, in_lat, in_alt).T

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
                first = True
                numcol_for_pbar = 1
                for varname in cols:
                    if show_progress:
                        pbar.set_description('Reading in GITM data')
                    darr = GITM.read_bin_to_nparrays(
                        gitm_dir=gitm_data_path,
                        cols=[varname],
                        progress_bar=False)
                    interpd = []
                    if show_progress:
                        pbar.set_description(
                            'interpolating %s (%i/%i)' %
                            (varname, numcol_for_pbar, len(cols)))
                    for t in range(numfiles):
                        interp = LinearNDInterpolator(
                            tri,
                            darr['gitmbins'][t, 0, :, :].T.flatten())
                        interpd.append(interp(out_lon_lat_alt.T))
                        if show_progress:
                            pbar.update()
                    if show_progress:
                        pbar.set_description('writing Dataset...')
                    ds = xr.Dataset(coords={
                        'time': (['time'], times),
                        'alt': (['alt'], altout),
                        'lon': (['lat'], lonout),
                        'lat': (['lon'], latout)},)
                    ds[varname] = (('time', 'lon', 'lat', 'alt'),
                                   np.array(interpd).reshape(
                        len(times),
                        len(lonout),
                        len(latout),
                        len(altout)))
                    ds.to_netcdf(os.path.join(
                        out_path,
                        '%sGITM-INTERP.nc' % (out_runname + '_' if
                                              out_runname != '' else '')),
                                 engine=engine,
                                 mode='w' if first else 'a',
                                 encoding={'time': {'dtype': float}})
                    first = False
                    numcol_for_pbar += 1
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
                        'lon': (['lon'], lonout),
                        'lat': (['lat'], latout)},)
                    for varnum, varname in enumerate(cols):
                        # print(darr['gitmbins'][t, varnum, :, :].shape,
                        #       darr['gitmbins'][
                        #           t, varnum, :, :].T.flatten().shape)
                        interp = LinearNDInterpolator(
                            tri,
                            darr['gitmbins'][0, varnum, :, :].T.flatten())

                        ds[varname] = (
                            ('time', 'lon', 'lat', 'alt'), np.array(
                                interp(out_lon_lat_alt.T)).reshape(
                                    1,  # single time value
                                    len(lonout),
                                    len(latout),
                                    len(altout)))

                        if show_progress:
                            pbar.update()

                    if show_progress:
                        pbar.set_description('writing Dataset...')
                    fname = make_ccmc_name('GITM-REGRID',
                                           times[t],
                                           out_runname if out_runname != ''
                                           else '')
                    ds.to_netcdf(os.path.join(out_path, fname),
                                 engine=engine,
                                 mode='w',
                                 encoding={'time': {'dtype': float}})
                    del ds, interpd, darr  # clean up memory
