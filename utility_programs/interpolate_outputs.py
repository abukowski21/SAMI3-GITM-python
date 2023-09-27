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

from tqdm import tqdm
import numpy as np

import os
from utility_programs.read_routines import SAMI, GITM
from utility_programs.utils import str_to_ut, make_ccmc_name
import pandas as pd
import glob
import pickle
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
from multiprocessing import Pool
from itertools import repeat

def interpolate_var(tri1, tri2, outpts1, outpts2, 
                    indata, ntime, out_shape, mask):
    """
    This is a refactor of the SAMI interps so that we can thread it.

    """
    np.seterr(all="ignore")

    interp1=LinearNDInterpolator(
                        tri1,
                        indata.T[mask]
                        .flatten())
    tmp1=interp1(outpts1).reshape(out_shape)
    interp2=LinearNDInterpolator(tri2, indata.T[mask].flatten())
    return np.nanmean([tmp1, interp2(outpts2).reshape(out_shape)],
                    axis=0)


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
    sami_mintime=0,
    save_delauney=False,
    max_alt=None,
    engine='h5netcdf',
    return_ds_too=False,
    num_workers=16,
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
    :param show_progress: (bool) Show progress bars? Default: False.
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
        and the interpolation is different from one already done.
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
        latout = np.arange(-90, 90, 1)
        lonout = np.arange(0, 360, 2)
        if sami_data_path is not None:
            altout = np.arange(100, 2600, 25)
        else:
            altout = np.arange(100, 725, 25)

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

    out_pts1 = np.array([out_lons.flatten(),
                                out_lats.flatten(),
                                out_alts.flatten()]).T
    out_pts2 = np.array([np.where(out_lons.flatten() > 180, 
                                     out_lons.flatten() - 360,
                                     out_lons.flatten()),
                         out_lats.flatten(),
                         out_alts.flatten()]).T

    # deal with sami first
    if sami_data_path is not None:
        if isinstance(dtime_sim_start, str):
            dtime_sim_start = str_to_ut(dtime_sim_start)

        nz, nf, nlt, nt = SAMI.get_grid_elems_from_parammod(sami_data_path)
        # old_shape = [nlt, nf, nz]
        grid = SAMI.get_sami_grid(sami_data_path, nlt, nf, nz)

        # specify max alt to build delauney at
        if max_alt is None:
            max_alt = 3000

        mask = np.where(grid['alt'] < max_alt)
        grid2 = {}
        for k in grid.keys():
            grid2[k] = grid[k][mask].flatten()
        del grid

        in_pts1 = np.array([grid2['glon'],
                             grid2['glat'],
                             grid2['alt']]).T
        in_pts2 = np.array([np.where(grid2['glon'] > 180, 
                                     grid2['glon'] - 360,
                                     grid2['glon']),
                            grid2['glat'],
                            grid2['alt']]).T

        if os.path.exists(os.path.join(sami_data_path,
                                       'delauney_max-%i-1.pkl' % max_alt)):
            if save_delauney:
                print('attempting to reuse existing triangulation file')
                with open(os.path.join(
                        sami_data_path, 'delauney_max-%i-1.pkl' % max_alt),
                        'rb') as f:
                    tri1 = pickle.load(f)
                with open(os.path.join(
                        sami_data_path, 'delauney_max-%i-2.pkl' % max_alt),
                        'rb') as f:
                    tri2 = pickle.load(f)
            else:
                print('Found existing triangulation file. Recalculating...',
                      '\n(Specify save_delauney=True to reuse)')
                tri1 = Delaunay(in_pts1)
                tri2 = Delaunay(in_pts2)
        else:
            print('Calculating Delauney Triangulation..')
            tri1 = Delaunay(in_pts1)
            tri2 = Delaunay(in_pts2)
            if save_delauney:
                print('Saving')
                with open(os.path.join(sami_data_path,
                                       'delauney_max-%i-1.pkl' % max_alt),
                          'wb') as f:
                    pickle.dump(tri1, f)
                with open(os.path.join(sami_data_path,
                                       'delauney_max-%i-2.pkl' % max_alt),
                          'wb') as f:
                    pickle.dump(tri2, f)

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
        out_shape = [len(lonout),len(latout), len(altout)]

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
                ds=xr.Dataset(coords={
                    'time': (['time'], times[sami_mintime:]),
                    'alt': (['alt'], altout),
                    'lat': (['lat'], latout),
                    'lon': (['lon'], lonout)})
                print('These warnings are OK! just means '
                      'the computer is thinking ;)')

            # If we have a list of points that aren't a grid, the index variable
            #   is a little different
            else:
                ds=xr.Dataset(
                    coords={
                        'sami_time': (['sami_time'], times),
                        'glat': (['sat_step'], out_lats),
                        'glon': (['sat_step'], out_lons),
                        'alt': (['sat_step'], out_alts),
                        'sat_time': (['sat_step'], [
                            pd.Timestamp(i) for i in sat_times])})
                ds[data_var]=(('sami_time', 'sat_step'),
                                np.zeros([len(times), len(sat_times)]))

            if is_grid:
                with Pool(num_workers) as pool:
                    ds[data_var]=(('time', 'lon', 'lat', 'alt'),
                                  pool.starmap(interpolate_var,
                                               zip(repeat(tri1), repeat(tri2),
                                               repeat(out_pts1), repeat(out_pts2),
                                               data['data'][data_var].T,
                                               range(len(times)-sami_mintime),
                                               repeat(out_shape),
                                               repeat(mask))))
                    if show_progress:
                        pbar.update(len(times)-sami_mintime)
            else:
                for t in range(len(times) - sami_mintime):
                    interp1=LinearNDInterpolator(
                        tri1,
                        data['data'][data_var][:, :, :, t + sami_mintime][mask]
                        .flatten())
                    interp2=LinearNDInterpolator(
                        tri2,
                        data['data'][data_var][:, :, :, t + sami_mintime][mask]
                        .flatten())
                    ds[data_var][t]=np.nanmean(np.array([
                                                    interp1(out_pts1),
                                                    interp2(out_pts2)]),
                                                axis=0)

                if show_progress:
                    pbar.update()
            if show_progress:
                pbar.set_description('writing Dataset...')

            ds.to_netcdf(
                os.path.join(
                    out_path,
                    out_runname +
                    '_SAMI-REGRID.nc' if is_grid else
                    out_runname + '_SAMI-INTERP.nc'),
                engine=engine,
                mode='w' if first else 'a')
            first=False
            numcol_for_pbar += 1
            if return_ds_too:
                return ds
            del ds, data  # clean up memory

    if gitm_data_path is not None:
        doraw=False
        do_recalculate_delauney=False

        if len(glob.glob(os.path.join(gitm_data_path, '*.bin'))) == 0:
            if len(glob.glob(os.path.join(gitm_data_path, '*.nc'))) != 0:
                raise NotImplementedError('NetCDF files not yet supported')

            raise ValueError(
                'No GITM files found in {}'.format(gitm_data_path),
                'Go run pGITM and rerun this with the .bin'
                ' files.')
        else:
            doraw=True

        if doraw:
            if cols == 'all':
                cols=['all']
            elif isinstance(cols, str):
                cols=[cols]
            else:
                cols=np.asarray(cols)

            # Double check varnames are correct.
            f0=GITM.read_bin_to_nparrays(
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
                cols=f0['gitmvars']

            # make/load Delauney triangulation
            if os.path.exists(os.path.join(gitm_data_path,
                                           'delauney-1.pkl')):
                if save_delauney:
                    print('attempting to reuse existing triangulation file')
                    with open(os.path.join(
                            gitm_data_path, 'delauney-1.pkl'),
                            'rb') as f:
                        tri1=pickle.load(f)
                    with open(os.path.join(
                            gitm_data_path, 'delauney-2.pkl'),
                            'rb') as f:
                        tri2=pickle.load(f)
                else:
                    print(
                        'Found existing triangulation file. Recalculating...',
                        '\n(Specify save_delauney=True to reuse)')
                    do_recalculate_delauney=True
                    
            else:
                do_recalculate_delauney=True

            if do_recalculate_delauney:
                print('Calculating Delauney Triangulation..')

                in_lat=f0['gitmgrid']['latitude'].flatten()
                in_lon=f0['gitmgrid']['longitude'].flatten()
                in_alt=f0['gitmgrid']['altitude'].flatten()

                in_pts1=np.array([in_lon, in_lat, in_alt]).T
                in_pts2=np.array([np.where(in_lon>180,in_lon-360, in_lon),
                                 in_lat, in_alt]).T

                tri1=Delaunay(in_pts1)
                tri2 = Delaunay(in_pts2)
                if save_delauney:
                    print('Saving')
                    with open(os.path.join(gitm_data_path,
                                           'delauney-1.pkl'), 'wb') as f:
                        pickle.dump(tri1, f)
                    with open(os.path.join(gitm_data_path,
                                           'delauney-2.pkl'), 'wb') as f:
                        pickle.dump(tri2, f)

            # Now start grabbing model outputs...
            files=glob.glob(os.path.join(gitm_data_path, '*.bin'))
            times=GITM.gitm_times_from_filelist(files)
            numfiles=len(files)

            if gitm_output_each_var:
                if show_progress:
                    pbar=tqdm(total=len(cols) * numfiles)
                first=True
                numcol_for_pbar=1
                for varname in cols:
                    if show_progress:
                        pbar.set_description('Reading in GITM data')
                    darr=GITM.read_bin_to_nparrays(
                        gitm_dir=gitm_data_path,
                        cols=[varname],
                        progress_bar=False)
                    interpd=[]
                    if show_progress:
                        pbar.set_description(
                            'interpolating %s (%i/%i)' %
                            (varname, numcol_for_pbar, len(cols)))
                    for t in range(numfiles):
                        interp1=LinearNDInterpolator(
                            tri1,
                            darr['gitmbins'][t, 0, :, :].T.flatten())
                        interp2=LinearNDInterpolator(
                            tri2,
                            darr['gitmbins'][t, 0, :, :].T.flatten())
                        interpd.append(np.nanmean(np.array([
                                                           interp1(out_pts1).T,
                                                           interp2(out_pts2).T],
                                                           axis=0)))
                        if show_progress:
                            pbar.update()
                    if show_progress:
                        pbar.set_description('writing Dataset...')
                    ds=xr.Dataset(coords={
                        'time': (['time'], times),
                        'alt': (['alt'], altout),
                        'lon': (['lat'], lonout),
                        'lat': (['lon'], latout)},)
                    ds[varname]=(('time', 'lon', 'lat', 'alt'),
                                   np.array(interpd).reshape(len(times),
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
                    first=False
                    numcol_for_pbar += 1
                    del ds, interpd, darr  # clean up memory

            elif gitm_output_each_time:
                for t in range(numfiles):
                    if show_progress:
                        pbar.set_description('Reading in GITM data')
                    darr=GITM.read_bin_to_nparrays(
                        gitm_dir=gitm_data_path,
                        cols=cols,
                        start_idx=t,
                        end_idx=t + 1,
                        return_vars=True,
                        progress_bar=False)
                    interpd=[]
                    ds=xr.Dataset(coords={
                        'time': (['time'], [times[t]]),
                        'alt': (['alt'], altout),
                        'lon': (['lon'], lonout),
                        'lat': (['lat'], latout)},)
                    for varnum, varname in enumerate(cols):
                        # print(darr['gitmbins'][t, varnum, :, :].shape,
                        #       darr['gitmbins'][
                        #           t, varnum, :, :].T.flatten().shape)
                        interp1=LinearNDInterpolator(
                            tri,
                            darr['gitmbins'][0, varnum, :, :].T.flatten())
                        interp2 = LinearNDInterpolator(
                            tri2,
                            darr['gitmbins'][0, varnum, :, :].T.flatten())

                        ds[varname]=(
                            ('time', 'lon', 'lat', 'alt'), 
                            np.nanmean(np.array([
                                interp1(out_lon_lat_alt.T).reshape(
                                    1,  # single time value
                                    len(lonout),
                                    len(latout),
                                    len(altout)),
                                interp2(out_pts2.T).reshape(1,  
                                # single time value
                                    len(lonout),
                                    len(latout),
                                    len(altout))]),
                                axis=0))

                        if show_progress:
                            pbar.update()

                    if show_progress:
                        pbar.set_description('writing Dataset...')
                    fname=make_ccmc_name('GITM-REGRID',
                                           times[t],
                                           out_runname if out_runname != ''
                                           else '')
                    ds.to_netcdf(os.path.join(out_path, fname),
                                 engine=engine,
                                 mode='w',
                                 encoding={'time': {'dtype': float}})
                    del ds, interpd, darr  # clean up memory