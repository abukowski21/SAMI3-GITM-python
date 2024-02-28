
"""read sami data.


"""


import datetime
import os

import numpy as np
import pandas as pd
from utility_programs.utils import make_ccmc_name
from tqdm import tqdm

global sami_og_vars
sami_og_vars = {
    'deneu.dat': 'edens',
    'deni1u.dat': 'h+dens',
    'deni2u.dat': 'o+dens',
    'deni3u.dat': 'n+odens',
    'deni4u.dat': 'o2+dens',
    'deni5u.dat': 'he+dens',
    'deni6u.dat': 'n2+dens',
    'deni7u.dat': 'n+dens',
    'denn1u.dat': 'hdens',
    'denn2u.dat': 'odens',
    'denn3u.dat': 'nodens',
    'denn4u.dat': 'o2dens',
    'denn5u.dat': 'hedens',
    'denn6u.dat': 'n2dens',
    'denn7u.dat': 'ndens',
    'teu.dat': 'etemp',
    'ti1u.dat': 'h+temp',
    'ti2u.dat': 'o+temp',
    'ti5u.dat': 'he+temp',
    'vsi1u.dat': 'h+vel_parallel',
    'vsi2u.dat': 'o+vel_parallel',
    'u1pu.dat': 'mer_exb',
    'u3hu.dat': 'zon_exb',
    'u1u.dat': 'zon_neut',
    'u2u.dat': 'mer_neut', }
# Needed for all reads.


def get_grid_elems_from_parammod(sami_data_path):
    """

    Go into sami data directory and get the grid elements
        from the parameter_mod.f90 file.

    Args:
        sami_data_path (str): data path for sami outputs

    Returns:
        nz: num. grid points along each field line
        nf: num. field lines along each magnetic longitude
        nlt: num. magnetic longitudes
        nt: num. time steps
    """

    # Make sure that we only grab the first instance of each var in the file.
    # SOmetimes they repeat and we don't want them
    found = [False, False, [False, False], False]

    nz = 0
    nf = 0
    numwork = 0
    nl = 0
    nt = 0

    try:
        with open(os.path.join(sami_data_path,
                               'parameter_mod.f90'), 'r') as fp:
            # read all lines in a list
            lines = fp.readlines()
            for line in lines:
                # check if string present on a current line

                if not found[0]:
                    if line.find('nz0') != -1:
                        nz0 = []
                        for char in line:
                            if char.isdigit():
                                nz0.append(char)
                        if len(nz0[1:4]) == 3:
                            nz = int(''.join(nz0[1:4]))
                            found[0] = True

                if not found[1]:
                    if line.find('nf') != -1:
                        nf = []
                        for char in line:
                            if char.isdigit():
                                nf.append(char)
                        nf = int(''.join(nf))
                        found[1] = True

                if not found[2][0]:
                    if line.find('nl ') != -1:
                        nl = []
                        for char in line:
                            if char.isdigit():
                                nl.append(char)
                        nl = int(''.join(nl))
                        found[2][0] = True

                if not found[2][1]:
                    if line.find('numwork ') != -1:
                        numwork = []
                        for char in line:
                            if char.isdigit():
                                numwork.append(char)
                        numwork = int(''.join(numwork))
                        found[2][1] = True
    except FileNotFoundError:
        print('Could not read file! Check that the path is correct.',
              'here is what I see in %s:' % sami_data_path,)
        print(os.listdir(sami_data_path))

    # time
    with open(os.path.join(sami_data_path, 'time.dat'), 'r') as fp:
        lines = fp.readlines()
        nt = len(lines) - 1
        found[3] = True

    if not all(found):
        raise ValueError('Could not find all of the required variables in '
                         'parameter_mod.f90 and time.dat. Please check that '
                         'the file is '
                         'formatted correctly.')
    else:
        return nz, nf, numwork * (nl - 2), nt


def get_postprocessed_grid(sami_data_path):
    """

    Go into sami data directory and get the grid elements
        from the parameter_mod.f90 file.

    Args:
        sami_data_path (str): data path for sami outputs

    Returns:
        nx: num. grid points along each field line
        ny: num. field lines along each magnetic longitude
    """

    found = [False, False]

    nx = 0
    ny = 0

    # try:
    with open(os.path.join(sami_data_path, 'postrun-utils', 'TEC',
                           'param_diag.inc'), 'r') as fp:
        # read all lines in a list
        lines = fp.readlines()
        for line in lines:
            # check if string present on a current line
            if not found[0]:
                if line.find('nx') != -1:
                    nx0 = []
                    for char in line:
                        if char.isdigit():
                            nx0.append(char)
                    nx = int(''.join(nx0))
                    found[0] = True

            if not found[1]:
                if line.find('ny') != -1:
                    ny0 = []
                    for char in line:
                        if char.isdigit():
                            ny0.append(char)
                    ny = int(''.join(ny0))
                    found[1] = True
    if not all(found):
        print('found nx? %s and ny? %s' % (found[0], found[1]))
        raise ValueError('Could not find all of the required variables in '
                         'parameter_mod.f90 and time.dat. Please check that '
                         'the file is '
                         'formatted correctly.')
    else:
        return nx, ny


def make_times(nt, sami_data_path, dtime_sim_start,
               dtime_storm_start=None,
               hrs_before_storm=None, hrs_after_storm=None,
               need_help=False, skip_time_check=False):
    """
    Make a list of datetime objects for each time step from the time.dat file.

    Args:
        nt (int): Number of time steps (from get_grid_elems_from_parammod)
        sami_data_path (str): Path to sami data
        dtime_storm_start (datetime.datetime): Datetime of the start
            of the storm
        hrs_before_storm (int, optional):Hours from the onset
            of the storm (or any event, really) to begin processing.
            Set to -1 to run for the whole entire simulation.
            Defaults to None.
        hrs_after_storm (int, optional): Hours from the end of the
            storm to stop processing.
            Set to -1 to run for the whole entire simulation.
            Defaults to None.
        help (bool, optional): If help is set to true,
            we will print the time list.
            (useful when getting acquainted with the run)

    Raises:
        ValueError: Sometimes SAMI outputs fake time steps.
        ValueError: You only set one of hrs_before_storm or hrs_after_storm.

    Returns:
        (tuple) tuple containing:

            times (list):
                List of datetime objects for each time step
            hrs_since_storm_start (list):
                List of (float) hours since the storm start
            start_idx (int):
                Start index for the times list,
                    calculated from hrs_before_storm (ONLY if hrs_before_storm)
            end_idx (int):
                End index for the times list,
                    calculated from hrs_after_storm (ONLY if hrs_after_storm)

    """

    times = []
    if dtime_storm_start is not None:
        hrs_since_storm_start = []

    if not isinstance(dtime_sim_start, pd.Timestamp):
        dtime_sim_start = pd.Timestamp(dtime_sim_start)

    for t in range(nt):
        time_here = dtime_sim_start + t * pd.Timedelta(5, 'minutes')
        times.append(time_here.to_pydatetime())

        if dtime_storm_start is not None:
            hrs = (time_here - dtime_storm_start) / pd.Timedelta(1, 'hour')
            hrs_since_storm_start.append(hrs)

    times_df = pd.read_fwf(os.path.join(sami_data_path, 'time.dat'),
                           names=['istep', 'hour', 'minute', 'second',
                                  'hrdelta'],
                           infer_nrows=115)
    times_df.pop('istep')

    times_list = []
    for hr in times_df['hrdelta']:
        times_list.append(dtime_sim_start + datetime.timedelta(hours=hr))

    truths = np.array([pd.Timestamp(times_list[t]).round(
        'min') == times[t] for t in range(len(times))])
    if truths.sum() != len(truths) and not skip_time_check:
        raise ValueError(
            '\nThe times are wrong! Somehow this needs to be fixed.'
            ' SAMI may be outputting fake files again... \n'
            ' (this is a known issue - run with --skip_time_check'
            ' to ignore this error)')
    elif skip_time_check:
        times = times_list

    # maybe chop the time lists, depending on if the plot start/end are given.
    # adjusted to allow for -1 in plot start/end deltas (plot all times)

    if dtime_storm_start is None:
        return times

    if hrs_before_storm is not None and hrs_after_storm is not None:
        if dtime_storm_start is None:
            raise ValueError('cant call hrs from storm without storm info')

        if hrs_before_storm is not None:
            start_idx = np.argmin(np.abs(
                np.array(times)
                - (dtime_storm_start
                   - pd.Timedelta(hrs_before_storm, 'hour'))))

        else:
            start_idx = 0

        if hrs_after_storm is not None:
            end_idx = np.argmin(np.abs(
                np.array(times)
                - (dtime_storm_start
                   + pd.Timedelta(hrs_after_storm, 'hour'))))
        else:
            end_idx = len(times)

        times = times[start_idx:end_idx]
        hrs_since_storm_start = hrs_since_storm_start[start_idx:end_idx]
        times_list = times_list[start_idx:end_idx]

    elif hrs_before_storm is not None or hrs_after_storm is not None:
        raise ValueError('You cannot specify one and not the other!')

    if need_help:
        print(times, '\n', hrs_since_storm_start)

    if dtime_storm_start is None:
        return times

    else:
        return times, hrs_since_storm_start, (start_idx, end_idx)


def get_sami_grid(sami_data_path, nlt, nf, nz):
    """Read in SAMI grid files.

    Args:
        sami_data_path (str): path to SAMI data
        nlt (int):
            Number of magnetic local times (lons)
        nf (int):
            Number of field lines along each longitude
        nz (int):
            Number of grid cells along each field line
        geo_grid_files (dict, optional):
            Files to use for getting the grid. Defaults to {
            'glat': 'glatu.dat', 'glon': 'glonu.dat',
            'alt': 'zaltu.dat', 'mlat': 'blatu.dat',
            'mlon': 'blonu.dat', 'malt': 'baltu.dat'}.

    Returns:
        dict:
            SAMI3 grid in a dictionary with keys:
            'glat', 'glon', 'alt', 'mlat', 'mlon', 'malt'
    """
    geo_grid_files = {
        'glat': 'glatu.dat', 'glon': 'glonu.dat', 'alt': 'zaltu.dat',
        'mlat': 'blatu.dat', 'mlon': 'blonu.dat', 'malt': 'baltu.dat'
    }

    grid = {}

    for f in geo_grid_files:
        file = open(os.path.join(sami_data_path, geo_grid_files[f]), 'rb')
        raw = np.fromfile(file, dtype='float32')[1:-1].copy()
        file.close()

        grid[f] = raw.reshape(nlt, nf, nz).copy()
    return grid


def read_to_nparray(sami_data_path, dtime_sim_start,
                    dtime_storm_start=None,
                    hrs_before_storm=None, hrs_after_storm=None, pbar=False,
                    cols='all', need_help=False,
                    skip_time_check=False):
    """Automatically read in SAMI data.

    Args:
        sami_data_path (str):
            Path to SAMI data
        dtime_storm_start (datetime.datetime):
            Datetime of the start of the storm
        dtime_sim_start (datetime.datetime-like or str):
            Datetime of the start of the simulation
        t_start_idx (int, optional):
            Time index of the start of the data return. Defaults to None.
        t_end_idx (int, optional):
            Time index of the end of the data return. Defaults to None.
        pbar (bool, optional):
            Do you want to show a progress bar? It is automatically set if
            tqdm is successfully imported. Defaults to False.
        cols (str or list-like, optional):
            List of columns to get data for. Defaults to 'all'.
        help (bool, optional):
            Prints time and variable info. Defaults to False.

    Raises:
        KeyError: If given cols is not valid
        FileNotFoundError: If the filepath is invalid

    Returns:
        dict:
            Dictionary of SAMI data with keys: ['grid', 'data']
            data is in np arrays with the shape [nlt,nf,nz]
        np.array:
            Times of the data

    """

    sami_data = {}

    # Check cols:
    if cols == 'all':
        cols = sami_og_vars.values()
    else:
        if isinstance(cols, str):
            cols = [cols]
    data_files = {}
    for ftype in sami_og_vars:
        if sami_og_vars[ftype] in cols:
            data_files[ftype] = sami_og_vars[ftype]

    if need_help or len(data_files.keys()) == 0:
        print('the available data files are: \n', sami_og_vars.keys())

        return

    # Get the grid
    nz, nf, nlt, nt = get_grid_elems_from_parammod(sami_data_path)

    grid = get_sami_grid(sami_data_path, nlt, nf, nz)
    sami_data['grid'] = grid

    if dtime_storm_start is not None:
        times, hrs_since_storm_start, (start_idx, end_idx) = make_times(
            nt, sami_data_path, dtime_sim_start, dtime_storm_start,
            hrs_before_storm, hrs_after_storm,
            skip_time_check=skip_time_check)
    else:
        times = make_times(nt, sami_data_path, dtime_sim_start,
                           skip_time_check=skip_time_check)
        start_idx = 0
        end_idx = len(times)

    ntimes = len(times)

    if need_help:
        print('nz, nlt, nf, ntimes = ', nz, nlt, nf, ntimes)
        return

    if pbar:
        progress = tqdm(total=len(cols) * ntimes, desc='reading SAMI data')

    sami_data['data'] = {}
    for f in data_files:
        sami_data['data'][data_files[f]] = np.zeros((nlt, nf, nz, ntimes))

        try:
            file = open(os.path.join(sami_data_path, f), 'rb')
        except KeyError:
            print('the column name you entered is not valid. \n',
                  'the available columns are: \n', data_files.keys())
            raise
        except FileNotFoundError:
            print('the file you entered is not valid. \n',
                  'the model results may not be in the path you specified:')
            raise

        for t in range(nt):
            raw = np.fromfile(file, dtype='float32',
                              count=(nz * nf * nlt) + 2)[1:-1]
            if t > start_idx and t < end_idx:
                sami_data['data'][data_files[f]][:, :, :, t - start_idx] = \
                    raw.reshape(nlt, nf, nz).copy()
            if pbar:
                progress.update(1)
        file.close()
    if pbar:
        progress.close()

    return sami_data, np.array(times)


def read_sami_dene_tec_MAG_GRID(sami_data_path,
                                dtime_sim_start=None,
                                reshape=True):
    """ Read in TEC (and interpolated dene) data!

    :param sami_data_path: path to SAMI data
    :type sami_data_path: str
    :param dtime_sim_start: datetime of the start of the simulation
    :type dtime_sim_start: datetime.datetime
    :param reshape: reshape the data to the correct shape, defaults to True.
        Otherwise, the data will be returned as a 1D array.
    :type reshape: bool, optional
    :return: SAMI data, times
    :rtype: dict, np.array

    """
    # TODO: Add in all of the data files. This is just a placeholder.
    # TODO: remove hard-coding shapes.
    data_files = {'edens': 'dene0B.dat', 'tec': 'tecuB.dat'}

    sami_data = {'grid': {}, 'data': {}}

    # Get the grid
    geo_grid_files = {
        'glat': 'glat0B.dat', 'glon': 'glon0B.dat', 'alt': 'zalt0B.dat',
        'mlat': 'blat0.dat', 'mlon': 'blon0.dat', 'malt': 'balt0.dat'}

    for f in geo_grid_files:
        try:
            file = open(os.path.join(sami_data_path, geo_grid_files[f]), 'rb')
            raw = np.fromfile(file, dtype='float32')[1:-1].copy()
            file.close()
        except FileNotFoundError:
            print("No TEC/BtoG files found. Make sure:")
            print(geo_grid_files.keys())
            print("exist in %s directory." % (str(sami_data_path)))
            print("hint: you may need to run post-processing scripts")

        sami_data['grid'][f] = raw

    for f in data_files:
        file = open(os.path.join(sami_data_path, data_files[f]), 'rb')
        raw = np.fromfile(file, dtype='float32')[1:-1].copy()
        file.close()

        sami_data['data'][f] = raw

    nz, nf, nlt, nt = get_grid_elems_from_parammod(sami_data_path)

    nt += 1

    # defaults
    nx, ny = get_postprocessed_grid(sami_data_path)

    # Reshape everything!
    if reshape:
        sami_data['data']['edens'] = sami_data['data']['edens'].reshape(
            nt, nlt, nx, ny)
        sami_data['data']['tec'] = sami_data['data']['tec'].reshape(
            nt, nlt, nx)
        sami_data['grid']['glat'] = sami_data['grid']['glat'].reshape(
            nlt, nx, ny)
        sami_data['grid']['glon'] = sami_data['grid']['glon'].reshape(
            nlt, nx, ny)

    with open(os.path.join(sami_data_path, 'time.dat'), 'r') as fp:
        lines = fp.readlines()
        nt = len(lines) - 1

    times = make_times(
        nt, sami_data_path, dtime_sim_start=dtime_sim_start)

    return sami_data, np.array(times)


def read_raw_to_xarray(sami_data_path, dtime_sim_start, cols='all',
                       hrs_before_storm_start=None,
                       hrs_after_storm_start=None,
                       dtime_storm_start=None,
                       start_dtime=None, end_dtime=None,
                       start_idx=None, end_idx=None,
                       progress_bar=False, skip_time_check=False):

    """
    Read in (raw) SAMI data and return an xarray dataset.
        ! This only works on raw, pre-processed SAMI data !
        (not TEC or anything like that)

    Args:
        sami_data_path (str- path-like): Directory of SAMI files.
        dtime_sim_start (datetime): Start time of simulation.
        cols (str or list-like, optional): Model outputs to read.
            Defaults to 'all'.
        hrs_before_storm_start (int, optional): Hours before storm onset
            to read data from. Need to set dtime_storm_start. Defaults to None.
        hrs_after_storm_start (int, optional): Hours after storm onset
            to read data from. Need to set dtime_storm_start. Defaults to None.
        dtime_storm_start (datetime, optional): storm/event start time.
            Only used if hrs_before/after is set. Defaults to None.
        start_dtime (datetime, optional): datetime to start reading data.
            Defaults to None.
        end_dtime (datetime, optional): datetime to stop reading data.
            Defaults to None.
        start_idx (int, optional): Index of time list to start.
            Defaults to None.
        end_idx (int, optional): Index of time list to stop.
            Defaults to None.
        progress_bar (bool, optional): Show progress bar. Defaults to False.
            (Requires tqdm)

    Raises:
        ValueError: Invalid inputs
        ValueError: Missing Files

    Returns:
        xarray.Dataset: Dataset of SAMI data.
    """

    import xarray as xr

    nz, nf, nlt, nt = get_grid_elems_from_parammod(sami_data_path)
    times = make_times(nt, sami_data_path, dtime_sim_start,
                       skip_time_check=skip_time_check)
    times = np.array(times)
    grid = get_sami_grid(sami_data_path, nlt, nf, nz)

    if start_idx is None and end_idx is None:
        if dtime_storm_start is not None:
            if hrs_before_storm_start is not None:
                start_idx = np.argmin(np.abs(
                    times - (dtime_storm_start - pd.Timedelta(
                        np.abs(hrs_before_storm_start), 'hours'))))
            if hrs_after_storm_start is not None:
                end_idx = np.argmin(
                    np.abs(times - (dtime_storm_start + pd.Timedelta(
                        np.abs(hrs_before_storm_start), 'hours'))))
            if start_idx is None and end_idx is None:
                raise ValueError(
                    'You must specify either start_idx or end_idx',
                    'to use dtime_storm_start')
        elif hrs_after_storm_start is not None or\
                hrs_before_storm_start is not None:
            raise ValueError(
                'why did you give storm start time but no time args?')
        elif start_dtime is not None or end_dtime is not None:
            if end_dtime is not None:
                start_idx = np.argmin(np.abs(times - start_dtime))
            if end_dtime is not None:
                end_idx = np.argmin(np.abs(times - end_dtime))

        if start_idx is None:
            start_idx = 0
        if end_idx is None:
            end_idx = nt
    nt = end_idx - start_idx

    ds = xr.Dataset(
        coords=dict(
            time=(('time'), times),
            mlat=(('nlt', 'nf', 'nz'), grid['mlat'].round(4)),
            mlon=(('nlt', 'nf', 'nz'), grid['mlon'].round(2)),
            malt=(('nlt', 'nf', 'nz'), grid['malt'].round(4)),
            glat=(('nlt', 'nf', 'nz'), grid['glat'].round(4)),
            glon=(('nlt', 'nf', 'nz'), grid['glon'].round(4)),
            alt=(('nlt', 'nf', 'nz'), grid['alt']),),
    )

    if dtime_storm_start is not None:
        ds = ds.assign_attrs(
            dtime_event_start=dtime_storm_start,)

    dimnames = ('time', 'nlt', 'nf', 'nz')

    if progress_bar:
        pbar1 = tqdm(total=len(cols), desc='Reading SAMI binaries')

    if cols != 'all':
        if isinstance(cols, str):
            cols = [cols]

    none_found = True
    for fname in sami_og_vars:
        if cols != 'all':
            if sami_og_vars[fname] not in cols:
                continue
            else:
                none_found = False
        curr_arr = np.zeros((len(times), nlt, nf, nz))
        try:
            with open(os.path.join(sami_data_path, fname), 'rb') as f:
                t_ins = 0
                for t in range(len(times)):
                    raw = np.fromfile(f, dtype='float32',
                                      count=(nz * nf * nlt) + 2)[1:-1]
                    if t > start_idx and t <= end_idx:
                        curr_arr[t_ins] = raw.reshape(nlt, nf, nz).copy()
                        t_ins += 1

                ds[sami_og_vars[fname]] = (dimnames, curr_arr)

                if progress_bar:
                    pbar1.update()
        except FileNotFoundError:
            print(f'File {fname} not found!\n',
                  'in directory {sami_data_path}')
    if none_found:
        raise KeyError('Key not found',
                       'Please double check key names.',
                       f'You gave: {cols} But only {sami_og_vars.values()} are valid')

    return ds


def process_all_to_cdf(sami_data_path,
                       dtime_sim_start,
                       dtime_storm_start=None,
                       progress_bar=False,
                       start_dtime=None,
                       end_dtime=None,
                       out_dir=None,
                       use_ccmc=True,
                       split_by_time=True,
                       split_by_var=False,
                       whole_run=False,
                       run_name=None,
                       OVERWRITE=False,
                       delete_raw=False,
                       append_files=False,
                       low_mem=True,
                       cols='all',
                       skip_time_check=False
                       ):
    """Process SAMI binary files to netcdf format.

    Args:
        sami_data_path (str): Path to SAMI data.
        dtime_sim_start (datetime): Simulation start time.
        progress_bar (bool, optional): Show progress bar. Defaults to False.
            Requires tqdm
        start_dtime (datetime, optional): datetime to start reading data.
            Defaults to None.
        end_dtime (datetime, optional): datetime to stop reading data.
            Defaults to None.
        out_dir (str, optional): Directory to save netcdf files.
            Defaults to sami_data_path.
        split_by_time (bool, optional): Split files by time. Defaults to False.
        split_by_var (bool, optional): Split files by variable.
            Defaults to False.
        whole_run (bool, optional): Save whole model run (in time range)
            as one netcdf. Defaults to False.
        OVERWRITE (bool, optional): Overwrite existing files.
            Defaults to False.
        append_files (bool, optional): Append to existing files.
        low_mem (bool, optional): Read data in chunks to save memory.
            Defaults to False.
        cols (list-like or str, optional): List of columns to read.
            Defaults to 'all'.

    Raises:
        ValueError: If incorrect time args are given.
        ValueError: If files exist and OVERWRITE is False.
        ValueError: If cols is not in available columns.

    returns:
        None

    """

    if out_dir is None:
        out_dir = sami_data_path

    if whole_run and run_name is None:
        raise ValueError('You must set the run name if outputting'
                         ' to a single file')

    if low_mem:
        if cols != 'all':
            if isinstance(cols, str):
                cols = [cols]
            for c in cols:
                if c not in sami_og_vars.values():
                    raise ValueError(
                        f'Column {c} not in available columns!\n',
                        ' available columns are:\n',
                        ' '.join(sami_og_vars.values()))
        else:
            cols = sami_og_vars.values()

        if progress_bar:
            total = len(cols) * np.sum(
                [split_by_time, split_by_var, whole_run])
            print(
                'Processing in lowmem mode. This will take substantially',
                'longer than "normal" mode.')
            pbar = tqdm(total=total,
                        desc='Variable loop')

        did_one = False

        if split_by_time:
            nz, nf, nlt, nt = get_grid_elems_from_parammod(
                sami_data_path)
            times = make_times(nt, sami_data_path, dtime_sim_start,
                               skip_time_check=skip_time_check)
            # So this gets complicated. First check the files.
            # Then go ahead and process...

            file_list = []

            for t in times:
                if use_ccmc:
                    fname = make_ccmc_name('SAMI', t, data_type='RAW')
                else:
                    fname = 't' + t.strftime('%y%m%d_%H%M%S') + '.nc'

                out_file = os.path.join(out_dir, fname)
                file_list.append(out_file)
                if os.path.exists(out_file):
                    if OVERWRITE:
                        os.remove(out_file)
                    else:
                        raise FileExistsError(
                            out_file,
                            ' already exists! Cannot rewrite!')

            if progress_bar:
                pbar2 = tqdm(total=nt * len(cols), desc='Vars & Times')

            for ftype in sami_og_vars:
                ds = read_raw_to_xarray(
                    sami_data_path=sami_data_path,
                    dtime_sim_start=dtime_sim_start,
                    cols=sami_og_vars[ftype],
                    start_dtime=start_dtime,
                    end_dtime=end_dtime,
                    skip_time_check=skip_time_check)
                for nfile in range(len(times)):
                    try:
                        ds.isel(time=nfile).to_netcdf(
                            file_list[nfile],
                            mode='a')
                    except FileNotFoundError:
                        ds.isel(time=nfile).to_netcdf(
                            file_list[nfile])
                    pbar2.update()
                pbar.update()
                did_one = True

        if split_by_var or whole_run:

            for ftype in sami_og_vars:
                did_var = False
                if sami_og_vars[ftype] in cols:

                    ds = read_raw_to_xarray(sami_data_path, dtime_sim_start,
                                            cols=sami_og_vars[ftype],
                                            skip_time_check=skip_time_check)

                    if start_dtime is not None or end_dtime is not None:
                        if start_dtime is not None:
                            start_idx = np.argmin(
                                np.abs(ds.time.values - start_dtime))
                        else:
                            start_idx = 0
                        if end_dtime is not None:
                            end_idx = np.argmin(
                                np.abs(ds.time.values - end_dtime))
                        else:
                            end_idx = len(ds.time)
                        ds = ds.isel(time=slice(start_idx, end_idx))

                    if split_by_var:
                        out_file = os.path.join(
                            out_dir, f'{sami_og_vars[ftype]}.nc')
                        ds.to_netcdf(out_file, mode='a')
                        did_one = True
                        did_var = True
                    if whole_run:
                        try:
                            ds.to_netcdf(
                                os.path.join(
                                    out_dir,
                                    run_name +
                                    '_SAMI-RAW.nc'),
                                mode='a',
                                engine='h5netcdf')
                        except FileNotFoundError:
                            ds.to_netcdf(
                                os.path.join(
                                    out_dir,
                                    run_name +
                                    '_SAMI-RAW.nc'))
                        did_one = True
                        did_var = True

                    if did_var and progress_bar:
                        pbar.update()

        if not did_one:
            raise ValueError('Your columns were not found in the data!')
        elif delete_raw:
            for ftype in sami_og_vars:
                if sami_og_vars[ftype] in cols:
                    os.remove(os.path.join(sami_data_path,
                                           sami_og_vars[ftype]))

    else:

        if delete_raw:
            raise NotImplementedError(
                'delete_raw not implemented for non-low_mem mode!')

        ds = read_raw_to_xarray(sami_data_path, dtime_sim_start,
                                progress_bar=progress_bar, cols=cols)

        if start_dtime is not None or end_dtime is not None:
            if start_dtime is not None:
                start_idx = np.argmin(np.abs(ds.time.values - start_dtime))
            else:
                start_idx = 0

            if end_dtime is not None:
                end_idx = np.argmin(np.abs(ds.time.values - end_dtime))
            else:
                end_idx = len(ds.time)

            ds = ds.isel(time=slice(start_idx, end_idx))

        write_mode = 'w'

        if split_by_time:
            if progress_bar:
                pbar = tqdm(total=len(ds.time.values),
                            desc='Saving to netcdf')
            for t in ds.time.values:

                if use_ccmc:
                    fname = make_ccmc_name('SAMI', t, data_type='RAW')
                else:
                    fname = 't' + \
                        pd.Timestamp(t).strftime('%y%m%d_%H%M%S') + '.nc'
                if os.path.exists(os.path.join(out_dir, fname)):
                    if OVERWRITE and not append_files:
                        os.remove(os.path.join(out_dir, fname))
                    elif append_files:
                        write_mode = 'a'
                    else:
                        raise FileExistsError('%s already exists!' % fname)

                ds.sel(time=t).to_netcdf(os.path.join(out_dir, fname),
                                         mode=write_mode)

                if progress_bar:
                    pbar.update()

        elif split_by_var:
            if progress_bar:
                pbar = tqdm(total=len(ds.time.values),
                            desc='Saving to netcdf')
            for var in ds.data_vars:
                if os.path.exists(os.path.join(out_dir, f'{var}.nc')):
                    if OVERWRITE and not append_files:
                        os.remove(os.path.join(out_dir, f'{var}.nc'))
                    elif append_files:
                        write_mode = 'a'
                    else:
                        raise FileExistsError(f'{var}.nc already exists!')

                ds[var].to_netcdf(os.path.join(out_dir, f'{var}.nc',
                                               mode=write_mode))
                if progress_bar:
                    pbar.update()

        else:
            ds.to_netcdf(os.path.join(out_dir,
                                      run_name + '_SAMI-RAW.nc'))


def auto_read(sami_dir,
              cols='all',
              split_by_time=False,
              split_by_var=False,
              whole_run=False,
              return_xarray=True,
              filetype='SAMI-REGRID',
              force_nparrays=False,
              dtime_sim_start=None,
              parallel=True,
              start_dtime=None,
              start_idx=None,
              end_dtime=None,
              end_idx=None,
              hrs_before_storm_start=None,
              hrs_after_storm_start=None,
              dtime_storm_start=None,
              progress_bar=False,
              use_dask=False,
              engine='h5netcdf',
              skip_time_check=False,
              ):
    """

    Automatically reads in SAMI data and returns it in a format of your choice.

    - Preference is to read/return xarray datasets, but can read
        and return numpy arrays.
    - Prefer whole files, fall back on time, variable, split.


    Args:
        sami_dir (str: path-like): Path to the directory containing the SAMI
            data
        cols (str or list-like, optional): Variables to return.
            Defaults to 'all'.
        split_by_time (bool, optional): If files are output by time
            (and whether to prefer those files). Defaults to False.
        split_by_var (bool, optional): If files are output by variable.
            And to prefer those files. Defaults to False.
        whole_run (bool, optional): If the whole run is in one file.
            Defaults to False.
        return_xarray (bool, optional): Return xarray dataset?
            Defaults to True.
        force_nparrays (bool, optional): Force program to return dicts of numpy
            arrays. Defaults to False.
        dtime_sim_start (datetime, optional): Datetime of the start of the
            simulation. Defaults to None. Required if netCDF files aren't made.
            If netCDF files are made, You had the option to add
            this as an attribute to the file.
        parallel (bool, optional): Force parallel reading of files.
            NetCDF files are read weird. Might be buggy. Defaults to True.
        start_dtime (datetime, optional): Datetime to start reading data.
            Defaults to None.
        start_idx (int, optional): Index of the first time to read.
            Defaults to None.
        end_dtime (datetime, optional): Datetime to stop reading data.
            Defaults to None
        end_idx (int, optional): Index of the last time to read.
            Defaults to None.
        hrs_before_storm_start (int, optional): Hours before the storm start
            to read data. Defaults to None. (dtime_storm_start must be set)
        hrs_after_storm_start (int, optional): Hours after the storm start
            to read data. Defaults to None. (dtime_storm_start must be set)
        dtime_storm_start (datetime, optional): Datetime of the storm start.
            Defaults to None. (hrs_before_storm_start and hrs_after_storm_start
            must be set)
        progress_bar (bool, optional): Show progress bar? Defaults to False.
            Requires tqdm.

    Returns:
        xarray Dataset: Dataset of the SAMI data
            (If return_xarray is True)
        dict: Dictionary of numpy arrays of the SAMI data
            (If force_nparrays is False)
    """

    from glob import glob

    ncfiles = glob(os.path.join(sami_dir, '*.nc'))
    if len(ncfiles) > 0:
        if len(glob(os.path.join(sami_dir, 'SAMI*T*.nc'))) > 0:
            split_by_time = True
        if len(glob(os.path.join(sami_dir, 'sami_data.nc'))) > 0:
            whole_run = True
        if cols != 'all':
            if isinstance(cols, str):
                cols = [cols]
            multivars = []
            for col in cols:
                p = os.path.join(sami_dir, f'{col}.nc')
                if os.path.exists(p):
                    split_by_var = True
                    multivars.append(p)

        if not force_nparrays:
            import xarray as xr

            if whole_run:
                ds = xr.open_dataset(os.path.join(sami_dir, 'sami_data.nc'),
                                     chunks='auto')

            elif split_by_var:
                ds = xr.open_mfdataset(os.path.join(sami_dir, '*.nc'),
                                       parallel=parallel,
                                       combine_attrs='drop_conflicts',
                                       data_vars='minimal',
                                       concat_dim="time", combine="nested",
                                       coords='minimal', compat='override')

            elif split_by_time:
                files = np.sort(
                    glob(
                        os.path.join(
                            sami_dir,
                            filetype +
                            '*.nc')))
                ret_early = False
                if start_idx is None:
                    start_idx = 0
                    ret_early = True
                if end_idx is None:
                    end_idx = -1
                    ret_early = True

                files = files[start_idx:end_idx]

                if len(files) > 1:
                    if use_dask:
                        ds = xr.open_mfdataset(
                            files,
                            parallel=parallel,
                            data_vars=cols,
                            combine_attrs='drop_conflicts',
                            engine=engine,
                            concat_dim="time",
                            combine="nested",
                            coords='minimal',
                            compat='override')
                    else:
                        drops = []
                        ds0 = xr.open_dataset(files[0])
                        for v in ds0.data_vars:
                            if v not in cols:
                                drops.append(v)
                        del ds0
                        dss = []
                        for f in files:
                            dss.append(xr.open_dataset(
                                f, drop_variables=drops, engine=engine))

                        ds = xr.concat(dss, dim='time')
                        del dss
                else:
                    ds = xr.open_dataset(files[0])
                if ret_early:
                    return ds

            else:
                print('Something went wrong with naming netcdf files.\n',
                      'Switching to nparray read')

            if ds is not None:
                # if cols != 'all':
                # if type(cols) is str:
                # cols = [cols]
                # ds = ds[cols]

                if start_dtime is not None or end_dtime is not None:
                    if start_dtime is not None:
                        start_idx = np.argmin(
                            np.abs(ds.time.values - start_dtime))
                    if end_dtime is not None:
                        end_idx = np.argmin(np.abs(ds.time.values - end_dtime))

                # if start_idx is not None or end_idx is not None:
                    # if start_idx is None:
                        # start_idx = 0
                    # if end_idx is None:
                        # end_idx = len(ds.time.values)
                    # ds = ds.isel(time=slice(start_idx, end_idx))
                return ds

        else:
            print('found netcdf files but forcing nparray read')

    if return_xarray and not force_nparrays:
        if dtime_sim_start is None:
            raise ValueError(
                'dtime_sim_start must be set if not reading netcdf files')
        ds = read_raw_to_xarray(sami_dir, dtime_sim_start,
                                progress_bar=progress_bar, cols=cols,
                                hrs_before_storm_start=hrs_before_storm_start,
                                hrs_after_storm_start=hrs_after_storm_start,
                                dtime_storm_start=dtime_storm_start,
                                start_dtime=start_dtime,
                                end_dtime=end_dtime, start_idx=start_idx,
                                end_idx=end_idx,
                                skip_time_check=skip_time_check)
        return ds

    else:
        sami_data = read_to_nparray(sami_dir, dtime_sim_start,
                                    dtime_storm_start=dtime_storm_start,
                                    t_end_idx=start_idx,
                                    t_start_idx=end_idx,
                                    pbar=progress_bar, cols=cols,
                                    skip_time_check=skip_time_check)

        return sami_data
