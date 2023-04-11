
"""read sami data.

call:

read_sami_data(sami_data_path, dtime_sim_start, dtime_storm_start,
                   t_start_idx=None, t_end_idx=None, pbar=False,
                   cols='all', help=False, chop_times=False

 Automatically read in SAMI data.

    Args:
        sami_data_path (str):
            Path to SAMI data
        dtime_storm_start (datetime.datetime):
            Datetime of the start of the storm
        dtime_sim_start (datetime.datetime):
            Datetime of the start of the simulation
        t_start_idx (int, optional):
            Time index of the start of the data return. Defaults to None.
        t_end_idx (int, optional):
            Time index of the end of the data return. Defaults to None.
        pbar (bool, optional):
            Do you want to show a progress bar? It is automatically set if
            tqdm is successfully imported. Defaults to False.
        cols (str, optional):
            List of columns to get data for. Empty is all. Defaults to 'all'.
        help (bool, optional):
            Prints time and variable info. Defaults to False.

    Raises:
        ValueError: if t_start_idx and t_end_idx are not both given

    Returns:
        dict:
            Dictionary of SAMI data with keys: ['grid', 'data']
            data is in np arrays with the shape [nlt,nf,nz]
        np.array:
            Times of the data
"""


import datetime
import os

import numpy as np
import pandas as pd

try:
    from tqdm import tqdm
    pbar = True
except ImportError:
    pbar = False


def get_grid_elems_from_parammod(sami_data_path):
    """Go into sami data directory and get the grid elements
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
        return nz, nf, numwork*(nl - 2), nt


def make_times(nt, sami_data_path,  dtime_sim_start, dtime_storm_start=None,
               plot_start_delta=None, plot_end_delta=None, help=False):
    """_summary_

    Args:
        nt (int):
            Number of time steps
        sami_data_path (str):
            Path to sami data
        dtime_storm_start (datetime.datetime):
            Datetime of the start of the storm
        plot_start_delta (int, optional):
            Hours from the onset of the storm to begin processing.
            Set to -1 to run for the whole entire simulation.
            Defaults to None.
        plot_end_delta (int, optional): Hours from the end of the
            storm to stop processing.
            Set to -1 to run for the whole entire simulation.
            Defaults to None.
        help (bool, optional):
            If help is set to true, we will print the time list.
            (useful when getting acquainted with the run)

    Raises:
        ValueError: Sometimes SAMI outputs fake time steps.
        ValueError: You only set one of plot_start_delta or plot_end_delta.

    Returns:
        times (list):
            List of datetime objects for each time step
        hrs_since_storm_start (list):
            List of (float) hours since the storm start
        start_idx (int, optional):
            Start index for the times list, calculated from plot_start_delta
        end_idx (int, optional):
            End index for the times list, calculated from plot_end_delta

    """

    times = []
    hrs_since_storm_start = []

    for t in range(nt):
        time_here = pd.Timestamp(dtime_sim_start) + \
            t * pd.Timedelta(5, 'minutes')
        times.append(time_here.to_pydatetime())
        if dtime_storm_start:
            hrs = (time_here - dtime_storm_start)/pd.Timedelta(1, 'hour')
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
        'T') == times[t] for t in range(len(times))])
    if truths.sum() != len(truths):
        raise ValueError(
            'The times are wrong! Somehow this needs to be fixed.'
            'probably outputting fake files again... \n')

    # maybe chop the time lists, depending on if the plot start/end are given.
    # adjusted to allow for -1 in plot start/end deltas (plot all times)

    if dtime_storm_start is None:
        return times

    if plot_start_delta and plot_end_delta:
        if plot_start_delta != -1:
            start_idx = np.argmin(np.abs(
                np.array(times)
                - (dtime_storm_start
                   - pd.Timedelta(plot_start_delta, 'hour'))))
        else:
            start_idx = 0

        if plot_end_delta != -1:
            end_idx = np.argmin(np.abs(
                np.array(times)
                - (dtime_storm_start
                   + pd.Timedelta(plot_end_delta, 'hour'))))
        elif plot_end_delta == -1:
            end_idx = len(times)
        else:
            end_idx = len(times)
        times = times[start_idx:end_idx]
        hrs_since_storm_start = hrs_since_storm_start[start_idx:end_idx]
        times_list = times_list[start_idx:end_idx]

    elif plot_start_delta != plot_end_delta:
        raise ValueError('You cannot specify one and not the other!')
    else:
        start_idx = 0
        end_idx = len(times)
    if help:
        print(times, '\n', hrs_since_storm_start)

    if plot_start_delta is None and plot_end_delta is None:
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
            TODO: Maybe move this to a regular definition?

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


def read_sami_data(sami_data_path, dtime_sim_start,
                   dtime_storm_start=None,
                   t_start_idx=None, t_end_idx=None, pbar=False,
                   cols='all', help=False):
    """Automatically read in SAMI data.

    Args:
        sami_data_path (str):
            Path to SAMI data
        dtime_storm_start (datetime.datetime):
            Datetime of the start of the storm
        dtime_sim_start (datetime.datetime):
            Datetime of the start of the simulation
        t_start_idx (int, optional):
            Time index of the start of the data return. Defaults to None.
        t_end_idx (int, optional):
            Time index of the end of the data return. Defaults to None.
        pbar (bool, optional):
            Do you want to show a progress bar? It is automatically set if
            tqdm is successfully imported. Defaults to False.
        cols (str, optional):
            List of columns to get data for. Empty is all. Defaults to 'all'.
        help (bool, optional):
            Prints time and variable info. Defaults to False.

    Raises:
        ValueError: if t_start_idx and t_end_idx are not both given

    Returns:
        dict:
            Dictionary of SAMI data with keys: ['grid', 'data']
            data is in np arrays with the shape [nlt,nf,nz]
        np.array:
            Times of the data

    """

    # TODO: Add in all of the data files. This is just a placeholder
    data_files = {'edens': 'deneu.dat', 'hplusdens': 'deni1u.dat',
                  'oplusdens': 'deni2u.dat', 'noplusdens': 'deni3u.dat',
                  'o2plusdens': 'deni4u.dat', 'heplusdens': 'deni5u.dat',
                  'n2plusdens': 'deni6u.dat', 'nplusdens': 'deni7u.dat',
                  'hdens': 'denn1u.dat', 'odens': 'denn2u.dat',
                  'nodens': 'denn3u.dat', 'o2dens': 'denn4u.dat',
                  'hedens': 'denn5u.dat', 'n2dens': 'denn6u.dat',
                  'ndens': 'denn7u.dat'}

    sami_data = {}

    # Check cols:
    if cols == 'all':
        cols = data_files.keys()
    if help:
        print('the available columns are: \n', data_files.keys())

    # Get the grid
    nz, nf, nlt, nt = get_grid_elems_from_parammod(sami_data_path)

    grid = get_sami_grid(sami_data_path, nlt, nf, nz)
    sami_data['grid'] = grid

    if dtime_storm_start is None:
        times = make_times(nt, sami_data_path)

    elif t_start_idx is not None or t_end_idx is not None:
        times, hrs_since_storm_start, (start_idx, end_idx) = make_times(
            nt, sami_data_path, dtime_sim_start, dtime_storm_start,
            t_start_idx, t_end_idx, help)
    else:
        times, hrs_since_storm_start = make_times(
            nt, sami_data_path, dtime_storm_start, dtime_sim_start)

    ntimes = len(times)

    if help:
        print('nz, nlt, nf, ntimes = ', nz, nlt, nf, ntimes)
        return

    if pbar:
        progress = tqdm(total=len(cols) * nt, desc='reading SAMI data')

    sami_data['data'] = {}
    for f in cols:
        sami_data['data'][f] = np.zeros((nlt, nf, nz, ntimes))

        try:
            file = open(os.path.join(sami_data_path, data_files[f]), 'rb')
        except KeyError:
            print('the column name you entered is not valid. \n',
                  'the available columns are: \n', data_files.keys())
            raise
        except FileNotFoundError:
            print('the file you entered is not valid. \n',
                  'the model results may not be in the path you specified:')
            raise

        for t in range(nt):
            raw = np.fromfile(file, dtype='float32', count=(nz*nf*nlt)+2)[1:-1]
            if t in range(start_idx, end_idx):
                sami_data['data'][f][:, :, :, t-start_idx] = raw.reshape(
                    nlt, nf, nz).copy()
            if pbar:
                progress.update(1)
        file.close()
    if pbar:
        progress.close()

    print(sami_data['data']['edens'].shape, len(times), start_idx, end_idx)

    return sami_data, np.array(times)


def read_sami_dene_tec(sami_data_path, dtime_sim_start, 
                       dtime_storm_start=None,
                       reshape=True):
    """ Read in TEC (and interpolated dene) data!

    """
    # TODO: Add in all of the data files. This is just a placeholder
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
            print("exist in %s directory." %(str(sami_data_path)))
            print("hint: you may need to run post-processing scripts")

        sami_data['grid'][f] = raw

    for f in data_files:
        file = open(os.path.join(sami_data_path, data_files[f]), 'rb')
        raw = np.fromfile(file, dtype='float32')[1:-1].copy()
        file.close()

        sami_data['data'][f] = raw

    # Reshape everything!
    if reshape:
        sami_data['data']['edens'] = sami_data['data']['edens'].reshape(
            625, 80, 100, 100)
        sami_data['data']['tec'] = sami_data['data']['tec'].reshape(
            625, 80, 100)
        sami_data['grid']['glat'] = sami_data['grid']['glat'].reshape(
            80, 100, 100)
        sami_data['grid']['glon'] = sami_data['grid']['glon'].reshape(
            80, 100, 100)

    with open(os.path.join(sami_data_path, 'time.dat'), 'r') as fp:
        lines = fp.readlines()
        nt = len(lines) - 1

    times = make_times(
        nt, sami_data_path, dtime_sim_start=dtime_sim_start)

    return sami_data, np.array(times)