# Misc Utility Functions
import numpy as np
from datetime import datetime
import pandas as pd
import xarray as xr


def str_to_ut(in_str):
    """Convert a string to a datetime object.

    :param in_str: String to convert
    :type in_str: str
    :return: Datetime object
    :rtype: datetime
    """

    return datetime.strptime(
        in_str.ljust(14, '0'), '%Y%m%d%H%M%S')


def make_ccmc_name(modelname, ut, data_type=None):
    """
    Make a CCMC-formatted filename.

    :param modelname: Name of model.
    :type modelname: str
    :param ut: Datetime object of the time of the file.
    :type ut: datetime
    :param data_type: Type of data in the file (3DALL/2DANC/REGRID/etc.).
    :type data_type: str, optional
    :return: CCMC-formatted filename.
    :rtype: str

    :raises TypeError: If ut is not a datetime object.

    :example:

    >>> from datetime import datetime
    >>> make_ccmc_name('model', datetime(2022, 1, 1), '3DALL')
    'MODEL_3DALL_2022-01-01T00-00-00.nc'
    """

    if not isinstance(ut, pd.Timestamp):
        raise TypeError('ut must be a datetime object')

    ut_str = ut.strftime('%Y-%m-%dT%H-%M-%S')

    if data_type is not None:
        return '{}_{}_{}.nc'.format(
            modelname.upper(), data_type.upper(), ut_str)
    else:
        return '{}_{}.nc'.format(modelname.upper(), ut_str)


def get_var_names(dir, models):
    """Print out a list of variable names from a NetCDF file.

    :param dir: Directory to look in.
    :type dir: str
    :param models: Model name(s) to look for.
    :type models: str or list-like
    :return: None
    """

    import xarray as xr
    import os
    from glob import glob

    if isinstance(models, str):
        models = [models]

    for model in models:
        print(model)
        files = glob(os.path.join(dir, model + '*.nc'))
        ds = xr.open_dataset(files[0])
        print(ds.data_vars)
        print('\n\n')
        ds.close()


def autoread(file_list,
             columns_to_return=None,
             concat_dim='time',
             ):
    """Automatically read in a list of files and concatenate them.

    :param file_list: List of files to read in.
    :type file_list: list-like
    :param columns_to_return: Columns (data_vars) to return, defaults to None
        (Read in all comunns)
    :type columns_to_return: str or list-like, optional
    :param concat_dim: concat along this dimension.
        No reason to ever change. Defaults to 'time'.
    :type concat_dim: str, optional
    :raises ValueError: If the requested columns are not in the file.
    :return: Concatenated dataset
    :rtype: xarray.Dataset
    """

    import xarray as xr

    if isinstance(file_list, str):
        file_list = [file_list]

    ds0 = xr.open_dataset(file_list[0])
    drops = []
    for v in ds0.data_vars:
        if v not in columns_to_return:
            drops.append(v)

    if len(drops) == len(ds0.data_vars):
        raise ValueError('No variables to return!\n '
                         'You gave: {}\n'
                         'and available variables are: {}'
                         .format(columns_to_return, ds0.data_vars))

    ds = []
    for filename in file_list:
        ds.append(xr.open_dataset(filename,
                                  drop_variables=drops))

    ds = xr.concat(ds, dim=concat_dim)
    return ds


def ut_to_lt(time_array, glon):
    """Compute local time from date and longitude.

    :param time_array: Array-like of datetime objects in universal time
    :type time_array: array-like
    :param glon: Float or array-like of floats containing geographic longitude
        in degrees. If single value or array of a different shape, all
        longitudes are applied to all times. If the shape is the same as
        `time_array`, the values are paired in the SLT calculation.
    :type glon: array-like or float
    :return: List of local times in hours
    :rtype: array of floats
    """

    time_array = np.asarray(time_array)
    glon = np.asarray(glon)

    # Get UT seconds of day
    try:  # if numpy timestamps
        utsec = [(ut.hour * 3600.0 + ut.minute * 60.0 + ut.second
                  + ut.microsecond * 1.0e-6) / 3600.0 for ut in time_array]
    except BaseException:
        utsec = []
        for ut in time_array:
            ut = pd.Timestamp(ut)
            utsec.append((ut.hour * 3600.0 + ut.minute * 60.0 + ut.second
                          + ut.microsecond * 1.0e-6) / 3600.0)
    # Determine if the calculation is paired or broadcasted
    if glon.shape == time_array.shape:
        lt = np.array([utime + glon[i] / 15.0 for i,
                      utime in enumerate(utsec)])
    else:
        lt = np.array([utime + glon / 15.0 for utime in utsec])

    # Adjust to ensure that 0.0 <= lt < 24.0
    while np.any(lt < 0.0):
        lt[lt < 0.0] += 24.0

    while np.any(lt >= 24.0):
        lt[lt >= 24.0] -= 24.0

    return lt


def add_lt_to_dataset(ds,
                      localtimes=None):
    """Add local time to a dataset.

    :param ds: Dataset to add local time to
    :type ds: xarray.Dataset
    :param localtimes: Localtimes, defaults to None
    :type localtimes: list or array, optional
    :return: Dataset with local time added
    :rtype: xarray.Dataset
    """
    if localtimes is None:
        localtimes = len(ds.lon)

    if isinstance(localtimes, int):
        # Make linspace of localtimes, 0-24, then chop at the ends.
        localtimes = np.linspace(0, 24, localtimes + 1,
                                 endpoint=False)[:-1]
    else:
        localtimes = np.asarray(localtimes)
    if 'lt' in ds.coords:
        ds = ds.rename('lt', 'localtime')

    elif 'localtime' not in ds.coords:
        ds['localtime'] = (('time', 'lon'),
                           ut_to_lt([pd.Timestamp(i) for i in ds.time.values],
                                    ds.lon))

    lt_dss = []

    for t in localtimes:
        ltds = []

        lons_to_iter = np.argmin(np.abs(ds.localtime.values - t), axis=1)

        for do_lon in lons_to_iter:
            b = ds.isel(lon=do_lon).swap_dims(time='localtime')
            ltds.append(b.interp(localtime=t, method='cubic'))

        lt_dss.append(xr.concat(ltds, dim=ds.time))
    return xr.concat(lt_dss, dim='localtime')


def hours_from_storm_onset_into_ds(ds, onset_ut):
    """Add hours from storm onset to a dataset.
    Column is added as "HoursFromStormOnset".

    :param ds: Dataset to add hours from storm onset to
    :type ds: xarray.Dataset
    :param onset_ut: Datetime of the storm onset
    :type onset_ut: datetime
    :raises ValueError: If the dataset spans multiple days
    :return: Dataset with hours from storm onset added
    :rtype: xarray.Dataset
    """

    if ds.day[0] != ds.day[-1]:
        raise ValueError(' Does not yet support multiple days!')
    ds['HoursFromStormOnset'] = ((ds.time.dt.hour - (onset_ut.hour)) +
                                 (ds.time.dt.minute - (onset_ut.minute)) / 60)

    return ds
