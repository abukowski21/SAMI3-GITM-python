# Misc Utility Functions
import numpy as np
from datetime import datetime
import pandas as pd
import xarray as xr


def str_to_ut(in_str):
    """Convert a string to a datetime object.

    Parameters
    ----------
    in_str : str
        String to convert to datetime object.

    Returns
    -------
    datetime
        Datetime object.

    """

    return datetime.strptime(
        in_str.ljust(14, '0'), '%Y%m%d%H%M%S')


def make_ccmc_name(
        modelname,
        ut,
        data_type=None):
    """Make a CCMC-formatted filename.
    Returns a string of the form:
        modelname_datatype_YYYY-MM-DDThh-mm-ss.nc

    Parameters
    ----------
    modelname : str
        Name of model.
    ut : datetime
        Datetime object of the time of the file.
    data_type : str
        Type of data in the file (3DALL/2DANC/ALL/etc.).

    Returns
    -------
    str
        CCMC-formatted filename.

    """

    # format ut str as YYYY-MM-DDThh-mm-ss
    if type(ut) is np.datetime64:
        ut = pd.Timestamp(ut)
    ut_str = ut.strftime('%Y-%m-%dT%H-%M-%S')
    # Make sure modelname & filt_type is all caps

    if data_type is not None:
        return '{}_{}_{}.nc'.format(modelname, data_type, ut_str)

    else:
        return '{}_{}.nc'.format(modelname, ut_str)


def get_var_names(dir, models):
    """Print out a list of variable names.

    Args:
        dir (str: path-like): directory of outputs
        models (str/list): name of model.
    """

    import xarray as xr
    import os
    from glob import glob

    if type(models) is str:
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

    Args:
        file_list (str or list of paths): List of files to read in.
        columns_to_return (str or list, optional):
            Columns (data_vars) to return. Defaults to None.
        concat_dim (str, optional): concat along this dimension.
                No reason to ever change. Defaults to 'time'.

    Raises:
        ValueError: Column not found in files.

    Returns:
        xarray.Dataset: Dataset holding requested variable(s).
    """

    import xarray as xr

    if type(file_list) is str:
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

    Parameters
    ----------
    time_array : array-like
        Array-like of datetime objects in universal time
    glon : array-like or float
        Float or array-like of floats containing geographic longitude in
        degrees. If single value or array of a different shape, all longitudes
        are applied to all times. If the shape is the same as `time_array`,
        the values are paired in the SLT calculation.

    Returns
    -------
    lt : array of floats
        List of local times in hours

    Raises
    ------
    TypeError
        For badly formatted input

    """

    time_array = np.asarray(time_array)
    glon = np.asarray(glon)

    # Get UT seconds of day
    utsec = [(ut.hour * 3600.0 + ut.minute * 60.0 + ut.second
              + ut.microsecond * 1.0e-6) / 3600.0 for ut in time_array]

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


def add_lt_to_dataset(ds,  # xarray.Dataset or xarray.Dataarray
                      localtimes):  # int (for number of localtimes)
    # or list-like for converting those localtimes

    if type(localtimes) == int:
        # Make linspace of localtimes, 0-24, then chop at the ends.
        localtimes = np.linspace(0, 24, localtimes+1,
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
            ltds.append(b.interp(localtime=t))

        lt_dss.append(xr.concat(ltds, dim=ds.time))
    return xr.concat(lt_dss, dim='localtime')