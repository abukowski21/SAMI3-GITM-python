# Misc Utility Functions
import numpy as np
from datetime import datetime
import pandas as pd
import xarray as xr
import glob
import os


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
    if isinstance(ut, np.datetime64):
        ut = pd.Timestamp(ut)
    ut_str = ut.strftime('%Y-%m-%dT%H-%M-%S')
    # Make sure modelname & filt_type is all caps

    if data_type is not None:
        return '{}-{}_{}.nc'.format(modelname, data_type, ut_str)

    else:
        return '{}_{}.nc'.format(modelname, ut_str)


def get_var_names(dir, models):
    """
    Print out a list of variable names.

    Parameters
    ----------
    dir : str
        Directory of outputs.
    models : str or list
        Name of model.

    Returns
    -------
    None

    Notes
    -----
    This function prints out a list of variable names for a given model or list of models. It does this by searching for
    netCDF files in the specified directory that match the model name(s), opening the first file found, and printing out
    the names of the data variables in the file.

    Examples
    --------
    >>> get_var_names('/path/to/outputs', 'model1')
    model1
    <xarray.Dataset>
    Dimensions:  (time: 10, x: 100, y: 100)
    Coordinates:
      * time     (time) datetime64[ns] 2000-01-01 2000-02-01 ... 2000-10-01
      * x        (x) float64 0.0 1.0 2.0 3.0 4.0 ... 96.0 97.0 98.0 99.0 100.0
      * y        (y) float64 0.0 1.0 2.0 3.0 4.0 ... 96.0 97.0 98.0 99.0 100.0
    Data variables:
        var1    (time, y, x) float64 ...
        var2    (time, y, x) float64 ...
        var3    (time, y, x) float64 ...
        var4    (time, y, x) float64 ...
        var5    (time, y, x) float64 ...
        var6    (time, y, x) float64 ...
        var7    (time, y, x) float64 ...
        var8    (time, y, x) float64 ...
        var9    (time, y, x) float64 ...
        var10   (time, y, x) float64 ...

    """

    if isinstance(models, str):
        models = [models]

    for model in models:
        print(model)
        files = glob(os.path.join(dir, model + '*.nc'))
        ds = xr.open_dataset(files[0])
        print(ds.data_vars)
        print('\n\n')
        ds.close()


def autoread(file_list, columns_to_return=None, concat_dim='time'):
    """
    Automatically read in a list of files and concatenate them.

    Parameters
    ----------
    file_list : str or list of paths
        List of files to read in.
    columns_to_return : str or list, optional
        Columns (data_vars) to return. Defaults to None.
    concat_dim : str, optional
        Concatenate along this dimension. No reason to ever change. Defaults to 'time'.

    Raises
    ------
    ValueError
        Column not found in files.

    Returns
    -------
    xarray.Dataset
        Dataset holding requested variable(s).


    Notes
    -----
    This is not well supported. Best to use xarray.open_mfdataset() instead.

    """

    if isinstance(file_list, str):
        file_list = [file_list]

    ds0 = xr.open_dataset(file_list[0])
    drops = []
    for v in ds0.data_vars:
        if v not in columns_to_return:
            drops.append(v)

    if len(drops) == len(ds0.data_vars):
        raise ValueError(
            "No variables to return!\nYou gave: %s\nand available"
            "variables are: %s " % (str(columns_to_return, str(ds0.data_vars))))

    ds = []
    for filename in file_list:
        ds.append(xr.open_dataset(filename, drop_variables=drops))

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
    try:
        utsec = [(ut.hour * 3600.0 + ut.minute * 60.0 + ut.second
                  + ut.microsecond * 1.0e-6) / 3600.0 for ut in time_array]
    except AttributeError:  # if datetime objects
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
                      localtimes=[2, 6, 10, 14, 18, 22],
                      pbar=False):
    """Add localtime as a coordinate to an existing dataset/dataarray

    Parameters
    ----------
    ds : xarray.dataset or xarray.dataarray
        xarray object to add the column to. Can be dataset (takes longer) 
        or dataarray. ** MUST HAVE `time` AS A DIMENSION**.
    localtimes : array-like or int
        Localtimes to be returned! If int, this is the number of 
        localtimes. They will be evenly spaced from 0-24.
        If list-like, these are the actual localtimes to return. Be careful 
        adding too many of these. Since we do some interpolating, the more 
        localtimes you want, the longer this will take. Optional to show a 
        progress bar is it's taking a long time.
    pbar : bool
        Set this to true to show a progress bar. 

    Returns
    -------
    ds : xarray.(dataset/dataarray)
        Same exact format as the input dataset, but with a new coordinate of 
        'localtime''

    """

    if isinstance(localtimes, int):
        # Make linspace of localtimes, 0-24, then chop at the ends.
        localtimes = np.linspace(0, 24 + (24/localtimes), localtimes + 1,
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

    # debug
    if pbar:
        from tqdm import tqdm
        pbar = tqdm(total=ds.time.size)

    ltdss = []

    for itime in range(ds.time.size):
        timeds = ds.isel(time=itime).swap_dims(lon='localtime')
        ltdss.append(timeds.interp(localtime=localtimes))
        if pbar:
            pbar.update()

    return xr.concat(ltdss, dim='time')


def hours_from_storm_onset_into_ds(ds, onset_ut):
    """
    Calculate hours from an event and add to dataset.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset to add hours from storm onset to
    onset_ut : datetime.datetime
        Datetime object of storm/event

    Raises
    ------
    ValueError
        If multiple days are in the dataset

    Returns
    -------
    xarray.Dataset
        Dataset with hours from storm onset added
    """

    if ds.day[0] != ds.day[-1]:
        raise ValueError(' Does not yet support multiple days!')
    ds['HoursFromStormOnset'] = ((ds.time.dt.hour - (onset_ut.hour)) +
                                 (ds.time.dt.minute - (onset_ut.minute)) / 60)

    return ds
