# Misc Utility Functions
import numpy as np
from datetime import datetime
import pandas as pd


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
