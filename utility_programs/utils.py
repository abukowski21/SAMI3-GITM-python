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
    modelname = modelname.upper()

    if data_type is not None:
        data_type = data_type.upper()
        return '{}_{}_{}.nc'.format(modelname, data_type, ut_str)

    else:
        return '{}_{}.nc'.format(modelname, ut_str)
