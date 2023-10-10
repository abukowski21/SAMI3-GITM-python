
import numpy as np
from scipy import signal


def make_fits(da,
              freq=5,
              lims=[40, 85],
              order=1,
              percent=True):
    """Calculate bandpass filter for all data previously read in.
    this is a wrapper for the scipy bandpass filter.

    :param da: DataArray to be filtered
    :type da: xarray DataArray
    :param freq: Output frequency, units must be same as lims, defaults to 5
    :type freq: int, optional
    :param lims: Limits of bandpass filter, defaults to [40, 85]
    :type lims: list, optional
    :param order: Order of the filter, defaults to 1
    :type order: int, optional
    :param percent: Return DataArray as percent over background?
        Set to False to return the fit. Defaults to True
    :type percent: bool, optional
    :return: Filtered data
    :rtype: _xarray.DataArray
    """

    # Define sampling frequency and limits in minutes
    lower_limit = min(lims)
    upper_limit = max(lims)

    # Convert limits to corresponding indices
    lower_index = int(lower_limit / freq)
    upper_index = int(upper_limit / freq)

    # Design the bandpass filter
    nyquist_freq = 0.5 * freq
    lower_cutoff = lower_index / nyquist_freq
    upper_cutoff = upper_index / nyquist_freq
    b, a = signal.butter(order, [1 / upper_cutoff, 1 / lower_cutoff],
                         btype='band', analog=False)

    # Apply the filter to the data
    filtd = signal.filtfilt(b, a, da, axis=0)
    # filtd = xr.apply_ufunc(filtfilt, b, a, da, dask='allowed')

    if percent:
        return (100 * (filtd) / da)

    else:
        da.values = filtd
        return da


def remove_outliers(array):
    """Remove outliers from an array by replacing them with the median value.

    :param array: Data to be filtered
    :type array: numpy array or xarray DataArray
    :return: Filtered data
    :rtype: numpy array or xarray DataArray
    """

    arr2 = array.copy()
    # calculate mean, standard deviation, and median over all elements
    mean, std, median = np.mean(arr2), np.std(arr2), np.median(arr2)
    # set outlier threshold (in terms of number of standard deviations)
    outlier_threshold = 5
    outliers = np.logical_or(
        arr2 < mean - outlier_threshold * std, arr2 > mean +
        outlier_threshold * std)  # find outliers
    arr2[outliers] = median  # set outliers to median
    return arr2


def filter_xarray_DA_diff(da, dim='time', order=2, percent=False,
                          label='lower'):
    """Calculate the difference of a DataArray along a given dimension.

    :param da: DataArray to be filtered
    :type da: xarray DataArray
    :param dim: Dimension over which to calculate the diff, defaults to 'time'
    :type dim: str, optional
    :param order: Order of the diff, defaults to 2
    :type order: int, optional
    :param percent: Return the percent over diff (False),
        or just return the fit (True), defaults to False
    :type percent: bool, optional
    :param label: Label values to the 'lower' or 'upper' bound,
        defaults to 'lower'
    :type label: str, optional
    :return: Filtered data
    :rtype: _xarray.DataArray
    """

    if percent:
        return 100 * (da.diff(dim, order, label) / da)
    else:
        return da.diff(dim, order, label)
