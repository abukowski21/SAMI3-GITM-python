
import numpy as np
from scipy import signal


def make_fits(da,
              freq=5,
              lims=[40, 85],
              order=1,
              percent=True):
    """
    calculate bandpass filter for all data previously read in.
    this is a wrapper for the scipy bandpass filter.

    Args:
        da (xarray DataArray): data to be filtered
        freq (int, optional): sampling frequency . Defaults to 5.
        lims (list, optional): limits of the bandpass filter.
            Defaults to [40, 85].
        order (int, optional): order of the filter. Defaults to 1.
        percent (bool, optional): whether to return the filtered data as
            a percentage of the original data. Defaults to True.

    Returns:
        xarray DataArray: filtered data


    """

    # Define sampling frequency and limits in minutes
    sampling_freq = freq
    lower_limit = min(lims)
    upper_limit = max(lims)

    # Convert limits to corresponding indices
    lower_index = int(lower_limit / sampling_freq)
    upper_index = int(upper_limit / sampling_freq)

    # Design the bandpass filter
    nyquist_freq = 0.5 * sampling_freq
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

    if percent:
        return 100 * (da.diff(dim, order, label) / da)
    else:
        return da.diff(dim, order, label)
