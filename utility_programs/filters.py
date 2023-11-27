
import numpy as np
from scipy import signal


def make_fits(da,
              freq=5,
              lims=[40, 85],
              order=1,
              percent=True):
    """Calculate bandpass filter for all data previously read in,
    this is a wrapper for the scipy bandpass filter.
    
    Parameters
    ----------
    da : xarray DataArray or numpy array
        DataArray to be filtered
    freq : int, optional
        Time between outputs, units must be same as lims, by default 5; for 
        data every 5 minutes. 
    lims : list, optional
        Period limits of bandpass filter, by default [40, 85]. Units must be 
        same as `freq`.
    order : int, optional
        Order of the filter, by default 1
    percent : bool, optional
        Return Array as percent over background?
        Set to False to return the absolute perturbation over background.
        Defaults to True.
        
    Returns
    -------
    xarray DataArray
        DataArray with filtered data.
        
    Notes
    -----
    The bandpass filter is applied to the data using the scipy.signal.filtfilt
    function. This function is a forward-backward filter, meaning that the
    data is filtered in both the forward and reverse directions. This results
    in zero phase shift in the filtered data.

    The bandpass filter is designed using the scipy.signal.butter function.
    This function designs a Butterworth filter, which is a type of infinite
    impulse response filter. The filter is designed using the specified
    order, and the cutoff frequencies are calculated from the specified
    limits and the sampling frequency.
    """

    # Define sampling frequency and limits in minutes
    lower_limit = min(lims)
    upper_limit = max(lims)

    # Convert limits to corresponding indices
    # lower_index = int(lower_limit / freq)
    # upper_index = int(upper_limit / freq)

    # Design the bandpass filter
    nyquist_freq = 0.5 / freq
    lower_cutoff = (1/lower_limit) / nyquist_freq
    upper_cutoff = (1/upper_limit) / nyquist_freq
    b, a = signal.butter(order, [upper_cutoff, lower_cutoff],
                         btype='band', analog=False)

    # Apply the filter to the data
    filtd = signal.filtfilt(b, a, da, axis=0)
    # filtd = xr.apply_ufunc(filtfilt, b, a, da, dask='allowed')

    if percent:
        return (100 * (filtd) / da)

    else:
        if isinstance(da, np.ndarray):
            return filtd
        else:
            da.values = filtd
            return da

def remove_outliers(array):
    """ Remove outliers from an array by replacing them with the median value.
    
    Parameters
    ----------
    array : numpy array or xarray DataArray
        Data to be filtered
        
    Returns
    -------
    numpy array or xarray DataArray
        Filtered data
        
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
    
    """ Calculate the difference of a DataArray along a given dimension.
    
    Parameters
    ----------
    da : xarray DataArray
        DataArray to be filtered
    dim : str, optional
        Dimension over which to calculate the diff, by default 'time'
    order : int, optional
        Order of the diff, by default 2
    percent : bool, optional
        Return the percent over diff (False), or just return the fit (True),
        by default False
    label : str, optional
        Label values to the 'lower' or 'upper' bound, by default 'lower'
        
    Returns
    -------
    xarray DataArray
        Filtered data
        
    """

    if percent:
        return 100 * (da.diff(dim, order, label) / da)
    else:
        return da.diff(dim, order, label)
