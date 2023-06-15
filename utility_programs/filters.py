
from scipy import signal
import numpy as np
import xarray as xr


def make_filter(lowcut=80, highcut=40, order=3):
    """_summary_

    Args:
        lowcut (int, optional): lowcut of the filter. Defaults to 100.
        highcut (int, optional): highcut of the filter. Defaults to 30.
        order (int, optional): Order of the filter (how steep the filter
            cuts at the cutoff frequencies). Defaults to 2.

    Returns:
        scipy butterworth filter: the filter with settings defined by the user.
    """
    # Define the cutoff frequencies
    lowcut_f = 1 / (lowcut / 60)  # 100 minutes in units of sample^-1
    highcut_f = 1 / (highcut / 60)  # 30 minutes in units of sample^-1

    # Define the Butterworth filter
    nyquist = 0.5 * 5  # 5 minutes is the sampling frequency
    low = lowcut_f / nyquist
    high = highcut_f / nyquist
    sos = signal.butter(order, [low, high], btype="bandstop", output="sos")
    return sos


def make_fits(gitm_bins, lowcut=80, highcut=40, order=3):
    """
    calculate bandpass filter for all data previously read in.

    inputs: nparray (or xarray.DataArray) of gitmdata

    returns:
    fits: np array indexed at fits[time][col][ilon][ilat][ialt]


    todo: you can thread this by splitting the alts into different threads.
    then just append the fits_full later.

    """
    sos = make_filter(lowcut, highcut, order)
    

    
    if type(gitm_bins) == xr.Dataset or type(gitm_bins) == xr.DataArray:
        filtered_arr = xr.apply_ufunc(signal.sosfiltfilt, 
                                      sos, gitm_bins,
                                     dask='allowed')
    else:
        filtered_arr = signal.sosfiltfilt(sos, gitm_bins, axis=0)
    return filtered_arr


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
        return 100*(da.diff(dim, order, label)/da)
    else:
        return da.diff(dim, order, label)

