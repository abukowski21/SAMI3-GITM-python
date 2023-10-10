import numpy as np
import pandas as pd


def UT_from_Storm_onset(itime, dtime_storm_start):
    """Calculate the UT from storm onset, returning a string.
    Before storm onset, the UT is negative and must be treated
        differently.

    :param itime: Datetime of the time step
    :type itime: Datetime
    :param dtime_storm_start: Datetime of the start of the storm
    :type dtime_storm_start: Datetime
    :return: UT from storm onset (as HH:MM)
    :rtype: str
    """

    minute_diff = ((pd.Timestamp(itime) - dtime_storm_start) /
                   pd.Timedelta('1 minute'))
    if minute_diff > 0:
        hrs = np.floor(minute_diff / 60)
        hrs = str(int(hrs)).rjust(2, '0')
    else:
        hrs = np.ceil(minute_diff / 60)
        hrs = '-' + str(np.abs(int(hrs))).rjust(2, '0')
    mins = str(int(np.abs(minute_diff) % 60)).rjust(2, '0')
    ut = hrs + ':' + mins
    return ut
