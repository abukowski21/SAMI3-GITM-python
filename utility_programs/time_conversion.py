#!/usr/bin/env python
# Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md

# TAKEN FROM https://github.com/aethermodel/aetherpy on April 10 2023

"""Routines to perform temporal calculations."""

import datetime as dt
import numpy as np
import pandas as pd


def epoch_to_datetime(epoch_time):
    """Convert from epoch seconds to datetime.

    :param epoch_time: Seconds since 1 Jan 1965
    :type epoch_time: float
    :return: Datetime object corresponding to `epoch_time`
    :rtype: datetime
    """

    dtime = dt.datetime(1965, 1, 1) + dt.timedelta(seconds=epoch_time)

    return dtime


def datetime_to_epoch(dtime):
    """Convert datetime to epoch seconds.

    :param dtime: Datetime object
    :type dtime: datetime
    :return: Seconds since 1 Jan 1965
    :rtype: float
    """

    epoch_time = (dtime - dt.datetime(1965, 1, 1)).total_seconds()

    return epoch_time


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


def lt_to_ut(lt, glon):
    """Compute universal time in hours from local time and longitude.

    :param lt: Local time(s) in hours
    :type lt: float or array-like
    :param glon: Geographic longitude(s) in degrees.
    :type glon: float or array-like
    :return: Universal time in hours
    :rtype: float
    """
    uth = np.asarray(lt) - np.asarray(glon) / 15.0

    # Ensure the values range from 0-24 h
    uth[uth <= 0] += 24.0
    uth[uth >= 24.0] -= 24.0

    return uth


def calc_time_shift(utime):
    """Calculate the time shift needed to orient a polar dial.

    :param utime: Datetime object of the time we're plotting
    :type utime: datetime
    :return: Time shift in degrees
    :rtype: float
    """

    uth = utime.hour + (utime.minute + (
        utime.second + utime.microsecond * 1.0e-6) / 60.0) / 60.0
    shift = uth * 15.0

    return shift
