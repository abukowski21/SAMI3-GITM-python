#!/usr/bin/env python
# Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md
"""Routines to read Aether files.

Taken from github.com/AetherModel/Aetherpy


"""

import datetime as dt
import numpy as np
import os
from struct import unpack


def parse_line_into_int_and_string(line, parse_string=True):
    """Parse a data string into integer and string components.

    Parameters
    ----------
    line : str
        Data line with a known format, described in Notes section
    parse_string : bool
        Construct a single white space separated string consisting of the
        remaining portion of `line`

    Returns
    -------
    line_num : int
        Integer corresponding to the first number in `line`
    line_str : str
        Whitespace separated string containing the rest of the `line` data if
        `parse_string` is True, otherwise it returns the original line

    Notes
    -----
    Splits data line using a single space.  The first element is
    expected to be an integer. If desired, the remaining data is returned as a
    single space-separated string.

    """

    split_line = line.split(" ")
    line_num = int(split_line[0])

    if parse_string:
        line_str = " ".join(split_line[1:])
    else:
        line_str = str(line)

    return line_num, line_str


def read_gitm_headers(filelist, finds=-1):
    """Read ancillary information from a GITM file.

    Parameters
    ----------
    filelist : array-like
        Array-like object of names for netCDF files
    finds : int, list, or slice
        Index(es) for file(s) from which the header information will be read.
        (default=-1)

    Returns
    -------
    header : dict
        A dictionary containing header information from the netCDF files,
        including:
        nlons - number of longitude grids
        nlats - number of latitude grids
        nalts - number of altitude grids
        vars - list of data variable names
        time - list of datetimes for the processed file start times
        filename - list of the processed filenames
        version - file version number

    Raises
    ------
    IOError
        If any one of the input files encounters an unexpected value or
        dimension.

    Notes
    -----
    This routine obtains the same info as `read_aether_headers`

    """

    # Initialize the output
    header = {"vars": [], "time": [], "filename": np.asarray(filelist)}

    # Ensure the filelist is array-like, allowing slicing of input
    if header['filename'].shape == ():
        header['filename'] = np.asarray([filelist])
    else:
        header['filename'] = list(header['filename'])

    # Ensure selected files are list-like
    hfiles = np.asarray(header['filename'][finds])
    if hfiles.shape == ():
        hfiles = np.asarray([hfiles])

    for filename in hfiles:
        # Read in the header from the binary file
        file_vars = list()

        with open(filename, 'rb') as fin:
            # Test to see if the correct endian is being used
            end_char = '>'
            raw_rec_len = fin.read(4)
            rec_len = (unpack(end_char + 'l', raw_rec_len))[0]
            if rec_len > 10000 or rec_len < 0:
                # Ridiculous record length implies wrong endian, fix here
                end_char = '<'
                rec_len = (unpack(end_char + 'l', raw_rec_len))[0]

            # Read version; read fortran footer+header.
            file_version = unpack(end_char + 'd', fin.read(rec_len))[0]
            _, rec_len = unpack(end_char + '2l', fin.read(8))

            # Test the version number
            if 'version' not in header.keys() == 0:
                header["version"] = file_version
            elif header['version'] != file_version:
                raise IOError(''.join(['unexpected version number in file ',
                                       filename]))

            # Read grid size information.
            nlons, nlats, nalts = unpack(end_char + 'lll', fin.read(rec_len))
            _, rec_len = unpack(end_char + '2l', fin.read(8))

            # Test the dimensions
            if np.any([dim_var not in header.keys()
                       for dim_var in ['nlons', 'nlats', 'nalts']]):
                header["nlons"] = nlons
                header["nlats"] = nlats
                header["nalts"] = nalts
            elif (header['nlons'] != nlons or header['nlats'] != nlats
                  or header['nalts'] != nalts):
                raise IOError(''.join(['unexpected dimensions in file ',
                                       filename]))

            # Read number of variables.
            num_vars = unpack(end_char + 'l', fin.read(rec_len))[0]
            _, rec_len = unpack(end_char + '2l', fin.read(8))

            # Collect variable names.
            for ivar in range(num_vars):
                vcode = unpack(end_char + '%is' % (rec_len),
                               fin.read(rec_len))[0]
                var = vcode.decode('utf-8').replace(" ", "")
                file_vars.append(var)
                _, rec_len = unpack(end_char + '2l', fin.read(8))

            # Test the variable names
            if len(header["vars"]) == 0:
                header["vars"] = list(file_vars)
            elif header["vars"] != file_vars:
                raise IOError(''.join(['unexpected number or name of ',
                                       'variables in file ', filename]))

            # Extract time
            out_time = np.array(
                unpack(
                    end_char +
                    'lllllll',
                    fin.read(rec_len)))
            out_time[-1] *= 1000  # Convert from millisec to microsec
            header["time"].append(dt.datetime(*out_time))

    return header


def read_gitm_file(filename, file_vars=None):
    """Read list of variables from one GITM file.

    Parameters
    ----------
    filename : str
        GITM file to read
    file_vars : list or NoneType
        List of desired variable names to read or None to read all
        (default=None)

    Returns
    -------
    data : dict
        Dict with keys 'time', which contains a datetime object specifying the
        time of the file and zero-offset indices, corresponding to the
        variable names in `file_vars` that holds arrays of the specified data.
        Also contains version number, dimensions, and a list of the variable
        names.

    """

    data = {"vars": []}

    if not os.path.isfile(filename):
        raise IOError('input file does not exist')

    with open(filename, 'rb') as fin:
        # Determine the correct endian
        end_char = '>'
        raw_rec_len = fin.read(4)
        rec_len = (unpack(end_char + 'l', raw_rec_len))[0]
        if rec_len > 10000 or rec_len < 0:
            # Ridiculous record length implies wrong endian.
            end_char = '<'
            rec_len = (unpack(end_char + 'l', raw_rec_len))[0]

        # Read version; read fortran footer+data.
        data["version"] = unpack(end_char + 'd', fin.read(rec_len))[0]

        _, rec_len = unpack(end_char + '2l', fin.read(8))

        # Read grid size information.
        data["nlons"], data["nlats"], data["nalts"] = unpack(
            end_char + 'lll', fin.read(rec_len))
        _, rec_len = unpack(end_char + '2l', fin.read(8))

        # Read number of variables.
        num_vars = unpack(end_char + 'l', fin.read(rec_len))[0]
        _, rec_len = unpack(end_char + '2l', fin.read(8))

        if file_vars is None:
            file_vars = np.arange(0, num_vars, 1)

        # Collect variable names in a list
        for ivar in range(num_vars):
            vcode = unpack(end_char + '%is' % (rec_len),
                           fin.read(rec_len))[0]
            var = vcode.decode('utf-8').replace(" ", "")
            data['vars'].append(var)
            dummy, rec_lec = unpack(end_char + '2l', fin.read(8))

        # Extract time
        rec_time = np.array(unpack(end_char + 'lllllll', fin.read(28)))
        rec_time[-1] *= 1000  # convert from millisec to microsec
        data["time"] = dt.datetime(*rec_time)

        # Header is this length:
        # Version + start/stop byte
        # nlons, nlats, nalts + start/stop byte
        # num_vars + start/stop byte
        # variable names + start/stop byte
        # time + start/stop byte

        iheader_length = 84 + num_vars * 48

        ntotal = data["nlons"] * data["nlats"] * data["nalts"]
        idata_length = ntotal * 8 + 8

        # Save the data for the desired variables
        for ivar in file_vars:
            fin.seek(iheader_length + ivar * idata_length)
            sdata = unpack(end_char + 'l', fin.read(4))[0]
            data[ivar] = np.array(
                unpack(end_char + '%id' % (ntotal), fin.read(sdata))).reshape(
                    (data["nlons"], data["nlats"], data["nalts"]), order="F")

    return data
