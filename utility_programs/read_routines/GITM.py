import datetime
import glob
import os

import numpy as np
from aetherpy.io import read_routines
from tqdm.auto import tqdm


def read_gitm_into_nparrays(gitm_dir, dtime_storm_start=None,
                            gitm_file_pattern='3DALL*.bin',
                            cols=['all'],
                            t_start_idx=0, t_end_idx=-1,
                            century_prefix='20',
                            return_vars=False):
    """reads in gitm data into a dictionary of numpy arrays

    Args:
        gitm_dir (str): path to gitm files
        dtime_storm_start (datetime.datetime): time of storm onset
        gitm_file_pattern (str, optional): file pattern to match.
            Defaults to '3DALL*.bin'.
        cols (list, optional): which columns to read (strs).
            Defaults to ['all'].
        t_start_idx (int, optional): hrs before storm onset to start.
            Defaults to 0 (all data).
        t_end_idx (int, optional): hrs after onset to end.
            Defaults to -1 (all data).
        century_prefix (str, optional): century. Defaults to '20'.

    Raises:
        ValueError: _description_
        ValueError: _description_
        ValueError: _description_

    Returns:
        gitm_dtimes (list): list of datetimes
        gitmgrid (dict): grid. keys of ['latitude', 'longitude', 'altitude']
        gitmbins (dict): data. keys are variables we opted to get.
            shape is [ntimes, nvars, nlons, nlats, nalts]
    """""""""

    flist = np.sort(glob.glob(os.path.join(gitm_dir, gitm_file_pattern)))
    if len(flist) == 0:
        raise ValueError("No %s files found in %s" %
                         (gitm_file_pattern, gitm_dir),
                         "\n \n instead there is: ",
                         glob.glob(os.path.join(gitm_dir, gitm_file_pattern)))

    gitm_dtimes = []
    for i in flist:
        yy, MM, dd = i[-17:-15], i[-15:-13], i[-13:-11]
        hr, mm, sec = i[-10:-8], i[-8:-6], i[-6:-4]
        try:
            gitm_dtimes.append(
                datetime.datetime(
                    int(century_prefix + yy), int(MM), int(dd),
                    int(hr), int(mm), int(sec)))
        except ValueError:
            raise ValueError(
                "GITM file name does not match expected format,",
                "filename %s cannot be parsed" % i)

    if dtime_storm_start is not None:
        if t_start_idx != 0:
            start_idx = gitm_dtimes.index(
                dtime_storm_start - datetime.timedelta(hours=t_start_idx))
            gitm_dtimes = gitm_dtimes[start_idx:]
            flist = flist[start_idx:]

        if t_end_idx != -1:
            end_idx = gitm_dtimes.index(
                dtime_storm_start + datetime.timedelta(hours=t_end_idx))
            gitm_dtimes = gitm_dtimes[:end_idx]
            flist = flist[:end_idx]

    f = read_routines.read_gitm_file(flist[0])
    if '3DALL' in gitm_file_pattern:
        gitmgrid = {f["vars"][k].lower(): f[k][2:-2, 2:-2, 2:-2]
                    for k in [0, 1, 2]}
        nlons, nlats, nalts = np.array(f[0].shape) - 4  # ghost cells
    elif '2DANC' in gitm_file_pattern:
        gitmgrid = {f["vars"][k].lower(): f[k]
                    for k in [0, 1, 2]}
        nlons, nlats, nalts = np.array(f[0].shape)  # NO ghost cells

    # Don't get these variables
    ignore_cols = ['Longitude', 'Latitude', 'Altitude']

    if 'all' in cols:  # get all variables
        gitmvars = [i for i in f["vars"] if i not in ignore_cols]
    else:  # get only the variables in cols
        gitmvars = [i for i in f["vars"] if i in cols and i not in ignore_cols]
    if len(gitmvars) == 0:  # if no variables found
        raise ValueError("No GITM variables found in file!",
                         "found these variables: ",
                         f["vars"], ' cols we want are: ', cols,
                         ' ignore cols: ', ignore_cols)
    try:
        gitmbins = np.zeros([len(flist), len(gitmvars), nlons, nlats, nalts])
    except ValueError:
        raise ValueError("GITM file has different dimensions than expected!",
                         len(flist), len(gitmvars), nlons, nlats, nalts,
                         '\n flist, gitmvars, nlons, nlats, nalts')

    for ifile, file_name in enumerate(tqdm(flist)):
        f = read_routines.read_gitm_file(file_name)

        for num_var, real_var in enumerate(gitmvars):
            num_v_src = f["vars"].index(real_var)
            if '3DALL' in gitm_file_pattern:
                gitmbins[ifile, num_var] = f[num_v_src][2:-2, 2:-2, 2:-2]
            elif '2DANC' in gitm_file_pattern:
                gitmbins[ifile, num_var] = f[num_v_src]

    gitmgrid["latitude"] = np.rad2deg(gitmgrid["latitude"])
    gitmgrid["longitude"] = np.rad2deg(gitmgrid["longitude"])

    # Fix the ordering of the longitudes and go from -180-180 not 0->360
    newlons_for_order = []
    for ilon in range(len(gitmgrid["longitude"])):
        oldlon = gitmgrid["longitude"][ilon, 0, 0]
        if oldlon <= 180:
            newlons_for_order.append(int(oldlon))

        else:
            newlons_for_order.append(int(oldlon) - 360)
            gitmgrid["longitude"][ilon] = gitmgrid["longitude"][ilon] - 360

    new_lons_sorted = np.sort(newlons_for_order)
    new_order = np.array(
        [newlons_for_order.index(new_lons_sorted[i])
         for i in range(len(new_lons_sorted))])

    gitmbins = gitmbins[:, :, new_order, :, :]
    gitmgrid["longitude"] = np.sort(gitmgrid["longitude"], axis=0)

    if return_vars:
        return gitm_dtimes, gitmgrid, gitmbins, gitmvars
    else:
        return gitm_dtimes, gitmgrid, gitmbins
