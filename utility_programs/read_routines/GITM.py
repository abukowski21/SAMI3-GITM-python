import datetime
import glob
import os

import numpy as np
from aetherpy.io import read_routines
from tqdm.auto import tqdm
from struct import unpack


def read_gitm_into_nparrays(gitm_dir, dtime_storm_start,
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
    
    

def read_gitm_bin_xarray(filename, 
                        add_time=True,
                        drop_ghost_cells=True,
                        cols=None):

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
        version = unpack(end_char + 'd', fin.read(rec_len))[0]

        _, rec_len = unpack(end_char + '2l', fin.read(8))

        # Read grid size information.
        nlons, nlats, nalts = unpack(
            end_char + 'lll', fin.read(rec_len))
        _, rec_len = unpack(end_char + '2l', fin.read(8))

        # Read number of variables.
        num_vars = unpack(end_char + 'l', fin.read(rec_len))[0]
        _, rec_len = unpack(end_char + '2l', fin.read(8))

        file_vars = np.arange(0, num_vars, 1)

        varnames = []

        # Collect variable names in a list
        for ivar in range(num_vars):
            vcode = unpack(end_char + '%is' % (rec_len),
                           fin.read(rec_len))[0]
            var = vcode.decode('utf-8').replace(" ", "")
            var=var.replace('!N','').replace('!U','').replace('!D','')\
                .replace('[','').replace('[','').replace(']','')
            varnames.append(var)
            dummy, rec_lec = unpack(end_char + '2l', fin.read(8))

        # Extract time
        rec_time = np.array(unpack(end_char + 'lllllll', fin.read(28)))
        rec_time[-1] *= 1000  # convert from millisec to microsec
        time_here = datetime(*rec_time)

        # Header is this length:
        # Version + start/stop byte
        # nlons, nlats, nalts + start/stop byte
        # num_vars + start/stop byte
        # variable names + start/stop byte
        # time + start/stop byte

        iheader_length = 84 + num_vars * 48

        ntotal = nlons * nlats * nalts
        idata_length = ntotal * 8 + 8

        data_vars={}

        # Save the data for the desired variables
        dimnames=['lon','lat','alt']
        
        for ivar in file_vars:
            fin.seek(iheader_length + ivar * idata_length)
            sdata = unpack(end_char + 'l', fin.read(4))[0]

            if ivar == 0:
                lons =  np.rad2deg(np.unique(np.array(
                        unpack(end_char + '%id' % (ntotal), fin.read(sdata))).reshape(
                        (nlons, nlats, nalts), order="F")))

            elif ivar == 1:
                lats =  np.rad2deg(np.unique(np.array(
                        unpack(end_char + '%id' % (ntotal), fin.read(sdata))).reshape(
                        (nlons, nlats, nalts), order="F")))

            elif ivar == 2:
                alts =  np.unique(np.array(
                        unpack(end_char + '%id' % (ntotal), fin.read(sdata))).reshape(
                        (nlons, nlats, nalts), order="F"))/1000


            else:
                data_vars[varnames[ivar]] = dimnames,np.array(
                    unpack(end_char + '%id' % (ntotal), fin.read(sdata))).reshape(
                        (nlons, nlats, nalts), order="F")
                # break
    ds = xr.Dataset(coords={'time':[time_here],'lon':lons,'lat':lats,'alt':alts},
                    data_vars=data_vars,
                    attrs={'version':version,
                          'dropped-ghost-cells':str(drop_ghost_cells),
                          'with_time':str(add_time)})

    if drop_ghost_cells:
        if nalts > 1:
            ds = ds.drop_isel(lat=[0,1,-2,-1],lon=[0,1,-1,-2],alt=[0,1,-1,-2])
        else:
            ds = ds.drop_isel(lat=[0,1,-2,-1],lon=[0,1,-1,-2])
    if not add_time:
        ds = ds.drop_vars('time')

    if cols is not None:
        ds = ds.get(cols)
                
    return ds


def gitm_times_from_filelist(file_list, century_prefix='20'):
    gitm_dtimes = []
    for i in file_list:
        yy, MM, dd = i[-17:-15], i[-15:-13], i[-13:-11]
        hr, mm, sec = i[-10:-8], i[-8:-6], i[-6:-4]
        try:
            gitm_dtimes.append(
                datetime(
                    int(century_prefix + yy), int(MM), int(dd),
                    int(hr), int(mm), int(sec)))
        except ValueError:
            raise ValueError(
                "GITM file name does not match expected format,",
                "filename %s cannot be parsed" % i)
    return gitm_dtimes



def read_gitm_multiple_bins(file_list,
                            start_dtime=None,
                            end_dtime=None,
                            start_idx=0,
                            end_idx=-1,
                            drop_ghost_cells=False,
                            cols=None):
    
    # Check inputs! Cannot specify start time & idx:
    if start_dtime is not None and start_idx is not None:
        raise ValueError("Cannot specify both Start idx & dtime")
    if end_dtime is not None and end_idx is not None:
        raise ValueError("Cannot specify both End idx & dtime")

    file_list=file_list[start_idx:end_idx]
    
    if start_dtime is not None or end_dtime is not None:
        times = gitm_times_from_filelist(file_list)
    if start_dtime is not None:
        time_mask = np.where(times>=start_dtime)
        file_list = file_list[time_mask]
        times = times[time_mask]
    if end_dtime is not None:
        time_mask = np.where(times<=end_dtime)
        file_list = file_list[time_mask]
    
    ds=[]
    for file in file_list:
        ds.append(read_gitm_bin_xarray(file,
                                       drop_ghost_cells=drop_ghost_cells,
                                       cols=cols))
        
    ds = xr.concat(ds,'time')
    
    return ds