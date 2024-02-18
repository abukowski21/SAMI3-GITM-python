import datetime
import glob
import os

import numpy as np
from tqdm import tqdm
from struct import unpack
import xarray as xr
from utility_programs.utils import make_ccmc_name
from dask.diagnostics import ProgressBar


def read_bin_to_nparrays(gitm_dir,
                         gitm_file_pattern='3DALL*.bin',
                         cols=['all'],
                         dtime_start=None,
                         dtime_end=None,
                         start_idx=0, end_idx=-1,
                         century_prefix='20',
                         return_vars=False,
                         progress_bar=False):
    """reads in gitm data into a dictionary of numpy arrays
    (deprecated, use read_to_xarray instead)

    Args:
        gitm_dir (str): path to gitm files
        gitm_file_pattern (str, optional): file pattern to match.
            Defaults to '3DALL*.bin'.
        cols (list, optional): which columns to read (strs).
            Defaults to ['all'].
        dtime_start (datetime, optional): datetime to start read.
            Defaults to 0 (all data).
        dtime_end (datetime, optional): datetime to end read.
            Defaults to -1 (all data).
        start_idx (int, optional): index to start reading at.
        end_idx (int, optional): index to end reading at.
        century_prefix (str, optional): century. Defaults to '20'.
        progress_bar (bool, optional):
            show progress bar. Defaults to False. (Requires tqdm)

    Raises:
        ValueError: If GITM files don't exist or are in a weird format.

    Returns:
        dict: dictionary of numpy arrays with keys:

            gitmdtimes:
                times of gitm outputs
            gitmbins:
                gitm data
            gitmgrid:
                dictionary of grid variables
            gitmvars (optionaln only returned if return_vars):
                list of variables

    """

    try:
        import read_from_aether as ather_read

    except ModuleNotFoundError:
        import sys
        sys.path.insert(0, 'utility_programs/read_routines/')
        # This does not work when rendering 'docs/source/plotting.ipynb',
        # or any of the notebooks in refs folder,so we need to
        # take those into account with:
        try:
            import read_from_aether as ather_read
        except ModuleNotFoundError:
            sys.path.insert(0, '../utility_programs/read_routines')
            sys.path.insert(0, '../../utility_programs/read_routines')
            import read_from_aether as ather_read

    flist = np.sort(glob.glob(os.path.join(gitm_dir, gitm_file_pattern)))
    if len(flist) == 0:
        raise ValueError("No %s files found in %s" %
                         (gitm_file_pattern, gitm_dir),
                         "\n \n instead there is: ",
                         glob.glob(os.path.join(gitm_dir, gitm_file_pattern)))

    gitm_dtimes = gitm_times_from_filelist(file_list=flist,
                                           century_prefix=century_prefix)

    if dtime_start is not None:
        start_idx_time = gitm_dtimes.index(dtime_start)
        gitm_dtimes = gitm_dtimes[start_idx_time:]
        flist = flist[start_idx_time:]

    if dtime_end is not None:
        end_idx_time = gitm_dtimes.index(dtime_end)
        gitm_dtimes = gitm_dtimes[:end_idx_time]
        flist = flist[:end_idx_time]

    if start_idx != 0 and end_idx != -1:
        gitm_dtimes = gitm_dtimes[start_idx:end_idx]
        flist = flist[start_idx:end_idx]
    elif start_idx != 0:
        gitm_dtimes = gitm_dtimes[start_idx:]
        flist = flist[start_idx:]
    elif end_idx != -1:
        gitm_dtimes = gitm_dtimes[:end_idx]
        flist = flist[:end_idx]

    f = ather_read.read_gitm_file(flist[0])
    if '3D' in gitm_file_pattern:
        gitmgrid = {f["vars"][k].lower(): f[k][2:-2, 2:-2, 2:-2]
                    for k in [0, 1, 2]}
        nlons, nlats, nalts = np.array(f[0].shape) - 4  # ghost cells
    elif '2D' in gitm_file_pattern:
        gitmgrid = {f["vars"][k].lower(): f[k]
                    for k in [0, 1, 2]}
        nlons, nlats, nalts = np.array(f[0].shape)  # NO ghost cells in 2d

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

    if progress_bar:
        pbar = tqdm(total=len(flist))
    for ifile, file_name in enumerate(flist):
        f = ather_read.read_gitm_file(file_name)

        for num_var, real_var in enumerate(gitmvars):
            num_v_src = f["vars"].index(real_var)
            if '3DALL' in gitm_file_pattern:
                gitmbins[ifile, num_var] = f[num_v_src][2:-2, 2:-2, 2:-2]
            elif '2DANC' in gitm_file_pattern:
                gitmbins[ifile, num_var] = f[num_v_src]
        if progress_bar:
            pbar.update()

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

    rets = {
        'gitmdtimes': gitm_dtimes,
        'gitmbins': gitmbins,
        'gitmgrid': gitmgrid}

    if return_vars:
        rets['gitmvars'] = gitmvars

    return rets


def read_bin_to_xarray(filename,
                       drop_ghost_cells=True,
                       cols='all'):
    """Reads GITM binary file into xarray
    Works for all GITM files, including 3DALL, 2DALL, 2DANC, etc.
    - (Taken and modified from aetherpy)

    Args:
        filename (str, path): Path to the file to read.
        drop_ghost_cells (bool, optional):
            Drop GITM ghost cells. See GITM manual for details
            on ghost cells. Defaults to True.
        cols (str/list-like, optional):
            Set which columns to read. On systems with limited memory
            this will make datasets too large to fit into memory.
            Defaults to 'all' (all columns).

    Raises:
        IOError: File does not exist

    Returns:
        xarray.Dataset: Dataset holding the data.
            Indexed with glat, glon, alt (converted to deg, deg, km)

    """

    if not os.path.isfile(filename):
        raise IOError('input file does not exist')

    if 'HME' in filename:
        raise ValueError('HME files not supported yet')

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
            var = var\
                .replace('!N', '').replace('!U', '').replace('!D', '')\
                .replace('[', '').replace('[', '').replace(']', '')\
                .replace('/', '-').replace('(', '_').replace(')', '')\
                .replace('+', '_plus').replace(',', '_').replace('__', '_')\
                .replace('-', '_')

            varnames.append(var)
            dummy, rec_lec = unpack(end_char + '2l', fin.read(8))

        # Extract time
        rec_time = np.array(unpack(end_char + 'lllllll', fin.read(28)))
        rec_time[-1] *= 1000  # convert from millisec to microsec
        time_here = datetime.datetime(*rec_time)

        # Header is this length:
        # Version + start/stop byte
        # nlons, nlats, nalts + start/stop byte
        # num_vars + start/stop byte
        # variable names + start/stop byte
        # time + start/stop byte

        iheader_length = 84 + num_vars * 48

        ntotal = nlons * nlats * nalts
        idata_length = ntotal * 8 + 8

        data_vars = {}

        # Save the data for the desired variables
        if '2D' in filename:
            dimnames = ['lon', 'lat']
            newshape = [nlons, nlats]
        else:
            dimnames = ['lon', 'lat', 'alt']
            newshape = [nlons, nlats, nalts]
        for ivar in file_vars:
            fin.seek(iheader_length + ivar * idata_length)
            sdata = unpack(end_char + 'l', fin.read(4))[0]

            if ivar == 0:
                lons = np.rad2deg(np.unique(np.array(
                    unpack(end_char + '%id' % (ntotal),
                           fin.read(sdata))).reshape(
                               newshape, order="F")))

            elif ivar == 1:
                lats = np.rad2deg(np.unique(np.array(
                    unpack(end_char + '%id' % (ntotal),
                           fin.read(sdata))).reshape(
                               newshape, order="F")))

            elif ivar == 2:
                alts = np.unique(np.array(
                    unpack(end_char + '%id' % (ntotal),
                           fin.read(sdata))).reshape(
                               newshape, order="F")) / 1000

            else:
                data_vars[varnames[ivar]] = dimnames, np.array(
                    unpack(end_char + '%id' % (ntotal),
                           fin.read(sdata))).reshape(
                               newshape, order="F")
                # break
    if '2D' in filename:
        ds = xr.Dataset(coords={
            'time': [time_here], 'lon': lons, 'lat': lats},
            data_vars=data_vars,
            attrs={'version': version,
                   'dropped-ghost-cells': str(drop_ghost_cells)})
    else:
        ds = xr.Dataset(coords={
            'time': [time_here], 'lon': lons, 'lat': lats, 'alt': alts},
            data_vars=data_vars,
            attrs={'version': version,
                   'dropped-ghost-cells': str(drop_ghost_cells)})
        if drop_ghost_cells:
            ds = ds.drop_isel(lat=[0, 1, -2, -1],
                              lon=[0, 1, -1, -2], alt=[0, 1, -1, -2])

    cols = [cols] if type(cols) is str else cols
    if cols != 'all' and 'all' not in cols:
        for c in cols:
            if c not in ds.data_vars:
                raise ValueError(
                    f'Column {c} not found. \n  Available columns: \n',
                    list(ds.data_vars))
        ds = ds.get(cols)

    return ds


def gitm_times_from_filelist(file_list, century_prefix='20'):
    """Generate datetimes from a list of GITM files.

    Args:
        file_list (list-like): list of gitm files to parse
        century_prefix (str, optional):
            Which century? Defaults to '20'.

    Raises:
        ValueError: Incorrect file format.

    Returns:
        list: List of datetimes in the same order as the filelist input.
    """
    gitm_dtimes = []
    for i in file_list:
        yy, MM, dd = i[-17:-15], i[-15:-13], i[-13:-11]
        hr, mm, sec = i[-10:-8], i[-8:-6], i[-6:-4]
        try:
            gitm_dtimes.append(
                datetime.datetime(
                    int(century_prefix + yy), int(MM), int(dd),
                    int(hr), int(mm), int(sec)))
        except ValueError:

            try:
                yy, MM, dd = i[-20:-18], i[-17:-15], i[-14:-12]
                hr, mm, sec = i[-11:-9], i[-8:-6], i[-5:-3]

                gitm_dtimes.append(
                    datetime.datetime(
                        int(century_prefix + yy), int(MM), int(dd),
                        int(hr), int(mm), int(sec)))

            except ValueError:
                raise ValueError(
                    "GITM file name does not match expected format,",
                    "filename %s cannot be parsed" % i)
    return gitm_dtimes


def read_multiple_bins_to_xarray(file_list,
                                 start_dtime=None,
                                 end_dtime=None,
                                 start_idx=0,
                                 end_idx=-1,
                                 drop_ghost_cells=True,
                                 cols='all',
                                 pbar=False):
    """Read a list-like of GITM files into an xarray Dataset.

    Args:
        file_list (list-like): files to pull from.
        start_dtime (datetime, optional):
            Time to start read at. Not necessary (especially if you have
             pre-filtered the file_list. Defaults to None.
        end_dtime (datetime, optional):
            Time to end reads at. See above. Can be used exclusively.
             Defaults to None.
        start_idx (int, optional): Index of file_list to start reading.
            Defaults to 0.
        end_idx (int, optional): Index of file_list to end reading at.
            Defaults to -1.
        drop_ghost_cells (bool, optional): Remove Ghost cells?
            Defaults to True.
        cols (str or list-like, optional): Specific columns to read.
            Defaults to 'all'.
        pbar (bool, optional): Whether or not to show progress bar.
            Requires tqdm. Defaults to False.

    Raises:
        ValueError: If start/end inputs are mixed up.

    Returns:
        xarray.Dataset:
            Dataset containing all variables in the file_list at the
             times specified.
    """

    # Check inputs! Cannot specify start time & idx:
    if start_dtime is not None and start_idx is not None:
        raise ValueError("Cannot specify both Start idx & dtime")
    if end_dtime is not None and end_idx is not None:
        raise ValueError("Cannot specify both End idx & dtime")

    file_list = file_list[start_idx:end_idx]

    if start_dtime is not None or end_dtime is not None:
        times = gitm_times_from_filelist(file_list)
    if start_dtime is not None:
        time_mask = np.where(times >= start_dtime)
        file_list = file_list[time_mask]
        times = times[time_mask]
    if end_dtime is not None:
        time_mask = np.where(times <= end_dtime)
        file_list = file_list[time_mask]

    if pbar:
        progress = tqdm(total=len(file_list))

    ds = []
    for file in file_list:
        ds.append(read_bin_to_xarray(file,
                                     drop_ghost_cells=drop_ghost_cells,
                                     cols=cols))
        if pbar:
            progress.update()

    ds = xr.merge(ds)
    return ds


def process_all_to_cdf(gitm_dir,
                       out_dir=None,
                       dtime_storm_start=None,
                       delete_bins=False,
                       replace_cdfs=False,
                       progress_bar=True,
                       drop_ghost_cells=True,
                       drop_before=None,
                       drop_after=None,
                       skip_existing=False,
                       file_types='all',
                       use_ccmc=True,
                       single_file=False,
                       run_name=None,
                       tmp_dir=None
                       ):
    """Process all GITM .bin files in a directory to .cdf files.

    Args:
        gitm_dir (str: path-like): Directory containing GITM .bin files.
        out_dir (str: path-like, optional): Directory to output .cdf files.
            If None, will go into the same directory as the .bin files.
            Defaults to None.
        dtime_storm_start (datetime, optional): Attribute added to the
            netCDF file. Defaults to None.
        delete_bins (bool, optional): Delete GITM bins after making Datasets?
            Defaults to False.
        replace_cdfs (bool, optional): Replace pre-existing netCDF files?
            Defaults to False.
        progress_bar (bool, optional):  Whether or not to show progress bar.
            Requires tqdm. Defaults to True. If outputting to a single file,
            a progress bar will be added when writing files to disk. This
            cannot be changed.
        drop_ghost_cells (bool, optional): Drop GITM ghost cells?
            Defaults to True.
        drop_before (datetime, optional): Similar to start_dtime.
            When to start processing files. Will delete files before this time.
            Defaults to None.
        drop_after (datetime, optional): Similar to start_dtime.
            When to start processing files. Will delete files before this time.
            Defaults to None.
        skip_existing (bool, optional): Skip existing netCDF files?
            Defaults to False. This will slow down the program significantly.
        file_types (str or list-like, optional): Which file types to process.
            Defaults to 'all'. Can be a list of strings or a single string.
            Example usage is ['3DALL', '2DALL'] or '3DALL'.
        use_ccmc (bool, optional): Write files with CCMC naming convention?
            Defaults to True. Recommended if not using single_file.
        single_file (bool, optional): Output to a single file?
            Defaults to False. If True, will output to a single netCDF file.
            If False, will output to multiple netCDF files, one for each time.
        run_name (str, optional): Name of the run. Only used if single_file.
            Defaults to None. '_GITM.nc' will be appended to this.
        tmp_dir (str, optional): Temporary directory to write files to.
            Only used if single_file. Defaults to None.
            *Some systems have a local temp directory that's much faster
            than the standard output_directory.*

    """

    if not drop_ghost_cells:
        print('Not dropping Ghost cells.',
              'This will cause issues if you are processing both 3D and 2D',
              'files. Not robust enough to deal with that, unfortunately.')
    else:
        print('dropping ghost cells')

    if out_dir is None:
        out_dir = gitm_dir

    if single_file and run_name is None:
        raise ValueError('You must set the run name if outputting'
                         ' to a single file')

    if file_types == 'all':
        files = np.sort(glob.glob(os.path.join(gitm_dir, '*.bin')))
    else:
        if isinstance(file_types, str):
            file_types = [file_types]
        files = []
        for f in file_types:
            files.extend(
                glob.glob(os.path.join(gitm_dir, '%s*.bin' % f)))
        files = np.sort(files)

    if drop_after is not None or drop_before is not None:
        times = np.array(gitm_times_from_filelist(files))
        b4 = len(files)
        if drop_after is not None:
            mask = np.where(times <= drop_after)
            times = times[mask]
            files = files[mask]
        if drop_before is not None:
            mask = np.where(times >= drop_before)
            times = times[mask]
            files = files[mask]
        print('filtered %i out of %i files by time' % (b4 - len(files), b4))

    indiv_ends = []

    for i in files:
        if i[-19:] not in indiv_ends:
            indiv_ends.append(i[-19:])

    num_existing_cdfs = len(glob.glob(os.path.join(out_dir, '*.nc')))
    if len(indiv_ends) == num_existing_cdfs and not replace_cdfs:
        print('looks like all files have already been processed.',
              'you can run again if you need to reprocess all netcdfs')
    elif num_existing_cdfs != 0 and not skip_existing and not single_file:
        import warnings
        warnings.warn(
            '\nThere are %i existing netcdfs in this directory,\n'
            ' but %i should be made. This may be because you\n'
            ' have already processed some of these files.\n'
            ' It is possible that some files have been deleted,\n'
            ' or that the writing of a file was interrupted.\n'
            ' If you want to skip existing files, set skip_existing=True\n'
            % (num_existing_cdfs, len(indiv_ends)))

    to_remove = []

    if single_file:
        # things we need to keep track of for single_file
        # existing vars cannot be replaced in a netcdf file so
        #   we have to write temp files, combine them, then delete them.

        if tmp_dir is None:
            tmp_dir = out_dir

        print('writing temp files to %s' % tmp_dir)

        files_written = []
        if not os.path.exists(os.path.join(tmp_dir, run_name + '_tmp')):
            os.makedirs(os.path.join(tmp_dir, run_name + '_tmp'))

    elif tmp_dir is not None:
        raise ValueError('tmp_dir is only used if single_file=True')

    if skip_existing:

        if tmp_dir is None:
            tmp_dir = out_dir
        
        files_written = glob.glob(os.path.join(tmp_dir, '*.nc'))
        indiv_ends = indiv_ends[len(files_written):]

    if progress_bar:
        if skip_existing and not single_file:
            pbar = tqdm(total=len(indiv_ends) - num_existing_cdfs)
        else:
            pbar = tqdm(total=len(indiv_ends), desc='Processing GITM')

    for fileend in indiv_ends:

        if file_types == 'all':
            files_here = glob.glob(gitm_dir + '/*' + fileend)
        else:
            files_here = []
            for filetype in file_types:
                files_here.extend(
                    glob.glob(
                        os.path.join(
                            gitm_dir,
                            filetype + '*' + fileend)))
        # make sure 3DALL files are first:
        if len(files_here) > 1:
            files_here = np.flip(np.sort(files_here))
        ds_now = []

        for f in files_here:
            ds_now.append(read_bin_to_xarray(
                filename=f,
                drop_ghost_cells=drop_ghost_cells,
                cols='all'))
            to_remove.append(f)

        ds_now = xr.merge(ds_now)

        if dtime_storm_start is not None:
            ds_now = ds_now.assign_attrs(
                dtime_event_start=dtime_storm_start,)

        if use_ccmc and not single_file:
            ut = ds_now.time.values[0]
            outfile = os.path.join(
                out_dir,
                make_ccmc_name('GITM', ut))

        elif single_file:
            outfile = os.path.join(tmp_dir, run_name + '_tmp',
                                   fileend[fileend.rfind('t'):
                                           ].replace('.bin', '.nc'))

        else:
            outfile = os.path.join(
                out_dir,
                fileend[fileend.rfind('t'):].replace('.bin', '.nc'))

        if skip_existing:
            if os.path.exists(outfile):
                if progress_bar:
                    pbar.update()
                if single_file:
                    files_written.append(outfile)
                continue

        if single_file:
            ds_now.to_netcdf(outfile, engine='h5netcdf')
            files_written.append(outfile)
        else:
            ds_now.to_netcdf(outfile, mode='w')

        if progress_bar:
            pbar.update()
    if progress_bar:
        pbar.close()

    if single_file:
        # read in all written files(with dask).
        # then write to a new netCDF file.
        # and then clean up.
        print('reading in temp files...')
        ds = xr.open_mfdataset(files_written, engine='h5netcdf',
                               concat_dim='time', combine='nested')

        print('writing... (this takes a while)')
        with ProgressBar():
            ds.to_netcdf(os.path.join(out_dir, run_name + '_GITM.nc'),
                         encoding={'time': {'dtype': float}})
        print('cleaning up temp files')
        for f in files_written:
            os.remove(f)
        os.removedirs(os.path.join(tmp_dir, run_name + '_tmp'))

    if delete_bins:
        print('FILES WILL BE DELETED. YOU HAVE BEEN WARNED.',
              ' 10 seconds to cancel.')
        import time
        time.sleep(10)
        for f in glob.glob(gitm_dir + '/*.bin'):
            os.remove(f)

    print('Done!')


def find_variable(gitm_dir, varname=None,
                  varhelp=False, nc=True,):
    """
    Help function. Finds a variable in a directory of GITM files.
        Return the filetype and/or all of the variables available.


    Args:
        gitm_dir (str: path-like): Directory of GITM files.
        varname (str, optional): Variable you're looking for. Not setting this
            will just print all variables.
            Defaults to None.
        varhelp (bool, optional): If True, will print out all vaiables
            available. Think of it as "just checking".
            Defaults to False.
        nc (bool, optional): Whether to only look at .nc files.
            Defaults to True.

    Raises:
        ValueError: If you don't specify either varhelp or varname.

    Returns:
        str (optional): The filetype holding the variable you're loooking for.

    """

    if not varhelp and varname is None:
        raise ValueError('Must specify either varhelp or varname')

    if nc:
        files = np.sort(glob.glob(os.path.join(gitm_dir, '*GITM*.nc')))
        if len(files) == 0:
            nc = False
            print('no netcdf files found, trying binary files')
    if nc:
        ftypes_checked = []
        ds = xr.open_dataset(files[0])
        if varname in list(ds.data_vars.keys()):
            print('Found %s in %s' % (varname, files[0]))

        else:
            print('Did not find %s in %s\n' % (varname, files[0]),
                  'Instead found: \n', list(ds.data_vars.keys()),
                  '\nmoving to bin files')

        ds.close()

    files = np.sort(glob.glob(os.path.join(gitm_dir, '*.bin')))
    if len(files) == 0:
        print('no binary files found, exiting')
    ftypes_checked = []
    for f in files:
        ftype = f.split('/')[-1][:5]
        if ftype not in ftypes_checked:
            ftypes_checked.append(ftype)
            binary = read_bin_to_nparrays(f)
            for col in binary['vars']:
                if col == varname:
                    if varhelp:
                        print('Found %s in %s' % (varname, f))
                    else:
                        return ftype
                else:
                    col = col.replace('!N', '').replace('!U', '')\
                        .replace('!D', '').replace('[', '')\
                        .replace('[', '').replace(']', '')\
                        .replace('/', '-')
                    if col == varname:
                        print('Found %s in %s' % (varname, f))


def auto_read(gitm_dir,
              single_file=False,
              start_dtime=None,
              start_idx=None,
              end_dtime=None,
              end_idx=None,
              cols='all',
              progress_bar=True,
              drop_ghost_cells=True,
              file_type=None,
              return_xarray=True,
              force_dict=False,
              parallel=True,
              engine='h5netcdf',
              use_dask=False):
    """
    Automatically reads in a directory of GITM files.

    :param gitm_dir: Directory of GITM files.
    :type gitm_dir: str
    :param single_file: Whether to read in a single file. Defaults to False.
    :type single_file: bool, optional
    :param start_dtime: Start time of the data you want. Defaults to None.
    :type start_dtime: datetime, optional
    :param start_idx: Start index of the data you want. Defaults to None.
    :type start_idx: int, optional
    :param end_dtime: End time of the data you want. Defaults to None.
    :type end_dtime: datetime, optional
    :param end_idx: End index of the data you want. Defaults to None.
    :type end_idx: int, optional
    :param cols: List of columns you want to read in. Defaults to 'all'.
    :type cols: list-like or str, optional
    :param progress_bar: Whether to show a progress bar. Defaults to True.
        Requires tqdm.
    :type progress_bar: bool, optional
    :param drop_ghost_cells: Whether to drop ghost cells. Defaults to True.
    :type drop_ghost_cells: bool, optional
    :param file_type: File type of the data you want to read in.
        Defaults to None.
    :type file_type: str, optional
    :param return_xarray: Whether to return an xarray. Defaults to True.
    :type return_xarray: bool, optional
    :param force_dict: Whether to force a dictionary return.
        Defaults to False.
    :type force_dict: bool, optional
    :param parallel: Whether to read in files in parallel. Defaults to True.
        This will use Dask, which can get hairy.
        If you're having issues, try setting this to False.
        Needs dask and dask.distributed
    :type parallel: bool, optional
    :param engine: The engine to use for reading in the data.
        Defaults to 'h5netcdf'.
    :type engine: str, optional
    :param use_dask: Whether to use Dask for reading in the data.
        Defaults to False.
    :type use_dask: bool, optional
    :return: The data read in from the GITM files.
    :rtype: xarray.Dataset or dict

    """

    if single_file:
        try:
            data = read_bin_to_xarray(
                filename=gitm_dir,
                drop_ghost_cells=drop_ghost_cells,
                cols=cols)
        except ValueError:
            data = read_bin_to_nparrays(
                filename=gitm_dir,
                drop_ghost_cells=drop_ghost_cells,
                cols=cols)
        return data

    files = np.sort(glob.glob(os.path.join(gitm_dir, '*GITM*.nc')))
    if len(files) == 0 and force_dict:
        if not force_dict:
            print("""No NetCDF files found, You should probably convert
                  from '.bin' to '.nc' first! (use process_all_to_cdf)\n
                  Continuing with your read...""")

        if cols != 'all' and file_type is None:
            file_type = find_variable(gitm_dir, varname=cols[0], nc=False)
        elif file_type is not None:
            files = np.sort(
                glob.glob(os.path.join(gitm_dir, file_type + '*.bin')))
        else:
            files = np.sort(
                glob.glob(os.path.join(gitm_dir, '3DALL' + '*.bin')))
            print('Defaulting to 3DALL files.')
        if return_xarray:
            ds = read_multiple_bins_to_xarray(
                file_list=files,
                start_dtime=start_dtime,
                start_idx=start_idx,
                end_dtime=end_dtime,
                end_idx=end_idx,
                drop_ghost_cells=drop_ghost_cells,
                cols=cols,
                pbar=progress_bar)
            return ds
        else:
            datadict = read_bin_to_nparrays(
                gitm_dir=gitm_dir,
                start_dtime=start_dtime,
                start_idx=start_idx,
                end_dtime=end_dtime,
                end_idx=end_idx,
                cols=cols,
                progress_bar=progress_bar)
            return datadict
    else:
        if file_type is not None:
            raise ValueError('Cannot specify file_type if using NetCDF files.')

        if start_idx is not None and end_idx is not None:
            files = files[start_idx:end_idx]
        elif start_idx is not None:
            files = files[start_idx:]
        elif end_idx is not None:
            files = files[:end_idx]

        if start_dtime is not None and end_dtime is not None:
            files = files[(np.array(gitm_times_from_filelist(files))
                           >= start_dtime) &
                          (np.array(gitm_times_from_filelist(files))
                          <= end_dtime)]

        if isinstance(cols, str):
            cols = [cols]

        if use_dask:

            ds = xr.open_mfdataset(files, parallel=parallel,
                                   combine_attrs='drop_conflicts',
                                   data_vars=cols,
                                   # concat_dim="time", combine="nested",
                                   coords='minimal', compat='override',
                                   engine=engine)

        else:
            drops = []
            ds0 = xr.open_dataset(files[0])
            for v in ds0.data_vars:
                if v not in cols:
                    drops.append(v)
            del ds0
            dss = []
            for f in files:
                dss.append(xr.open_dataset(
                    f, drop_variables=drops, engine=engine))
            ds = xr.concat(dss, dim='time')
            del dss

        return ds
