"""
Postprocess SAMI data to be read into xarrays.

- This is just an entrypoint to the interpolation functions
    in the utility_programs folder (utility_programs/interpolate_outputs.py)
- These functions can be called from the command line or from other scripts.

"""


import os
import numpy as np
import pandas as pd
from test_interpolate_outputs import do_interpolations
from utility_programs.utils import str_to_ut
from utility_programs.interpolate_outputs import do_interpolations
import argparse


def main(
        sami_data_path,
        out_path=None,
        save_weights=True,
        cols='all',
        dtime_sim_start=None,
        lat_step=1,
        lon_step=4,
        alt_step=50,
        minmax_alt=[150, 2200],
        out_coord_file=None,
        sami_mintime=0,
        run_name=None,
        skip_time_check=False,
        progress_bar=True,):
    """Interpolate SAMI3 outputs to a new grid.

    Args:
        sami_data_path (str, path-like):
            Path to SAMI3 output files.
        out_path (str; path-like, optional):
            Location to save output files. Defaults to "./"
        save_weights (bool, optional):
            Whether or not to save weight and index files. Defaults to True.
        use_saved_weights (bool, optional):
            Read in existing weight & index files (if they exist). Will not
                error if files don't exist, just create them.
                Defaults to True.
        apply_weights (bool, optional):
            Whether or not to apply weights. Otherwise weight files are made
                and the program exits. Defaults to True.
        cols (str or list-like, optional):
            Columns to read data from. Can be string or list-like.
                Defaults to 'all'.
        dtime_sim_start (datetime.datetime, optional):
            Datetime of simulation start. Required to read raw SAMI outputs.
                Defaults to None.
        out_coord_file (_type_, optional):
            Output coordinates from a file instead of using the default grid.
            Cannot be used with finerinterps.
            Defaults to None (no coordinate file).
            ** MUST HAVE "time, lat, lon, alt" COLUMNS **
        lat_step (int, optional):
            Integer step to use in output grid (deg). Defaults to 2.
        lon_step (int, optional):
            Integer step to use in output grid (deg). Defaults to 5.
        alt_step (int, optional):
            Integer step to use in output grid (in km). Defaults to 50.
        minmax_alt (list, optional):
            Min & Max altitude to output in km. Defaults to [100, 2200].
        lat_finerinterps (int, optional):
            Finer interpolation factor to use in latitude.
            Will by coarsened to the desired output resolution before writing.
            Defaults to 3.
        lon_finerinterps (int, optional):
            Finer interpolation factor to use in longitude.
            Defaults to 3.
        alt_finerinterps (int, optional):
            Finer interpolation factor to use in altitude. Defaults to 2.
        sami_mintime (int, optional):
            Minimum time to read in SAMI data. Defaults to 0.
            Use this to skip the first few hours of SAMI data and save
            time & memory.
        split_by_var (bool, optional):
            Write each variable to a separate file. Defaults to False.
        single_file (bool, optional):
            Write all variables to a single file. Defaults to False.
        split_by_time (bool, optional):
            Write each time to a separate file. Defaults to True.
        use_ccmc (bool, optional):
            Whether or not to use CCMC naming conventions for outputs. Must be
                used with split_by_time. Defaults to False.
        numba (bool, optional):
            Use numba to speed up applying weights. Defaults to False.

    Raises:
        ValueError: If out_coord_file does not have the reuired variables.
        ValueError: If split_by_var and single_file are both True.
    """

    if out_path is None:
        out_path = sami_data_path

    if not os.path.exists(out_path):
        print('making outpath:  %s' % out_path)
        os.makedirs(out_path)

    if out_coord_file is None:
        latout = np.arange(-90, 90, lat_step)
        lonout = np.arange(0, 360, lon_step)
        altout = np.arange(minmax_alt[0], minmax_alt[1] + 1, alt_step, )

        out_lats = []
        out_lons = []
        out_alts = []

        out_lats, out_lons, out_alts = np.meshgrid(latout, lonout, altout)

        do_interpolations(sami_data_path=sami_data_path,
                          dtime_sim_start=dtime_sim_start,
                          out_path=out_path,
                          save_delauney=save_weights,
                          out_lat_lon_alt=np.array(
                              [out_lats, out_lons, out_alts]),
                          is_grid=True,
                          out_runname=run_name,
                          cols=cols,
                          sami_mintime=sami_mintime,
                          skip_time_check=skip_time_check,
                          show_progress=progress_bar,
                          )
    else:
        # read in the file
        # try:
        coord_df = pd.read_csv(out_coord_file, sep=',')
        out_time = pd.to_datetime(coord_df['time']).values
        out_latlonalt = coord_df[['lat', 'lon', 'alt']].values.T

        do_interpolations(
            sami_data_path=sami_data_path,
            dtime_sim_start=dtime_sim_start,
            out_path=out_path,
            save_delauney=save_weights,
            out_lat_lon_alt=out_latlonalt,
            sat_times=out_time,
            cols=cols,
            sami_mintime=sami_mintime,
            skip_time_check=skip_time_check,
            show_progress=progress_bar,
            is_grid=False,
            out_runname=run_name if run_name is not None else 'custom',
        )


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('sami_data_path', type=str, help='path to sami data')

    parser.add_argument('-d', '--dtime_sim_start', type=str, default=None,
                        help='datetime string of simulation start time.'
                        'Required to read raw SAMI files')

    parser.add_argument('-o', '--out_path', type=str, default=None,
                        help='path to save regridded data,'
                        'defaults to sami_data_path')

    parser.add_argument('-r', '--run_name', type=str, default=None,
                        help='name of run to save outputs as. Defaults to'
                        ' SAMI-REGRID or SAT-INTERP, depending if the output'
                        ' is a grid or not')

    parser.add_argument('--save_weights', action='store_true',
                        help='Save/reuse Delauney Triangulation? Saves time.')

    parser.add_argument('--cols', type=str, default='all', nargs='*',
                        help='columns to read from sami data. Defaults to all'
                        'input as a list with spaces between.')

    parser.add_argument('--skip_time_check', action='store_true',
                        help='Skip verifying accuracy of times. Useful when'
                        ' SAMI has been configured to skip some outputs '
                        '(hrpr != 0)')

    parser.add_argument('--custom_grid', action='store_true', default=False,
                        help='Launches interactive script to set custom grid.'
                        'Default grid is 4 x 1 deg lonxlat, 50 alts.'
                        'minimum alt is 100, max is 2200, global lon x lat')

    parser.add_argument('--input_coord_file', type=str, default=None,
                        help='path to file containing input coordinates.'
                        ' Input coordinates must be in a CSV file with columns'
                        ' [lat, lon, alt] (in degrees and km.)')

    args = parser.parse_args()

    if args.dtime_sim_start is None:
        raise ValueError('You must specify a simulation start time'
                         ' when using these programs to read data')
    else:
        dtime_sim_start = str_to_ut(args.dtime_sim_start)

    if args.custom_grid:
        print('We are going to set a custom grid for you. ')
        print('Press enter to accept the default values in parentheses')
        latstep = int(input('latitude step size in degrees (1):', default=1))
        lonstep = int(input('longitude step size in degrees: (4):', default=4))
        altstep = int(input('altitude step size in km (50):', default=50))
        minalt = int(input('minimum altitude in km (100):', default=100))
        maxalt = int(input('maximum altitude in km (2200):', default=2200))
        # print('Now for the options to interpolate at a finer resolution'
        #       ' and then coarsen afterwards. If you dont know what this'
        #       ' means you can run with 1s and it will be faster. if you'
        #       ' see weird artifacts in your outputs you can try '
        #       ' adjusting this. Number given multiplies the step size')
        # latfiner = int(input(
        #     'interpolate a finer resolution in latitude? (1):',
        #     default=1))
        # lonfiner = int(input(
        #     'interpolate a finer resolution in longitude? (1):',
        #     default=1))
        # altfiner = int(input(
        #     'interpolate a finer resolution in altitude? (1):',
        #     default=1))

        main(args.sami_data_path,
             out_path=args.out_path,
             save_weights=args.save_weights,
             cols=args.cols,
             run_name=args.run_name,
             dtime_sim_start=dtime_sim_start,
             lat_step=latstep,
             lon_step=lonstep,
             alt_step=altstep,
             minmax_alt=[minalt, maxalt],
             skip_time_check=args.skip_time_check)

    elif args.input_coord_file is None:

        main(args.sami_data_path,
             out_path=args.out_path,
             save_weights=args.save_weights,
             cols=args.cols,
             run_name=args.run_name,
             dtime_sim_start=dtime_sim_start,
             skip_time_check=args.skip_time_check)

    else:

        main(args.sami_data_path,
             out_path=args.out_path,
             save_weights=args.save_weights,
             cols=args.cols,
             run_name=args.run_name,
             dtime_sim_start=dtime_sim_start,
             out_coord_file=args.input_coord_file,
             skip_time_check=args.skip_time_check)
