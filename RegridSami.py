"""
Postprocess SAMI data to be read into xarrays.

- This is just an entrypoint to the interpolation functions
  in the utility_programs folder (utility_programs/interpolate_outputs.py)
- These functions can be called from the command line or from other scripts.

* Sometimes the interpolations can get a little wonky, especially if the
  grid is too coarse. If you see weird artifacts in the output, try
  increasing the grid resolution, or runnnig the interpolations manually in
  utility_programs/interpolate_outputs.py

"""


import os
import numpy as np
import pandas as pd
from utility_programs.interpolate_outputs import do_interpolations
from utility_programs.utils import str_to_ut
import argparse


def main(
        sami_data_path,
        out_path=None,
        save_weights=True,
        cols='all',
        dtime_sim_start=None,
        lat_step=1,
        lon_step=4,
        alt_step=25,
        minmax_alt=[100, 2600],
        out_coord_file=None,
        sami_mintime=0,
        run_name=None,
        skip_time_check=False,
        progress_bar=True,
        num_workers=16):
    """
    Interpolate SAMI3 outputs to a new grid.

    Parameters
    ----------
    sami_data_path : str or os.pathLike
        Path to SAMI3 output files.
    out_path : str or os.pathLike, optional 
        Location to save output files. Defaults to sami_data_path.
    save_weights : bool, optional
        Whether or not to save Delauney Triangulation. Defaults to True.
    cols : str or list-like, optional
        Columns to read data from. Can be string or list-like. 
        Defaults to 'all'.
    dtime_sim_start : datetime, optional
        Datetime of simulation start. Required to read raw SAMI outputs.
        Defaults to None.
    lat_step : int, optional
        Integer step to use in output grid (deg). Defaults to 2.
    lon_step : int, optional
        Integer step to use in output grid (deg). Defaults to 4.
    alt_step : int, optional
        Integer step to use in output grid (in km). Defaults to 50.
    minmax_alt : list, optional
        Min & Max altitude to output in km. Defaults to [100, 2200].
    out_coord_file : str or os.pathLike, optional
        Output coordinates from a file instead of using the default grid.
        Defaults to None (no coordinate file).
        ** MUST HAVE "time, lat, lon, alt" COLUMNS **
        Use this to specify a user-defined grid or to interpolate to a set of
        satellite coordinates.
    sami_mintime : int, optional
        Minimum time to read in SAMI data. Defaults to 0. Use this to skip the
        first few hours of SAMI data and save time & memory.
    run_name : str optional
        Name of the run, used to name the output file.
        run_name='test' will generate a file called 'test_SAMI_REGRID.nc'.
        Defaults to None.
    skip_time_check : bool, optional
        Skip checking if the time range is valid. SAMI can sometimes output
        fake timesteps, or be configured to start outputting data after
        several hours. This will skip checking that the times are valid.
        Defaults to False.
    progress_bar : bool, optional
        Show progress bar? Defaults to True.
    num_workers : int, optional
        Number of workers to use when interpolating. Defaults to 16.
        (16 workers => 1.3 GB of RAM/10 time-steps of SAMI at 80/72/256
        resolution)

    Raises
    ------
    ValueError
        If out_coord_file does not have the required variables.

    Returns
    -------
    None

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
                          num_workers=num_workers,
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
            num_workers=num_workers,
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

    parser.add_argument('--sami_mintime', type=int, default=0,
                        help='Set this to help with memory management (especially '
                        'when using aarons_mod). Setting this to 5, for example, '
                        'will skip the first 5 output times of the SAMI run. '
                        'This is relevant even without aarons_mod as the first 12-24 '
                        'hours of SAMI outputs are not to be trusted.')

    parser.add_argument('--custom_grid', action='store_true', default=False,
                        help='Launches interactive script to set custom grid.'
                        'Default grid is 4 x 1 deg lonxlat, 50 alts.'
                        'minimum alt is 100, max is 2200, global lon x lat')

    parser.add_argument('--input_coord_file', type=str, default=None,
                        help='path to file containing input coordinates.'
                        ' Input coordinates must be in a CSV file with columns'
                        ' [lat, lon, alt] (in degrees and km.)')

    parser.add_argument('--num_workers', type=int, default=16,
                        help='When doing a regrid of the SAMI data, we need '
                        'to do a lot of calculations. By default this will use '
                        '16 workers, but you can change it if you want. '
                        "(higher workers = faster, to a point... "
                        "The number concurrent procs can be limited for "
                        "*reasons*. Email me if this is an issue & we'll chat"
                        # '16 workers => 1.3 GB of RAM/10 time-steps of SAMI '
                        # 'at 80/72/256 resolution)'
                        )

# parser.add_argument('--aarons_mod', action='store_true',
#                     help='Interpolating SAMI data is hard. Through a lot of '
#                     'testing, I found that you can interpolate to 2x '
#                     'resolution and then coarsen the results and things look '
#                     'normal. If you notice bad stuff in the regridded files, '
#                     'use this option. \n'
#                     'This will take a lot memory and time than default runs, '
#                     'so be careful!')

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
             skip_time_check=args.skip_time_check,
             sami_mintime=args.sami_mintime,
             num_workers=args.num_workers)

    elif args.input_coord_file is None:

        main(args.sami_data_path,
             out_path=args.out_path,
             save_weights=args.save_weights,
             cols=args.cols,
             run_name=args.run_name,
             dtime_sim_start=dtime_sim_start,
             sami_mintime=args.sami_mintime,
             skip_time_check=args.skip_time_check,
             num_workers=args.num_workers)

    else:

        main(args.sami_data_path,
             out_path=args.out_path,
             save_weights=args.save_weights,
             cols=args.cols,
             run_name=args.run_name,
             dtime_sim_start=dtime_sim_start,
             out_coord_file=args.input_coord_file,
             sami_mintime=args.sami_mintime,
             skip_time_check=args.skip_time_check,
             num_workers=args.num_workers)
