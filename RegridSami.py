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
        alt_step=50,
        minmax_alt=[150, 2200],
        out_coord_file=None,
        sami_mintime=0,
        run_name=None,
        skip_time_check=False,
        progress_bar=True,):
    """
    Interpolate SAMI3 outputs to a new grid.

    :param sami_data_path: Path to SAMI3 output files.
    :type sami_data_path: str or os.pathLike
    :param out_path: Location to save output files. Defaults to
        sami_data_path.
    :type out_path: str or os.pathLike, optional
    :param save_weights: Whether or not to save Delauney Triangulation.
        Defaults to True.
    :type save_weights: bool, optional
    :param cols: Columns to read data from. Can be string or list-like.
        Defaults to 'all'.
    :type cols: str or list-like, optional
    :param dtime_sim_start: Datetime of simulation start. Required to read raw
        SAMI outputs. Defaults to None.
    :type dtime_sim_start: datetime
    :param lat_step: Integer step to use in output grid (deg). Defaults to 2.
    :type lat_step: int, optional
    :param lon_step: Integer step to use in output grid (deg). Defaults to 4.
    :type lon_step: int, optional
    :param alt_step: Integer step to use in output grid (in km).
        Defaults to 50.
    :type alt_step: int, optional
    :param minmax_alt: Min & Max altitude to output in km. Defaults to
        [100, 2200].
    :type minmax_alt: list, optional
    :param out_coord_file: Output coordinates from a file instead of using the
        default grid. Cannot be used with finerinterps. Defaults to None (no
        coordinate file). ** MUST HAVE "time, lat, lon, alt" COLUMNS **
        Use this to specify a user-defined grid or to interpolate to a set of
        satellite coordinates.
    :type out_coord_file: str or os.pathLike, optional
    :param sami_mintime: Minimum time to read in SAMI data. Defaults to 0. Use
        this to skip the first few hours of SAMI data and save time & memory.
    :type sami_mintime: int, optional
    :param run_name: Name of the run, used to name the output file.
        run_name='test' will generate a file called 'test_SAMI_REGRID.nc'.
        Defaults to None.
    :type run_name: Optional[str]
    :param skip_time_check: Skip checking if the time range is valid. SAMI can
        sometimes output fake timesteps, or be configured to start outputting
        data after several hours. This will skip checking that the times are
        valid. Defaults to False.
    :type skip_time_check: bool, optional
    :param progress_bar: Show progress bar? Defaults to True.
    :type progress_bar: bool, optional
    :raises ValueError: If out_coord_file does not have the required variables.
    :return: None
    :rtype: None
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
