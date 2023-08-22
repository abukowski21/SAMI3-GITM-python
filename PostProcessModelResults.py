"""
Process GITM & SAMI data to NetCDF format.

- Can use just one model output, if preferred.
- More functionality is available in the individual model modules.
- This program will process every column into netCDF files by time.
- SAMI is regridded according to the default values in RegridSami.main()
  - This can be adjusted with a custom grid, or use the individual model's
    post-processing routines for more fine control.

## TODO:
- Automatically process all model outputs from a given directory.

"""

import os
import glob
import subprocess
from utility_programs.read_routines import GITM, SAMI
from utility_programs.utils import str_to_ut
import RegridSami
import argparse


def main(args):

    psami = False
    pgitm = False

    if args.sami_dir != './sami_dir':
        if len(glob.glob(os.path.join(args.sami_dir, '*.dat'))) != 0:
            psami = True
    if args.gitm_dir != './gitm_dir':
        if len(glob.glob(os.path.join(args.gitm_dir, '*.bin'))) != 0:
            pgitm = True

    if args.dtime_sim_start is not None:
        args.dtime_sim_start = str_to_ut(args.dtime_sim_start)

    if args.dtime_event_start is not None:
        args.dtime_event_start = str_to_ut(args.dtime_event_start)

    if pgitm is None and psami is None:
        raise ValueError(
            'You must specify at least one model output directory.')

    if args.output_dir is None:
        output_dir = os.path.join(os.getcwd(), 'OUTPUTS')
    else:
        output_dir = args.output_dir

    if not os.path.exists(output_dir):
        print('making output directory: {}'.format(output_dir))
        os.makedirs(output_dir)

    if pgitm:

        header_files = glob.glob(os.path.join(args.gitm_dir, '*.header'))
        if len(header_files) > 0:
            print('GITM headers found in {}'.format(args.gitm_dir))
            print('Attempting to postprocess...')
            print('(This is not very robust.)')
            gitm_parent_dir = args.gitm_dir[:args.gitm_dir.rfind('/data')]

            cmd = os.path.join('.', gitm_parent_dir, 'pGITM')
            print('Running: {}'.format(cmd))
            if args.verbose:
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
                for line in p.stdout:
                    print(line)
                p.wait()
                print(p.returncode)
            else:
                subprocess.Popen(cmd, shell=True).wait()

            if len(glob.glob(os.path.join(args.gitm_dir, '*.bin'))) == 0:
                raise ValueError(
                    'No GITM files found in {}'.format(args.gitm_dir),
                    'Double check the directory or go run pGITM.')

        if len(glob.glob(os.path.join(output_dir,
                                      'GITM*.nc'))) > 0 or\
            len(glob.glob(os.path.join(output_dir,
                                       args.single_file + '*GITM.nc'))) > 0:
            if args.replace:
                print('Replacing existing netCDF files...')

            else:
                raise ValueError(
                    "Postprocessed GITM files already exist in: {}\n"
                    "Run with --replace to overwrite.\n"
                    "   Contents: \n    {}".format(
                        output_dir, os.listdir(output_dir)))

        GITM.process_all_to_cdf(
            gitm_dir=args.gitm_dir,
            out_dir=output_dir,
            delete_bins=args.delete_bins,
            replace_cdfs=args.replace,
            drop_ghost_cells=args.ghost_cells,
            progress_bar=args.progress,
            use_ccmc=args.ccmc,
            file_types=args.gitm_types,
            single_file=True if args.single_file else False,
            run_name=args.single_file if args.single_file else None)

    if psami:

        existing_sami_files = glob.glob(os.path.join(output_dir, 'SAMI*.nc'))

        do_write_raw = False
        do_write_regrid = False

        if args.sami_type == 'all' or args.sami_type == 'regrid':
            do_write_regrid = True
        if args.sami_type == 'all' or args.sami_type == 'raw':
            do_write_raw = True

        if len(existing_sami_files) > 0:
            if args.replace:
                print('Replacing existing netCDF files...')
                print(' not implemented yet. Go clear directory manually ')
            else:
                if 'RAW' in str(existing_sami_files) and do_write_raw:
                    raise ValueError(
                        'RAW files already exist in output_dir. You may want'
                        ' to set --replace, or delete them.')
                if 'REGRID' in str(existing_sami_files) and do_write_regrid:
                    raise ValueError(
                        'REGRID files already exist in output_dir. You may'
                        ' want to set --replace, or delete them.')

        if args.ccmc and not args.single_file:
            use_ccmc = True
            split_by_time = True
            split_by_var = False
        else:
            use_ccmc = False
            split_by_time = False
            if not args.single_file:
                split_by_var = True
            else:
                split_by_var = False

        if do_write_raw:

            print('Attempting to convert raw -> netCDF...')

            if args.dtime_sim_start is None:
                raise ValueError(
                    'You must specify a start time for the SAMI simulation,',
                    'in the format YYYYMMDDHHMMSS')

            SAMI.process_all_to_cdf(
                sami_data_path=args.sami_dir,
                out_dir=output_dir,
                use_ccmc=use_ccmc,
                split_by_time=split_by_time,
                split_by_var=split_by_var,
                dtime_sim_start=args.dtime_sim_start,
                progress_bar=args.progress,
                OVERWRITE=args.replace,
                low_mem=args.low_mem,
                delete_raw=args.delete_bins,
                dtime_storm_start=args.dtime_event_start,
                skip_time_check=args.skip_time_check,
                whole_run=True if args.single_file else False,
                run_name=args.single_file if args.single_file else None)

        if do_write_regrid:
            print('attempting to regrid!')
            try:
                # ignore flake8 linting on the following line.
                import numba  # noqa: F401
                numba_installed = True
                print('numba installed! Will speed up regridding.')
            except ImportError:
                numba_installed = False
                print('Module numba not found.')

            weights_exist = False
            if args.reuse_weights:
                if os.path.exists(os.path.join(output_dir, 'weights')):
                    weights_exist = True
                    print('found weights to reuse')
                else:
                    print('No existing weight file found')

            if args.set_custom_grid:
                RegridSami.main(sami_data_path=args.sami_dir,
                                out_path=output_dir,
                                save_weights=weights_exist,
                                use_saved_weights=weights_exist,
                                dtime_sim_start=args.dtime_sim_start,
                                lat_step=latstep,
                                lon_step=lonstep,
                                alt_step=altstep,
                                minmax_alt=[minalt, maxalt],
                                lat_finerinterps=latfiner,
                                lon_finerinterps=lonfiner,
                                alt_finerinterps=altfiner,
                                use_ccmc=args.ccmc,
                                split_by_time=args.ccmc,
                                split_by_var=not (
                                    args.ccmc and args.single_file),
                                numba=numba_installed and not args.low_mem,
                                skip_time_check=args.skip_time_check,
                                whole_run=True if args.single_file else False,
                                run_name=args.single_file if args.single_file
                                else None)

            else:
                RegridSami.main(
                    sami_data_path=args.sami_dir,
                    out_path=output_dir,
                    save_weights=weights_exist,
                    use_saved_weights=weights_exist,
                    dtime_sim_start=args.dtime_sim_start,
                    use_ccmc=args.ccmc,
                    split_by_time=args.ccmc,
                    split_by_var=not (
                        args.ccmc and args.single_file),
                    numba=numba_installed and not args.low_mem,
                    skip_time_check=args.skip_time_check,
                    whole_run=True if args.single_file else False,
                    run_name=args.single_file if args.single_file else None)

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-gitm', '--gitm_dir', type=str, default='./gitm_data',
                        help='GITM Directory. Defaults to ./gitm_data')
    parser.add_argument('-sami', '--sami_dir', type=str, default='./dami_dir',
                        help='SAMI directory. Defaults to ./sami_data')
    parser.add_argument('-out', '--output_dir', type=str, default=None,
                        help='If you want to save the files to another'
                        ' directory, specify it here. Defaults to a new'
                        ' folder at "./OUTPUTS/".')
    parser.add_argument('--sami_type', type=str, default='all',
                        help='Which SAMI data to process? (Default: all)'
                        '(Options: "all", "raw", "regrid")')
    parser.add_argument('--gitm_types', type=str, default='all',
                        nargs='*', help='Which GITM data to process?'
                        ' (EX: 3DALL, 3DNEU, etc.) (Default: all)')
    parser.add_argument('--single_file', type=str, default=False,
                        help='Set this to the run name to output the entire'
                        ' model run data to a single netCDF file.'
                        ' Note: model name will be added automatically')
    parser.add_argument('--set_custom_grid', type=bool, default=False,
                        help='Set a custom grid for SAMI regridding?'
                        ' Default: False')
    parser.add_argument('--low_mem', type=bool, default=True,
                        help='Process SAMI files in low memory mode?'
                        ' (NOTE: Memory usage is still 30GB+, without lowmem'
                        ' the entire run is read in at once.) Default: True')
    parser.add_argument('-c', '--ccmc', action='store', type=bool,
                        default=True,
                        help='Use CCMC naming conventions? (Default: True)')
    parser.add_argument('-r', '--replace', action='store_true',
                        help='Replace existing netCDF files? (Default: False)'
                        ' (not implemented yet)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print out more information? (Default: False)')
    parser.add_argument('--reuse_weights', action='store_true',
                        default=False,
                        help='Reuse weight & index file? (default False \n'
                        ' If True, weights will be re-used if they exist.'
                        ' If the file is not found, a new one will be written'
                        ' to the --output_dir')
    parser.add_argument('-d', '--delete_bins', action='store_true',
                        help='Delete binary files after processing? '
                        'Caution: This is irreversible! (Default: False)')
    parser.add_argument('-g', '--ghost_cells', action='store_false',
                        help='Retain GITM Ghost Cells? (Default: True).'
                        ' Not fully implemented yet.')
    parser.add_argument('--skip_time_check', action='store_true',
                        help='Skip verifying accuracy of times. Useful when'
                        ' SAMI has been configured to skip some outputs '
                        '(hrpr != 0)')
    parser.add_argument('--no_progress', action='store_false', dest='progress',
                        help='Show progress bar? (Default: True)'
                        ' - recommended since things can take a LONG time.'
                        ' Requires tqdm')
    parser.add_argument('--dtime_sim_start', type=str, default=None,
                        help='Start time of the simulation, in the format: '
                        'YYYYMMDDHHmmSS Required to process SAMI from *.dat')
    parser.add_argument('--dtime_event_start', type=str, default=None,
                        help='Event (storm) start time, in the format: '
                        'YYYYMMDDHHmmSS (Optional. added as attr to netCDFs)')

    args = parser.parse_args()

    opts = ['all', 'raw', 'regrid']
    # make sure args.sami_type is one of opts
    if args.sami_type not in opts:
        raise ValueError('sami_type must be one of {}'.format(opts))\

    if args.set_custom_grid and (args.sami_type == 'regrid'
                                 or args.sami_type == 'all'):

        global latstep, lonstep, altstep, minalt, maxalt
        global latfiner, lonfiner, altfiner
        print('We are going to set a custom grid for your sami regridding. ')
        latstep = input('latitude step size in degrees (1):')
        lonstep = input('longitude step size in degrees: (4):')
        altstep = input('altitude step size in km (50):')
        minalt = input('minimum altitude in km (100):')
        maxalt = input('maximum altitude in km (2200):')
        print('Now for the options to interpolate at a finer resolution'
              ' and then coarsen afterwards. If you dont know what this'
              ' means you can run with 1s and it will be faster. if you'
              ' see weird artifacts in your outputs you can try '
              ' adjusting this. Number given multiplies the step size')
        latfiner = input('interpolate a finer resolution in latitude? (1):')
        lonfiner = input('interpolate a finer resolution in longitude? (1):')
        altfiner = input('interpolate a finer resolution in altitude? (1):')

    elif args.set_custom_grid and args.sami_type == 'raw':
        raise ValueError('You cannot set a custom grid for raw SAMI files,'
                         ' since nothing is being regridded.')

    main(args)
