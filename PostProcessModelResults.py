"""
Process GITM & SAMI data to NetCDF format.

- Can use just one model output, if preferred.
- More functionality is available in the individual model modules.
- This program will process every column into netCDF files by time.

"""

import os
import glob
import subprocess
from utility_programs.read_routines import GITM, SAMI
from utility_programs.utils import str_to_ut
import argparse
import time


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

    need_output_dir = False
    if args.output_dir is not None:
        if pgitm:
            d = os.path.join(args.output_dir, 'gitm')
            if not os.path.exists(d):
                print('making output directory: {}'.format(d))
                os.makedirs(d)
        if psami:
            d = os.path.join(args.output_dir, 'sami')
            if not os.path.exists(d):
                print('making output directory: {}'.format(d))
                os.makedirs(d)
    else:
        need_output_dir = True


    if pgitm:
        if need_output_dir:
            output_dir = args.gitm_dir
        else:
            output_dir = os.path.join(args.output_dir, 'gitm')

        gitm_files = glob.glob(os.path.join(args.gitm_dir, '*.bin'))
        if len(gitm_files) == 0:
            print('No GITM files found in {}'.format(args.gitm_dir))
            print('Attempting to postprocess...')
            gitm_parent_dir = args.gitm_dir[:args.gitm_dir.rfind('/')]

            cmd = './'+gitm_parent_dir+'/pGITM'
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
                    'And I could not postprocess myself. Go run pGITM.')

        if len(glob.glob(os.path.join(args.output_dir, '*.nc'))) > 0:
            if args.replace:
                print('Replacing existing netCDF files...')

            else:
                raise ValueError(
                    'netCDF files already exist in {}'.format(args.gitm_dir),)

        GITM.process_all_to_cdf(
            gitm_dir=args.gitm_dir,
            out_dir=output_dir,
            delete_bins=args.delete_bins,
            replace_cdfs=args.replace,
            drop_ghost_cells=args.ghost_cells,
            progress_bar=args.progress,
        )

    if psami:
        if need_output_dir:
            output_dir = args.sami_dir
        else:
            output_dir = os.path.join(args.output_dir, 'sami')

        sami_files_nc = glob.glob(os.path.join(output_dir, '*.nc'))

        process_from_scratch = True
        regrid = True
        if len(sami_files_nc) != 0:

            if not args.replace:
                process_from_scratch = False
            for f in sami_files_nc:
                if 'grid' in f and args.dont_regrid:
                    regrid = False

            if not regrid and not process_from_scratch:
                raise ValueError(
                    'Looks like Im not doing anything here. Exiting.')

        if process_from_scratch:
            print(
                'No netCDF files found in samidir: {}'.format(args.sami_dir))
            print('Attempting to convert raw -> netCDF...')

            if args.dtime_sim_start is None:
                raise ValueError(
                    'You must specify a start time for the SAMI simulation,',
                    'in the format YYYYMMDDHHMMSS')

            SAMI.process_all_to_cdf(
                sami_data_path=args.sami_dir,
                out_dir=output_dir,
                dtime_sim_start=args.dtime_sim_start,
                progress_bar=args.progress,
                OVERWRITE=args.replace,
                low_mem=args.low_mem,
                delete_raw=args.delete_bins,
                dtime_storm_start=args.dtime_event_start,)

        if regrid:
            print('attempting to regrid!')
            calc_weight = True
            if os.path.exists(os.path.join(args.sami_dir, 'weight.nc')):
                calc_weight = False
            if os.path.exists(os.path.join(output_dir, 'weight.nc')):
                calc_weight = False

            if calc_weight:
                SAMI.make_input_to_esmf(sami_data_path=args.sami_dir,)
                print('and return here when done.')
                return

        else:
            raise NotImplementedError("Haven't gotten here yet.")

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-gitm', '--gitm_dir', type=str, default='./gitm_data',
                        help='GITM Directory. Defaults to ./gitm_data')
    parser.add_argument('-sami', '--sami_dir', type=str, default='./dami_dir',
                        help='SAMI directory. Defaults to ./sami_data')
    parser.add_argument('--dont_regrid', action='store_false',
                        help='Regrid SAMI data? (Default: True)')
    parser.add_argument('--low_mem', type=bool, default = True,
                        help='Process SAMI files in low memory mode?'
                        'Default: True')
    parser.add_argument('-out', '--output_dir', type=str, default=None,
                        help='If you want to save the files to another'
                        'directory, specify it here.')
    parser.add_argument('-r', '--replace', action='store_true',
                        help='Replace existing netCDF files? (Default: False)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print out more information? (Default: False)')
    parser.add_argument('-d', '--delete_bins', action='store_true',
                        help='Delete binary files after processing?'
                        'Caution: This is irreversible! (Default: False)')
    parser.add_argument('-g', '--ghost_cells', action='store_false',
                        help='Drop Ghost Cells? (Default: True)')
    parser.add_argument('-p', '--progress', action='store_true',
                        help='Show progress bar? (Default: False)')
    parser.add_argument('--dtime_sim_start', type=str, default=None,
                        help='Start time of the simulation, in the format:'
                        'YYYYMMDDHHmmSS Required to process SAMI from *.dat')
    parser.add_argument('--dtime_event_start', type=str, default=None,
                        help='Event (storm) start time, in the format:'
                        'YYYYMMDDHHmmSS (Optional. added as attr to netCDFs)')

    args = parser.parse_args()
    main(args)
