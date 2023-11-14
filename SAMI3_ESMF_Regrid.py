"""
This module contains functions for processing SAMI raw data for use in ESMF.
Input and output grids are calculated and written to file as a 3D UGRID mesh
and the ESMF weight file is computed and then applied.

All functions use spherical coordinates with `lon` & `lat` in degrees (with
`lon` from 0-360), and `alt` in km from Earth's surface.

"""

import os
import subprocess
from tqdm import tqdm
import xarray as xr
import pandas as pd
import numpy as np
from scipy.io import netcdf_file
from utility_programs.read_routines import SAMI
from utility_programs.utils import str_to_ut


def generate_interior_points_sami_raw(old_shape, progress=False):
    """Generates mesh points from SAMI raw data for use in ESMF.

    Args:
        in_cart (numpy.ndarray): 3xN array of SAMI points (any coord system)
        old_shape (list): (nlt, nf, nz) shape of original sami outputs.
        proogress (bool): Whether or not to show `tqdm` progress bar

    Returns:
        list: list of 1D indices of the corners of the cuboids
            to be used as the grid corners in ESMF. These are used for
            the "vertids" variable in the UGRID mesh ESMF input.


    Notes:
        - This is not well documented. Contact me with questions
        - It should "just work", but may not...
        - The mesh generated does follow ESMF conventions, however
            more points are thrown out than is probably necessary.
            With the size of the SAMI grid this is probably not an issue.
    """

    nlt, nf, nz = old_shape

    idxs = []
    if progress:
        pbar = tqdm(total=np.prod(old_shape),
                    desc='Generating interior points from SAMI outputs')
    badbadbad = []  # used for debugging, not returned.

    for lt in range(nlt):
        for f in range(nf):
            for z in range(nz):

                if lt == old_shape[0] - 1:
                    l2 = 0
                else:
                    l2 = lt + 1

                if z == 0:
                    z2 = -1
                else:
                    z2 = z - 1

                f2 = f + 1
                if f == old_shape[1] - 1:
                    badbadbad.append([lt, f, z])
                    continue

                cs = [[lt, f, z],
                      [l2, f, z],
                      [l2, f, z2],
                      [lt, f, z2],
                      [lt, f2, z],
                      [l2, f2, z],
                      [l2, f2, z2],
                      [lt, f2, z2]]

                id_pt = []

                for c in cs:
                    try:
                        index = np.ravel_multi_index(c, old_shape)
                    except ValueError:
                        break
                    id_pt.append(index)

                idxs.append(id_pt)

                if progress:
                    pbar.update()

    if progress:
        pbar.close()
    return idxs


def generate_interior_points_output_grid(longitudes, latitudes, altitudes,
                                         progress=False):
    """ Generate a 3D mesh from arrays of desired output coordinate values.
    The input points are 1D arrays of various lengths, denoting the
    desired output grid (in deg and km above Earth's surface).
    The output is a 1D array of the corner points of the cuboids to be
    used as the grid corners in ESMF. These are used for
    the "vertids" variable in the UGRID mesh ESMF input.

    Args:
        longitudes (numpy.ndarray or list): 1D array of longitudes
        latitudes (numpy.ndarray or list): 1D array of latitudes
        altitudes (numpy.ndarray or list): 1D array of altitudes
        progress (bool): Whether or not to show `tqdm` progress bar

    Returns:
        tuple: tuple of 1D arrays of the corner points of the cuboids
            to be used as the grid corners in ESMF. These are used for
            the "vertids" variable in the UGRID mesh ESMF input.
            [lons, lats, alts, connections]

    """
    lons = []
    lats = []
    alts = []

    node_conns = []
    outcoords = []  # for debugging - not returned.

    lon_dim = len(longitudes)
    lat_dim = len(latitudes)
    alt_dim = len(altitudes)

    out_shape = [lon_dim, lat_dim, alt_dim]

    if progress:
        pbar = tqdm(total=np.prod([len(longitudes),
                                   len(latitudes),
                                   len(altitudes)]),
                    desc='Generating interior points for output grid')

    for a, lon in enumerate(longitudes):
        for b, lat in enumerate(latitudes):
            for c, alt in enumerate(altitudes):
                lons.append(lon)
                lats.append(lat)
                alts.append(alt)

                a1 = a+1 if a < lon_dim else 0
                # if we're near the poles or alt lims, throw points out:
                b1 = b+1 if b < lat_dim else 'a'
                c1 = c+1 if c < alt_dim else 'a'

                cs = [[a, b, c],
                      [a1, b, c],
                      [a1, b1, c],
                      [a, b1, c],
                      [a, b, c1],
                      [a1, b, c1],
                      [a1, b1, c1],
                      [a, b1, c1]]

                id_pt = []
                yes = True
                for c in cs:
                    try:
                        index = np.ravel_multi_index(c, out_shape)
                    except ValueError:
                        yes = False
                        break
                    id_pt.append(index)

                if yes:
                    outcoords.append(cs)

                    node_conns.append(id_pt)

                if progress:
                    pbar.update()
    if progress:
        pbar.close()

    return lons, lats, alts, node_conns


def generate_interior_points_custom_grid(lats, lons, alts,
                                         cell_radius=0.5,
                                         progress=False):
    """Generate a 3D mesh from given 3D coordinates and write it to a UGRID
    file for use in ESMF. The mesh has overlap (degeneracy) which *should*
    be fine. The lats, lons, and alts arrays should contain the same number of
    points.

    Args:
        lats (numpy.ndarray): 1D array of latitudes
        lons (numpy.ndarray): 1D array of longitudes
        alts (numpy.ndarray): 1D array of altitudes
        cell_radius (float): Distance from point to corner of cell,
            in degrees, (altitude is 10* cell_radius). Default is 1 degree.
        progress (bool): Whether or not to show `tqdm` progress bar

    Returns:
        tuple: tuple of 1D arrays of the corner points of the cuboids
            to be used as the grid corners in ESMF. These are used for
            the "vertids" variable in the UGRID mesh ESMF input.
            [lons, lats, alts, connections]
    """

    assert len(lats) == len(lons) == len(alts), \
        'lats, lons, and alts must be the same length!'

    lat_corners = []
    lon_corners = []
    alt_corners = []

    node_conns = []

    for centerpt in range(len(lats)):
        lon_below = lons[centerpt] - cell_radius
        lon_above = lons[centerpt] + cell_radius

        lat_below = lats[centerpt] - cell_radius
        lat_above = lats[centerpt] + cell_radius

        alt_below = alts[centerpt] - 10*cell_radius
        alt_above = alts[centerpt] + 10*cell_radius

        # Make sure we don't go over the poles:
        if lat_below < -90:
            lat_below = -90
        if lat_above > 90:
            lat_above = 90

        lat_corners.append([lat_below, lat_above,
                            lat_above, lat_below,
                            lat_below, lat_above,
                            lat_above, lat_below])
        lon_corners.append([lon_below, lon_below,
                            lon_above, lon_above,
                            lon_below, lon_below,
                            lon_above, lon_above])
        alt_corners.append([alt_below, alt_below,
                            alt_below, alt_below,
                            alt_above, alt_above,
                            alt_above, alt_above])

        node_conns.append([8*centerpt, 8*centerpt+1,
                           8*centerpt+2, 8*centerpt+3,
                           8*centerpt+4, 8*centerpt+5,
                           8*centerpt+6, 8*centerpt+7])

    return np.array(lon_corners).flatten(), \
        np.array(lat_corners).flatten(), \
        np.array(alt_corners).flatten(), \
        np.array(node_conns)


def write_UGRID_mesh(lon, lat, alt, indices, fname):
    """ Write a UGRID mesh file for given coordinates and
    indices (interconnects between nodes). See the ESMF documentation
    (UGRID section) for more info.

    Args:
        lon (numpy.ndarray): 1D array of longitudes
        lat (numpy.ndarray): 1D array of latitudes
        alt (numpy.ndarray): 1D array of altitudes
        indices (numpy.ndarray): 2D array of indices, shape (N, 8)
        fname (str): Output filename

    Returns:
        None

    Notes:
        - The mesh must be created elsewhere, this just writes the file.
        - This function resides in a file with other functions to help
            generate the mesh interconnects.
        - The provided indices must be in the correct order, see the
            ESMF documentation for more info.
        - This function is used for both the input and output meshes,
            and can be used for anyything else, provided the arguments are
            correct.

    """

    f = netcdf_file(fname, "w")  # make file
    num_cells = np.sum([len(i) == 8 for i in indices])

    # add dimensions
    f.createDimension("eight", 8)
    f.createDimension("nnodes", len(lon))
    f.createDimension("ncells", num_cells)

    # add variables
    # mesh, empty variable to describe topology:
    mesh = f.createVariable("mesh", "i", {})
    mesh.cf_role = "mesh_topology"
    mesh.topology_dimension = 3
    mesh.node_coordinates = "nodelon nodelat height"
    mesh.volume_node_connectivity = "vertids"
    mesh.volume_shape_type = "meshtype"

    nodelon = f.createVariable("nodelon", float, ("nnodes",))
    nodelon.standard_name = "longitude"
    nodelon.units = "degrees_east"

    nodelat = f.createVariable("nodelat", float, ("nnodes",))
    nodelat.standard_name = "latitude"
    nodelat.units = "degrees_north"

    height = f.createVariable("height", float, ("nnodes",))
    height.standard_name = "altitude"
    height.units = "kilometers"

    vertids = f.createVariable("vertids", "i", ("ncells", "eight"))
    vertids.cf_role = "volume_node_connectivity"
    vertids.start_index = 0

    meshtype = f.createVariable("meshtype", "i", ("ncells",))
    meshtype.cf_role = "volume_shape_type"
    meshtype.flag_range = 1.0
    meshtype.flag_values = 1.0
    meshtype.flag_meanings = "hexahedron"

    # and add location data...
    nodelon[:] = lon
    nodelat[:] = lat
    height[:] = alt

    # this is the most important part - the connection between nodes.
    vertids[:] = np.array(
        [*np.array(indices, dtype="object")
         [np.where([len(i) == 8 for i in indices])[0]]]
    )

    meshtype[:] = [1 for i in range(num_cells)]

    f.close()

    return


def apply_weight_file(sami_data_path,
                      dtime_sim_start,
                      out_dir,
                      cols='all',
                      progress=True,
                      output_filename=None
                      ):
    """Apply the ESMF weight file to the SAMI raw data.

    Args:
        sami_data_path (str): Path to SAMI3 raw data
        dtime_sim_start (str): Simulation start date (YYYYMMDD)
        out_dir (str): Output directory, Default is same as sami_data_path
        cols (str or list): Columns to interpolate. Default is 'all'
        progress (bool): Whether or not to show `tqdm` progress bar
        output_filename (str): Output filename, Default is to make a new file
            for each variable.

    Returns:
        None

    Notes:
        - This has been sped up and takes ~1 hour for a full SAMI run
        - Multithreading this should be straighforward, but it is not a
            priority since the speed is good enough (and Xarray can lock
            files which will cause issues).
        - Memory usage is ~10GB, so it is likely too heavy to run on a
            login node.

    """

    #  Start with loading weight file and the output grid file...
    weights = xr.open_dataset(os.path.join(
        sami_data_path, 'esmf_weightfile.nc'))
    dst_ds = xr.open_dataset(os.path.join(sami_data_path, 'dst_ugrid.nc'))

    # Get the output dimensions outside the loop the loop:
    outlon = np.unique(dst_ds.nodelon.values)
    outlat = np.unique(dst_ds.nodelat.values)
    outalt = np.unique(dst_ds.height.values)
    # newshape = [len(outlon), len(outlat), len(outalt)]

    if cols != 'all':
        if type(cols) is str:
            cols = [cols]
    else:
        cols = SAMI.sami_og_vars.values()

    # Put weight file into np arrays for speed:
    new_s = weights.S.values
    new_col = weights.col.values - 1
    new_row = weights.row.values - 1

    first = True
    for nvar, datavar in enumerate(cols):

        sds = SAMI.read_raw_to_xarray(sami_data_path,
                                      dtime_sim_start,
                                      cols=datavar,
                                      progress_bar=False)
        # Progress bar:
        if progress:
            if first:
                pbar = tqdm(total=len(sds.time) * len(cols))

            pbar.set_description('interpolating %s \t vars:(%i/%i)'
                                 % (datavar, nvar, len(cols)))

        # Xarray dataset for output data:
        out_ds = xr.Dataset()
        out_ds['lon'] = outlon
        out_ds['lat'] = outlat
        out_ds['alt'] = outalt
        out_ds['time'] = sds.time

        # loop thru time, esmf does not hold any of our time info:
        dstpts = np.zeros((sds.time.size, weights.n_b.size))
        unique_rows, row_indices, row_counts = np.unique(
            new_row[weights.n_s.values],
            return_inverse=True,
            return_counts=True)

        for t in range(sds.time.size):
            srcData = sds[datavar].isel(time=t).values.flatten()
            # Accumulate values for each unique row index using np.add.at
            np.add.at(dstpts[t],
                      unique_rows,
                      np.bincount(row_indices,
                                  weights=new_s[weights.n_s.values]
                                  * srcData[new_col[weights.n_s.values]],
                                  minlength=len(unique_rows)))

            if progress:
                pbar.update()

        if progress:
            pbar.set_description('Writing %s \t vars:(%i/%i)'
                                 % (datavar, nvar, len(cols)))

        datavar = datavar.replace('+', '_plus_').replace('-', '_')
        out_ds[datavar] = (('time', 'lon', 'lat', 'alt'),
                           dstpts.reshape([sds.time.size,
                                           len(outlon),
                                          len(outlat),
                                           len(outalt)]))

        if output_filename:
            out_ds.to_netcdf(os.path.join(out_dir,
                                          output_filename+'_SAMI_REGRID.nc'),
                             mode='a' if os.path.exists(
                                 os.path.join(sami_data_path,
                                              output_filename
                                              + '_SAMI_REGRID.nc')) else 'w',
                             engine='h5netcdf')
        else:
            out_ds.to_netcdf(os.path.join(out_dir, datavar+'_SAMI_REGRID.nc'),
                             engine='h5netcdf')

        del out_ds

        first = False


def main(sami_data_path,
         dtime_sim_start,
         num_lons=90,
         num_lats=180,
         num_alts=100,
         alt_step=None,  # use either this or num_alts!
         min_alt=100,
         max_alt=2400,
         custom_input_file=None,
         # User-defined input file (.csv with comma sep and a header)
         # (i.e. sat track w/ `glon, glat, alt` columns.)
         custom_grid_size=0.5,
         cols='all',
         progress=False,
         remake_files=False,
         out_dir=None,
         # if outname is none, output new files for each var.
         output_filename=None):
    """ Main function for processing SAMI raw data for use in ESMF.

    Args:
        sami_data_path (str): Path to SAMI3 raw data
        dtime_sim_start (str): Simulation start date (YYYYMMDD)
        num_lons (int): Number of longitudes in output grid
        num_lats (int): Number of latitudes in output grid
        num_alts (int): Number of altitudes in output grid
        alt_step (int): Altitude step size in output grid. Use this or num_alts
            (Default is None)
        min_alt (int): Minimum altitude in output grid
        max_alt (int): Maximum altitude in output grid
        custom_input_file (str): User-defined input file (.csv with comma sep
            and a header) (i.e. sat track w/ `glon, glat, alt` columns.)
        custom_grid_size (float): Size of grid cells to use when regridding to
            a satellite file. Measured in degrees from center (so it is a
            radius), and altitude is 10*custom_grid_size. Default is 0.5.
        cols (str or list): Columns to interpolate (comma sep)
        progress (bool): Whether or not to show `tqdm` progress bar
        remake_files (bool): Remake ESMF input files? Default is False
            (reuse existing files)
        out_dir (str): Output directory, Default is same as sami_data_path
        output_filename (str): Output filename, Default is to make a new file
            for each variable.

    Returns:
        None

    Notes:
        - This is not *super* well documented. Contact me with questions
        - It should "just work", but may not...
        - The mesh generated does follow ESMF conventions, however
            more points are thrown out than is probably necessary.
            With the size of the SAMI grid this is probably not an issue.

    """

    if type(dtime_sim_start) is str:
        dtime_sim_start = str_to_ut(dtime_sim_start)

    # read in the SAMI grid:
    nz, nf, nlt, nt = SAMI.get_grid_elems_from_parammod(sami_data_path)
    sami_grid = SAMI.get_sami_grid(sami_data_path, nlt, nf, nz)

    # Make the SAMI interior points...

    # Will need these later:
    raw_glon = sami_grid['glon'].flatten()
    raw_glat = sami_grid['glat'].flatten()
    raw_alt = sami_grid['alt'].flatten()

    make_esmf_inputs = False
    if os.path.exists(os.path.join(sami_data_path, 'src_ugrid.nc')) and\
            os.path.exists(os.path.join(sami_data_path, 'dst_ugrid.nc')):
        if remake_files:
            make_esmf_inputs = True
        else:
            print('ESMF inputs already found. Will interpolate',
                  'using the same input/output files.\n',
                  'Run with remake_files=True to make the input/output files',
                  'again.')
    else:
        make_esmf_inputs = True

    # only redo if files don't exist (and user doesn't want them remade...)
    if make_esmf_inputs:
        # Make connection indices:
        sami_idxs = generate_interior_points_sami_raw(
            np.array([raw_glat, raw_glon, raw_alt]),
            sami_grid['mlat'].shape,  # Shape of the SAMI grid (nlt, nf, nz)
            progress=progress)

        # write this to a NetCDF file:
        write_UGRID_mesh(raw_glon, raw_glat, raw_alt, sami_idxs,
                         os.path.join(sami_data_path, 'src_ugrid.nc'))

    if make_esmf_inputs:
        if custom_input_file:

            custom_grid_df = pd.read_csv(custom_input_file)

            # Now make the outputs:
            out_lat = custom_grid_df['glat'].values
            out_lon = custom_grid_df['glon'].values
            out_alt = custom_grid_df['alt'].values

            # And generate interior points
            flat_lon, flat_lat, flat_alt, output_idxs = \
                generate_interior_points_custom_grid(
                    out_lon, out_lat, out_alt,
                    cell_radius=custom_grid_size,
                    progress=progress)

            write_UGRID_mesh(flat_lon, flat_lat, flat_alt, output_idxs,
                             os.path.join(sami_data_path, 'dst_ugrid.nc'))

        else:
            # Now make the outputs:
            out_lat = np.linspace(-90, 90, num_lats)
            out_lon = np.linspace(0, 360, num_lons, endpoint=False)

            if alt_step:
                out_alt = np.arange(min_alt, max_alt+1, alt_step)
            else:
                out_alt = np.linspace(min_alt, max_alt, num_alts)

            # And generate interior points
            flat_lon, flat_lat, flat_alt, output_idxs = \
                generate_interior_points_output_grid(out_lon, out_lat, out_alt,
                                                     progress=progress)

            # write this to a NetCDF file, just like before:
            write_UGRID_mesh(flat_lon, flat_lat, flat_alt, output_idxs,
                             os.path.join(sami_data_path, 'dst_ugrid.nc'))

    make_esmf_weights = False
    if os.path.exists(os.path.join(sami_data_path, 'esmf_weightfile.nc')):
        if remake_files:
            print('Found ESMF weight file, making it again...')
            make_esmf_weights = True
        else:
            print('Reusing existing ESMF weight file.')
    else:
        make_esmf_weights = True

    if make_esmf_weights:
        # And now we can call ESMF!
        print('calling ESMF...')
        esmf_command = ['ESMF_RegridWeightGen -s',
                        os.path.join(sami_data_path, 'src_ugrid.nc'),
                        '-d', os.path.join(sami_data_path, 'dst_ugrid.nc'),
                        '--src_loc corner --dst_loc',
                        'center' if custom_input_file else 'corner',
                        '-l greatcircle -i -w',
                        os.path.join(sami_data_path, 'esmf_weightfile.nc'),
                        ]

        esmf_result = subprocess.run(
            ' '.join(esmf_command), shell=True, check=True)

        # Check the result
        if esmf_result.returncode == 0:
            print("Command was successful")
        else:
            print(f"Command failed with return code {esmf_result.returncode}")
            print('\n\n\nError in using ESMF. Make sure it is loaded!\n',
                  'If you do have ESMF loaded, try to run:\n\n\t',
                  ' '.join(esmf_command),
                  '\n\nCheck output logs in PET*.Log for more info ',
                  '(if the file does not exist, ESMF was unable to run).\n',
                  'After ESMF is finished, come back and apply the weights.\n')
            print("Error output:", esmf_result.stderr)
            raise ValueError('ESMF did not work!')

    # Now we can apply weights!
    apply_weight_file(sami_data_path,
                      dtime_sim_start,
                      out_dir,
                      cols=cols,
                      progress=progress,
                      output_filename=output_filename)

    return


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Process SAMI3 raw data for use in ESMF.'
        ' The output grid can be user specified.')

    parser.add_argument('sami_data_path', type=str,
                        help='Path to SAMI3 raw data')
    parser.add_argument('dtime_sim_start', type=str,
                        help='Simulation start date (YYYYMMDD)')
    parser.add_argument('--no_pbar', action='store_true',
                        help='Disable progress bar')
    parser.add_argument('--cols', type=str, default='all',
                        help='Columns to interpolate (comma sep)')
    parser.add_argument('--out_dir', type=str, default=None,
                        help='Output directory, '
                        'Default is same as sami_data_path')
    parser.add_argument('--output_filename', type=str, default=None,
                        help='Output filename, '
                        'Default is to make a new file for each variable.')
    parser.add_argument('--remake_files', action='store_true',
                        help='Remake ESMF input files? '
                        'Default is False (reuse existing files)',
                        default=False)
    parser.add_argument('--num_lons', type=int, default=90,
                        help='Number of longitudes in output grid')
    parser.add_argument('--num_lats', type=int, default=180,
                        help='Number of latitudes in output grid')
    parser.add_argument('--num_alts', type=int, default=100,
                        help='Number of altitudes in output grid')
    parser.add_argument('--alt_step', type=int, default=None,
                        help='Altitude step size in output grid. '
                        'Use this or num_alts! (Default is None)')
    parser.add_argument('--min_alt', type=int, default=100,
                        help='Minimum altitude in output grid')
    parser.add_argument('--max_alt', type=int, default=2400,
                        help='Maximum altitude in output grid')
    parser.add_argument('--custom_input_file', type=str, default=None,
                        help='User-defined input file (.csv with comma sep '
                        'and a header) (i.e. sat ephemeris file w/ '
                        '`glon, glat, alt` columns.)')
    parser.add_argument('--custom_grid_size', type=float, default=1,
                        help='Size of grid cells to use when regridding to a '
                        'satellite file. Measured in degrees from center (so '
                        'it is a radius), and altitude is '
                        '10*custom_grid_size. Default is 0.5.')

    args = parser.parse_args()

    if args.out_dir is None:
        args.out_dir = args.sami_data_path

    main(args.sami_data_path,
         args.dtime_sim_start,
         args.num_lons,
         args.num_lats,
         args.num_alts,
         args.alt_step,
         args.min_alt,
         args.max_alt,
         args.custom_input_file,
         args.custom_grid_size,
         args.cols,
         not args.no_pbar,
         args.remake_files,
         args.out_dir,
         args.output_filename
         )
