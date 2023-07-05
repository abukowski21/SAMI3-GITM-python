"""
Code to interpolate any model data to either a:
- user-defined grid
- collection of satellite points


It's easiest to interface with the  `interpolate` function,
    or use the `main` function from the command line...

Can read output coords from cdf/csv or by user-defined grid.
> Default output grid is 5deg lon x 2deg lat x 50km alt

"""


import xarray as xr

from tqdm.auto import tqdm
import numpy as np

import os
from utility_programs.read_routines import SAMI
from utility_programs.utils import str_to_ut, make_ccmc_name
import argparse
import pickle
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator



def latlonalt_to_cart(lat, lon, radius):
    """Convert lat, lon, alt to cartesian coordinates.

    Args:
        lat (numpy.ndarray): latitude (in degrees)
        lon (numpy.ndarray): longitude (in degrees)
        radius (numpy.ndarray): radius in Re from center of Earth

    Returns:
        numpy.ndarray: 3xN array of cartesian coordinates
    """
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)
    x = radius * np.cos(lat) * np.cos(lon)
    y = radius * np.cos(lat) * np.sin(lon)
    z = radius * np.sin(lat)
    return np.array([x, y, z])


    
def do_interpolations(
    sami_data_path=None,
    dtime_sim_start=None,
    out_latlonalt=None,
    out_path=None,
    out_runname='',
    save_delauney=False,
    max_alt=None,
    cols='all',
    show_progress=False,
    engine='h5netcdf',
    ):
    

    # deal with sami first
    if sami_data_path is not None:
        
        if type(dtime_sim_start) == str:
            dtime_sim_start = str_to_ut(dtime_sim_start)
        
        nz, nf, nlt, nt=SAMI.get_grid_elems_from_parammod(sami_data_path)
        old_shape=[nlt, nf, nz]
        grid=SAMI.get_sami_grid(sami_data_path, nlt, nf, nz)
        
        # specify max alt to build delauney at
        if max_alt is None:
            if out_latlonalt is not None:
                #alt will be the biggest coord:
                max_alt = np.max(out_lat_lon_alt) + 300
                # add 300 to make sure we have enough points above top
            else:
                max_alt = 2500
        
        mask = np.where(grid['alt'] < max_alt)
        grid2={}
        for k in grid.keys():
            grid2[k]=grid[k][mask].flatten()
        del grid
        
        in_cart=latlonalt_to_cart(grid2['glat'],
                                  grid2['glon'],
                                  grid2['malt']).T
        
        
        if out_latlonalt is None:
            latout=np.arange(-90, 90, 2)
            lonout=np.arange(0, 360, 5)
            altout=np.arange(200, 2200, 50)
            out_lats=[]
            out_lons=[]
            out_alts=[]

            for a in latout:
                for o in lonout:
                    for l1 in altout:
                        out_lats.append(a)
                        out_lons.append(o)
                        out_alts.append(l1)

            out_latlonalt=latlonalt_to_cart(
                out_lats, out_lons, np.array(out_alts) + 6371)
        else:
            out_latlonalt=latlonalt_to_cart(
                out_latlonalt[0], out_latlonalt[1],
                np.array(out_latlonalt[2]) + 6371)
        
        if os.path.exists(os.path.join(sami_data_path, 
                                       'delauney_max-%i.pkl' %max_alt)):
            if save_delauney:
                print('attempting to reuse existing triangulation file')
                with open(os.path.join(
                    sami_data_path, 'delauney_max-%i.pkl' %max_alt), 'rb') as f:
                    tri = pickle.load(f)
            else:
                print('Found existing triangulation file. Recalculating...',
                     '\n(Specify save_delauney=True to reuse)')
                tri = Delaunay(in_cart)
        else:
            print('Calculating Delauney Triangulation..')
            tri = Delaunay(in_cart)
            if save_delauney:
                print('Saving')
                with open(os.path.join(sami_data_path,
                                       'delauney_max-%i.pkl' %max_alt),
                          'wb') as f:
                    pickle.dump(tri, f)
        
        #format 'cols' variable
        if cols == 'all':
            cols = SAMI.sami_og_vars.items()
        else:
            if type(cols) == str:
                cols = [cols]
            else:
                cols = np.asarray(cols)
                
        if show_progress:
            pbar = tqdm(total=len(cols)*nt,
                       desc='Reading in SAMI data')
            
        first=True #for choosing which mode to write
        for data_var in cols:
            interpd = []
            data, times=SAMI.read_to_nparray(
                sami_data_path, dtime_sim_start,
                cols=data_var,
                skip_time_check=True)
            
            if show_progress:
                pbar.set_description('interpolating')
            for t in tqdm(range(len(times))):
                interp = LinearNDInterpolator(
                    tri,
                    data['data'][data_var][:,:,:,t][mask].flatten())
                interpd.append(interp(out_latlonalt.T))
                if show_progress:
                    pbar.update()
            if show_progress:
                pbar.set_description('writing Dataset...')
            ds=xr.Dataset(coords={
                'time': (['time'], times),
                'alt': (['alt'], altout),
                'lat': (['lat'], latout),
                'lon': (['lon'], lonout)},)
            ds[data_var]=(('time', 'lat', 'lon', 'alt'),
                         np.array(interpd).reshape(
                             len(times),
                             len(latout),
                             len(lonout),
                             len(altout)))
            ds.to_netcdf(os.path.join(
                out_path,'SAMI_REGRID'+out_runname),
                         engine=engine,
                         mode='w' if first else 'a',
                         encoding={'time':{'dtype':float}})
            first=False
            del ds, interpd, data #clean up memory

    
# if __name__ == '__main__':