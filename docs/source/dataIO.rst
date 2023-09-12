.. _rw_data:
Reading Data
============

While I prefer to use NetCDF files for storing data, there is support to read & interface with the raw model outputs as well. I will assume you are trying to do this from a python script. This process is very similar (though *slightly* different for SAMI3 and GITM).

Both these methods will assume you are running a script from either the root directory of this repository, or have followed the instructions in  :ref:`postinstall`

GITM
----

Two options exist (and are supported) for reading GITM raw GITM data.

Numpy Arrays
^^^^^^^^^^^^

We can read any number (and format) of GITM outputs into a dictionary of numpy arrays. This is not especially fast, memory efficient, or user-friendly, but hey, better than nothing! The code used here has been adapted from the GITM reads in the (develop branch of the) Aetherpy_ repository.

.. _Aetherpy: https://github.com/AetherModel/aetherpy


By just specifying a directory, the ``utility_programs.read_routines.GITM.read_bin_to_nparrays`` module will read all the GITM outputs and return a dictionary with the data. They keys of the dictionary are ``gitmdtimes``, corresponding to the time of midel output; ``gitmbins``, corresponding to the actual data, ``gitmgrid``, which contains the grid information, and optionally ``gitmvars``, which contains the names of the columns of data output.


For example, to read the ``Rho`` outputs from time-steps 250 - 650, and plot a keogram at altitude 250 km and 120 degrees longitude, we can do the following:

.. code-block:: python
    
    from utility_programs.read_routines.GITM import read_bin_to_nparrays
    import matplotlib.pyplot as plt
    import numpy as np

    gitmdata = read_bin_to_nparrays('path/to/GITM/output',
                                    gitm_file_pattern='3DALL*.bin',
                                    cols=['all'],
                                    start_idx=250, end_idx=650)

    lonidx = np.argmin(np.abs(gitmdata['gitmgrid']['longitude'] - 120))
    altidx = np.argmin(np.abs(gitmdata['gitmgrid']['altitude'] - 250))

    plt.imshow(gitmdata['gitmbins'][:, lonidx, :, altidx].T, aspect='auto',
                extent=[gitmdata['gitmdtimes'][0], 
                        gitmdata['gitmdtimes'][-1],
                        gitmdata['gitmgrid']['latitude'][0],
                        gitmdata['gitmgrid']['latitude'][-1]])

    plt.title('GITM Rho at 250 km and 120 degrees longitude')

    plt.show()

This is a little clunky and highlights why we prefer NetCDF files. The numpy arrays and dictionaries are good for starting out, but do not scale well to large datasets or more complicated analysis.


NetCDF Files
^^^^^^^^^^^^

In the backend, this code-base uses Xarray_ for all handling of NetCDF files. The reading from binary to Xarray DataSet is very similar to the methodology for reading to numpy arrays. There are separate scripts to read one & multiple files. To accomplish the same as above, we just need to run:

.. _Xarray: https://docs.xarray.dev/en/stable/user-guide/data-structures.html

.. code-block:: python

    from utility_programs.read_routines import GITM 

    # For a single file (one time)
    gitmdata = GITM.read_bin_to_xarray('path/to/GITM/output/file.bin', 
                                        cols='Rho')

    # And for multiple files:
    gitm_files = GITM.read_multiple_bins_to_xarray('path/to/GITM/output',
                                                    cols='Rho',
                                                    start_idx=250, end_idx=650)

    gitm_files.Rho.sel(lon=120, alt=250, method='nearest').plot(x='time')


Approximately the same amount of work to read the files, but SO much easier to plot & analyze! 

There is also a script to automatically read all of the GITM data in a directory, with preference to reading in to xarray from NetCDF, though it can also read to numpy arrays, or to xarray from binaries. This is the recommended method for reading GITM data.

.. automodule:: utility_programs.read_routines.GITM
    :members: auto_read
    :noindex:


SAMI3
------

This is very similar to reading GITM data, but the format of the non-xarray reads is very different. SAMI3 runs on a magnetic grid, so the data has to be read in a little differently. We cannot make a plot at a single altitude, for example, because the altitude is not constant. Instead, we can make a plot at a single magnetic longitude, or at a range of altitudes.

To read in SAMI3 data to numpy arrays, we need to specify both the path to the data, as well as the start time of the simulation. All other parameters are automatically read from the settings files.

.. code-block:: python

    from utility_programs.read_routines import SAMI
    from datetime import datetime

    sami_data, times = SAMI.read_to_nparray('/path/to/sami/data', 
                    dtime_sim_start=datetime(2011,5,16),
                    cols='edens')

This ``sami_data`` is a python dictionary with keys ['grid', 'data'], where ``sami_data['data']`` is indexed with [varname][nt, nlt, nf, nz].

It's complicated, so there's also a script to read the data to xarray, which is much easier to use.

.. code-block:: python

    from utility_programs.read_routines import SAMI
    from datetime import datetime

    sami_ds = read_raw_to_xarray('/path/to/sami/data',
                     dtime_sim_start=datetime(2011,5,16),
                     cols='edens')


And this sami_ds is slightly easier to deal with! It is indexed with [nt, nlt, nf, nz], and has the same variables as the numpy array read, but also has the grid information as coordinates.

Of course, we also have and auto_read module:

.. automodule:: utility_programs.read_routines.SAMI
    :members: auto_read
    :noindex:
