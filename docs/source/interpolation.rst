.. _Interpolation:

Interpolation
============

After *so* much experimentation, limitations were found with the ``LinearNDInterpolator`` in SciPy. Mainly, it is slow (even when re-using a Delauney Triangulation) and the outputs are heavily influenced by the distance from source to destination points. After trying to write several interpolators, I have decided to implement the Earth System Modeling Framework (ESMF) to interpolate SAMI3 results. Modules to interpolate GITM are currently half-baked but will be finalized soon.

Overview
*********

ESMF is a super powerful framework designed to link various models together. Of the entire functionality, we are only using a tiny bit of it. Of the various interpolation routines I have tried, ESMF has been by far the fastest and most accurate. 

ESMPy is extremely hard to debug, at least for me. For this reason, we are calling ESMF from the command line within Python code. Files are autmomatically written with the source and destination "grids", and then the weights are automatically applied. ESMF must be installed with NetCDF support, so the ``conda`` or ``pip`` versions will not work. See the :ref:`Installing ESMF` page for more details.

Some more details on why ``LinearNDInterpolator`` is not used, as well as the details on how ESMF is used are available on the indepth-interpolation page. #TODO: (add link and finish that, haha). This is a placeholder right now...

Usage
*****

Two options on how to use the ESMF functionality currently. From the SAMI3 outputs you can interpolate to a grid or to a collection of user-specified points.

The entrypoint into the ESMF functionality is the ``SAMI3_ESMF_Regrid.py`` program. Run it with ``-h`` to see helpful information.

Interpolating to a Grid
=======================

If you are trying to work with global data, this is the choice for you. You can tell the program your required output grid spacing and the rest will be done for you. Here is an example:

.. code-block:: bash
	python SAMI3_ESMF_Regrid.py /path/to/sami/outputs [simulation start date as a YYYYMMDD str] --cols edens --out_dir /path/to/output/location --num_lons 90 --num_lats 90 --num_alts 100 --min_alt 100 --max_alt 2500

Interpolating to User-Specified Points
======================================

Often, one wishes to simulate satellite measurements, or not interpolate an entire grid. If that sounds like something you want to do, this is the mode for you!

As ESMF is tricky, the requires a specific format for the locations to be output. If you are having errors, you probably did not format this file correctly. 
The code will look for columns titled **glon**, **glat**, and **alt**. If these are not found, the program will look for columns titled **lon**, **lat**, and **alt**. The file must be a ``.csv`` file with these column names. Time does not matter, as the grid will be interpolated to these points at all simulation output times. The units are degrees, degrees, and km above Earth's surface, respectively. 

.. note:: I have not noticed any differences between using longitudes from (0,360) and (-180,180), but SAMI3 gives longitude in (0,360) so if (-180,180) is not working for you, that is why. Let me know and I will update this to reflect your findings.

Here is an example:

.. code-block:: bash
	python SAMI3_ESMF_Regrid.py /path/to/sami/outputs [simulation start date as a YYYYMMDD str] --cols edens --out_dir /path/to/output/location --custom_input_file satfiletmp.csv


Under the hood, a ESMF Mesh object is created for each output location. the ``----custom_grid_size`` variable specifies the size of this mesh, in degrees. If you see a lot of missing points, are working at high altitudes (above ~3,000 km), or notice the interpolated data are not as smooth as you wopuld like, you can try to adjust this. The default should be fine in almost all cases though.


Tips
****

If the weight generation is taking a while, you can use MPI to give the work to multiple processors. This does not affect results. Use the ``--use_mpi`` flag with the number of processors you would like to use if the weight generation is taking a while. Note: MPI is only used in the weight calculation, not applying the weights.

Only run with the variables you want interpolated. Interpolating all variables takes a long time. The weight application is vectorized but it is still slow.

Interpolation to a custom grid makes a 2D NetCDF file. The dimensions are for the simulation output time step and the index of the location. If you only want 1D files, let me know and we can try to make this happen, or you can use some code I've written that may or may not be included #TODO: (make sure it's here and add this functionality.)

ESMF is complicated and I've tried to build this out in an approachable way. Please don't hesitate to fill out a GitHub issue or contact me with questions or issues.

