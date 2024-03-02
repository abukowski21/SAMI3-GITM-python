.. _Interpolation:

Interpolation
#############

After *so* much experimentation, limitations were found with the ``LinearNDInterpolator`` in SciPy. Mainly, it is slow (even when re-using a Delauney Triangulation) and the outputs are heavily influenced by the distance from source to destination points. After trying to write several interpolators, I have decided to implement the Earth System Modeling Framework (ESMF) to interpolate SAMI3 results. Modules to interpolate GITM are currently half-baked but will be finalized soon.

Overview
*********

ESMF is a super powerful framework designed to link various models together. Of the entire functionality, we are only using a tiny bit of it. ESMF has far out-performed the various interpolation routines I have tried, from standard Python routines to exotic machine learning routines. ESMF has been fast, accurate, reliable (when set up right), and reproducible.

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



Known errors
************

This section provides a brief overview of known errors that may occur when running the ESMF interpolation, as well as how to fix them. If you encounter an error that is not listed here, please let me know and I will add it to the list.

A disclaimer first. As best stated `here <https://xesmf.readthedocs.io/en/latest/other_tools.html#other-geospatial-regridding-tools>`_,

    "However, ESMPy is a complicated Python API that controls a huge Fortran beast hidden underneath. It is not as intuitive as native Python packages, and even a simple regridding task requires more than 10 lines of arcane code."

Ideally we would use ESMPy for interpolations since it does not rely on the user installing ESMF, managing the modules with that process, and calling command line functions within Python. But ESMPy is a bit difficult to work with. One day the functionality to use ESMPy will be added and this section can be deleted, but that day will need to wait until my dissertation is done.

The following are known errors that may occur when running the ESMF interpolation, as well as how to fix them.


1. ``Error in system call pthread_mutex_destroy: Device or resource busy
    src/mpi/init/initthread.c``, ``[system_name:mpi_rank_0]``, or ``application called MPI_Abort(comm=...`` errors.
	- ESMF cannot set up the MPI interface.
	- You are likely trying to run MPI programs on a login node. This is bad practice and system administrators have put in place measures to prevent this. You will achieve higher throughput and not take up resources from other users by allocating yourself a compute node (or using a dedicated analysiis node) and running things there.
2. Error code 127 can be from several things. First, ``Error loading shared library: lib[...].so: cannot open shared object file: No such file or directory``
	- This error occurs when the Fortran and C compilers used during the ESMF install are not loaded. ESMF is looking to use library files that it cannot find. 
	- To fix this, you need to load the same modules used during the install before running any ESMF modules. To save time, this module stack can be saved to a collection with ``module save [name]`` and then reloaded with ``module restore [name]``.
	- On a system without ``modules``, you will need to add the libraries used to install ESMF to your ``LD_LIBRARY_PATH`` or ``$PATH``.
	- Either can be placed into your startup scripts (.bashrc/.bashprofile/.zshrc/etc.) to be configured automatically when you log iinto the system. Setting default modules is also a good idea, but deprecated so not advised in case you screw things up.
3. Second, error 127 and: ``forrtl: severe (174): SIGSEGV, segmentation fault occurred``.
	- This error occurs when the MPI modules used during the install are not loaded. This one takes a while to show up so oyu might feel like you got lucky and then it will crash.
	- See previous for how to fix. Just change your MPI modules to the ones used during the install.
4. Third, error code 127 and ``Command [...] not found``.
	- The subprocess call could not find the ESMF_RegridWeightGen executable.
	- To fix this, set the ``ESMF_DIR`` flag (unfortunately named since it's the same name as the variable used during ESMF install) to the path to the ESMF executables. From the $ESMF_DIR used inn the install, go to apps/[...]/ and you will see the executables. Get into the apps directory and hit tab till you find some programs. The directory you found is what ``ESMF_DIR`` should be set to.
	- This error could also be caused by the ESMF executables not being listed correctly in $PATH. If you *did* add them to $PATH and installed something like ESMPy into your Python environment, they could be in the wrong order. Run ``echo $PATH$`` or ``which ESMF_RegridWeightGen`` to see where the executable is being called from. If it is not the same as the one you installed, you need to fix your $PATH (or use the ``ESMF_DIR`` flag).
4. Error code 1. Take a look in the logs (generated in the folder you ran the code from). If ``ESMF_PIO  : disabled`` and ``ESMF_PNETCDF  : disabled``, it could be two things.
	- Easy: Check the output of ``which ESMF_PrintInfo``. If it points to a conda environment, ESMF cannot find the executable you installed. You need to fix your $PATH or use the ``ESMF_DIR`` flag when running the Python code (it does not look at your environment variables on its own).
	- Less easy: If that command is not found, ESMF cannot find the executables. Hopefully you just typed the ESMF_DIR wrong or something. Try running the ``ESMF_PrintInfo`` executable (same directory as the other installed scripts) and see if PIO & pnetcdf are enabled. If not, you need to reinstall ESMF with these enabled. Sorry! 

If these don't work or you find new errors, let me know and I'll update this page. Good luck!

