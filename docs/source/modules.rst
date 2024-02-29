API REFERENCE
=============


Processing Model Results
------------------------

To post-process model results, use PostProcessModelResults.py. The help information can be accessed by running:

``python PostProcessMOdelResults.py --help``

.. automodule:: PostProcessModelResults
    :members:
    :undoc-members:
    :show-inheritance:

This will rewrite both GITM and SAMI3 output files to NetCDF files. SAMI3 can optionally be interpolated:

Interpolation
-------------

SAMI3_ESMF_Regrid
^^^^^^^^^^^^^^^^^

.. automodule:: SAMI3_ESMF_Regrid
    :members:
    :undoc-members:
    :show-inheritance:



Reading Data
------------

GITM
^^^^

.. automodule:: utility_programs.read_routines.GITM
    :members:
    :undoc-members:
    :show-inheritance:

SAMI3
^^^^^

.. automodule:: utility_programs.read_routines.SAMI
    :members:
    :undoc-members:
    :show-inheritance:



Plotting
--------


After converting files to netCDF, you can plot them using the following:

.. automodule:: basic_plots_from_netcdf
    :members:

More plotting routines can be accessed through the utility_programs.plotting_routines module:

.. automodule:: utility_programs.plotting_routines
    :members:
    :undoc-members:
    :show-inheritance:

Utilities
---------

A number of useful utilities are available. 

These are not very well organized, but they are all available in the utility_programs module:

.. automodule:: utility_programs.filters
    :members:
    :undoc-members:
    :show-inheritance:

.. automodule:: utility_programs.utils
    :members:
    :undoc-members:
    :show-inheritance:

.. automodule:: utility_programs.time_conversion
    :members:
    :undoc-members:
    :show-inheritance:


