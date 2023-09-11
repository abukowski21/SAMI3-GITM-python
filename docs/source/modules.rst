API REFERENCE
=============


Processing Model Results
------------------------

To post-process model results, use PostProcessModelResults.py. The help information can be accessed by running:

``python PostProcessMOdelResults.py --help``


More functionality can be unlocked specifically through running RegridSami...

RegridSami
^^^^^^^^^^

.. automodule:: RegridSami
    :members:
    :undoc-members:


Under the hood, all interpolations are performed with `LinearNDInterpolator <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.LinearNDInterpolator.html>`_. This is a wrapper around `Qhull <http://www.qhull.org/>`_.


For more information (and to unlock even more functionality), use the utility_programs.interpolate_outputs module:

utility_programs.interpolate_outputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: utility_programs.interpolate_outputs
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


