Plotting
########


This section still needs to be polished, A LOT!


The Basics
==========

After postprocessing data, we can plot it using the `basic_plots_from_netcdf.py` script. This script is a very basic plotter and serves to demonstrate how to plot data from a netcdf file. It is not meant to be a comprehensive plotting script. It's a great way to check that your data is being postprocessed correctly.





More Advanced Plotting
======================

Plotting can be very personal. The `basic_plots_from_netcdf.py` script is a good starting point, but you will probably want to make more detailed plots. 


This package provides a lot of functions to help make plots. From the command line, polar_dial_plots.py proivides a nice introduction to additional plotting methids. This will plot any requested column on polar dials, both for the Northern and Southern hemispheres simultaneously.




Even more plotting routines are available in utility_programs/plotting_routines.py. The functions there should be mostly well documented, though are not accessible from the command line. They must be imported into a python script and run from there. See API reference for more information.