Introduction
============

``SAMI3-GITM-Python`` provides python scripts to allow easy manipulation of SAMI3 and GITM model outputs.

This project will not run any models for you, only manipulate the outputs into a more user-friendly format. 

For example, GITM outputs can be converted from several ``.bin`` files (one for each time-step) to a single (or multiple) netCDF files. SAMI3 outputs can be converted to netCDF files indexed with the native grid ``(nlt, nf, nz)`` or interpolated to a standard geographic grid.




Project Scope
*************

Again, these scripts do not interface with the models at all. They will only read outputs and convert files.

Some basic plotting scripts have been provided, though more complicated plotting is left to the user. 

This package has no plans be ported to a user-installable python package. If you wish to use anything in this project in your own python scripts, add the package to your `$PATH` (see usage for more info).


I've done my best to capture as many use-cases as possible. Please fill out an `Issue <https://github.com/abukowski21/SAMI3-GITM-python/issues>`_ if you notice any problems.
