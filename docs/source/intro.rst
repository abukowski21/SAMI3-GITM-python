Introduction
============

``SAMI3-GITM-Python`` provides python scripts to allow easy manipulation of SAMI3 and GITM model outputs.

This project will not run any models for you, only manipulate the outputs into a more user-friendly format. 

For example, GITM outputs can be converted from several ``.bin`` files (one for each time-step) to a single (or multiple) netCDF files. SAMI3 outputs can be converted to netCDF files indexed with the native grid (nlt, nf, nz) or interpolated to a standard geographic grid.




Project Scope
***********

Again, these scripts do not interface with the models at all. They will only read outputs and convert files.

Some basic plotting scripts have been provided, though more complicated and/or involved plotting is left to the user. 

Thsi package will not be ported to a user-installable python package. If you wish to use anything in this project ion your own python scripts, add the package to your $PATH (see usage for more info).


I've done my best to capture as many use-cases as possible. Please fill out an Issue if you notice any problems.
