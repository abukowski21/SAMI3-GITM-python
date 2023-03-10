# SAMI3-GITM-Python

A collection of Python scripts for reading and plotting both SAMI3 and GITM results.

SAMI3 has been modified to accept GITM neutral atmosphere data, but the code should work for any SAMI3 or GITM run.


All work is being done on the `develop` branch. Here's how to get to it:

`git clone` the repository

`cd` into it

`git checkout develop`

----

Feel free to fork this, make changes you like to the existing files & functions, or make your own stuff! Some ideas on changes that need to be made:




## TODO:
- Move the TODO list to GitHub issues page.
- Turn the GITM & SAMI plotting routines into regular python scripts that can be imported into another script. (read data automatically, return data in a nice way, etc.)
 - Use argparse to call specific scripts
 - read data from linked and/or set directories.
- Make more plot types:
  - GITM longitude plots (to see the altitudinal distribution of features)
  - SAMI3 TEC maps & keograms
  - Latitudinal keograms (to see longitudinal distribution of features)
  - GITM Auroral Joule heating, (3DHME), E-field plots.
- Handle & plot satellite outputs
- Make things a little more user-friendly, and/or document things a bit better. 
 - Adding argparse to the python scripts will make "help" sections available
 - Documentation pages or README's or something similar...
