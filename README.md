# SAMI3-GITM-python

Here is all the current work being done. See the main branch for info on what needs to be done.


## Please feel free to make any changes you would like!!


---

## USAGE:

Git clone, get on this branch. 

`git clone git@github.com:abukowski21/SAMI3-GITM-python`

`cd SAMI3-GITM-python/`

`git checkout satellite_interpolations`

To ensure compatibility, an Anaconda environment is available. Install it with:

`conda env create -f python-env.yml && conda activate SAMI3-GITM`
> To create the environment with another name, edit the first line of `python-env.yml`

> Anaconda installation information can be found at [this link](https://conda.io/projects/conda/en/latest/index.html)

Data can be read directly from binary format and plotted, though there are scripts to postprocess
these data into netCDF format. This will also interpolate SAMI3 model outputs to a "regular" grid
in geographic coordinates.

To postprocess model outputs: `python PostProcessModelResults.py [args]`

To generate plots with these postprocessed outputs: `python basic_plots_from_netcdf.py [args]`


## Notes:
- Run any python script with the `-h` flag to see available arguments
- All routines (including those in `utility_programs`) can be called in other scripts for further analysis (in a Jupyter Notebook, for example)
- These scripts can handle both GITM and SAMI model outputs.
- Modifications to existing functions should be fairly easy (adding more models, different types of plots, etc.)
- Contact the author with any questions, suggestions, issues, etc.
