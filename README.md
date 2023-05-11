# SAMI3-GITM-python

Here is all the current work being done. See the main branch for info on what needs to be done.


## Please feel free to make any changes you would like!!


Scripts are available for both SAMI3 and GITM data analysis. 
SAMI3 data needs to be post-processed before it can be manipulated with these programs. 
> 

---

## USAGE:

Git clone, get on this branch. 

`git clone git@github.com:abukowski21/SAMI3-GITM-python`

`cd SAMI3-GITM-python/`

`git checkout develop`

To ensure compatibility, an Anaconda environment is available. Install it with:

`conda env create -f python-env.yml && conda activate SAMI3-GITM`
> To create the environment with another name, edit the first line of `python-env.yml`

> Anaconda installation information can be found at [this link](https://conda.io/projects/conda/en/latest/index.html)

Data can be read directly from binary format and plotted, though there are scripts to postprocess
these data into netCDF format. This will also interpolate SAMI3 model outputs to a "regular" grid
in geographic coordinates.

To Postprocess model outputs: `python PostProcessModelResults.py [args]`

To generate plots with these postprocessed outputs: `python basic_plots_from_netcdf.py [args]`

Run any python script with the `-h` flag to see available arguments


