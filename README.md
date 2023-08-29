[![Documentation Status](https://readthedocs.org/projects/sami3-gitm-python/badge/?version=latest)](https://sami3-gitm-python.readthedocs.io/en/latest/?badge=latest)
    

# SAMI3-GITM-python

This repository contains many scripts to deal with both SAMI3 and GITM outputs. Primarily, they are converted to NetCDF format and the remaining analysis is left to the user. However, many examples of analyses and other processing options are available. 


Please contact the author or fill out a GitHub issue if you notice any problems. 


Further documentation is available at https://sami3-gitm-python.readthedocs.io/en/latest/

---

## Installation

Git clone, get on this branch. 

```
git clone git@github.com:abukowski21/SAMI3-GITM-python.git
cd SAMI3-GITM-python/
git switch -c develop origin/develop
```

> Older versions of Git may require different commands... If the above doesn't work, try `git checkout develop`


To ensure compatibility, an Anaconda environment is available. Install it with:

`conda env create -f python-env.yml && conda activate SAMI3-GITM`

> To create the environment with another name, edit the first line of `python-env.yml`, or use the `-n` flag.
> Anaconda installation information can be found at [this link](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file)


## USAGE

Data can be read directly from binary format and plotted, though there are scripts to postprocess
these data into netCDF format. This will also interpolate SAMI3 model outputs to a "regular" grid
in geographic coordinates.

To postprocess model outputs: `python PostProcessModelResults.py [args]`

To generate plots with these postprocessed outputs: `python gitm_basic_plots.py -gitm_data_path gitm_data 201105211340 -k -m --cols Rho [e-] --plot_start_delta=3 --plot_end_delta=8 --save_or_show=show -lat_lim 65`

`python sami-fieldline-plots.py 2011052112 20110520 ~/scratch/GITM-simstorm-run1/sami-gitm-coupled -out_path . --plot_start_delta 2 --plot_end_delta 6 --cols "edens" --interpolate --plot_type diff --fpeak`


## Notes:
- Run any python script with the `-h` flag to see available arguments
- All routines (including those in `utility_programs`) can be called in other scripts for further analysis (in a Jupyter Notebook, for example)
- These scripts can handle both GITM and SAMI model outputs.
- Modifications to existing functions should be fairly easy (adding more models, different types of plots, etc.)
- Contact the author with any questions, suggestions, issues, etc.
