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


For users that do not have Anaconda, the required dependencies can also be installed with `pip`. While supported, this is not necessarily encouraged. Nonetheless, to set up the correct python environment with pip, run:

`python -m pip install -r requirements.txt`


## USAGE

Data can be read directly from binary format and plotted, though there are scripts to postprocess
these data into netCDF format. Users can also interpolate SAMI3 model outputs to a "regular" grid
in geographic coordinates.

To postprocess model outputs into NetCDF format: `python PostProcessModelResults.py [args]`

To generate plots with these postprocessed outputs: `python basic_plots_from_netcdf.py [args]`

Run any script with the `-h` flag to see the available arguments.

### EXAMPLE:

Here is an example of a workflow to look at SAMI outputs. Here we will just create the regridded data & then plot out the values at all times on a map.

```
python PostProcessModelResults.py -sami /path/to/sami_data/ -out /path/to/outputs/ --dtime_sim_start 20110521 --sami_type regrid

python basic_plots_from_netcdf.py /path/to/outputs/ -col edens -out_dir /path/to/save/plots/ -m --alt_cut 650 --loop_var time

```

A ton of extra functionality is available. Run the scripts with the `-h` or `--help` flags to see available arguments, or check the documentation to read more.


Good luck!


## Notes:

- All routines (including those in `utility_programs`) can be called in other scripts for further analysis (in a Jupyter Notebook, for example)
- These scripts can handle both GITM and SAMI model outputs.
- Modifications to existing functions should be fairly easy (adding more models, different types of plots, etc.)
- Contact the author with any questions, suggestions, issues, etc.
- A number of examples of analysis scripts & other misc. usage is available in REFERENCE-examplenotebooks & ACTIVE_ANALYSIS. These will likely need to be moved to the root directory of the repository (and have their paths changed) to run. They are primarily posted to just show the capabilities & usage of the codebase, not for users to run & build on... but nothing is preventing you from doing that (other than missing the source data).

## PULBICATIONS:

Source code for publications can be found in the src_PUBLICATIONS folder. A README in that filder will direct you where to look for the source scripts for each paper which has used this repository.



## Contributing/Licensing


Users are open to collaborate on this project. I do not intend for it to be closed-source or private ever. Feel free to edit all work posted to your heart's content.

Officially, this code is licensed under the GNU General Public License v3.0 - meaning that if you take this source code and modify it, the modifications must be made public, and the original work cited, when sharing with others. What this means for you is that if you do use this for a publication, your modifications must be shared somehow. I would appreciate your contributions to be made into a Pull Request so the source can be (hopefully) improved and to be able to add your publication to src_PUBLICATIONS, though that is not necessary.




