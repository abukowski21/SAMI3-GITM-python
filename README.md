# SAMI3-GITM-python

Here is all the current work being done. See the main branch for info on what needs to be done.


## Please feel free to make any changes you would like!!


Scripts are available for both SAMI3 and GITM data analysis. 
SAMI3 data needs to be post-processed before it can be manipulated with these programs. 
> (See utility_programs for information.)

SAMI3 data can be processed further using InterpSAMI3toGrid.ipynb, allowing it to be plot 
with grographic coordinate systems.

---

## USAGE:

Clone, get on this branch. Best practice is to leave this in the home folder that is backed up. Then link your GITM and SAMI directories:

From inside this directory, create links to GITM and SAMI data directories:

`ln -sfn [path-to-gitm-data] gitm_data`

`ln -sfn [path-to-sami-data] sami_dir`

I like to put the plots into the same directory ad the model outputs. For now, since only GITM plots work:

`mkdir gitm_data/out_plots_gitm` (or whatever you want)

`ln -sfn gitm_data/out_plots_gitm out_plots_gitm` (make sure that first directory is pointing to what you made before)

Then you can call some gitm plots. For help:

`python gitm_basic_plots.py -h`

This will show you all of the available arguments and settings and what not.

To make some beautiful maps and keograps (assuming you set the directories like I have):

`python gitm_basic_plots.py -gitm_data_path gitm_data 201105211340 -k -m --cols Rho [e-] --plot_start_delta=3 --plot_end_delta=8 --save_or_show=show -lat_lim 65`

`python sami-fieldline-plots.py 2011052112 20110520 ~/scratch/GITM-simstorm-run1/sami-gitm-coupled -out_path . --plot_start_delta 2 --plot_end_delta 6 --cols "edens" --interpolate --plot_type diff --fpeak`

These takes a litle while, though progress bars will be displayed.
