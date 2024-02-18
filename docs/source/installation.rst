Installation
============

Instructions
************

This package is based on GitHub and Conda [1]_. To begin, go to a clean folder on your computer where you want to do the install.

First clone the package:

``git clone git@github.com:abukowski21/SAMI3-GITM-python.git``

Checkout the ESMF branch:

``git checkout -b esmf``

Enter the directory:

``cd SAMI3-GITM-python``

Create a new conda[1] environment according to the requirements specified by this project:

``conda env create -f python-env.yml && conda activate SAMI3-GITM``

.. note::
    Conda is not required, but it makes things easier. If you do not want to use conda, you can install the packages with pip or manually. To use pip:

    ``pip install -r requirements.txt``

    This usage is not supported, so you may run into dependency issues doing this.


.. [1] Additional information on ``conda`` can be found on their website at the Installation_ and Environment_ pages.

.. _Installation: https://conda.io/projects/conda/en/latest/user-guide/install/index.html

.. _Environment: https://conda.io/projects/conda/en/latest/user-guide/install/index.html

.. _Install-esmf

Installing ESMF
---------------

ESMF is required to run the interpolation code on this branch. This section of the documentation is a work in progress, so contact me with issues & questions.

ESMF (Earth System Modeling Framework) is a robust software suite designed to facilitate the communication between models. We are only using a tiny aspect so can cut some corners in the installation, but the software package is massive and care must be taken to ensure everything iscorrectly installed.

For complete details on how to install ESMF, see their documentation_. I will provide a bit of information on the crucial steps, but this is by no means complete.
The necessary pages in the documentation are the the quick start page_ and this page on more details on the installation specifics_. 

.. _page: https://earthsystemmodeling.org/docs/release/latest/ESMF_usrdoc/node6.html

.. _specifics: https://earthsystemmodeling.org/docs/release/latest/ESMF_usrdoc/node10.html

We need to install ESMF with just the command line tools, howerver we need support for PIO, and NetCDF. Here are general instructions along with some tips. We will be following these instructions_.

.. _instructions: https://earthsystemmodeling.org/docs/release/latest/ESMF_usrdoc/node6.html#SECTION00063000000000000000

#. Ensure that a suitable version of ESMF is downloaded and uncompressed. Not all versions support 3D spherical interpolations in offline mode. I recommend version 8.4.0, but it is not required.
#. Set the required environmental variables. These are: **ESMF_DIR**, **ESMF_COMM**, **ESMF_COMPILER**, **ESMF_NETCDF**, **ESMF_PNETCDF**, and **ESMF_PIO**. Other variables in the documentation can be set, but my testing has shown that these will work.
	* **ESMF_DIR** is the path to the installation location
	* **ESMF_COMM** denotes if ESMF should be installed with support for MPI. Set this to the correct implementation of MPI that your system is using.
	* **ESMF_COMPILER** denotes the Fortran/C++ compiler that ESMF will be built with. typing ``make info`` shows which compiler & MPI impletation are currently planned on being used.
	* **ESMF_NETCDF** enables the ESMF command line tools to read NetCDF files from disk. This is required. See this link_ for more information, but usually this can be set to "nc-config" if the commands ``nc-config`` and ``nf-config`` run successfully on your system. If not, either install (or load) the C and Fortran NetCDF modules.
	* **ESMF_PNETCDF** can be set to "pnetcdf-config" if the command ``pnetcdf-config`` runs successfully. If not, either install or load the modules.
	* **ESMF_PIO** enables MPI executables to read and write NetCDF files. This can be set to "internal" if you would like for ESMF to install this, or "external" if you have this software installed already.
#. Double check everything looks right by running ``make info``. Take a close look at the directories and include paths.
#. Install the command line tools with ``make build_apps``. This will take a while. You can double check that it is correctly installed by running ``./$(ESMF_DIR)/apps/..../ESMF_PrintInfo`` and the information on the EMSF installation will be displayed. Ensure PIO and NetCDF support is enabled.

.. _link: https://earthsystemmodeling.org/docs/release/latest/ESMF_usrdoc/node10.html#sec:netcdf


.. _postinstall:

Post-Installation
-----------------

To run the scripts provided with this package, make sure you are in the correct directory and call them from the command line.

To call any of the scripts in you own code, you will need to add the package to your ``$PATH``. This can be done by adding the following to the top of your code:

.. code-block:: python
    
    import sys
    sys.path.append('/path/to/SAMI3-GITM-python')


For more information on how to use the provided code, continue reading! 


Issues
******

- Older GitHub versions may require a slightly different command to access the correct branch.
- If you do not have ``conda`` installed and don't want to install it, you can install the required packages with ``pip``. 
- If you encounter more problems, please either fill out a GitHub issue or contact the author directly.
