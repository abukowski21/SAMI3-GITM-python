Installation
=====

Quickstart
-----


This package is based on GitHub and Conda [1]_. To begin, go to a clean folder on your computer where you want to do the install.

First clone the package:

``git clone git@github.com:abukowski21/SAMI3-GITM-python.git``

Checkout the correct branch:

``git checkout -b develop``

Enter the directory:

``cd SAMI3-GITM-python```

Create a new conda[1] environment accoring to the specs listed:
``conda env create -f python-env.yml``

And activate:

``conda activate SAMI3-GITM```


.. [1] Additional information on ``conda`` can be found on their website at the Installation_ and Environment_ pages.

.. _Installation: https://conda.io/projects/conda/en/latest/user-guide/install/index.html

.. _Environment: https://conda.io/projects/conda/en/latest/user-guide/install/index.html


Known Issues
-----

- Older GitHub versions may require a slightly different command to access the correct branch.
- If you do not have ``conda`` installed and don't want to install it, you can install the required packages with ``pip``. For help on this, contact the author as the documentation has not been built yet.

