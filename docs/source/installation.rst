Installation
============

Instructions
************

This package is based on GitHub and Conda [1]_. To begin, go to a clean folder on your computer where you want to do the install.

First clone the package:

``git clone git@github.com:abukowski21/SAMI3-GITM-python.git``

Checkout the develop branch:

``git checkout -b develop``

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
