"""
Here we have some more complicated GITM plotting routines.

Currently, the data we can plot are:
    - difference between two storms
    - difference between two bandpass filtered storms

And the plot types are:
    - keograms
    - maps
"""

import gc
import argparse
import geopandas
from scipy import signal
from aetherpy.io import read_routines
from utility_programs.plot_help import UT_from_Storm_onset, plotting_routines
import datetime
import time
import numpy as np
from multiprocessing import Pool
import os
from tqdm.auto import tqdm
import matplotlib
import matplotlib.pyplot as plt
import glob

matplotlib.use("Agg")






if __name__ == "__main__":
    
