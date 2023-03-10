import datetime
import gc
import glob
import os
import time
from multiprocessing import Pool

import geopandas
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from aetherpy.io import read_routines
from scipy import signal
from tqdm.auto import tqdm

from utility_programs.plot_help import UT_from_Storm_onset


def get_gitm_data(filepath:str, 
                  t_start_idx:int, 
                  i_end_idx:int, 
                  cols:list, 
                  dtime_storm_start:datetime.datetime):
    