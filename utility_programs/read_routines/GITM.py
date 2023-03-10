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
                  dtime_storm_start:datetime.datetime
                  file_pattern:str = '3DALL*.bin',
                  t_start_idx:int = 0, 
                  t_end_idx:int = -1, 
                  cols:list = ['all'], 
                  century_prefix:str = '20',
                  ):
    
    
    
    gitm_files = np.sort(
    glob.glob(os.path.join(file_path, file_pattern)))
    if len(gitm_files) == 0:
        raise ValueError("No GITM Binary files found!",
                         "found these files: ",
                         os.listdir(file_path))
    
    gitm_dtimes = []
    for i in gitm_files:
        yy, MM, dd = i[-17:-15], i[-15:-13], i[-13:-11]
        hr, mm, sec = i[-10:-8], i[-8:-6], i[-6:-4]
        gitm_dtimes.append(
            try:
                datetime.datetime(
                    int(century_prefix + yy), int(MM), int(dd),
                    int(hr), int(mm), int(sec)))
            except ValueError:
                raise ValueError(
                    "GITM file name does not match expected format, "
                                 "filename %s cannot be parsed" % i)
    
    if t_start_idx != 0:
        start_idx = gitm_dtimes.index(
            dtime_storm_start - datetime.timedelta(hours=t_start_idx))
        gitm_dtimes = gitm_dtimes[t_start_idx:]
        gitm_files = gitm_files[t_start_idx:]
    
    if t_end_idx != -1:
        end_idx = gitm_dtimes.index(
            dtime_storm_start + datetime.timedelta(hours=t_end_idx))
        gitm_dtimes = gitm_dtimes[:t_end_idx]
        gitm_files = gitm_files[:t_end_idx]
        
    gitm_file_0 = read_routines.read_GITM(gitm_files[0])
    vars_to_read = []
    for nvar, var_name in enumerate(gitm_file_0.keys()):
        if var_name in cols or 'all' in cols:
            vars_to_read.append(var_name)
