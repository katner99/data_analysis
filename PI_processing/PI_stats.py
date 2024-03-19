import os
import numpy as np
import xarray as xr
from funcs import find_nearest
from stats import stat_timeseries
from directories_and_paths import *
from config_options import *
from mitgcm_python.grid import Grid
from mitgcm_python.utils import add_time_dim

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

def stats_vs_timeseries():
    var = "theta"
    ensemble = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    experiments = ["CTRL", "TEMP"]
    
    data = read_timeseries(experiments, ensemble, var)
    experiment_names =  ["pre-industrial control", "thermodynamic forcing"]
    
    stat_timeseries(data, ensemble, experiments, experiment_names)
    
if __name__ == '__main__':
    stats_vs_timeseries()