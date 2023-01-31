import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np
import sys
import datetime
import netCDF4 as nc
from plots import make_interannual_timeseries
import numpy as np


if __name__ == "__main__":
    filepath = "/data/oceans_output/shelf/katner33/PIctrl_output/PAS_ctrl07/output/"
    filename = "output.nc"
    start_year = 2070
    n_years = 29
    var = "THETA"
    make_interannual_timeseries(filepath, filename, start_year, n_years, var)

