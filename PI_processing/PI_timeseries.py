import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np
import sys
import datetime
import netCDF4 as nc
from plots import make_interannual_timeseries, compare_timeseries_TSM
import numpy as np


if __name__ == "__main__":
    filepath_ens07 = "/data/oceans_output/shelf/katner33/PIctrl_output/PAS_ctrl07/output/"
    filepath_ens08 = "/data/oceans_output/shelf/katner33/PIctrl_output/PAS_ctrl08/output/"
    filename = "output.nc"
    start_year = 2070
    n_years = 29
    var = "THETA"
    [theta07,time]=make_interannual_timeseries(filepath_ens07, filename, start_year, n_years, var, plot = False, save = False, show = False)
    [theta08,time]=make_interannual_timeseries(filepath_ens08, filename, start_year, n_years, var, plot = False, save = False, show = False)
    var = "SALT"
    [salt07,time]=make_interannual_timeseries(filepath_ens07, filename, start_year, n_years, var, plot = False, save = False, show = False)
    [salt08,time]=make_interannual_timeseries(filepath_ens08, filename, start_year, n_years, var, plot = False, save = False, show = False)
    var = "SIfwmelt"
    [melt07,time]=make_interannual_timeseries(filepath_ens07, filename, start_year, n_years, var, plot = False, save = False, show = False)
    [melt08,time]=make_interannual_timeseries(filepath_ens08, filename, start_year, n_years, var, plot = False, save = False, show = False)

    compare_timeseries_TSM(time, theta07, theta08, salt07, salt08, melt07, melt08, "ens_07", "ens08")
