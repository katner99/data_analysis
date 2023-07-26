import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt

import numpy as np
import xarray as xr
import sys

from directories_and_paths import *
from funcs import find_nearest, make_timeseries
from plots_1D import comparison
from scipy.ndimage.filters import gaussian_filter1d
from mitgcm_python.grid import Grid

def main():
    var_name = [ "wind 07", "wind 08", "wind 09"]
    title = "SEA ICE"
    
    filepath_1 = output_path + "CUR_wind07_slice.nc"
    filepath_2 = output_path + "CUR_wind08_slice.nc" 
    filepath_3 = output_path + "CUR_wind09_slice.nc"

    data_1 = xr.open_dataset(filepath_1)
    data_2 = xr.open_dataset(filepath_2)
    data_3 = xr.open_dataset(filepath_3)
    # set up grid
    #grid = Grid(grid_filepath)
    grid_file = xr.open_dataset(grid_filepath)
    
    # set lat, lon, and depth range
    depth_range = [find_nearest(grid_file.Z.values, -200), find_nearest(grid_file.Z.values, -700)]
    lon_range = [find_nearest(grid_file.XC.values, 250), find_nearest(grid_file.XC.values, 260)]
    lat_range = [find_nearest(grid_file.YC.values, -76), find_nearest(grid_file.YC.values, -72)]

    #timeseries_new = make_timeseries("SIheff", data_new, grid, lat_range, lon_range, depth_range, time = len(data_new.time.values))
    #timeseries_old = make_timeseries("SIheff", data_old, grid, lat_range, lon_range, depth_range, time = len(data_new.time.values))
    
    timeseries_1 = np.nanmax(data_1.SIheff.values, axis = (-1, -2))
    timeseries_2 = np.nanmax(data_2.SIheff.values, axis = (-1, -2))   
    timeseries_3 = np.nanmax(data_3.SIheff.values, axis = (-1, -2))  

    #print(timeseries_new[30], timeseries_old[30])
    data = [timeseries_1, timeseries_2, timeseries_3]
    comparison(data, "SIheff", var_name, title, "current_wind_ice.png", ylabel = "Ses-ice thickness" )
    
if __name__ == "__main__":
    main()
    
