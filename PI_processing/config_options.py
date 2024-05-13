"""
Module Docstring:

This module provides configurations and utilities for analyzing oceanographic data using 
MITgcm model outputs.

Dependencies:
- xarray (imported as xr)
- numpy (imported as np)
- Grid class from mitgcm_python.grid module
- find_nearest function from funcs module
- grid_filepath from directories_and_paths module

Variables:
- grid (Grid): Instance of the Grid class initialized with the grid file path.
- grid_file (xarray.Dataset): Opened dataset of the grid file with time decoding disabled.
- depth_range (list): Range of depths for analysis, specified as [min_depth, max_depth].
- sv (float): Scaling factor for velocity.
- lat_bin (int): Number of latitude bins.
- lon_bin (int): Number of longitude bins.
- time_bin (int): Number of time bins.
- days_in_month (list): Number of days in each month.
- slice_ranges (dict): range of areas to analyse
"""

import xarray as xr
import numpy as np
from mitgcm_python.grid import Grid
from funcs import find_nearest, read_variable
from directories_and_paths import grid_filepath

grid = Grid(grid_filepath)
grid_file = xr.open_dataset(grid_filepath, decode_times=False)

depth_range = [
    find_nearest(grid_file.Z.values, -200),
    find_nearest(grid_file.Z.values, -700),
]

sv = 10 ** (-6)

lat_bin, lon_bin, time_bin = 192, 288, 365

lon_slices = [238, 250, 243, 255]
lat_slices = [-74, -69.5]

days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

slice_ranges = {
    "lat_range_73": find_nearest(grid_file.YC.values, -73),
    "lon_range_73": [
        find_nearest(grid_file.XC.values, 252.8),
        find_nearest(grid_file.XC.values, 255),
    ],
    "lat_range_u1": [
        find_nearest(grid_file.YC.values, -73.5),
        find_nearest(grid_file.YC.values, -71.8),
    ],
    "lon_range_u1": find_nearest(grid_file.XC.values, 238),
    "lat_range_u2": [
        find_nearest(grid_file.YC.values, -73),
        find_nearest(grid_file.YC.values, -69.5),
    ],
    "lon_range_u2": find_nearest(grid_file.XC.values, 250),
    "lat_range_u3": [
        find_nearest(grid_file.YC.values, -74),
        find_nearest(grid_file.YC.values, -70),
    ],
    "lon_range_u3": find_nearest(grid_file.XC.values, 243),
    "lat_range_u4": [
        find_nearest(grid_file.YC.values, -72),
        find_nearest(grid_file.YC.values, -70),
    ],
    "lon_range_u4": find_nearest(grid_file.XC.values, 255),
    "lon_range_cont": [
        find_nearest(grid_file.XC.values, 250), 
        find_nearest(grid_file.XC.values, 260)],
    "lat_range_cont": [
        find_nearest(grid_file.YC.values, -75), 
        find_nearest(grid_file.YC.values, -71.5)],
}


def config_timeseries(var, data, experiments, ensemble):
    """
    Configure plotting information for time series data.

    Parameters:
    - var (str): Variable to be plotted. Either "melt" or "theta".
    - data (array_like): Time series data to be plotted.
    - experiments (str or list): Names of experiments.
    - ensemble (int): Number of ensemble members.

    Returns:
    - dict: Plotting information dictionary containing the following keys:
        - "var" (str): Variable name.
        - "data" (array_like): Time series data.
        - "experiments" (str or list): Names of experiments.
        - "ensemble_members" (int): Number of ensemble members.
        - "ylabel" (str): Label for the y-axis.
        - "time" (int): Total number of months.
        - "xlabel" (array_like): Labels for the x-axis.
        - "x_lim" (list): Limits for the x-axis.
        - "smooth" (int): Smoothing factor.
        - "shade_range" (bool): Indicates whether to shade the range.
        - "experiment_full" (list): Full names of experiments.
        - "linearity" (bool, optional): Indicates whether linearity is considered.
        - "warming" (bool, optional): Indicates whether warming is considered.
    """
    if var == "melt":
        plot_info = {
            "var": var,
            "data": data,
            "experiments": experiments,
            "ensemble_members": ensemble,
            "ylabel": "melt (Gt/yr)",
            "time": 181 * 12,
            "xlabel": np.arange(1921, 2100, 10),
            "x_lim": [12, (2099 - 1920) * 12],
            "smooth": 2,
            "shade_range": True,
            "experiment_full": [
                "pre-industrial control",
                "high emissions",
                "thermodynamic forcing",
                "wind forcing",
            ],
        }

    elif var == "theta":
        plot_info = {
            "var": var,
            "data": data,
            "experiments": experiments,
            "ensemble_members": ensemble,
            "ylabel": "Temperature (degC)",
            "time": 181 * 12,
            "xlabel": np.arange(1921, 2100, 10),
            "x_lim": [12, (2099 - 1920) * 12],
            "smooth": 2,
            "shade_range": True,
            "experiment_full": [
                "pre-industrial control",
                "high emissions",
                "thermodynamic forcing",
                "wind forcing",
            ],
            "linearity": True,
            "warming": True,
        }

    return plot_info


def config_comparison(var, input_data, grid, period = None, var_name = None):
    """
    Configures plotting parameters and data for a comparison plot of various oceanographic variables.

    Parameters:
    - var (str): Variable to be plotted (e.g., "THETA", "SIheff", "SALT", "SHIfwFlx").
    - input_data (list): List of input datasets.
    - grid (object): Grid object for the model grid.
    - period (str): Time period for the comparison plot.

    Returns:
    - data (list): List of data arrays to be plotted.
    - graph_params (dict): Dictionary containing parameters for the main plot.
    - graph_params_anom (dict): Dictionary containing parameters for the anomaly plot.
    """
    if var == "trend":
        if var_name == "oceFWflx":
            data = [-input[var].values*((3600*24*365)**2)/10 for input in input_data]
            pval = [input["pvalue"].values for input in input_data]
            color_scheme = "PRGn_r"
            anom = 100000
            min_val = -250000
            max_val = 250000
            print(min_val, max_val)
            title = f"Freshwater Fluxes (m/yr/century)"
            interval = 50000
            interval_anom = 50000
            color_sceme_anom = "PRGn_r"
        elif var_name == "SHIfwFlx":
            data = [input[var].values*((3600*24*365)**2)/10 for input in input_data]
            pval = [input["pvalue"].values for input in input_data]
            color_scheme = "PRGn_r"
            anom = 100000
            min_val = -250000
            max_val = 250000
            print(min_val, max_val)
            title = f"melt (m/yr/century)"
            interval = 50000
            interval_anom = 50000
            color_sceme_anom = "PRGn_r"

        elif var_name == "SALT":
            data = [input[var].values[:,:,0]*3600*24*365*100 for input in input_data]
            color_scheme = "PRGn_r"
            anom = 0.25
            min_val = int(np.min(data[0]))
            max_val = int(np.max(data[0])) + 1
            title = f"Salinity trend per century"
            interval = 1
            interval_anom = 0.25
            color_sceme_anom = "PRGn_r"

        graph_params = {
            "font_size": 12,
            "low_val": min_val,
            "high_val": max_val,
            "interval": interval,
            "color_scheme": color_scheme,
            "step": 20,
            "title": title,
            "pvalue": pval
        }

    else:    
        if var == "THETA":
            depth_range = [
                find_nearest(input_data[0]["Z"].values, -200),
                find_nearest(input_data[0]["Z"].values, -700),
            ]
            data = [read_variable(input, var, grid, depth_range) for input in input_data]
            color_scheme = "coolwarm"
            color_sceme_anom = color_scheme
            anom = 1.5
            min_val = -2
            max_val = 2.1
            title = f"Average potential temperature between 200 and 700m {period}"
            interval = 1
            interval_anom = 0.5
        # sea ice thickness
        elif var == "SIheff":
            data = [read_variable(input, var, grid) for input in input_data]
            color_scheme = "jet"
            anom = 1
            min_val = 0
            max_val = 2
            interval = 0.25
            interval_anom = 0.25
            title = f"Sea ice thickness (m) {period}"
            color_sceme_anom = "PRGn_r"
        # sea ice tickness
        elif var in ["SIheff", "oceFWflx", "SIfwmelt", "SIfwfrz", "EXFvwind", "oceQnet"]:
            data = [
                read_variable(input, var, grid) * 3600 * 24 * 365 / 1000
                for input in input_data
            ]
            color_scheme = "PRGn_r"
            anom = 1.5
            min_val = -5
            max_val = 5
            title = f"Freshwater fluxes m/yr {period}"
            color_sceme_anom = "PRGn_r"

        # salinity
        elif var == "SALT":
            data = [read_variable(input, var, grid) for input in input_data]
            color_scheme = "PRGn_r"
            anom = 0.25
            min_val = int(np.min(data[0]))
            max_val = int(np.max(data[0])) + 1
            title = f"Surface Salinity {period}"
            interval = 1
            interval_anom = 0.25
            color_sceme_anom = "PRGn_r"

        elif var == "SHIfwFlx":
            data = [
                -read_variable(input, var, grid) * 3600 * 24 * 30 * 10 ** (-3)
                for input in input_data
            ]
            color_scheme = "rainbow"
            color_sceme_anom = "PRGn_r"
            anom = 25
            min_val = 0
            max_val = 76
            title = f"ice shelf basal melt rate ({period}) (m.w.e./yr)"
            interval = 20
            interval_anom = 10

        graph_params = {
            "font_size": 12,
            "low_val": min_val,
            "high_val": max_val,
            "interval": interval,
            "color_scheme": color_scheme,
            "step": 20,
            "title": title,
        }
            

    # graph parameters 1.333112e-05 -1.0214189e-06
    

    graph_params_anom = {
        "font_size": 12,
        "low_val": -anom,
        "high_val": anom,
        "interval": interval_anom,
        "color_scheme": color_sceme_anom,
        "step": 20,
    }
    return data, graph_params, graph_params_anom
