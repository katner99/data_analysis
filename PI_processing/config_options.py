"""
Module Docstring:

This module provides configurations and utilities for analyzing oceanographic data using 
MITgcm model outputs.
"""

import xarray as xr
import numpy as np
from mitgcm_python.grid import Grid
import matplotlib
from funcs import find_nearest, read_variable
from directories_and_paths import grid_filepath

### HANDY VARIABLES ###
sv = 10 ** (-6) # sverdrups
g = 9.81        # acceleration due to gravity

days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

### GRID DATA AND INFO ###
grid = Grid(grid_filepath)
grid_file = xr.open_dataset(grid_filepath, decode_times=False)

lat_bin, lon_bin, time_bin = 192, 288, 365

lon_slices = [238, 250, 243, 255]
lat_slices = [-74, -69.5]

depth_range = [
    find_nearest(grid_file.Z.values, -200),
    find_nearest(grid_file.Z.values, -700),
]

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

### FUNCTIONS CONTAININT CONFIG INFO FOR PLOTTING  CODES ###
def config_timeseries(var, data, experiments, ensemble):
    """
    Configure plotting information for time series data.

    Parameters:
    - var (str): Variable to be plotted. Either "melt" or "theta".
    - data (array_like): Time series data to be plotted.
    - experiments (str or list): Names of experiments.
    - ensemble (int): Number of ensemble members.

    Returns:
    - dict: Plotting information dictionary.
    """
    plot_info = {
        "var": var,
        "data": data,
        "experiments": experiments,
        "ensemble_members": ensemble,
        "time": 181 * 12,
        "xlabel": np.arange(1921, 2100, 10),
        "x_lim": [12, (2099 - 1920) * 12],
        "smooth": 2,
        "shade_range": True,
        "experiment_full": [
            "NONE",
            "ALL",
            "WIND",
            "THERMO",
        ],
    }

    if var == "melt":
        plot_info.update({
            "ylabel": "melt (Gt/yr)",
        })
    elif var == "theta":
        plot_info.update({
            "ylabel": "Temperature on continental shelf, 200-700m (°C)",
            "linearity": False,
            "warming": True,
        })

    return plot_info

def config_comparison(var, input_data, grid, period=None, var_name=None):
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
    common_params = {
        "font_size": 12,
        "step": 20,
    }

    if var == "trend":
        data, pval = None, None
        if var_name == "SALT":
            data = [input[var].values[:, :, 0] * 3600 * 24 * 365 * 100 for input in input_data]
            min_val, max_val, interval, interval_anom, anom, color_scheme = int(np.min(data[0])), int(np.max(data[0])) + 1, 1, 0.25, 0.25, "PRGn_r"
            title = "Salinity trend per century"
        elif var_name in ["oceFWflx", "SIfwfrz", "SIfwmelt"]:
            data = [input[var].mean(dim="ensemble_member").values * 100 for input in input_data]
            pval = [np.where(np.sum(p["pvalue"] < 0.05, axis=0) >= 5, 0, 1) for p in input_data]
            min_val, max_val, interval, interval_anom, anom, color_scheme = -2, 2, 0.5, 0.25, 1, "PRGn_r"
            title = "Freshwater Fluxes (m/yr/century)"
        elif var_name == "SIheff":
            data = [input[var].mean(dim="ensemble_member").values * 100 for input in input_data]
            pval = [input["pvalue"].max(dim="ensemble_member").values for input in input_data]
            min_val, max_val, interval, interval_anom, anom, color_scheme = -0.5, 0.5, 0.5, 0.25, 0.25, "nipy_spectral"
            title = "Sea Ice thickness trend"
        elif var_name == "SHIfwFlx":
            data = [-input[var].mean(dim="ensemble_member").values * 100 for input in input_data]
            pval = [input["pvalue"].max(dim="ensemble_member").values for input in input_data]
            c_purp, c_red, c_blue, c_super, c_white_trans = (matplotlib.colors.colorConverter.to_rgba(c) for c in ['purple', 'orchid', 'skyblue', 'steelblue', 'white'])
            color_scheme = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap', [c_super, c_blue, c_white_trans, c_red, c_purp], 512)
            min_val, max_val, interval, interval_anom, anom = -6, 6, 2.5, 5, 10
            title = "melt (m/yr/century)"
        
        graph_params = {**common_params, "low_val": min_val, "high_val": max_val, "interval": interval, "color_scheme": color_scheme, "title": title, "pvalue": pval}
    elif var == "mean":
        if var_name == "oceFWflx":
            data = [np.nanmean(input[var_name].values[(2070 - 1920) * 12:(2075 - 1920) * 12, ...] * 3600 * 24 * 365 / 1000, axis=0) for input in input_data]
            min_val, max_val, interval, interval_anom, anom, color_scheme = -4, 4, 1, 0.25, 1, "PRGn_r"
            title = "Freshwater fluxes m/yr 2070-2075"
        
        graph_params = {**common_params, "low_val": min_val, "high_val": max_val, "interval": interval, "color_scheme": color_scheme, "title": title}
    else:
        if var == "THETA":
            depth_range = [find_nearest(input_data[0]["Z"].values, -200), find_nearest(input_data[0]["Z"].values, -700)]
            data = [read_variable(input, var, grid, depth_range) for input in input_data]
            min_val, max_val, interval, interval_anom, anom, color_scheme = -2, 2.1, 1, 0.5, 1.5, "coolwarm"
            title = f"Average potential temperature between 200 and 700m {period}"
        elif var == "SIheff":
            data = [read_variable(input, var, grid) for input in input_data]
            min_val, max_val, interval, interval_anom, anom, color_scheme = 0, 2, 0.25, 0.25, 1, "jet"
            title = f"Sea ice thickness (m) {period}"
        elif var in ["SIheff", "oceFWflx", "SIfwmelt", "SIfwfrz", "EXFvwind", "oceQnet"]:
            data = [read_variable(input, var, grid) * 3600 * 24 * 365 / 1000 for input in input_data]
            min_val, max_val, interval, interval_anom, anom, color_scheme = -5, 5, 0.25, 0.25, 1.5, "PRGn_r"
            title = f"Freshwater fluxes m/yr {period}"
        elif var == "SALT":
            data = [read_variable(input, var, grid) for input in input_data]
            min_val, max_val, interval, interval_anom, anom, color_scheme = int(np.min(data[0])), int(np.max(data[0])) + 1, 1, 0.25, 0.25, "PRGn_r"
            title = f"Surface Salinity {period}"
        elif var == "SHIfwFlx":
            data = [-read_variable(input, var, grid) * 3600 * 24 * 30 * 10 ** -3 for input in input_data]
            min_val, max_val, interval, interval_anom, anom, color_scheme = 0, 76, 20, 10, 25, "rainbow"
            title = f"ice shelf basal melt rate ({period}) (m.w.e./yr)"
        
        graph_params = {**common_params, "low_val": min_val, "high_val": max_val, "interval": interval, "color_scheme": color_scheme, "title": title}

    graph_params_anom = {**common_params, "low_val": -anom, "high_val": anom, "interval": interval_anom, "color_scheme": color_scheme}

    return data, graph_params, graph_params_anom
