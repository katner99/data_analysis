"""
Module Docstring:

This module provides configurations and utilities for analyzing oceanographic data using 
MITgcm model outputs.
"""

import xarray as xr
import numpy as np
from mitgcm_python.grid import Grid
import matplotlib
from tools.funcs import find_nearest, ttest_pval
from tools.directories_and_paths import grid_filepath

## Constants
SV = 1e-6  # Sverdrups
G = 9.81  # Gravity (m/s^2)
DAYS_IN_MONTH = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

## Grid data
grid = Grid(grid_filepath)
grid_file = xr.open_dataset(grid_filepath, decode_times=False)
LAT_BIN, LON_BIN, TIME_BIN = 192, 288, 365

## Slice parameters
LON_SLICES = [238, 250, 243, 255]
LAT_SLICES = [-74, -69.5]

DEPTH_RANGE = [find_nearest(grid_file.Z.values, d) for d in (-200, -700)]

## Slice ranges
SLICE_RANGES = {
    "lat_range_73": find_nearest(grid_file.YC.values, -73),
    "lon_range_73": [find_nearest(grid_file.XC.values, lon) for lon in (252.8, 255)],
    "lat_range_u1": [find_nearest(grid_file.YC.values, lat) for lat in (-73.5, -71.8)],
    "lon_range_u1": find_nearest(grid_file.XC.values, 238),
    "lat_range_u2": [find_nearest(grid_file.YC.values, lat) for lat in (-73, -69.5)],
    "lon_range_u2": find_nearest(grid_file.XC.values, 250),
    "lat_range_u3": [find_nearest(grid_file.YC.values, lat) for lat in (-74, -70)],
    "lon_range_u3": find_nearest(grid_file.XC.values, 243),
    "lat_range_u4": [find_nearest(grid_file.YC.values, lat) for lat in (-72, -70)],
    "lon_range_u4": find_nearest(grid_file.XC.values, 255),
    "lon_range_cont": [find_nearest(grid_file.XC.values, lon) for lon in (250, 260)],
    "lat_range_cont": [find_nearest(grid_file.YC.values, lat) for lat in (-75, -71.5)],
}

## Plotting configurations functions
def config_timeseries(var, data, experiments, ensemble):
    """
    Configure plotting info for time series data.

    Parameters:
    - var (str): Variable ("melt" or "theta").
    - data (array_like): Time series data.
    - experiments (str or list): Experiment names.
    - ensemble (int): Number of ensemble members.

    Returns:
    - dict: Plotting configuration.
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
        "experiment_full": ["NONE", "ALL", "WIND", "THERMO"],
        "ylabel": "melt (Gt/yr)" if var == "melt" else "Temperature on continental shelf, 200-700m (°C)",
        **({"linearity": False, "warming": True} if var == "theta" else {}),
    }
    return plot_info


def config_comparison(var, input_data, loc=None, period=None, year=1920, calc_pval="trend"):
    """
    Configures parameters and data for a comparison plot of oceanographic variables.

    Parameters:
    - var (str): Variable to be plotted.
    - input_data (list): List of datasets.
    - loc (str): Location (optional).
    - period (str): Time period (optional).
    - year (int): Base year for mean calculation (optional).
    - calc_pval (str): P-value calculation method (optional).

    Returns:
    - data (list): Data arrays for plotting.
    - graph_params (dict): Parameters for main plot.
    - graph_params_anom (dict): Parameters for anomaly plot.
    """
    def configure_params(data, min_val, max_val, interval, color_scheme, title):
        common_params = {"font_size": 12, "step": 20}
        return {
            **common_params,
            "low_val": min_val,
            "high_val": max_val,
            "interval": interval,
            "color_scheme": color_scheme,
            "title": title,
        }

    def get_pval(data):
        return [np.where(np.sum(p["pvalue"] < 0.05, axis=0) >= 5, 0, 1) for p in data] if calc_pval == "trend" else ttest_pval(data)

    if loc == "trend":
         # the followinf makes it possible to display SHIfwflx with the ocean freshwater fluxes

        white_trans = matplotlib.colors.colorConverter.to_rgba("white", alpha=0.0)
        trends = {
            "SALT": {"factor": 3600 * 24 * 365 * 100, "min_val": -2, "max_val": 2, "interval": 1, "color_scheme": "PRGn_r", "title": "Salinity trend per century"},
            "oceFWflx": {"factor": 100, "min_val": -1, "max_val": 1, "interval": 0.25, "color_scheme": "PRGn_r", "title": "Freshwater Fluxes (m/yr/century)"},
            "SIheff": {"factor": 100, "min_val": -0.5, "max_val": 0.5, "interval": 0.5, "color_scheme": "nipy_spectral", "title": "Sea Ice thickness trend"},
            "SHIfwFlx": {"factor": -100, "min_val": -6, "max_val": 6, "interval": 2.5, "color_scheme": matplotlib.colors.LinearSegmentedColormap.from_list("rb_cmap", ["steelblue", "skyblue", white_trans, "orchid", "purple"], 512), "title": "Melt (m/yr/century)"},
        }
        trend = trends.get(var, {})
        data = [input.trend.mean(dim="ensemble_member").values * trend.get("factor") for input in input_data]
        graph_params = configure_params(data, trend.get("min_val"), trend.get("max_val"), trend.get("interval"), trend.get("color_scheme"), trend.get("title"))
        graph_params["pvalue"] = get_pval(input_data)
        
    elif loc == "mean":
        means = {
            "oceFWflx": {"factor": 3600 * 24 * 365 / 1000, "min_val": -4, "max_val": 4, "interval": 1, "color_scheme": "PRGn_r", "title": f"Freshwater fluxes ({var}) m/yr"},
            "SIfwmelt": {"factor": 3600 * 24 * 365 / 1000, "min_val": -4, "max_val": 4, "interval": 1, "color_scheme": "PRGn_r", "title": f"Freshwater fluxes ({var}) m/yr"},
            "SIfwfrz": {"factor": 3600 * 24 * 365 / 1000, "min_val": -4, "max_val": 4, "interval": 1, "color_scheme": "PRGn_r", "title": f"Freshwater fluxes ({var}) m/yr"},
            "SIarea": {"factor": 1, "min_val": 0, "max_val": 1, "interval": 0.25, "color_scheme": "seismic_r", "title": f"Sea Ice {var}"},
            "SIheff": {"factor": 1, "min_val": 0, "max_val": 1, "interval": 0.25, "color_scheme": "seismic_r", "title": f"Sea Ice {var}"},
            "THETA": {"factor": 1, "min_val": 0, "max_val": 1, "interval": 0.25, "color_scheme": "coolwarm", "title": "Surface temperature"},
        }
        mean = means.get(var, {})
        if var == "oceFWflx":
            data = [np.nanmean(input[var].mean(dim="ensemble_member").values[(year - 1920) * 12:(year + 10 - 1920) * 12,...] * mean.get("factor"), axis=0) for input in input_data]
        else:
            data = [np.nanmean(input[var].values[(year - 1920) * 12:(year + 10 - 1920) * 12,...] * mean.get("factor"), axis=0) for input in input_data]
        graph_params = configure_params(data, mean.get("min_val"), mean.get("max_val"), mean.get("interval"), mean.get("color_scheme"), mean.get("title"))
    
    else:
        default = {
            "THETA": {"depth_range": [find_nearest(input_data[0]["Z"].values, -200), find_nearest(input_data[0]["Z"].values, -700)], "factor": 1, "min_val": -2, "max_val": 2.1, "interval": 1, "color_scheme": "coolwarm", "title": f"Average potential temperature between 200 and 700m {period}"},
            "SIheff": {"factor": 1, "min_val": 0, "max_val": 2, "interval": 0.25, "color_scheme": "jet", "title": f"Sea ice thickness (m) {period}"},
            "oceFWflx": {"factor": 3600 * 24 * 365 / 1000, "min_val": -4, "max_val": 4, "interval": 1, "color_scheme": "PRGn_r", "title": f"Freshwater fluxes ({var}) m/yr"},
        }
        config = default.get(var, {})
        data = [np.nanmean(input[var].values[(year - 1920) * 12:(year + 10 - 1920) * 12] * config.get("factor"), axis=0) for input in input_data]
        graph_params = configure_params(data, config.get("min_val"), config.get("max_val"), config.get("interval"), config.get("color_scheme"), config.get("title"))
    
    return data, graph_params, graph_params