from plots import compare_contour_plots_LATLON, create_mask
from funcs import read_variable, find_nearest
from directories_and_paths import *

from mitgcm_python.grid import Grid

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import ticker

import numpy as np
import xarray as xr


def main():
    var = "SIfwfrz"
    
    # filepath_current = output_path + "PAS_wind01/output/207001/MITgcm/output.nc"
    filepath_future = output_path + "slice_averages/" + "WIND_ensemble_mean_2090.nc"
    filepath_current = output_path + "slice_averages/" + "CTRL_ensemble_mean_2090.nc"

    input_data_current = xr.open_dataset(filepath_current)
    input_data_future = xr.open_dataset(filepath_future)
    
    lat = input_data_current.YC.values[:300]
    lon = input_data_current.XC.values
    ice_mask = input_data_current.maskC.values[0,:300,:]
    depth = input_data_current.Depth.values[:300]
    grid = Grid(grid_filepath)

    depth_range = [find_nearest(input_data_current["Z"].values, -200), find_nearest(input_data_current["Z"].values, -700)]
    
    data_current = read_variable(input_data_current, var, grid, depth_range)[:300,...]
    data_current = data_current*(3600*24*31*365)/1000
    data_future = read_variable(input_data_future, var, grid, depth_range)[:300,...]
    data_future = data_future*(3600*24*31*360)/1000

    print(np.shape(data_current), np.shape(lat))

    # set mask
    [land_mask, mask, colors] = create_mask(depth, ice_mask)

    # set up the grid
    [X, Y] = np.meshgrid(lon, lat)
        
    color_scheme = "gnuplot2"
        
    # graph parameters
    graph_params = {
        "figsize": (18, 5),
        "font_size": 12,
        "low_val": -0.01,
        "high_val": 0,
        "step": 15,
        "low_val_anom": -0.5,
        "high_val_anom": 0.5,
        "ticks_anom": np.arange(-1.5, 2, 0.5)
    }

    # create subplots with each variable on a new line
    fig, axs = plt.subplots(nrows=1, ncols=3, gridspec_kw={"hspace": 0.5, "wspace": 0.4}, figsize=graph_params["figsize"])
        
    axs = axs.flatten()

    # , ticks=np.arange(graph_params["low_val"], graph_params["high_val"], 0.5) , levels=np.linspace(graph_params["low_val"], graph_params["high_val"], graph_params["step"]
    cs = axs[0].contourf(X, Y, data_current, cmap=color_scheme, levels=np.linspace(graph_params["low_val"], graph_params["high_val"], graph_params["step"]), extend="min")
    axs[0].contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
    axs[0].contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
    fig.colorbar(cs, ax=axs[0])
    axs[0].set_title("Pre-industrial", fontsize=graph_params["font_size"], weight="bold")
    axs[0].set_ylabel("Latitude", fontsize=graph_params["font_size"])
    axs[0].set_xlabel("Longitude", fontsize=graph_params["font_size"])

    cs = axs[1].contourf(X, Y, data_future, cmap=color_scheme, levels=np.linspace(graph_params["low_val"], graph_params["high_val"], graph_params["step"]), extend="min")
    axs[1].contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
    axs[1].contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
    fig.colorbar(cs, ax=axs[1])
    axs[1].set_title("High Winds", fontsize=graph_params["font_size"], weight="bold")
    axs[1].set_ylabel("Latitude", fontsize=graph_params["font_size"])
    axs[1].set_xlabel("Longitude", fontsize=graph_params["font_size"]) 
        
    cs = axs[2].contourf(X, Y, data_future - data_current, cmap="seismic", levels=np.linspace(graph_params["low_val_anom"], graph_params["high_val_anom"], graph_params["step"]), extend="both")
    axs[2].contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
    axs[2].contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
    fig.colorbar(cs, ax=axs[2])
    axs[2].set_title("pre-industrial - wind", fontsize=graph_params["font_size"])
    axs[2].set_ylabel("Latitude", fontsize=graph_params["font_size"])
    axs[2].set_xlabel("Longitude", fontsize=graph_params["font_size"]) 

    fig.suptitle(var, fontsize=16)
    
    file_out ="comparison"+var+".png"
    fig.savefig(file_out)

    plt.show()

if __name__ == '__main__':
    main() # run the program
