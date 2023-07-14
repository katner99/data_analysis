from plots import compare_contour_plots_LATLON, create_mask
from funcs import read_variable, find_nearest
from directories_and_paths import *
from mitgcm_python.grid import Grid

import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import sys
import numpy as np
import xarray as xr

def main():
    var = "EXFuwind"
    
    filepath_current = output_path + "PAS_temp01/output/207001/MITgcm/output.nc"
    filepath_future = output_path+ "PAS_ctrl01/output/207001/MITgcm/output.nc"


    input_data_current = xr.open_dataset(filepath_current)
    input_data_future = xr.open_dataset(filepath_future)
    
    lat = input_data_current.YC.values
    lon = input_data_current.XC.values
    ice_mask = input_data_current.maskC.values[0,:,:]
    depth = input_data_current.Depth.values
    grid = Grid(grid_filepath)

    depth_range = [find_nearest(input_data_current["Z"].values, -200), find_nearest(input_data_current["Z"].values, -700)]
    
    data_current = read_variable(input_data_current, var, grid, depth_range)
    data_future = read_variable(input_data_future, var, grid, depth_range)

    # set mask
    [land_mask, mask, colors] = create_mask(depth, ice_mask)

    # set up the grid
    [X, Y] = np.meshgrid(lon, lat)
        
    color_scheme = "coolwarm"
        
    # graph parameters
    graph_params = {
        "figsize": (18, 5),
        "font_size": 12,
        "low_val": -2.5,
        "high_val": 2.5,
        "step": 15,
        "low_val_anom": -0.05,
        "high_val_anom": 0.05,
        "ticks_anom": np.arange(-0.1, 0.1, 0.02)
    }

    # create subplots with each variable on a new line
    fig, axs = plt.subplots(nrows=1, ncols=3, gridspec_kw={"hspace": 0.5, "wspace": 0.4}, figsize=graph_params["figsize"])
        
    axs = axs.flatten()

    cs = axs[0].contourf(X, Y, data_current, cmap=color_scheme, extend="both", levels=np.linspace(graph_params["low_val"], graph_params["high_val"], graph_params["step"]))
    axs[0].contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
    axs[0].contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
    fig.colorbar(cs, ax=axs[0], ticks=np.arange(int(graph_params["low_val"]), int(graph_params["high_val"]), 0.5))
    axs[0].set_title("Thermodynamic", fontsize=graph_params["font_size"], weight="bold")
    axs[0].set_ylabel("Latitude", fontsize=graph_params["font_size"])
    axs[0].set_xlabel("Longitude", fontsize=graph_params["font_size"])

    cs = axs[1].contourf(X, Y, data_future, cmap=color_scheme, extend="both", levels=np.linspace(graph_params["low_val"], graph_params["high_val"], graph_params["step"]))
    axs[1].contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
    axs[1].contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
    fig.colorbar(cs, ax=axs[1], ticks=np.arange(int(graph_params["low_val"]), int(graph_params["high_val"]), 0.5))
    axs[1].set_title("Pre-industrial control", fontsize=graph_params["font_size"], weight="bold")
    axs[1].set_ylabel("Latitude", fontsize=graph_params["font_size"])
    axs[1].set_xlabel("Longitude", fontsize=graph_params["font_size"]) 
        
    cs = axs[2].contourf(X, Y, data_future - data_current, cmap="seismic", extend="both", levels=np.linspace(graph_params["low_val_anom"], graph_params["high_val_anom"], graph_params["step"]))
    axs[2].contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
    axs[2].contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
    fig.colorbar(cs, ax=axs[2], ticks=graph_params["ticks_anom"])
    axs[2].set_title("difference in the zonal winds", fontsize=graph_params["font_size"])
    axs[2].set_ylabel("Latitude", fontsize=graph_params["font_size"])
    axs[2].set_xlabel("Longitude", fontsize=graph_params["font_size"]) 

    fig.suptitle("Zonal winds", fontsize=16)
    
    file_out ="comparison"+var+".png"
    fig.savefig(file_out)

    plt.show()

if __name__ == '__main__':
    main() # run the program
