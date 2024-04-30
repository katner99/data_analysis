from funcs import read_variable, find_nearest, create_profile
from mitgcm_python.grid import Grid
from directories_and_paths import output_path, grid_filepath
from config_options import slice_ranges
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import sys
import numpy as np
import xarray as xr

def plot_period(ax, period, var):
    experiments = ["CTRL", "LENS", "TEMP", "WIND"]
    colors = ["forestgreen", "orangered", "purple", "dodgerblue"]
    filepaths = [
        output_path + "average_" + exp + "_" + period + ".nc"
        for exp in experiments
    ]
    for filepath in filepaths:
        try:
            open(filepath)
        except FileNotFoundError:
            sys.exit(f"Stopped - Could not find input file {filepath}")

    input_data = [xr.open_dataset(filepath, decode_times=False) for filepath in filepaths]
    grid = Grid(grid_filepath)
    lon_range = slice_ranges["lon_range_cont"]
    lat_range = slice_ranges["lat_range_cont"]
    
    
    profiles = [create_profile(data, var, grid, lat_range, lon_range, timeseries = False) for data in input_data]
    z = input_data[0].Z.values
      
    for i, profile in enumerate(profiles):
        ax.plot(profile, z, color=colors[i])
    
    if var == "THETA":
        ax.set_xlim([-2, 2])
        ax.set_xticks(np.arange(-2, 2, 1))
        ax.set_ylim([-1000, 0])
    else:
        ax.set_xlim([33, 35])
        ax.set_xticks(np.arange(33.5, 35, 0.5))
        ax.set_ylim([-600, 0])
    ax.set_title(period)
    ax.grid(alpha=0.8)
    if period != "1920-1950":
        ax.get_yaxis().set_ticklabels([])
    
def main():
    experiment_names = ["pre-industrial", "RCP8.5", "thermodynamic f.", "wind f."]
    fig, axs = plt.subplots(nrows=1, ncols=3, gridspec_kw={"hspace": 0.05, "wspace": 0.05})
    #axs = axs.flatten()
    var = "SALT"
    plot_period(axs[0], "1920-1950", var)
    plot_period(axs[1], "2000-2030", var)
    plot_period(axs[2], "2070-2100", var)
    plt.subplots_adjust(bottom=0.3, wspace=0.33)
    axs[0].set_ylabel("Depth (m)")
    axs[1].legend(labels=experiment_names,loc='upper center', 
             bbox_to_anchor=(0.5, -0.2),fancybox=False, shadow=False, ncol=4)
    if var == "THETA":
        plt.suptitle("Potential temperature (degC)")
    else:
        plt.suptitle("Salinity")
    fig.savefig(f"profile_{var}.png")
    plt.show()
        

if __name__ == "__main__":
    main()
    
