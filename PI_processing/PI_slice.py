import sys

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import numpy as np
import xarray as xr

from funcs import find_nearest
from directories_and_paths import *
from mitgcm_python.grid import Grid

def main():
    period = "2070-2100"
    var = "salt"
    save = True
    show = True
    
    # load up the file paths for the monster, needed 4
    filepaths = [output_path + "average_" + ens + "_"+period+".nc" for ens in ["CTRL", "LENS", "WIND", "TEMP"]]
    # check if the input files exist
    for filepath in filepaths:
        try:
            open(filepath)
        except FileNotFoundError:
            sys.exit(f"Stopped - Could not find input file {filepath}")

    experiment = ["pre-industrial", "high-emissions", "wind forcing", "thermodynamic f."]

    # load up the input data
    input_data = [xr.open_dataset(filepath, decode_times = False) for filepath in filepaths]
    
    # read in the general variables (these should be the same between the ensembles
    #[lat, lon, ice_mask_temp, depth] = [input_data[0][param].values for param in ["YC", "XC", "maskC", "Depth"]]
    #ice_mask = ice_mask_temp[0,:,:]

    grid = Grid(grid_filepath)

    lat_range = [find_nearest(input_data[0].YC.values, -71.5), find_nearest(input_data[0].YC.values, -70)]
    lon_range = find_nearest(input_data[0].XC.values, 255)
    print(lon_range)
    print(len(input_data[0].XC.values))
    z = input_data[0].Z.values
    lat = input_data[0]["YC"][lat_range[0]:lat_range[1]].values
    ice_mask = input_data[0].maskC.values[:,lat_range[0]:lat_range[1],lon_range]
    color_scheme = "PRGn"
       
    # load the variables and calculate the average over the year
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        
    ctrl_vvel = np.average(input_data[0].UVEL.values[:,:,lat_range[0]:lat_range[1],lon_range], axis = 0, weights = days_in_month)
    lens_vvel = np.average(input_data[1].UVEL.values[:,:,lat_range[0]:lat_range[1],lon_range], axis = 0, weights = days_in_month)
    wind_vvel = np.average(input_data[2].UVEL.values[:,:,lat_range[0]:lat_range[1],lon_range], axis = 0, weights = days_in_month)
    temp_vvel = np.average(input_data[3].UVEL.values[:,:,lat_range[0]:lat_range[1],lon_range], axis = 0, weights = days_in_month)
    
    ctrl_vvel[ice_mask == 0] = np.nan
    lens_vvel[ice_mask == 0] = np.nan
    wind_vvel[ice_mask == 0] = np.nan
    temp_vvel[ice_mask == 0] = np.nan
        
    low_val = -0.05
    high_val = -low_val
    step = 15
    font_size = 16
        
    # create subplots with each variable on a new line
    fig, axs = plt.subplots(nrows=1, ncols=4, gridspec_kw={"hspace": 0.5, "wspace": 0.4}, figsize=(15,5))
        
    axs = axs.flatten()

    # PI_ctrl
    cs = axs[0].contourf(lat, z, ctrl_vvel, cmap=color_scheme, extend="both", levels=np.linspace(low_val, high_val, step))
    axs[0].set_title("Pre-industrial", weight="bold", fontsize = font_size)
    axs[0].set_ylabel("Depth (m)", fontsize = font_size)
    axs[0].set_xlabel("Latitude", fontsize = font_size)
    axs[0].set_ylim([-2000, 0])
        
    # LENS
    cs = axs[1].contourf(lat, z, lens_vvel, cmap=color_scheme, extend="both", levels=np.linspace(low_val, high_val, step))
    axs[1].set_title("high emissions", weight="bold", fontsize = font_size)
    axs[1].set_xlabel("Latitude", fontsize = font_size)
    axs[1].set_ylim([-2000, 0])
        
    # WIND
    cs = axs[2].contourf(lat, z, wind_vvel, cmap=color_scheme, extend="both", levels=np.linspace(low_val, high_val, step))
    axs[2].set_title("Winds", weight="bold", fontsize = font_size)
    axs[2].set_xlabel("Latitude", fontsize = font_size)
    axs[2].set_ylim([-2000, 0])

    # TEMP
    cs = axs[3].contourf(lat, z, temp_vvel, cmap=color_scheme, extend="both", levels=np.linspace(low_val, high_val, step))
    axs[3].set_title("Thermodynamic", weight="bold", fontsize = font_size)
    axs[3].set_xlabel("Latitude", fontsize = font_size)
    axs[3].set_ylim([-2000, 0])
        
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(cs, cax=cbar_ax, ticks = np.arange(-0.5, 0.5, 0.025))
        
    plt.suptitle("Undercurrent 4 "+period, fontsize = font_size, weight = "bold")
        
    fig.savefig("underway_slice"+period+"_4.png")
    plt.show()
        
    
if __name__ == '__main__':
    main() # run the program
