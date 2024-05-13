import sys

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import numpy as np
import xarray as xr

from funcs import find_nearest
from directories_and_paths import *
from mitgcm_python.grid import Grid

from config_options import lat_slices, lon_slices

def main():
    var = "SALT"
    save = True
    show = True
    filepaths = [
        f"{output_path}{exp}_files_temp/{var}_trend.nc"
        for exp in ["CTRL", "LENS", "WIND", "TEMP"]
    ]
    
    for filepath in filepaths:
        try:
            open(filepath)
        except FileNotFoundError:
            sys.exit(f"Stopped - Could not find input file {filepath}")

    experiment = ["pre-industrial", "high-emissions", "wind forcing", "thermodynamic f."]

    # load up the input data
    input_data = [xr.open_dataset(filepath, decode_times = False) for filepath in filepaths]
    temp_data = xr.open_dataset(f"{output_path}average_CTRL_1920-1950.nc")
    
    sliced_data = temp_data["maskC"].sel(XC = lon_slices[3], method = "nearest")
    mask = sliced_data.sel(YC = slice(lat_slices[0], lat_slices[1]))

    # read in the general variables (these should be the same between the ensembles
    #[lat, lon, ice_mask_temp, depth] = [input_data[0][param].values for param in ["YC", "XC", "maskC", "Depth"]]
    #ice_mask = ice_mask_temp[0,:,:]

    grid = Grid(grid_filepath)

    z = input_data[0].depth.values
    lat = input_data[0]["lat"]
    color_scheme = "PRGn"
       
    # load the variables and calculate the average over the year
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    data_um = [input.trend.values[:,:,3] for input in input_data] 
    data = [np.ma.masked_where(mask == 0, data_idx) for data_idx in data_um]  
    print(np.shape(data[0]), np.shape(mask))

    data_um = [input.pvalue.values[:,:,3] for input in input_data] 
    pvalue = [np.ma.masked_where(mask == 0, data_idx) for data_idx in data_um]  
      
    low_val = -5e-06
    high_val = 5e-06
    print(low_val, high_val)
    step = 15
    font_size = 16
        
    # create subplots with each variable on a new line
    fig, axs = plt.subplots(nrows=1, ncols=4, gridspec_kw={"hspace": 0.05, "wspace": 0.05}, figsize=(25,8))
        
    axs = axs.flatten()

    for i in range(4):
        cs = axs[i].contourf(lat, z, data[i], cmap=color_scheme, extend="both", levels=np.linspace(low_val, high_val, step))
        axs[i].contourf(lat, z, pvalue[i], levels=[-np.inf, 0.01], colors='none', hatches=['..'], alpha=0)
        axs[i].set_title(experiment[i], weight="bold", fontsize = font_size)
        axs[i].set_xlabel("Latitude", fontsize = font_size)
        axs[i].set_ylim([-1000, 0])
        if i > 0:
            axs[i].get_yaxis().set_visible(False)
        axs[i].set_aspect("auto", adjustable="box")
    axs[0].set_ylabel("Depth (m)", fontsize = font_size)   
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(cs, cax=cbar_ax, ticks = np.arange(low_val, high_val, 2.5e-6))
        
    plt.suptitle("Salinity trend per cetury across 118W", fontsize = font_size, weight = "bold")
        
    fig.savefig("salt_trend_undercurrent_4.png")
    plt.show()
        
    
if __name__ == '__main__':
    main() # run the program
