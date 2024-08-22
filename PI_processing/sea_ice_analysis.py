from config_options import config_comparison
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import numpy as np
from plots import zoom_shelf, pretty_labels
from plots import read_mask
from plots_2d import contour_func
from directories_and_paths import output_path, grid_filepath
from mitgcm_python.grid import Grid

import sys
import xarray as xr

def read_data(var, var_name = None):
    """ read in analysed data """
    if var in ["trend", "pvalue"]:
        filename = f"{var_name}_spread.nc"
        file_out = f"mega_comparison_{var}_{var_name}.png"
    elif var == "mean":
        filename = f"{var_name}.nc"
        file_out = f"mega_comparison_{var}_{var_name}_mean.png"
    filepaths = [
        f"{output_path}{exp}_files_temp/{filename}"
        for exp in ["CTRL", "LENS", "TEMP", "WIND"]
    ]
    return filepaths

def main():
    var = "mean"             # either set to the parameter you will be using or set to mean or trend (tells it where to read the values in)
    save = True              # save figure
    show = False             # print the figure on script completion
    year = 2085

    filepaths_SIheff = read_data("mean", "SIheff")
    filepaths_SIarea = read_data("mean", "SIarea")
    filepaths_oceFWflx = read_data("mean", "oceFWflx")
    filepaths_SIfwmelt = read_data("mean", "SIfwmelt")
    filepaths_SIfwfrz = read_data("mean", "SIfwfrz")

    file_out = f"sea_ice_analysis.png"

    # experiment neames to add to the plot
    experiment = ["NONE", "ALL", "THERMO", "WIND"]

    input_data_SIheff = [
        xr.open_dataset(filepath, decode_times=False) for filepath in filepaths_SIheff
    ]
    input_data_SIarea = [
        xr.open_dataset(filepath, decode_times=False) for filepath in filepaths_SIarea
    ]
    input_data_oceFWflx = [
        xr.open_dataset(filepath, decode_times=False) for filepath in filepaths_oceFWflx
    ]
    input_data_SIfwmelt = [
        xr.open_dataset(filepath, decode_times=False) for filepath in filepaths_SIfwmelt
    ]
    input_data_SIfwfrz = [
       xr.open_dataset(filepath, decode_times=False) for filepath in filepaths_SIfwfrz
    ]

    grid = Grid(grid_filepath)

    set_up_data = xr.open_dataset(f"{output_path}average_CTRL_1920-1950.nc", decode_times=False)
    set_up = read_mask(set_up_data)
    [data_SIheff, graph_params, graph_params_seaice] = config_comparison(var, input_data_SIheff, grid, var_name = "SIheff", year = year)
    [data_SIarea, graph_params, graph_params_seaice] = config_comparison(var, input_data_SIarea, grid, var_name = "SIarea", year = year)
    [data_oceFWflx, graph_params, graph_params_fw] = config_comparison(var, input_data_oceFWflx, grid, var_name = "oceFWflx", year = year)
    [data_SIfwmelt, graph_params, graph_params_fw] = config_comparison(var, input_data_SIfwmelt, grid, var_name = "SIfwmelt", year = year)
    [data_SIfwfrz, graph_params, graph_params_fw] = config_comparison(var, input_data_SIfwfrz, grid, var_name = "SIfwfrz", year = year)
    
    fig, axs = plt.subplots(nrows=3, ncols=5, gridspec_kw={"hspace": 0.05, "wspace": 0.04}, figsize=(18, 10))
    
    titles = ["Sea-ice thickness (m)", "Sea-ice concentration (0-1)", "total freshwater flux \n(m yr$^{-1}$)", "Fresh water flux from \n melting (m yr$^{-1}$)", "Fresh water flux from \n freezing (m yr$^{-1}$)"]
    for i in range(3):
        for j, data in enumerate([data_SIheff, data_SIarea, data_oceFWflx, data_SIfwmelt, data_SIfwfrz]):
            if j == 0:
                axs[i, j].set_ylabel(f"{experiment[i+1]}", weight = "bold")
                hide_ticks_y = False
            else:
                hide_ticks_y = True

            if i == 0:
                axs[i,j].set_title(titles[j])
                hide_ticks_x = True
            elif i == 2:
                hide_ticks_x = False
            else:
                hide_ticks_x = True

            if j < 2:
                cs_seaice = contour_func(axs[i][j], data[i+1] - data[0], set_up, graph_params_seaice, hide_ticks_x, hide_ticks_y)
            else:
                cs_fw = contour_func(axs[i][j], data[i+1] - data[0], set_up, graph_params_fw, hide_ticks_x, hide_ticks_y)

            zoom_shelf(axs[i][j], "cont_shelf")
            pretty_labels(axs[i][j])
    ticks=np.arange(graph_params_seaice["low_val"], graph_params_seaice["high_val"]+0.1, graph_params_seaice["interval"])
    cbar_ax = fig.add_axes([0.13, 0.02, 0.3, 0.02]) 
    cbar = plt.colorbar(cs_seaice, cax=cbar_ax, orientation='horizontal')
    cbar.set_ticks(ticks)
    #cbar.ax.yaxis.set_ticks_position('left')
    cbar.set_label

    ticks=np.arange(graph_params_fw["low_val"], graph_params_fw["high_val"]+0.1, graph_params_fw["interval"])  
    cbar_ax = fig.add_axes([0.44, 0.02, 0.48, 0.02]) 
    cbar = plt.colorbar(cs_fw, cax=cbar_ax, orientation='horizontal')
    cbar.set_ticks(ticks)
    #cbar.ax.yaxis.set_ticks_position('left')

    fig.suptitle(f"Anomaly with respect to NONE {year} - {year + 10}", fontsize=16)

    # save figure
    if save == True:
        fig.savefig(file_out, bbox_inches='tight')

    # show figure
    if show == True:
        plt.show()

def sea_ice_total():
    var = "mean"             # either set to the parameter you will be using or set to mean or trend (tells it where to read the values in)
    save = True              # save figure
    show = False             # print the figure on script completion
    year = 2090

    filepaths_SIheff = read_data("mean", "SIheff")
    filepaths_SIarea = read_data("mean", "SIarea")
    filepaths_oceFWflx = read_data("mean", "oceFWflx")
    filepaths_SIfwmelt = read_data("mean", "SIfwmelt")
    filepaths_SIfwfrz = read_data("mean", "SIfwfrz")

    file_out = f"sea_ice_analysis_total.png"

    # experiment neames to add to the plot
    experiment = ["NONE", "ALL", "THERMO", "WIND"]

    input_data_SIheff = [
        xr.open_dataset(filepath, decode_times=False) for filepath in filepaths_SIheff
    ]
    input_data_SIarea = [
        xr.open_dataset(filepath, decode_times=False) for filepath in filepaths_SIarea
    ]
    input_data_oceFWflx = [
        xr.open_dataset(filepath, decode_times=False) for filepath in filepaths_oceFWflx
    ]
    input_data_SIfwmelt = [
        xr.open_dataset(filepath, decode_times=False) for filepath in filepaths_SIfwmelt
    ]
    input_data_SIfwfrz= [
       xr.open_dataset(filepath, decode_times=False) for filepath in filepaths_SIfwfrz
    ]

    grid = Grid(grid_filepath)

    set_up_data = xr.open_dataset(f"{output_path}average_CTRL_1920-1950.nc", decode_times=False)
    set_up = read_mask(set_up_data)
    [data_SIheff, graph_params_seaice, graph_params] = config_comparison(var, input_data_SIheff, grid, var_name = "SIheff", year = year)
    [data_SIarea, graph_params_seaice, graph_params] = config_comparison(var, input_data_SIarea, grid, var_name = "SIarea", year = year)
    [data_oceFWflx, graph_params_fw, graph_params] = config_comparison(var, input_data_oceFWflx, grid, var_name = "oceFWflx", year = year)
    [data_SIfwmelt, graph_params_fw, graph_params] = config_comparison(var, input_data_SIfwmelt, grid, var_name = "SIfwmelt", year = year)
    [data_SIfwfrz, graph_params_fw, graph_params] = config_comparison(var, input_data_SIfwfrz, grid, var_name = "SIfwfrz", year = year)
    
    fig, axs = plt.subplots(nrows=4, ncols=5, gridspec_kw={"hspace": 0.05, "wspace": 0.04}, figsize=(18, 10))
    
    titles = ["Sea-ice thickness (m)", "Sea-ice concentration (0-1)", "total freshwater flux \n(m yr$^{-1}$)", "Fresh water flux from \n melting(m yr$^{-1}$)",  "Fresh water flux from \n freezing (m yr$^{-1}$)"]
    for i in range(4):
        for j, data in enumerate([data_SIheff, data_SIarea, data_oceFWflx, data_SIfwmelt, data_SIfwfrz]):
            if j == 0:
                axs[i, j].set_ylabel(f"{experiment[i]}", weight = "bold")
                hide_ticks_y = False
            else:
                hide_ticks_y = True

            if i == 0:
                axs[i,j].set_title(titles[j])
                hide_ticks_x = True
            elif i == 3:
                hide_ticks_x = False
            else:
                hide_ticks_x = True

            if j < 2:
                cs_seaice = contour_func(axs[i][j], data[i], set_up, graph_params_seaice, hide_ticks_x, hide_ticks_y)
            else:
                cs_fw = contour_func(axs[i][j], data[i], set_up, graph_params_fw, hide_ticks_x, hide_ticks_y)

            zoom_shelf(axs[i][j], "cont_shelf")
            pretty_labels(axs[i][j])
    
    ticks=np.arange(graph_params_seaice["low_val"], graph_params_seaice["high_val"]+0.1, graph_params_seaice["interval"])
    cbar_ax = fig.add_axes([0.13, 0.02, 0.3, 0.02]) 
    cbar = plt.colorbar(cs_seaice, cax=cbar_ax, orientation='horizontal')
    cbar.set_ticks(ticks)
    #cbar.ax.yaxis.set_ticks_position('left')
    cbar.set_label

    ticks=np.arange(graph_params_fw["low_val"], graph_params_fw["high_val"]+0.1, graph_params_fw["interval"])  
    cbar_ax = fig.add_axes([0.44, 0.02, 0.48, 0.02]) 
    cbar = plt.colorbar(cs_fw, cax=cbar_ax, orientation='horizontal')
    cbar.set_ticks(ticks)
    #cbar.ax.yaxis.set_ticks_position('left')

    #ticks=np.arange(graph_params_theta["low_val"], graph_params_theta["high_val"]+0.1, graph_params_theta["interval"])
    #cbar_ax = fig.add_axes([0.755, 0.02, 0.135, 0.02])
    #cbar = plt.colorbar(cs_theta, cax=cbar_ax, orientation='horizontal')
    #cbar.set_ticks(ticks)
    #cbar.ax.yaxis.set_ticks_position('left')

    fig.suptitle(f"Sea ice producation and freshwater fluxes between {year} - {year + 10}", fontsize=16)

    # save figure
    if save == True:
        fig.savefig(file_out, bbox_inches='tight')

    # show figure
    if show == True:
        plt.show()

if __name__ == "__main__":
    main() 
