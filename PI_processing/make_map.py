import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.patches import Polygon

import numpy as np
import xarray as xr
import os

from tools.funcs import read_variable, find_nearest
from tools.plots import create_mask
from tools.directories_and_paths import *
from mitgcm_python.grid import Grid
#from PI_timeseries import plot_comparison

def load_data(experiments, ensemble, var):
    data_load = {exp: {ens: [] for ens in ensemble} for exp in experiments}
    for exp in experiments:
        timeseries = []
        for ens in ensemble:
            # read the filename of the experiment
            if exp == "LENS":
                if ens < 6:
                    filepath = output_path+"timeseries2101_experiment"+str(ens)+".nc"
                else:
                    path = output_path+"PAS_LENS00" + str(ens) + "_noOBC/"
                    valid_file = [filename for filename in os.listdir(path) if filename.startswith("timeseries")]
                    filepath = path+valid_file[0]
            else:
                valid_file = []
                path = output_path+exp+"_ens0"+str(ens)+"_noOBC/"
                valid_file = [filename for filename in os.listdir(path) if filename.startswith("timeseries")]
                filepath = path+valid_file[0]

            data_member = xr.open_dataset(filepath)
            timeseries = data_member[var].values
            data_load[exp][ens].append(timeseries)
    return data_load


def main():

    #variables = ["theta", "theta_c", "theta_s", "theta_a"]
    #experiments = ["WIND", "TEMP", "CTRL", "LENS"]
    #experiments = ["CTRL"]
    #ensemble = [1,2,3,4,5]
    
    # Define the output structure
    #data = {exp: {ens: [] for ens in ensemble} for exp in experiments}
    #data_c = {exp: {ens: [] for ens in ensemble} for exp in experiments}
    #data_s = {exp: {ens: [] for ens in ensemble} for exp in experiments}
    #data_a = {exp: {ens: [] for ens in ensemble} for exp in experiments}

    #data = load_data(experiments, ensemble, variables[0])
    #data_c = load_data(experiments, ensemble, variables[1])
    #data_s = load_data(experiments, ensemble, variables[2])
    #data_a = load_data(experiments, ensemble, variables[3])

    filepath = output_path + "average_WIND_2070-2100.nc"
    input_data = xr.open_dataset(filepath, decode_times=False)

    # read in the general variables (these should be the same between the ensembles
    [lat, lon, ice_mask_temp, depth] = [input_data[param].values for param in ["YC", "XC", "maskC", "Depth"]]
    ice_mask = ice_mask_temp[0,:,:]
    grid = Grid(grid_filepath)
    [land_mask, mask, colors] = create_mask(depth, ice_mask)

    # set up the grid
    [X, Y] = np.meshgrid(lon, lat)

    # pine island coastal box;
    lon_box_1 = [find_nearest(lon, 95+180), find_nearest(lon, 115+180)]
    lat_box_1 = [find_nearest(lat, -75), find_nearest(lat, -73)]
    #box = Polygon([[250, -75], [260, -75], [260, -71.5], [250, -71.5]], fill = False, closed = True, edgecolor = "black", linewidth = 3, label = "timeseries average")
    #box = 
    # box_1 = Polygon([[250, -75], [260, -75], [260, -71.5], [250, -71.5]], fill = False, closed = True, edgecolor = "magenta", linewidth = 3, label = "profile average")
    #box_2 = Polygon([[252.8, -73], [255, -73], [255, -73], [252.8, -73]], fill = False, closed = True, edgecolor = "mediumorchid", linewidth = 3, label = "trough")
    #box_3 = Polygon([[250, -72], [250, -72], [250, -70], [250, -70]], fill = False, closed = True, edgecolor = "hotpink", linewidth = 3, label = "undercurrent 2")
    #box_4 = Polygon([[243, -72], [243, -72], [243, -70.5], [243, -70.5]], fill = False, closed = True, edgecolor = "crimson", linewidth = 3, label = "undercurrent 3")
    #box_5 = Polygon([[255, -71.5], [255, -71.5], [255, -70], [255, -70]], fill = False, closed = True, edgecolor = "purple", linewidth = 3, label = "undercurrent 4")

    ylabel = "Temperature (°C)"
    xlabel = np.linspace(1920, 2100,21, dtype = int)
    time = 2172
    smooth = True
    linearity = True
    percentage = False
    #print(lat_box_1, lon_box_1)

    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(10, 7))
    
    #ax = ax.flatten()

    cs = ax.contourf(X, Y, depth, extend="max", cmap="bone_r", levels=np.linspace(0, 4000, 20))
    ax.contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
    ax.contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
    ax.vlines(x=238, ymin=-74, ymax=-71.5, color='magenta', linewidth = 3) # if plotting the velocity profiles, indicate where on the map
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    #ax.add_patch(box_1)
    #ax.add_patch(box_2)
    #ax.add_patch(box_3)
    #ax.add_patch(box_4)
    #ax.add_patch(box_5)
    #ax.legend()
    #ticks=np.arange(0, 4000, 500)
    #cbar = fig.colorbar(cs)
    #cbar.set_label('Depth (m)', rotation=270, labelpad=15)
    #cbar.set_ticks(ticks)

    #plot_comparison(ax[1], variables[0], data, experiments, ensemble, ylabel, xlabel, time, title = "general timeseries", file_out = None, smooth = smooth, linearity = linearity)
    #plot_comparison(ax[2], variables[1], data_c, experiments, ensemble, ylabel, xlabel, time, title = "coastal timeseries", file_out = None, smooth = smooth, linearity = linearity)
    #plot_comparison(ax[3], variables[2], data_s, experiments, ensemble, ylabel, xlabel, time, title = "shelf_break timeseries", file_out = None, smooth = smooth, linearity = linearity)
    #plot_comparison(ax[4], variables[3], data_a, experiments, ensemble, ylabel, xlabel, time, title = "Abbot timeseries", file_out = None, smooth = smooth, linearity = linearity)
    fig.savefig("amundsen_map_profile.png", transparent = True)
    plt.show()
            
if __name__ == '__main__':
    main() # run the program
