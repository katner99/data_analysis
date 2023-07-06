import matplotlib
matplotlib.use('TkAgg')

import numpy as np
import xarray as xr
from scipy.ndimage.filters import gaussian_filter1d

from directories_and_paths import *
from funcs import find_nearest, make_timeseries
from mitgcm_python.grid import Grid


def plot_timeseries(ax, time, var_ensemble, color, label, linestyle, plot_range = False, min_values = None, max_values = None):
    ax.plot(time, var_ensemble, linestyle, color=color, label=label)
    if plot_range:
        ax.fill_between(time, min_values, max_values, color=color, alpha=0.5)

def plot_comparison(timeseries_data, var, ylabel, xlabel, x_values, var_name, title, file_out, smooth = False, plot_range = False, min_values = None, max_values = None):
    var_ensemble = []
    linestyle = '-'
    
    for ensemble_member in timeseries_data:
        if smooth:
            var_smooth = gaussian_filter1d(ensemble_member, sigma=5)
            var_ensemble.append(var_smooth)
        else:
            var_ensemble.append(ensemble_member)
    
    fig, ax = plt.subplots(figsize=(12, 5))
    plt.rcParams.update({'font.size': 14})
    
    colors = ['seagreen', 'orchid', 'dodgerblue', 'orange']
    
    for i, var_member in enumerate(var_ensemble):
        plot_timeseries(ax, x_values, var_member, colors[i], var_name[i],
                        linestyle, plot_range, min_values, max_values)
    
    ax.set_ylabel(ylabel, fontsize=14)
    ax.legend(fancybox=True)
    step = int(len(x_values)/len(xlabel))
    ax.set_xticks(np.arange(0, len(x_values), step))
    ax.set_xticklabels(xlabel.astype(str), rotation=45)
    ax.set_xlim([0, len(x_values)])
    ax.grid(alpha=0.8)
    ax.set_title(title)
    
    fig.savefig(file_out)
    plt.show()


def main():

    var = "SIheff"
    var_name = ["wind forcing ensemble mean", "wind forcing new compiler"]
    title = "wind forcing"
    filepaths = [ensemble_mean_path + "wind_ensemble_mean.nc",
                 output_path + "PAS_wind09_slice.nc"]

    # use below for the comparison of four
    # var_name = ["pre-industrial forcing", "RCP 8.5 forcing", "wind forcing", "thermodynamic forcing"]
    # title = "Continental shelf mean temperature between 200m and 700m"
    # exp = "_ensemble_mean"
    # filepaths = [ensemble_mean_path + "ctrl" + ensemble_mean_file,
    #              ensemble_mean_path + "lens" + ensemble_mean_file,
    #              ensemble_mean_path + "wind" + ensemble_mean_file,
    #              ensemble_mean_path + "temp" + ensemble_mean_file]

    ylabel = "Temperature (°C)"
    xlabel = np.arange(2070, 2100)
    data = []
                  
    data_1 = xr.open_dataset(filepaths[0])
    data_2 = xr.open_dataset(filepaths[1])
    
    # set up grid
    grid = Grid(grid_filepath)
    grid_file = xr.open_dataset(grid_filepath)
    
    # set lat, lon, and depth range
    depth_range = [find_nearest(grid_file.Z.values, -200), find_nearest(grid_file.Z.values, -700)]
    lon_range = [find_nearest(grid_file.XC.values, 250), find_nearest(grid_file.XC.values, 260)]
    lat_range = [find_nearest(grid_file.YC.values, -76), find_nearest(grid_file.YC.values, -72)]

    time = data_1.time.values

    timeseries_1 = make_timeseries(var, data_1, grid, lat_range, lon_range, depth_range, time = len(time))
    timeseries_2 = make_timeseries(var, data_2, grid, lat_range, lon_range, depth_range, time = len(time))

    data = [timeseries_1, timeseries_2]    

    file_out = "comparison_loc_timeseries.png"

    plot_comparison(data, var, ylabel, xlabel, time, var_name, title, file_out)
    
if __name__ == "__main__":
    main()
    
    
