import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt

import numpy as np
import xarray as xr
import sys

from directories_and_paths import *
from funcs import moving_average
from scipy.ndimage.filters import gaussian_filter1d

def four_comparison(var, unit, var_name, title, exp, maxmin = False):
    simulation_names = ["ctrl", "lens", "wind", "temp"]
    data = []
    var_smoothed = []
    
    for simulation in simulation_names:
        file_path = ensemble_mean_path + simulation + ensemble_mean_file
        data_ensemble = xr.open_dataset(file_path)
        var_ensemble = gaussian_filter1d(data_ensemble[var].values, sigma = 5)
        data.append(data_ensemble)
        var_smoothed.append(var_ensemble)
    
    time = range(len(data[0].time.values))
    
    fig, ax = plt.subplots(figsize=(12, 5))
    plt.rcParams.update({'font.size': 14})
    
    colors = ['seagreen', 'orchid', 'dodgerblue', 'orange']
    
    for i, var_ensemble in enumerate(var_smoothed):
        if i>1:
            ax.plot(time, var_ensemble, '--', color=colors[i], label=var_name[i])
        else:
            ax.plot(time, var_ensemble, color=colors[i], label=var_name[i])
        #ax.plot(data[i].theta_timeseries.values, color=colors[i], label=var_name[i])
        ax.fill_between(time, data[i].theta_min.values, data[i].theta_max.values, color=colors[i], alpha = 0.5)
    
    ax.set_ylabel("Temperature (°C)", fontsize = 14)
    ax.legend(fancybox=True)
    ax.set_xticks(np.arange(0, len(data[0].time.values), 12))
    ax.set_xticklabels(np.arange(0, 31))
    ax.set_xlim([0, len(data[0].time.values)])
    ax.grid(alpha=0.8)
    ax.set_title(title)
    
    file_out = "poster_timeseries.png"
    fig.savefig(file_out)
    plt.show()

def comparison(timeseries_data, var, var_name, title, file_out, smooth = False, ylabel = None, xticks = "annual"):
    data_to_plot = []
    
    if smooth:
        timeseries = timeseries_data[0]
        var_smoothed = gaussian_filter1d(timeseries, sigma = 5)
        data_to_plot.append(var_smoothed)
        timeseries = timeseries_data[1]
        var_smoothed = gaussian_filter1d(timeseries, sigma = 5)
        data_to_plot.append(var_smoothed)
    else:
        data_to_plot = timeseries_data
    
    time = range(len(data_to_plot[0]))
    
    fig, ax = plt.subplots(figsize=(12, 5))
    plt.rcParams.update({'font.size': 14})
    
    colors = ['seagreen', 'orchid', 'dodgerblue', 'orange']
    
    for i, var_ensemble in enumerate(data_to_plot):
        ax.plot(time, var_ensemble, color=colors[i], label=var_name[i])
        print(var_name[i], colors[i])
        

    ax.set_ylabel(None, fontsize = 14)
    ax.legend(fancybox=True)
    if xticks == "annual":
        ax.set_xticks(np.arange(0, len(data_to_plot[0]), 12))
        ax.set_xticklabels(np.arange(31))
        ax.set_xlim([0, len(data_to_plot[0])])
    ax.grid(alpha=0.8)
    ax.set_title(title)
    
    fig.savefig(file_out)
    plt.show()
