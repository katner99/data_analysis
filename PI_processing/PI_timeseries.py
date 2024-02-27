import matplotlib

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

import numpy as np
import xarray as xr
from funcs import moving_average
import os
import sys
import scipy.stats

from directories_and_paths import *


def plot_timeseries(
    ax,
    time,
    var_ensemble,
    color,
    linestyle,
    plot_range=False,
    min_values=None,
    max_values=None,
    alpha=0.5,
    label=None,
):
    line = ax.plot(var_ensemble, linestyle, color=color, alpha=alpha, label=label)
    if plot_range:
        ax.fill_between(time, min_values, max_values, color=color, alpha=0.5)
    if label:
        line[0].set_label(label)  # Set label to the line, not to the ax

    return line  # Return the Line2D object


def plot_comparison(
    var,
    data,
    experiments,
    ensemble_members,
    ylabel,
    xlabel,
    time,
    file_out,
    smooth=0,
    linearity=False,
    warming=False,
    percentage=False,
):
    fig, ax = plt.subplots(figsize=(12, 5))
    colors = ["dodgerblue", "purple", "forestgreen", "orangered"]
    experiment_info = {}  # Dictionary to store experiment information
    all_means = []

    experiment_full = [
        "wind forcing",
        "thermodynamic f.",
        "pre-industrial control",
        "RCP 8.5",
    ]
    
    if warming:
        experiment_mean = np.nanmean([data["CTRL"][ens][ts_idx][:] for ens in ensemble_members for ts_idx in range(len(data["CTRL"][ens]))], axis=0)
        warming_diff = - np.nanmean(experiment_mean) 
    
    for exp_idx, exp in enumerate(experiments):

        experiment_mean = np.nanmean([data[exp][ens][ts_idx][:] for ens in ensemble_members for ts_idx in range(len(data[exp][ens]))], axis=0)
        experiment_max = np.nanmax([data[exp][ens][ts_idx][:] for ens in ensemble_members for ts_idx in range(len(data[exp][ens]))], axis=0)
        experiment_min = np.nanmin([data[exp][ens][ts_idx][:] for ens in ensemble_members for ts_idx in range(len(data[exp][ens]))], axis=0)
        
        if smooth > 0:
            smoothed_mean = moving_average(experiment_mean, smooth * 12)
            smoothed_max = moving_average(experiment_max, smooth * 12)
            smoothed_min = moving_average(experiment_min, smooth * 12)
        else:
            smoothed_mean = experiment_mean
            smoothed_max = experiment_max
            smoothed_min = experiment_min
        
        if warming == True:
            smoothed_mean = smoothed_mean + warming_diff
            smoothed_max = smoothed_max + warming_diff
            smoothed_min = smoothed_min + warming_diff
            print(warming_diff, len(smoothed_mean))

        x_values = time
        ax.plot(
            smoothed_mean,
            color=colors[exp_idx],
            label=experiment_full[exp_idx],
            alpha=1,
        )
        ax.fill_between(
            range(len(smoothed_min)),
            smoothed_min,
            smoothed_max,
            color=colors[exp_idx],
            alpha=0.3,
        )
        all_means.append(experiment_mean[:])

    if linearity == True:
        linearity = moving_average(
           all_means[3] + all_means[2] - all_means[0] - all_means[1], 12 * smooth
        )
        ax.plot(linearity, color="black", label="linearity")
        ax.axhline(0, color='grey')

    elif percentage == True:
        linearity = moving_average(all_means[3] + all_means[2] - all_means[0] - all_means[1], 12 * smooth)
        lens_warming = moving_average(all_means[3], 12 * smooth)
        ctrl_warming = moving_average(all_means[2], 12 * smooth)
        print(np.shape(linearity), np.shape(lens_warming), np.shape(ctrl_warming))
        perc_warming = (linearity / (lens_warming - ctrl_warming)) * 100
        ax2 = ax.twinx()
        ax2.set_ylabel("percentage of warming (linearity/total warming)")
        ax2.plot(perc_warming, color="black", label="linearity")
        ax2.set_ylim([-100, 100])
        ax2.axhline(0, color='grey')
    
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_xticks(np.arange(6*smooth, x_values-(6*smooth), 120))
    ax.set_xticklabels(xlabel.astype(str), rotation=45)
    ax.set_xlim([6*smooth, x_values-(6*smooth)])
    ax.set_title(f"{var} between 200 - 700m")
    if percentage:
        ax.legend()
        ax2.grid(alpha=0.8)
    else:
        if warming:
            ax.set_ylim([-1.6, 1.6])
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        ax.grid(alpha=0.8)

    fig.savefig(file_out, transparent=True)
    plt.show()


def plot_ensemble(
    var,
    data,
    experiments,
    ensemble_members,
    ylabel,
    xlabel,
    time,
    title,
    file_out,
    smooth=False,
    plot_range=False,
    min_values=None,
    max_values=None,
    linearity=False,
    percentage=False,
):

    colors = [
        "deepskyblue",
        "hotpink",
        "forestgreen",
        "orangered",
        "goldenrod",
        "purple",
        "blue",
        "dodgerblue",
        "red",
    ]

    fig, ax = plt.subplots(figsize=(12, 5))
    plt.rcParams.update({"font.size": 14})

    # run through all the experiments
    for ens in ensemble_members:
        print(ens)
        timeseries = data[experiments][ens]
        # print(max(timeseries))
        var_ensemble = np.transpose(timeseries)
        print(len(var_ensemble))
        plot_timeseries(
            ax,
            time,
            var_ensemble,
            colors[ens - 1],
            "-",
            plot_range,
            min_values,
            max_values,
            alpha=0.3,
            label="ens_" + str(ens),
        )
        x_values = time

    exp = experiments
    # shortest_length = np.min([len(data[exp][ens][ts_idx]) for ens in ensemble_members for ts_idx in range(len(data[exp][ens]))])
    # experiment_mean = np.nanmean([data[exp][ens][ts_idx][:] for ens in ensemble_members for ts_idx in range(len(data[exp][ens]))], axis=0)
    # smoothed_mean = moving_average(experiment_mean, 2*12)
    # print(type(smoothed_mean))
    x_values = range(2172)
    # ax.plot(smoothed_mean, color='k', label="mean", alpha = 1)

    # res = scipy.stats.linregress(range(len(smoothed_mean[2*12:])), smoothed_mean[2*12:])
    # print(res.slope, res.pvalue)
    # ax.plot(res.intercept + res.slope*smoothed_mean,'r', label='fitted line')
    # ax.legend()
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_xticks(np.arange(0, len(x_values), 120))
    ax.set_xticklabels(xlabel.astype(str), rotation=45)
    ax.set_xlim([0, len(x_values)])
    ax.grid(alpha=0.8)
    # ax.set_title("slope: "+str(np.round(res.slope,2))+" p-value: "+str( np.round(res.pvalue)))
    ax.set_title(title)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

    #fig.savefig(file_out)
    plt.show()

def comparison():
    # set variables for the function
    var = "transport_south"
    if var in ["transport", "transport_south"]:
        origin = "transport"
    else:
        origin = "timeseries"

    experiments = ["WIND", "TEMP", "CTRL", "LENS"]
    ensemble = [1, 2, 3, 4, 5, 6, 7, 8, 9]

    # Define the output structure
    data = {exp: {ens: [] for ens in ensemble} for exp in experiments}
    max_len = 181 * 12

    for exp in experiments:
        timeseries = []
        for ens in ensemble:
            if exp == "LENS":
                if ens < 6:
                    filepath = (
                        f"{output_path}{origin}2101_experiment" + str(ens) + ".nc"
                    )
                else:
                    path = output_path + "PAS_LENS00" + str(ens) + "_noOBC/"
                    valid_file = [
                        filename
                        for filename in os.listdir(path)
                        if filename.startswith(origin)
                    ]
                    filepath = path + valid_file[0]
            else:
                valid_file = []
                path = output_path + exp + "_ens0" + str(ens) + "_noOBC/"
                valid_file = [
                    filename
                    for filename in os.listdir(path)
                    if filename.startswith(origin)
                ]
                filepath = path + valid_file[0]

            data_member = xr.open_dataset(filepath)
            timeseries = data_member[var].values
            if len(timeseries) < max_len:
                timeseries = np.pad(
                    timeseries, (0, max_len - len(timeseries)), constant_values=np.nan
                )
            data[exp][ens].append(timeseries)

    if var == "sea_ice":
        title = exp + " max " + var
        ylabel = "Sea ice thickness (m)"
        smooth = 0
        linearity = False
    elif var == "salt":
        title = var
        ylabel = "salinity (psu)"
        smooth = 30
        linearity = False
    elif origin == "transport":
        title = var
        ylabel = "Transport (Sv)"
        smooth = 2
        linearity = True
        percentage = False
        warming = False
    else:
        title = var
        ylabel = "Temperature (C)"
        smooth = 2
        linearity = True
        percentage = False
        warming = True

    start_year = 1920 + (smooth/2)
    end_year = 2101 - (smooth/2)
    xlabel = np.arange(start_year, end_year, 10)
    time = 12 * 181

    if warming:
        file_out = "timeseries_warming_" + var + "_"+str(smooth)+".png"
    elif percentage:
        file_out = "timeseries_percentage_" + var + "_"+str(smooth)+".png"
    else:
        file_out = "timeseries_" + exp + "_" + var + ".png"

    # function to plot the timeseries
    plot_comparison(
        var,
        data,
        experiments,
        ensemble,
        ylabel,
        xlabel,
        time,
        file_out,
        smooth=smooth,
        linearity=linearity,
        warming=warming,
        percentage=percentage,
    )

    # plot_ensemble(var, data, exp, ensemble, ylabel, xlabel, time, title, file_out, smooth=smooth, linearity = linearity, percentage = percentage)

def one_on_one():
    filepath_OBC = f"{output_path}TEMP_ens01_OBC/"
    data_OBC = xr.open_dataset(f"{filepath_OBC}timeseries2006.nc", decode_times = False)

    filepath_noOBC = f"{output_path}TEMP_ens01_noOBC/"
    data_noOBC = xr.open_dataset(f"{filepath_noOBC}timeseries2101.nc", decode_times = False)

    variables = ["theta", "theta_pig", "theta_abbot", "theta_dotson", "theta_shelf_edge", "salt"]

    fig, axs = plt.subplots(nrows=6, ncols=1, gridspec_kw={"hspace": 0.5, "wspace": 0.4}, figsize=(10, 20))

    for var_idx, var in enumerate(variables):
        axs[var_idx].plot(data_OBC[var].values, color = "blue", label = "transient")
        axs[var_idx].plot(data_noOBC[var].values, color = "red", label = "climatology")
        axs[var_idx].set_title(f"{var}")
        axs[var_idx].legend()

    fig.savefig("initial_comparison_timeseries.png")
    plt.show()

if __name__ == "__main__":
    one_on_one()
