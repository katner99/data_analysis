import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import numpy as np
import xarray as xr
from funcs import moving_average
import os
import sys

from directories_and_paths import *

def plot_timeseries(ax, time, var_ensemble, color, linestyle, plot_range = False, min_values = None, max_values = None, alpha = 0.5, label = None):
    line = ax.plot(var_ensemble, linestyle, color=color, alpha=alpha, label = label)
    if plot_range:
        ax.fill_between(time, min_values, max_values, color=color, alpha=0.5)
    if label:
        line[0].set_label(label)  # Set label to the line, not to the ax
    
    return line  # Return the Line2D object


def plot_comparison(var, data, experiments, ensemble_members, ylabel, xlabel, time, title, file_out, smooth = False, plot_range = False, min_values = None, max_values = None, linearity = False):
    
    #linestyle = ['-', '--', '-.', '-', '--', '-.', '-', '--', '-.', '-', '--', '-.']
    colors = ["lightblue", "sandybrown", "lightgreen", "lightpink"]
    dark_colors = ["deepskyblue", "darkgoldenrod", "forestgreen", "orchid"]

    fig, ax = plt.subplots(figsize=(12, 5))
    plt.rcParams.update({'font.size': 14})
    experiment_info = {}  # Dictionary to store experiment information
    means = []

    for exp_idx, exp in enumerate(experiments):
        # Plot individual ensemble members
        shortest_length = np.min([len(data[exp][ens][ts_idx]) for ens in ensemble_members for ts_idx in range(len(data[exp][ens]))])
        print(exp, shortest_length)
        for ens in ensemble_members:
            timeseries=data[exp][ens]
            if smooth:
                var_smooth = moving_average(timeseries, 30*12)
                var_ensemble = var_smooth
            else:
                var_ensemble = np.transpose(timeseries)

            if exp not in experiment_info:
                experiment_info[exp] = {
                    "color": colors[exp_idx],
                    "label": f"{exp}",
                    "lines": [],  # Store Line2D objects for each experiment
                }

            line = plot_timeseries(ax, time, var_ensemble, experiment_info[exp]["color"], "-", plot_range, min_values, max_values, alpha = 0.5,label="_" + experiment_info[exp]["label"])
            experiment_info[exp]["lines"].append(line[0]) 
        # Calculate and plot the experiment mean 
        if var == "theta":
            experiment_mean = moving_average(np.mean([data[exp][ens][ts_idx][:shortest_length] for ens in ensemble_members for ts_idx in range(len(data[exp][ens]))], axis=0), 30*12)
            x_values = len(experiment_mean)
            ax.plot(experiment_mean, color=dark_colors[exp_idx], label=f"{exp} mean", alpha = 1)
            #if linearity == True:
             #   means.append(experiment_mean[:1189])
        else:
            x_values = len(time)
            
    if linearity == True:
        #print(np.shape(data["CTRL"][1][:1548]))
        linearity_members = []
        ctrls = []
        lenss = []
        
        for ens in ensemble_members:
            linearity_members.append(np.array(data["CTRL"][ens])+np.array(data["LENS"][ens][:][:1548])-np.array(data["WIND"][ens][:][:1548])-np.array(data["TEMP"][ens][:][:1548]))
            ctrls.append(np.array(data["CTRL"][ens]))
            lenss.append(np.array(data["LENS"][ens]))
        linearity = np.mean(linearity_members, axis = 0)
        ctrl_mean = np.mean(ctrls, axis = 0)
        lens_mean = np.mean(lenss, axis = 0)
        print(np.shape(linearity))
        percentage = moving_average(linearity*100/(lens_mean-ctrl_mean), 12*30)
        
        ax2 = ax.twinx()
        ax2.plot(percentage, color = 'black', label = "linearity")
        ax2.set_ylabel("Linearity/total change", fontsize=14)
        ax2.legend()

    for exp, info in experiment_info.items():
        if info["label"]:
            # Make sure the label is set only once
            info["lines"][0].set_label(info["label"])

    ax.legend()
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_xticks(np.arange(0, x_values, 100))
    ax.set_xticklabels(xlabel.astype(str), rotation=45)
    ax.set_xlim([0, x_values])
    ax.grid(alpha=0.8)
    ax.set_title(title)
    
    fig.savefig(file_out)
    plt.show()



def main():
    if len(sys.argv) != 2:
        sys.exit("Stopped - Incorrect number of arguements. Use python PI_timeseries.py <var>")
    
    # set up the variables you need
    var = str(sys.argv[1])
    experiments = ["WIND", "TEMP", "CTRL", "LENS"]
    ensemble = [1,2,3]
    
    # Define the output structure
    data = {exp: {ens: [] for ens in ensemble} for exp in experiments}
    time = []

    for exp in experiments:
        timeseries = []
        for ens in ensemble:
            # read the filename of the experiment
            if exp == "LENS":
                filepath = output_path+"timeseries2101_experiment"+str(ens)+".nc"
            else:
                valid_file = []
                path = output_path+exp+"_ens0"+str(ens)+"_noOBC/"
                valid_file = [filename for filename in os.listdir(path) if filename.startswith("timeseries")]
                filepath = path+valid_file[0]

            data_member = xr.open_dataset(filepath)
            timeseries = data_member[var].values
            #print(exp, ens, len(timeseries))
            data[exp][ens].append(timeseries[:1548])
    time = data_member.time.values

    xlabel = np.around(np.linspace(1920, 2101, 12),-1)
    xlabel.astype(int)

    if var == "sea_ice":
        title = "max "+var
        ylabel = "Sea ice thickness (m)"
        smooth = False
        linearity = True
    else:
        title = var
        ylabel = "Temperature (*C)"
        smooth = True
        linearity = True

    file_out = "timeseries_comparison_linearity_"+var+".png"

    plot_comparison(var, data, experiments, ensemble, ylabel, xlabel, time, title, file_out, smooth=smooth, linearity = linearity)
    
if __name__ == "__main__":
    main()
    
    
