import numpy as np
from scipy import stats
from plots import plot_timeseries_comparison
from directories_and_paths import *
from config_options import *
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

def calc_pvalue(data, experiments, ensemble):
    """
    The calc_pvalue function efficiently computes p-values for statistical comparison between 
    two specified experiments using ensemble data. It iterates over time intervals, extracting 
    temperature values for each experiment and ensemble member to perform independent t-tests. 
    The resulting p-values are collected and returned as a comprehensive list for further analysis 
    or visualization.
    Parameters:
        data: Dictionary containing experimental data.
        experiments: List containing identifiers for the two experiments to compare.
        ensemble: List of ensemble members.
    Returns:
        p_value: A list of p-values corresponding to each time interval, indicating the significance 
        of differences between the specified experiments.
    """
    p_value = []
    
    for time in range(181*12):
        values_1 = np.array([data[experiments[0]][ens][0] for ens in ensemble])[:, time]
        values_2 = np.array([data[experiments[1]][ens][0] for ens in ensemble])[:, time]
        
        test = stats.ttest_ind(values_1, values_2)
        p_value.append(test.pvalue)
    return p_value
    
def stat_timeserie(var, data, ensemble, experiments, experiment_names, file_out = None):
    """
    Parameters:
        data: Dictionary containing experimental data.
        ensemble: List of ensemble members.
        experiments: List containing two experiment identifiers for comparison.
        experiment_names: List of full experiment names.
        file_out (optional): Output file path to save generated plots.
    """
    
    p_value = calc_pvalue(data, experiments, ensemble)
    
    for time in range(181*12):
        values_1 = np.array([data[experiments[0]][ens][0] for ens in ensemble])[:, time]
        values_2 = np.array([data[experiments[1]][ens][0] for ens in ensemble])[:, time]
        
        test = stats.ttest_ind(values_1, values_2)
        p_value.append(test.pvalue)
    
    plot_info = {
        "var": var,
        "data": data,
        "experiments": experiments,
        "ensemble_members": ensemble,
        "ylabel": "Temperature (degC)",
        "time": 181*12,
        "xlabel": np.arange(1921, 2100, 10),
        "x_lim": [12, (2099-1920)*12],
        "smooth": 2,
        "shade_range": True,
        "experiment_full": experiment_names,
    }
    
    fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Plot p-value against time
    axs[0].plot(p_value, color='black')
    axs[0].set_xlim([12, 180*12])
    axs[0].set_ylabel('P-value', fontsize=14)
    axs[0].set_title(f"P-value {experiments[0]} vs. {experiments[1]}", fontsize=14)
    axs[0].axhline(0.05, color='red')

    # Plot time series comparison
    plot_timeseries_comparison(axs[1], var, data, experiments, ensemble, plot_info)

    plt.tight_layout()
    plt.show()
    
    if file_out:
        fig.savefig(file_out)
    else:
        fig.savefig(f"pvalue_{experiments[0]}vs{experiments[1]}.png")
    

