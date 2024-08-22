"""
This module generates a time series comparison plot for a given variable (e.g., temperature) 
across different experiments and ensemble members.

Dependencies:
- matplotlib (imported as plt)
- numpy (imported as np)
- read_timeseries function from funcs module
- plot_timeseries_comparison function from plots module
- config_timeseries function from config_options module

Main Function:
- main(): Main function that orchestrates the generation of the time series comparison plot. 
It defines the variable of interest, experiments, and ensemble members, reads the time series 
data, configures plotting information, generates the plot, saves it as a PNG file, and displays 
the plot.

Usage:
- Run the script to generate the time series comparison plot.
"""
import matplotlib.pyplot as plt

from .tools.funcs import read_timeseries
from .tools.plots import plot_timeseries_comparison
from config_options import config_timeseries


def main():
    var = "theta"

    experiments = ["CTRL", "LENS", "WIND", "TEMP"]
    #experiments = ["CTRL", "LENS", "TEMP", "WIND"]
    ensemble = [1, 2, 3, 4, 5, 6, 7, 8, 9]

    data = read_timeseries(experiments, ensemble, var)

    plot_info = config_timeseries(var, data, experiments, ensemble)

    file_out = "timeseries_" + var + ".png"

    fig, ax = plt.subplots(figsize=(15, 10), sharex=True)

    plot_timeseries_comparison(ax, data, experiments, ensemble, plot_info)
    ax.set_title("Potential temperature anomaly from pre-industrial mean", fontsize = 16, fontweight = 'bold')
    # ax.axvline(935, color = 'black', linewidth=2, alpha=0.8)
    # ax.axvline(967, color = 'black', linewidth=2, alpha=0.8)
    # ax.axvline(1169, color = 'black', linewidth=2, alpha=0.8)
    # ax.axvline(1827, color = 'black', linewidth=2, alpha=0.8)
    # plt.text((1996-1920)*12, -1.4, "THERMO\n diverges\n NONE", horizontalalignment='right', color = "black", fontsize = 12, fontweight = "bold")
    # plt.text((2002-1920)*12, -1.4, "ALL \ndiverges \nWIND", color="black", fontsize = 12, fontweight = "bold")
    # plt.text((2019-1920)*12, -1.4, "ALL\ndiverges \nNONE", color="black", fontsize = 12, fontweight = "bold")
    # plt.text((2070-1920)*12, -1.4, "ALL\n diverges\n THERMO", horizontalalignment='right', color="black", fontsize = 12, fontweight = "bold")
    
    fig.savefig(file_out, bbox_inches='tight', transparent=True)
    #plt.show()


if __name__ == "__main__":
    main()
