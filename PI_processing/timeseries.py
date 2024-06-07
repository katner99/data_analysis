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
import matplotlib

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

from funcs import read_timeseries
from plots import plot_timeseries_comparison
from config_options import config_timeseries


def main():
    var = "theta"

    experiments = ["CTRL", "LENS", "TEMP", "WIND"]
    ensemble = [1, 2, 3, 4, 5, 6, 7, 8, 9]

    data = read_timeseries(experiments, ensemble, var)

    plot_info = config_timeseries(var, data, experiments, ensemble)

    file_out = "timeseries_" + var + ".png"

    fig, ax = plt.subplots(figsize=(15, 10), sharex=True)

    plot_timeseries_comparison(ax, data, experiments, ensemble, plot_info)
    ax.set_title("Temperature normalised to the pre-industrial mean", fontsize = 16, fontweight = 'bold')
    ax.axvline((2000-1920)*12, color = 'dodgerblue', linewidth=2, alpha=0.5)
    ax.axvline((1997-1920)*12, color = 'orange', linewidth=2, alpha=0.5)
    plt.text((1996-1920)*12, -1.3, "thermodynamic scenario \n diverges from pre-industrial scenario ", horizontalalignment='right', color = "orange", fontsize = 14, fontweight = "bold")
    plt.text((2001-1920)*12, -1.3, " wind scenario diverges from \n pre-industrial scenario", color="dodgerblue", fontsize = 14, fontweight = "bold")
    fig.savefig(file_out, bbox_inches='tight', transparent=True)
    #plt.show()


if __name__ == "__main__":
    main()
