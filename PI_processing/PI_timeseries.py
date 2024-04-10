import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

import numpy as np
from funcs import read_timeseries
from plots import plot_timeseries_comparison

def main():
    # set variables for the function
    var = "melt"
    
    experiments = ["CTRL", "LENS", "TEMP", "WIND"]
    ensemble = [1, 2, 3, 4, 5, 6, 7, 8, 9]

    data = read_timeseries(experiments, ensemble, var)
    
    plot_info = {
        "var": var,
        "data": data,
        "experiments": experiments,
        "ensemble_members": ensemble,
        "ylabel": "melt (Gt/yr)",
        "time": 181*12,
        "xlabel": np.arange(1921, 2100, 10),
        "x_lim": [12, (2099-1920)*12],
        "smooth": 2,
        "shade_range": True,
        "experiment_full": experiments,
    }

    file_out = "timeseries_" + var + ".png"

    fig, ax = plt.subplots(figsize=(15, 10), sharex=True)

    plot_timeseries_comparison(ax, data, experiments, ensemble, plot_info)

    fig.savefig(file_out)
    plt.show()
   
if __name__ == "__main__":
    main()
