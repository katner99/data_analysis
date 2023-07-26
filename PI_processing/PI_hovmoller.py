import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt

import numpy as np
import xarray as xr

from directories_and_paths import *
    

def main():
    var = "theta_hovmoller"

    data = xr.open_dataset(ensemble_mean_path + "wind_ensemble_mean_hov.nc")
    dates = np.arange(2070, 2101, 2)
    
    # set up the grid
    [T, Z] = np.meshgrid(range(len(data.time.values)), data.depth.values)
    
    fig, ax = plt.subplots(figsize=(15, 5))

    plt.rc('font', size=16) #controls default text size
    plt.rc('axes', titlesize=16) #fontsize of the title
    plt.rc('axes', labelsize=16) #fontsize of the x and y labels
    plt.rc('xtick', labelsize=16) #fontsize of the x tick labels
    plt.rc('ytick', labelsize=16) #fontsize of the y tick labels

    # , ticks = np.arange(-1.5, 1.7, 0.25)
    cs = plt.contourf(T, Z, np.transpose(data.theta_std_hovmoller.values), levels = np.linspace(0.25, 1.5, 20), extend = "both", cmap = "RdPu")
    plt.ylim(-2000,0)
    plt.title("Temperature over depth in the wind simulation (°C)")
    plt.xticks(np.arange(0, len(data.time.values), 24), dates, rotation = 45)
    plt.ylabel("Depth (m)")
    plt.colorbar(cs)
    fig.savefig("std_hovmoller.png", transparent = True )
    plt.show()    
    
if __name__ == "__main__":
    main()
    
