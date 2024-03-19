from funcs import read_variable, find_nearest
from plots import create_mask
from directories_and_paths import *
from mitgcm_python.grid import Grid

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import ticker

import sys
import numpy as np
import xarray as xr

def main():
    """
    Script to look at differences in winds between runs 
    """
    
    filepaths = [output_path + "wind_" + ens + "_2080-2100.nc" for ens in ["old", "new"]]

    for filepath in filepaths:
        try:
            open(filepath)
        except FileNotFoundError:
            sys.exit(f"Stopped - Could not find input file {filepath}")

    input_data = [xr.open_dataset(filepath, decode_times = False) for filepath in filepaths]

    [lat, lon, ice_mask_temp, depth] = [input_data[0][param].values for param in ["YC", "XC", "maskC", "Depth"]]
    ice_mask = ice_mask_temp[0,:,:]

    grid = Grid(grid_filepath)
    
    var = "THETA"
    wind_u = [read_variable(input, var, grid) for input in input_data]
    var = "EXFuwind"
    wind_v = [read_variable(input, var, grid) for input in input_data]
    speed_old = np.sqrt(wind_u[0]**2 + wind_v[0]**2)
    speed_new = np.sqrt(wind_u[1]**2 + wind_v[1]**2)

    u_old = wind_u[0]
    v_old = wind_v[0]

    u_new = wind_u[1]
    v_new = wind_v[1]
    
    max_speed = np.max(speed_old)
    min_speed = 0
    # set mask
    [land_mask, mask, colors] = create_mask(depth, ice_mask)

    # set up the grid
    [X, Y] = np.meshgrid(lon, lat)
    [LON_short, LAT_short] = (lon[0:-1:50], lat[0:-1:50])

    fig, axs = plt.subplots(nrows=1, ncols=3, gridspec_kw={"hspace": 0.5, "wspace": 0.4}, figsize=(20, 7))
    axs = axs.flatten()

    # Plotting for the first subplot
    cs = axs[0].contourf(X, Y, speed_old, cmap="seismic", levels=np.linspace(0, max_speed, 15))
    axs[0].contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
    axs[0].contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
    quiv = axs[0].quiver(LON_short, LAT_short, u_old[0:-1:50, 0:-1:50], v_old[0:-1:50, 0:-1:50], scale = 50, color="white")
    axs[0].set_title("old", weight="bold")

    # Plotting for the second subplot
    cs = axs[1].contourf(X, Y, speed_new, cmap="seismic", levels=np.linspace(0, max_speed, 15))
    axs[1].contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
    axs[1].contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
    quiv = axs[1].quiver(LON_short, LAT_short, u_new[0:-1:50, 0:-1:50], v_new[0:-1:50, 0:-1:50], scale = 50, color="white")
    cbar = fig.colorbar(cs, ax=axs[1], ticks=np.arange(0, 13, 2))
    axs[1].set_title("new", weight="bold")

    # Plotting for the third subplot
    cs = axs[2].contourf(X, Y, speed_new - speed_old, cmap="seismic", extend = "both", levels=np.linspace(-max_speed, max_speed, 15))
    axs[2].contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
    axs[2].contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
    quiv = axs[2].quiver(LON_short, LAT_short, u_new[0:-1:50, 0:-1:50] - u_old[0:-1:50, 0:-1:50],
                         v_new[0:-1:50, 0:-1:50] - v_old[0:-1:50, 0:-1:50], scale = 50, color="black")
    qk = axs[2].quiverkey(quiv, 0.9, 0.9, 20, r'$20 \frac{m}{s}$', labelpos='E', coordinates='figure')  # Use the appropriate scale value
    cbar = fig.colorbar(cs, ax=axs[2], ticks=np.arange(-12, 13, 4))
    axs[2].set_title("anomaly", weight="bold")


    fig.savefig("relative_wind_diff.png")

    plt.show()
    
if __name__ == '__main__':
    main()
