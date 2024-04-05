import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import xarray as xr

from plots_2d import contour_func, quiver_func
from directories_and_paths import output_path, grid_filepath
from plots import read_mask, read_u_and_v
from funcs import find_nearest, read_variable
from mitgcm_python.grid import Grid

matplotlib.use("TkAgg")

def read_temp(filename, grid_filepath):
    var = "THETA"
    filepath = output_path + filename
    input_data = xr.open_dataset(filepath, decode_times=False)
    depth_range = [find_nearest(input_data["Z"].values, -200), find_nearest(input_data["Z"].values, -700)]
    grid = Grid(grid_filepath)
    temp = read_variable(input_data, var, grid, depth_range)
    return temp

def read_bottom_velocity(filepath):
    input_data = xr.open_dataset(filepath, decode_times=False)
    [speed, u, v] = read_u_and_v(input_data, "bottom")
    set_up = read_mask(input_data)
    graph_params = {
        "font_size": 12,
        "low_val": np.min(speed),
        "high_val": 0.15,
        "interval": 0.05,
        "step": 15,
        "color_scheme": "cool",
    }
    chunk = 10
    return speed, u, v, set_up, graph_params, chunk

def plot_contour(fig, ax, data, set_up, graph_params, title):
    cs = contour_func(fig, ax, data, set_up, graph_params)
    ax.set_xlim([230, 260])
    ax.set_ylim([-75.5, -68])
    ax.text(
        0.5,
        0.9,
        title,
        horizontalalignment="center",
        transform=ax.transAxes,
        bbox=dict(facecolor="white", alpha=0.9),
    )
    return cs

def main():
    ctrl_temp = read_temp("average_CTRL_1920-1950.nc", grid_filepath)
    lens_temp = read_temp("average_LENS_2070-2100.nc", grid_filepath)
    [speed, u, v, set_up, graph_params_vel, chunk] = read_bottom_velocity(output_path + "average_CTRL_1920-1950.nc")

    graph_params_temp = {
        "font_size": 12,
        "low_val": -2,
        "high_val": 2.1,
        "interval": 0.5,
        "step": 15,
        "color_scheme": "coolwarm",
    }

    fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True, gridspec_kw={'width_ratios': [1, 1]})
    plt.subplots_adjust(wspace=0.1, hspace=0.1)
    axs = axs.flatten()

    plot_contour(fig, axs[0], ctrl_temp, set_up, graph_params_temp, "Pre-industrial temperature (degC)", hide_ticks_y=False)
    cs_vel = plot_contour(fig, axs[1], speed, set_up, graph_params_vel, "Bottom Velocity (m/s)")
    quiver_func(axs[1], u, v, set_up["lat"], set_up["lon"], chunk)
    cs_temp = plot_contour(fig, axs[2], lens_temp, set_up, graph_params_temp, "End of century temperature (degC)", hide_ticks_x=False, hide_ticks_y=False)

    ticks = np.arange(graph_params_vel["low_val"], graph_params_vel["high_val"] + 0.1, graph_params_vel["interval"])
    cbar_ax = fig.add_axes([0.92, 0.50, 0.02, 0.4])
    cbar = plt.colorbar(cs_vel, cax=cbar_ax)
    cbar.set_ticks(ticks)

    ticks = np.arange(graph_params_temp["low_val"], graph_params_temp["high_val"] + 0.1, graph_params_temp["interval"])
    cbar_ax = fig.add_axes([0.05, 0.15, 0.02, 0.7])
    cbar = plt.colorbar(cs_temp, cax=cbar_ax, orientation='vertical')
    cbar.set_ticks(ticks)
    cbar.ax.yaxis.set_ticks_position('left') # Adjust padding here

    plt.show()

if __name__ == "__main__":
    main()

