import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import xarray as xr

from plots_2d import contour_func, quiver_func, trend_quiver_func
from directories_and_paths import output_path, grid_filepath
from plots import read_mask, read_u_and_v
from funcs import find_nearest, read_variable
from mitgcm_python.grid import Grid

matplotlib.use("TkAgg")

def read_temp(filename, grid_filepath, create_graph = False):
    var = "THETA"
    filepath = output_path + filename
    input_data = xr.open_dataset(filepath, decode_times=False)
    depth_range = [find_nearest(input_data["Z"].values, -200), find_nearest(input_data["Z"].values, -700)]
    grid = Grid(grid_filepath)
    temp = read_variable(input_data, var, grid, depth_range)
    if create_graph:
        graph_params = {
            "font_size": 12,
            "low_val": -2,
            "high_val": 2.1,
            "interval": 0.5,
            "step": 15,
            "color_scheme": "coolwarm",
        }
        return temp, graph_params
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

def read_winds(filename):
    data = xr.open_dataset(filename, decode_times=False)
    period_start = 2000
    start_index = (period_start - 1920) * 12

    uwinds = data.u_wind.values[start_index:,...]
    vwinds = data.v_wind.values[start_index:,...]
    time = data.time.values[start_index:,...]
    graph_params = {
        "font_size": 12,
        "low_val": 0,
        "high_val": 4000,
        "interval": 1000,
        "step": 15,
        "color_scheme": "bone_r",
    }
    return uwinds, vwinds, time, graph_params

def plot_contour(fig, ax, data, set_up, graph_params, title, hide_ticks_x = True, hide_ticks_y = True):
    cs = contour_func(fig, ax, data, set_up, graph_params, hide_ticks_x, hide_ticks_y)
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

def plot_intro(ctrl_temp, lens_temp, graph_params_temp, speed, u, v, set_up, graph_params_vel, chunk, graph_params_bathy, uwind, vwind, time):
    fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True, gridspec_kw={'width_ratios': [1, 1]})
    plt.subplots_adjust(wspace=0.05, hspace=0.05)
    axs = axs.flatten()

    plot_contour(fig, axs[0], ctrl_temp, set_up, graph_params_temp, "Pre-industrial temperature (degC)", hide_ticks_y=False)
    cs_vel = plot_contour(fig, axs[1], speed, set_up, graph_params_vel, "Bottom Velocity (m/s)")
    cs_temp = plot_contour(fig, axs[2], lens_temp, set_up, graph_params_temp, "End of century temperature (degC)", hide_ticks_x=False, hide_ticks_y=False)
    cs_bathy = plot_contour(fig, axs[3], set_up["depth"], set_up, graph_params_bathy, "Projected wind trends (m/s/century)", hide_ticks_x=False)
    
    quiver_func(axs[1], u, v, set_up["lat"], set_up["lon"], chunk)
    trend_quiver_func(axs[3], uwind, vwind, time, set_up)

    ticks = np.arange(graph_params_vel["low_val"], graph_params_vel["high_val"] + 0.1, graph_params_vel["interval"])
    cbar_ax = fig.add_axes([0.92, 0.55, 0.02, 0.35])
    cbar = plt.colorbar(cs_vel, cax=cbar_ax)
    cbar.set_ticks(ticks)

    ticks = np.arange(graph_params_temp["low_val"], graph_params_temp["high_val"] + 0.1, graph_params_temp["interval"])
    cbar_ax = fig.add_axes([0.05, 0.15, 0.02, 0.7])
    cbar = plt.colorbar(cs_temp, cax=cbar_ax, orientation='vertical')
    cbar.set_ticks(ticks)
    cbar.ax.yaxis.set_ticks_position('left') 

    ticks = np.arange(graph_params_bathy["low_val"], graph_params_bathy["high_val"] + 0.1, graph_params_bathy["interval"])
    cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.35])
    cbar = plt.colorbar(cs_bathy, cax=cbar_ax)
    cbar.set_ticks(ticks)

    fig.savefig("introduction_plot.png", transparent=True)
    
    plt.show()

def main():
    ctrl_temp = read_temp("average_CTRL_1920-1950.nc", grid_filepath)
    [lens_temp, graph_params_temp] = read_temp("average_LENS_2070-2100.nc", grid_filepath, True)
    [speed, u, v, set_up, graph_params_vel, chunk] = read_bottom_velocity(output_path + "average_CTRL_1920-1950.nc")
    [uwind, vwind, time, graph_params_bathy] = read_winds(f"{output_path}high_emissions_winds.nc")
    plot_intro(ctrl_temp, lens_temp, graph_params_temp, speed, u, v, set_up, graph_params_vel, chunk, graph_params_bathy, uwind, vwind, time)
    
if __name__ == "__main__":
    main()

