import matplotlib.pyplot as plt

import numpy as np
import xarray as xr

from plots_2d import plot_contour, quiver_func, trend_quiver_func
from directories_and_paths import output_path, grid_filepath
from plots import read_mask, read_u_and_v
from funcs import find_nearest, read_variable
from mitgcm_python.grid import Grid



def read_temp(filename, grid_filepath, create_graph=False):
    """
    Reads temperature data from a specified file and calculates the average potential temperature 
    between 200 and 700 meters depth. Optionally, returns parameters for creating a temperature plot.

    Parameters:
    filename (str): The name of the file containing temperature data.
    grid_filepath (str): The path to the grid file.
    create_graph (bool, optional): If True, returns parameters for creating a temperature plot. Default is False.

    Returns:
    temp (xr.DataArray): The average potential temperature between 200 and 700 meters depth.
    graph_params (dict, optional): If create_graph is True, returns a dictionary containing parameters 
                                   for creating a plot, including font size, value range, interval, 
                                   step, and color scheme.
    """
    var = "THETA"
    filepath = output_path + filename
    input_data = xr.open_dataset(filepath, decode_times=False)
    depth_range = [
        find_nearest(input_data["Z"].values, -200),
        find_nearest(input_data["Z"].values, -700),
    ]
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
    """
    Reads bottom velocity data from a specified file and prepares parameters for visualizing the data.

    Parameters:
    filepath (str): The path to the file containing the bottom velocity data.

    Returns:
    speed (xr.DataArray): The speed of the bottom current.
    u (xr.DataArray): The u-component (east-west) of the bottom velocity.
    v (xr.DataArray): The v-component (north-south) of the bottom velocity.
    set_up (xr.DataArray): The setup mask for the bottom velocity data.
    graph_params (dict): A dictionary containing parameters for creating a plot, including font size, 
                         value range, interval, step, and color scheme.
    chunk (int): The chunk size for plotting or processing the data.
    """
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
    """
    Reads wind data from a specified file and extracts the u and v wind components starting from the year 2000.
    Additionally, prepares parameters for visualizing the wind data.

    Parameters:
    filename (str): The name of the file containing wind data.

    Returns:
    uwinds (np.ndarray): The u-component (east-west) of the wind.
    vwinds (np.ndarray): The v-component (north-south) of the wind.
    time (np.ndarray): The time values corresponding to the wind data.
    graph_params (dict): A dictionary containing parameters for creating a plot, including font size, 
                         value range, interval, step, and color scheme.
    """
    data = xr.open_dataset(filename, decode_times=False)
    period_start = 2000
    start_index = (period_start - 1920) * 12

    uwinds = data.u_wind.values[start_index:, ...]
    vwinds = data.v_wind.values[start_index:, ...]
    time = data.time.values[start_index:, ...]
    graph_params = {
        "font_size": 12,
        "low_val": 0,
        "high_val": 4000,
        "interval": 1000,
        "step": 15,
        "color_scheme": "bone_r",
    }
    return uwinds, vwinds, time, graph_params

def plot_intro(
    ctrl_temp,
    lens_temp,
    graph_params_temp,
    speed,
    u,
    v,
    set_up,
    graph_params_vel,
    chunk,
    graph_params_bathy,
    uwind,
    vwind,
    time,
):
    """
    Creates a figure with multiple subplots to visualize temperature, bottom velocity, 
    and projected wind trends. 

    Parameters:
    ctrl_temp (xr.DataArray or np.ndarray): The pre-industrial averaged temperature data.
    lens_temp (xr.DataArray or np.ndarray): The end-of-century averaged temperature data.
    graph_params_temp (dict): Parameters for plotting temperature data, including font size, 
                              value range, interval, step, and color scheme.
    speed (xr.DataArray or np.ndarray): The speed of the bottom current.
    u (xr.DataArray or np.ndarray): The u-component (east-west) of the bottom velocity.
    v (xr.DataArray or np.ndarray): The v-component (north-south) of the bottom velocity.
    set_up (dict): A dictionary containing setup data, including depth, latitude, and longitude masks.
    graph_params_vel (dict): Parameters for plotting velocity data, including font size, 
                             value range, interval, step, and color scheme.
    chunk (int): The chunk size for plotting or processing the velocity data.
    graph_params_bathy (dict): Parameters for plotting bathymetry data, including font size, 
                               value range, interval, step, and color scheme.
    uwind (xr.DataArray or np.ndarray): The u-component (east-west) of the wind.
    vwind (xr.DataArray or np.ndarray): The v-component (north-south) of the wind.
    time (np.ndarray): The time values corresponding to the wind data.

    Returns:
    None
    """
    fig, axs = plt.subplots(
        2, 2, figsize=(10, 8), sharex=True, gridspec_kw={"width_ratios": [1, 1]}
    )
    plt.subplots_adjust(wspace=0.05, hspace=0.05)
    axs = axs.flatten()

    plot_contour(
        axs[0],
        ctrl_temp,
        set_up,
        graph_params_temp,
        "a. Pre-industrial averaged temperature (degC)",
        hide_ticks_y=False,
    )
    cs_vel = plot_contour(
        axs[1], speed, set_up, graph_params_vel, "b. Bottom Velocity (m/s)"
    )
    cs_temp = plot_contour(
        axs[2],
        lens_temp,
        set_up,
        graph_params_temp,
        "c. End of century averaged temperature (degC)",
        hide_ticks_x=False,
        hide_ticks_y=False,
    )
    cs_bathy = plot_contour(
        axs[3],
        set_up["depth"],
        set_up,
        graph_params_bathy,
        "d. Projected wind trends (m/s/century)",
        hide_ticks_x=False,
    )

    quiver_func(axs[1], u, v, set_up["lat"], set_up["lon"], chunk)
    trend_quiver_func(axs[3], uwind, vwind, time, set_up)

    map = plt.imread("maploc.png")
    map_ax = fig.add_axes([0.02, 0.825, 0.12, 0.12])
    map_ax.imshow(map)
    map_ax.axis("off")

    ticks = np.arange(
        graph_params_vel["low_val"],
        graph_params_vel["high_val"] + 0.1,
        graph_params_vel["interval"],
    )
    cbar_ax = fig.add_axes([0.92, 0.55, 0.02, 0.35])
    cbar = plt.colorbar(cs_vel, cax=cbar_ax)
    cbar.set_ticks(ticks)

    ticks = np.arange(
        graph_params_temp["low_val"],
        graph_params_temp["high_val"] + 0.1,
        graph_params_temp["interval"],
    )
    cbar_ax = fig.add_axes([0.05, 0.1, 0.02, 0.7])
    cbar = plt.colorbar(cs_temp, cax=cbar_ax, orientation="vertical")
    cbar.set_ticks(ticks)
    cbar.ax.yaxis.set_ticks_position("left")

    ticks = np.arange(
        graph_params_bathy["low_val"],
        graph_params_bathy["high_val"] + 0.1,
        graph_params_bathy["interval"],
    )
    cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.35])
    cbar = plt.colorbar(cs_bathy, cax=cbar_ax)
    cbar.set_ticks(ticks)

    fig.savefig("introduction_plot.png", transparent=False)

    plt.show()


def main():
    ctrl_temp = read_temp("average_CTRL_1920-1950.nc", grid_filepath)
    [lens_temp, graph_params_temp] = read_temp(
        "average_LENS_2070-2100.nc", grid_filepath, True
    )
    [speed, u, v, set_up, graph_params_vel, chunk] = read_bottom_velocity(
        output_path + "average_CTRL_1920-1950.nc"
    )
    [uwind, vwind, time, graph_params_bathy] = read_winds(
        f"{output_path}high_emissions_winds.nc"
    )
    plot_intro(
        ctrl_temp,
        lens_temp,
        graph_params_temp,
        speed,
        u,
        v,
        set_up,
        graph_params_vel,
        chunk,
        graph_params_bathy,
        uwind,
        vwind,
        time,
    )


if __name__ == "__main__":
    main()
