import matplotlib.pyplot as plt
from directories_and_paths import output_path
import sys
import xarray as xr
from matplotlib import animation

from PIL import Image, ImageFilter
import numpy as np
from calcs import moving_average
from mitgcm_python.grid import Grid

def read_var_fluxes(var_name):
    filepaths = [
        f"{output_path}{exp}_files_temp/{var_name}_spread.nc"
        for exp in ["LENS", "WIND", "TEMP"]
    ]
    
    for filepath in filepaths:
        try:
            open(filepath)
        except FileNotFoundError:
            sys.exit(f"Stopped - Could not find input file {filepath}")

    input_data = [
        xr.open_dataset(filepath, decode_times=False) for filepath in filepaths
    ]

    return input_data 

def read_var_profile():
    period = "2070-2100"
    experiments = ["CTRL", "LENS", "WIND", "TEMP"]

    filepaths = [
        output_path + "average_" + exp + "_" + period + ".nc"
        for exp in experiments
    ]

    for filepath in filepaths:
        try:
            open(filepath)
        except FileNotFoundError:
            sys.exit(f"Stopped - Could not find input file {filepath}")

    input_data = [xr.open_dataset(filepath, decode_times=False) for filepath in filepaths]

    return input_data

def pretty_labels(ax, both="all"):
    from mitgcm_python.plot_utils.labels import lat_label, lon_label
    """
    Adjusts the labels on the given matplotlib axis `ax` to display longitude and latitude
    values in a more readable format.

    Parameters:
    - ax (matplotlib.axis): The axis for which to adjust the labels.

    Returns:
    - None
    """
    if both == "all" or both =="lon":
        lon_ticks = ax.get_xticks() - 360
        lon_labels = []
        for x in lon_ticks:
            lon_labels.append(lon_label(x, 2))
        ax.set_xticklabels(lon_labels[:-1], size = 12)
        ax.tick_params(axis="x", labelrotation=45)
    if both == "all" or both =="lat":
        ax.locator_params(axis='y', nbins=6)
        lat_ticks = ax.get_yticks()
        lat_labels = []
        for y in lat_ticks:
            lat_labels.append(lat_label(y, 2))
        ax.set_yticklabels(lat_labels, size = 12)


def create_mask(depth, ice_mask):
    """
    Creates masks for land, ice shelf, and continental shelf based on depth and ice coverage.

    Parameters:
        depth (numpy.ndarray): Array representing depth information.
        ice_mask (numpy.ndarray): Array representing ice coverage.

    Returns:
        tuple: A tuple containing:
            - land_mask (numpy.ndarray): Mask for land areas.
            - mask (numpy.ndarray): Combined mask for ice shelf and continental shelf.
            - colors (list): List of RGBA colors for plotting.
    """
    land_mask = np.zeros(np.shape(depth))
    land_mask[depth == 0] = 1

    # apply mask over the ice shelf (determiend by the ice-mask) and the continental shelf (roughly where the depth is less than 1500m)
    mask = np.zeros(np.shape(depth))
    mask[depth < 1500] = 1
    mask[ice_mask == 0] = 2

    # set the colors to block the continent (set to grey)
    colors = [(1.0, 1.0, 1.0, 0), (0.7, 0.7, 0.7, 1), (0.6, 0.6, 0.6, 1)]
    return land_mask, mask, colors


# # ANIMATE CONTOUR: animate contour plot
# def animate_contour (x, y, data, year, var, depth, exp, color_scheme = "ocean", mask = None, lonlatplot = True, font_size = 12, low_val = None, high_val = None, xlabel = None, ylabel = False, save_as = None, show=False, save=True):
#     # prepare grid
#     [X, Y] = np.meshgrid(x, y)

#     # prepare values so all months have the same parameters
#     if low_val == None:
#         low_val = np.nanmin(data)
#     if high_val == None:
#         high_val = np.nanmax(data)
#     step = (high_val-low_val)/15

#     # prepares the title
#     labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

#     # First set up the figure, the axis, and the plot element we want to animate
#     fig = plt.figure(figsize=(8,6))
#     ax = plt.axes(xlim=(min(x), max(y)), ylim=(min(y), max(y)))
#     plt.xlabel(xlabel)
#     plt.ylabel(ylabel)

#     # animation function
#     def animate(i):
#         z = data[i,:,:]

#         #cont = contour_plots(x = x, y = y, data = z, year = year, var = var, depth = depth, exp = exp, color_scheme = color_scheme, mask = mask, lonlatplot = lonlatplot, font_size = font_size, low_val = low_val, high_val = high_val, xlabel = xlabel, ylabel = ylabel, title = None, save_as = None, show = False, save = False)

#         plt.title(var+" "+labels[i]+" "+year)
#         #return cont

#     plt.colorbar(animate(0))

#     # animate
#     anim = animation.FuncAnimation(fig, animate, frames=12, interval = 200)

#     if save == True:
#         if save_as is not None:
#             anim.save(save_as)
#         else:
#             anim.save(exp +"_cont_"+year+"_"+var+".gif", fps = 2)

#     # show figure
#     if show == True:
#         plt.show()


def read_u_and_v(input_data, option="avg"):
    """
    Reads and processes u and v velocity components from input data.

    Parameters:
        input_data (xarray.Dataset): Dataset containing velocity data.
        option (str, optional): Velocity option for plotting. Defaults to 'avg'.

    Returns:
        tuple: A tuple containing:
            - speed (numpy.ndarray): Array representing speed of velocity vectors.
            - u_plot (numpy.ndarray): Array representing processed u velocity component for plotting.
            - v_plot (numpy.ndarray): Array representing processed v velocity component for plotting.
    """
    from mitgcm_python.utils import mask_3d
    from mitgcm_python.plot_utils.latlon import prepare_vel
    from config_options import days_in_month
    from directories_and_paths import grid_filepath

    vvel = input_data.VVEL.values
    uvel = input_data.UVEL.values
    grid = Grid(grid_filepath)

    uvel = mask_3d(uvel, grid, gtype="u", time_dependent=True)
    vvel = mask_3d(vvel, grid, gtype="v", time_dependent=True)

    u = np.average(uvel, axis=0, weights=days_in_month)
    v = np.average(vvel, axis=0, weights=days_in_month)

    speed, u_plot, v_plot = prepare_vel(u, v, grid, vel_option=option)
    return speed, u_plot, v_plot


def zoom_shelf(ax, zoom):
    if zoom == "ice_shelf":
        ax.set_ylim([-75.6, -73])
        ax.set_xlim([245, 262])
    elif zoom == "cont_shelf":
        ax.set_xlim([230, 265])
        ax.set_ylim([-75.5, -68])


def read_mask(input_data=None, cut=None, lat_range=None, lon_range=None):
    """
    Reads and processes mask data from input data.

    Parameters:
        input_data (xarray.Dataset): Dataset containing mask data.

    Returns:
        dict: Dictionary containing setup information for plotting including latitude, longitude, depth, ice mask, land mask, combined mask, colors, X, and Y coordinates.
    """
    import xarray as xr
    from directories_and_paths import output_path

    if input_data is None:
        input_data = xr.open_dataset(
            f"{output_path}average_CTRL_1920-1950.nc", decode_times=False
        )

    if cut is None:
        [lat, lon, ice_mask_temp, depth] = [
            input_data[param].values for param in ["YC", "XC", "maskC", "Depth"]
        ]
        ice_mask = ice_mask_temp[0, :, :]
        [land_mask, mask, colors] = create_mask(depth, ice_mask)
        [X, Y] = np.meshgrid(lon, lat)

        set_up = {
            "lat": lat,
            "lon": lon,
            "depth": depth,
            "ice_mask": ice_mask,
            "land_mask": land_mask,
            "mask": mask,
            "colors": colors,
            "X": X,
            "Y": Y,
        }

    elif cut == "lat":
        z = input_data.Z.values
        lat = input_data["YC"][lat_range[0] : lat_range[1]].values
        lon = input_data["XC"][lon_range]
        ice_mask = input_data.maskC.values[:, lat_range[0] : lat_range[1], lon_range]

        set_up = {
            "lat": lat,
            "lon": lon,
            "z": z,
            "ice_mask": ice_mask,
        }

    return set_up


def plot_timeseries_comparison(
    ax, data, experiments, ensemble_members, plot_info, standard_deviation=False
):
    """
    This function generates a plot comparing different experiments based on the provided data and plot information.

    Parameters:
    - data: Dictionary containing the data for different experiments
    - experiments: List of experiment names to be compared
    - ensemble_members: List of ensemble members
    - plot_info: Dictionary containing plot configuration parameters including:
      - ylabel: Label for the y-axis
      - xlabel: Label for the x-axis
      - time: Time information
      - file_out: Output file name for saving the plot
      - smooth: Smoothing factor for the data (optional, default is 0)
      - linearity: Flag indicating whether to plot linearity (optional, default is False)
      - warming: Flag indicating whether to adjust for warming (optional, default is False)
      - percentage: Flag indicating whether to calculate percentage (optional, default is False)
      - shade_range: Flag indicating whether to shade the range (optional, default is True)
      - title: Title for the plot (optional)
      - experiment_full: List containing full names of experiments (optional)

    Returns:
    - ax: Axis object representing the main plot
    - all_means: List containing the means of all experiments' data
    """
    colors = ["slategrey", "deeppink", "dodgerblue", "orange"]
    experiment_full = plot_info["experiment_full"]
    all_means = []
    size = 16

    if plot_info.get("warming", False):
        ctrl_mean = np.nanmean(
            [
                data["CTRL"][ens][ts_idx][:]
                for ens in ensemble_members
                for ts_idx in range(len(data["CTRL"][ens]))
            ],
            axis=0,
        )
        warming_diff = -np.nanmean(ctrl_mean)

    for exp_idx, exp in enumerate(experiments):
        experiment_mean = np.nanmean(
            [
                data[exp][ens][ts_idx][:]
                for ens in ensemble_members
                for ts_idx in range(len(data[exp][ens]))
            ],
            axis=0,
        )
        if standard_deviation == False:
            experiment_max = np.nanmax(
                [
                    data[exp][ens][ts_idx][:]
                    for ens in ensemble_members
                    for ts_idx in range(len(data[exp][ens]))
                ],
                axis=0,
            )
            experiment_min = np.nanmin(
                [
                    data[exp][ens][ts_idx][:]
                    for ens in ensemble_members
                    for ts_idx in range(len(data[exp][ens]))
                ],
                axis=0,
            )
        else:
            experiment_max = experiment_mean + 2 * (
                np.nanstd(
                    [
                        data[exp][ens][ts_idx][:]
                        for ens in ensemble_members
                        for ts_idx in range(len(data[exp][ens]))
                    ],
                    axis=0,
                )
            )
            experiment_min = experiment_mean - 2 * (
                np.nanstd(
                    [
                        data[exp][ens][ts_idx][:]
                        for ens in ensemble_members
                        for ts_idx in range(len(data[exp][ens]))
                    ],
                    axis=0,
                )
            )

        if plot_info.get("smooth", 0) > 0:
            smoothed_mean = moving_average(experiment_mean, plot_info["smooth"] * 12)
            if plot_info.get("shade_range", True):
                smoothed_max = moving_average(experiment_max, plot_info["smooth"] * 12)
                smoothed_min = moving_average(experiment_min, plot_info["smooth"] * 12)
        else:
            smoothed_mean, smoothed_max, smoothed_min = (
                experiment_mean,
                experiment_max,
                experiment_min,
            )

        if plot_info.get("warming", False):
            smoothed_mean += warming_diff
            smoothed_max += warming_diff
            smoothed_min += warming_diff

        ax.plot(
            smoothed_mean,
            color=colors[exp_idx],
            label=experiment_full[exp_idx],
            alpha=1,
        )

        if plot_info.get("shade_range", True):
            ax.fill_between(
                range(len(smoothed_min)),
                smoothed_min,
                smoothed_max,
                color=colors[exp_idx],
                alpha=0.175,
            )

        all_means.append(experiment_mean[:])

    if plot_info.get("linearity", False) or plot_info.get("percentage", False):
        year_of_divergence = 1169
        linearity = moving_average(
            all_means[0] + all_means[1] - all_means[2] - all_means[3],
            12 * plot_info.get("smooth", 0),
        )
        linearity[:935] = np.nan
        ax.plot(linearity, color="black", LineStyle="dashed", label="non-linearity")
        ax.axvline(year_of_divergence, color="black", linewidth=2)
        ax.axhline(0, color="grey")

    ax.set_ylabel(plot_info["ylabel"], fontsize = size)
    ax.set_xticks(
        np.arange(
            6 * plot_info.get("smooth", 0),
            plot_info["time"] - (6 * plot_info.get("smooth", 0)),
            120,
        )
    )
    ax.set_xticklabels(plot_info["xlabel"].astype(str), rotation=45, fontsize = size)
    ax.set_xlim(plot_info["x_lim"])
    ax.set_title(plot_info.get("title"), fontsize = size)
    ax.tick_params(axis="y", labelsize=size)
    if plot_info.get("warming", False):
        ax.set_ylim([-1.6, 1.6])

    ax.legend(frameon=False, fontsize = size)
    ax.grid(alpha=0.8)

    return ax, all_means
