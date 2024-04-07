import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr 
from directories_and_paths import output_path

def contour_func(ax, data, set_up, graph_params, hide_ticks_x=True, hide_ticks_y=True):
    """
    Generates a contour plot with specified data and settings.

    Parameters:
        fig (matplotlib.figure.Figure): The figure object.
        ax (matplotlib.axes.Axes): The axes object to plot on.
        data (numpy.ndarray): The data to be plotted.
        set_up (dict): Dictionary containing setup information including X 
        and Y coordinates, land mask, and color scheme.
        graph_params (dict): Dictionary containing parameters for the contour 
        plot such as color scheme, levels, and interval.
        hide_ticks (bool, optional): Whether to hide tick labels on the plot. Defaults to True.

    Returns:
        matplotlib.axes.Axes: The axes object containing the contour plot.
    """

    cs = ax.contourf(
        set_up["X"],
        set_up["Y"],
        data,
        cmap=graph_params["color_scheme"],
        extend="both",
        levels=np.linspace(
            graph_params["low_val"], graph_params["high_val"], graph_params["step"]
        ),
    )
    ax.contourf(
        set_up["X"],
        set_up["Y"],
        set_up["land_mask"],
        cmap=matplotlib.colors.ListedColormap(set_up["colors"]),
    )
    ax.contour(
        set_up["X"], set_up["Y"], set_up["mask"], 2, cmap="Greys", linestyles="dashed"
    )
    
    if hide_ticks_x:
        ax.get_xaxis().set_visible(False)
    if hide_ticks_y:
        ax.get_yaxis().set_visible(False)
    ax.set_aspect("auto", adjustable="box")
    return cs


def quiver_func(ax, u, v, lat, lon, chunk, key = True):
    """
    Plots a quiver plot representing vector fields over specified coordinates.

    Parameters:
        ax (matplotlib.axes.Axes): The axes object to plot on.
        u (numpy.ndarray): Array representing the u component of vector field.
        v (numpy.ndarray): Array representing the v component of vector field.
        lat (numpy.ndarray): Array of latitude coordinates.
        lon (numpy.ndarray): Array of longitude coordinates.
        chunk (int): Spacing factor for subsampling the data for plotting.

    Returns:
        matplotlib.quiver.Quiver: The quiver object representing the plotted vector field.
    """
    q = ax.quiver(
        lon[0:-1:chunk],
        lat[0:-1:chunk],
        u[0:-1:chunk, 0:-1:chunk],
        v[0:-1:chunk, 0:-1:chunk],
        scale=1.1,
        color="indigo",
    )
    if key:
        ax.quiverkey(
            q, 0.94, 0.95, 0.1, r"$0.1 \frac{m}{s}$", labelpos="E", coordinates="figure"
        )

def trend_quiver_func(ax, u, v, time, set_up, timescale = 1200, key = True):
    import scipy.stats
    slope_u = np.empty(np.shape(u[0,...]))
    slope_v = np.empty(np.shape(v[0,...]))
    
    for lat in range(0, u.shape[1], 30):
        for lon in range(0, u.shape[2], 30):
            slope_u[lat, lon] = scipy.stats.linregress(u[:, lat, lon], time).slope / timescale
            slope_v[lat, lon] = scipy.stats.linregress(v[:, lat, lon], time).slope / timescale
            if scipy.stats.linregress(u[:, lat, lon], time).pvalue < 0.05 and scipy.stats.linregress(v[:, lat, lon], time).pvalue < 0.05:
                q = ax.quiver(set_up["X"][lat, lon], set_up["Y"][lat, lon], slope_u[lat, lon], slope_v[lat, lon], color='red', width=0.005)
            elif scipy.stats.linregress(u[:, lat, lon], time).pvalue < 0.05 or scipy.stats.linregress(v[:, lat, lon], time).pvalue < 0.05:
                q = ax.quiver(set_up["X"][lat, lon], set_up["Y"][lat, lon], slope_u[lat, lon], slope_v[lat, lon], color='coral', width=0.005)
            else:
                q = ax.quiver(set_up["X"][lat, lon], set_up["Y"][lat, lon],slope_u[lat, lon], slope_v[lat, lon], color='lightgray', width=0.005)
    if key:
        ax.quiverkey(q, 0.935, 0.5, 3, r'$3 \frac{m}{s}cent.$', labelpos='E', coordinates='figure')
