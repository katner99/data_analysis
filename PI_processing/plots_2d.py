import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import numpy as np


def contour_func(fig, ax, data, set_up, graph_params, hide_ticks=True):
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
    ticks = np.arange(
        graph_params["low_val"],
        graph_params["high_val"] + 0.1,
        graph_params["interval"],
    )
    cbar = fig.colorbar(cs, ax=ax, aspect=40)
    cbar.set_ticks(ticks)
    if hide_ticks:
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
    ax.set_aspect("auto", adjustable="box")
    return ax


def quiver_func(ax, u, v, lat, lon, chunk):
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
    ax.quiverkey(
        q, 0.94, 0.98, 0.1, r"$0.1 \frac{m}{s}$", labelpos="E", coordinates="figure"
    )
