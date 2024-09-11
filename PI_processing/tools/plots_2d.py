import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import numpy as np
from .plots import zoom_shelf, pretty_labels
from mitgcm_python.plot_utils.labels import lat_label, lon_label

def comparison(data, set_up, graph_params, graph_params_anom, experiment, file_out, save=True, show=False, linearity=False, zoom = None):
    """
    Generate a grid of subplots for visual comparison of multiple datasets.

    Parameters:
        data (list): A list of datasets to be plotted.
        set_up: (object): An object containing setup information.
        graph_params (dict): Parameters for main graphs.
        graph_params_anom (dict): Parameters for anomaly graphs.
        experiment (list): A list of experiment names corresponding to each dataset.
        title (str): Title for the entire figure.
        file_out (str): File path for saving the figure.
        save (bool, optional): Whether to save the figure. Defaults to True.
        show (bool, optional): Whether to display the figure. Defaults to False.
        linearity (bool, optional): Whether to include a residual plot. Defaults to False.
        residual: (object): An object containing residual plot data.
        zoom: (object): An object containing zoom parameters.

    Returns:
        None
    """
    from mitgcm_python.plot_utils.labels import lat_label, lon_label

    total = len(data)
    
    fig, axs = plt.subplots(nrows=total, ncols=total, gridspec_kw={"hspace": 0.05, "wspace": 0.04}, figsize=(15, 10))
    axs = axs.flatten()

    if linearity:
        residual = data[0] + data[1] - data[2] - data[3]
        position = total-1
        cs = contour_func(axs[position], residual, set_up, graph_params_anom)
        zoom_shelf(axs[position], zoom)

        fig.colorbar(cs, ax=axs[position], ticks=np.arange(-1, 1.1, 1))
        axs[position].set_title("Residual", fontsize=graph_params["font_size"], weight="bold")
        pretty_labels(axs[position])
    
    for i in range(total):
        # MAIN GRAPH
        diagonal = (total+1)*i
        if diagonal == 0:
            hide_ticks_y = False
        else:
            hide_ticks_y = True
        
        if diagonal == ((total*total)-1):
            hide_ticks_x = False
        else:
            hide_ticks_x = True

        cs_diag = contour_func(axs[diagonal], data[i], set_up, graph_params, hide_ticks_x, hide_ticks_y)
        if graph_params.get("pvalue", None):
            axs[diagonal].contourf(set_up["X"], set_up["Y"], graph_params["pvalue"][i], levels=[-np.inf, 0.05], colors='none', hatches=['....'], alpha=0)
        
        zoom_shelf(axs[diagonal], zoom)
        axs[diagonal].set_title(experiment[i], fontsize=graph_params["font_size"], weight="bold")
        pretty_labels(axs[diagonal])

        # ANOMALY
        for j in range(i+1, total):
            anomaly = (total*j)+i

            if anomaly % total == 0:
                hide_ticks_y = False
            else:
                hide_ticks_y = True

            if anomaly >= (total * (total - 1)):
                hide_ticks_x = False
            else:
                hide_ticks_x = True

            cs_anom = contour_func(axs[anomaly], data[j] - data[i], set_up, graph_params_anom, hide_ticks_x, hide_ticks_y)
            
            zoom_shelf(axs[anomaly], zoom)
            pretty_labels(axs[anomaly])
            
            axs[anomaly].text(
                0.5,
                0.85,
                experiment[j]+" - "+experiment[i],
                horizontalalignment="center",
                transform=axs[anomaly].transAxes,
                bbox=dict(facecolor="white", alpha=0.9),
            )
            
        # discarded graphs
        for k in range(i):
            axs[i + (total*k)].axis("off")
    
    ticks=np.arange(graph_params["low_val"], graph_params["high_val"]+0.1, graph_params["interval"])
    cbar_ax = fig.add_axes([0.05, 0.525, 0.02, 0.4])
    cbar = plt.colorbar(cs_diag, cax=cbar_ax, orientation='vertical')
    cbar.set_ticks(ticks)
    #cbar.set_label("Potential Temperature (°C)")
    cbar.ax.yaxis.set_ticks_position('left')

    ticks=np.arange(graph_params_anom["low_val"], graph_params_anom["high_val"]+0.1, graph_params_anom["interval"])
    cbar_ax = fig.add_axes([0.05, 0.1, 0.02, 0.4])
    cbar = plt.colorbar(cs_anom, cax=cbar_ax, orientation='vertical')
    cbar.set_ticks(ticks)
    #cbar.set_label("Potential Temperature Anomaly (°C)")
    cbar.ax.yaxis.set_ticks_position('left')

    fig.suptitle(graph_params["title"], fontsize=16)

    # save figure
    if save == True:
        fig.savefig(file_out, bbox_inches='tight')

    # show figure
    if show == True:
        plt.show()


def contour_func(ax, data, set_up, graph_params, hide_ticks_x=True, hide_ticks_y=True, option = "Greys"):
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
        set_up["X"], set_up["Y"], set_up["mask"], 2, cmap=option, linestyles="dashed"
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
        scale=1,
        width=0.008,
        headwidth=2.5,
        #headlength=3,
        color="indigo",
    )
    if key:
        ax.quiverkey(
            q, 0.94, 0.935, 0.1, r"$0.1 m/s$", labelpos="E", coordinates="figure"
        )

def trend_quiver_func(ax, u, v, time, set_up, key = True):
    import scipy.stats
    slope_u = np.full(np.shape(u[0,...]), np.nan)
    slope_v = np.full(np.shape(v[0,...]), np.nan)
    
    for lat in range(0, u.shape[1], 30):
        for lon in range(0, u.shape[2], 30):
            slope_u[lat, lon] = scipy.stats.linregress(time, u[:, lat, lon]).slope 
            slope_v[lat, lon] = scipy.stats.linregress(time, v[:, lat, lon]).slope
            print(slope_u[lat,lon], slope_v[lat, lon])
            if abs(slope_v[lat, lon]) or abs(slope_u[lat, lon]) > 10:
                slope_v[lat, lon] == 0
                slope_u[lat, lon] == 0
            if scipy.stats.linregress(u[:, lat, lon], time).pvalue < 0.05 and scipy.stats.linregress(v[:, lat, lon], time).pvalue < 0.05:
                q = ax.quiver(set_up["X"][lat, lon], set_up["Y"][lat, lon], slope_u[lat, lon], slope_v[lat, lon], color='white', edgecolor='k', linewidth = 1, scale = 10, width=0.008)
            elif scipy.stats.linregress(u[:, lat, lon], time).pvalue < 0.05 or scipy.stats.linregress(v[:, lat, lon], time).pvalue < 0.05:
                q = ax.quiver(set_up["X"][lat, lon], set_up["Y"][lat, lon], slope_u[lat, lon], slope_v[lat, lon], color='grey', edgecolor='k', linewidth = 1, scale = 10, width=0.008)
            else:
                q = ax.quiver(set_up["X"][lat, lon], set_up["Y"][lat, lon],slope_u[lat, lon], slope_v[lat, lon], color='black', edgecolor='k', linewidth = 1, scale = 10, width=0.001)
    print(np.nanmax(slope_v), np.nanmin(slope_v))
    if key:
        ax.quiverkey(q, 0.938, 0.495, 2, r'$2 \frac{m/s}{cent.}$', color='k', labelpos='E', coordinates='figure')

def plot_contour(
    ax, data, set_up, graph_params, title, hide_ticks_x=True, hide_ticks_y=True, option="Greys"
):
    """
    Plots a contour map on a given axis with specified data and parameters.

    Parameters:
    ax (matplotlib.axes.Axes): The axis on which to plot the contour map.
    data (np.ndarray or xr.DataArray): The data to be contoured.
    set_up (np.ndarray or xr.DataArray): The setup mask for the data.
    graph_params (dict): A dictionary containing parameters for creating the contour plot, 
                         including font size, value range, interval, step, and color scheme.
    title (str): The title of the plot.
    hide_ticks_x (bool, optional): If True, hides the x-axis ticks. Default is True.
    hide_ticks_y (bool, optional): If True, hides the y-axis ticks. Default is True.

    Returns:
    cs (QuadContourSet): The contour set created by the plot.
    """
    cs = contour_func(ax, data, set_up, graph_params, hide_ticks_x, hide_ticks_y, option=option)
    ax.set_xlim([230, 262])
    ax.set_ylim([-75.6, -68])
    ax.text(
        0.5,
        0.9,
        title,
        horizontalalignment="center",
        transform=ax.transAxes,
        bbox=dict(facecolor="white", alpha=0.9),
    )
    if hide_ticks_x == False:
        ax.locator_params(axis="x", nbins=6)
        lon_ticks = ax.get_xticks() - 360
        lon_labels = []
        for x in lon_ticks:
            lon_labels.append(lon_label(x, 2))
        ax.set_xticklabels(lon_labels[:-1] + [""])
        ax.tick_params(axis="x", labelrotation=45)

    if hide_ticks_y == False:
        ax.locator_params(axis="y", nbins=5)
        lat_ticks = ax.get_yticks()
        lat_labels = []
        for y in lat_ticks:
            lat_labels.append(lat_label(y, 2))
        ax.set_yticklabels(lat_labels)

    return cs
