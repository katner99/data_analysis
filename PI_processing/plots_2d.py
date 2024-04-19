import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import numpy as np
from plots import zoom_shelf

def comparison(data, set_up, graph_params, graph_params_anom, experiment, title, file_out, save=True, show=False, linearity=False, residual = None, zoom = None):
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
    
    fig, axs = plt.subplots(nrows=total, ncols=total, gridspec_kw={"hspace": 0.05, "wspace": 0.04}, figsize=(15, 15))
    axs = axs.flatten()


    if linearity:
        position = total-1
        cs = contour_func(axs[position], residual, set_up, graph_params_anom)
        zoom_shelf(axs[position], zoom)

        fig.colorbar(cs, ax=axs[position], ticks=np.arange(-10, 10.1, 10))
        axs[position].set_title("Residual", fontsize=graph_params["font_size"], weight="bold")

        lon_ticks = axs[position].get_xticks() - 360
        lon_labels = []
        for x in lon_ticks:
            lon_labels.append(lon_label(x,2))
        axs[position].set_xticklabels(lon_labels)
        axs[position].tick_params(axis='x', labelrotation=45)
        
        lat_ticks = axs[position].get_yticks()
        lat_labels = []
        for y in lat_ticks:
            lat_labels.append(lat_label(y,2))
        axs[position].set_yticklabels(lat_labels)
    
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
        zoom_shelf(axs[diagonal], zoom)
        axs[diagonal].set_title(experiment[i], fontsize=graph_params["font_size"], weight="bold")
        
        if hide_ticks_y == False:
            for y in lat_ticks:
                lat_labels.append(lat_label(y,2))
            axs[diagonal].set_yticklabels(lat_labels)

        if hide_ticks_x == False:
            lon_ticks = axs[diagonal].get_xticks() - 360
            lon_labels = []
            for x in lon_ticks:
                lon_labels.append(lon_label(x,2))
            axs[diagonal].set_xticklabels(lon_labels)
            axs[diagonal].tick_params(axis='x', labelrotation=45)
            axs[diagonal].set_xticklabels(lon_labels[:-1] + [""])

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
            
            if anomaly % total == 0:
                lat_ticks = axs[anomaly].get_yticks()
                lat_labels = []
                for y in lat_ticks:
                    lat_labels.append(lat_label(y,2))
                axs[anomaly].set_yticklabels(lat_labels)
            
            if hide_ticks_x == False:
                lon_ticks = axs[anomaly].get_xticks() - 360
                lon_labels = []
                for x in lon_ticks:
                    lon_labels.append(lon_label(x,2))
                axs[anomaly].set_xticklabels(lon_labels[:-1] + [""])
                axs[anomaly].tick_params(axis='x', labelrotation=45)

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
    cbar.ax.yaxis.set_ticks_position('left')

    ticks=np.arange(graph_params_anom["low_val"], graph_params_anom["high_val"]+0.1, graph_params_anom["interval"])
    cbar_ax = fig.add_axes([0.05, 0.1, 0.02, 0.4])
    cbar = plt.colorbar(cs_anom, cax=cbar_ax, orientation='vertical')
    cbar.set_ticks(ticks)
    cbar.ax.yaxis.set_ticks_position('left')

    fig.suptitle(title, fontsize=16)

    # save figure
    if save == True:
        fig.savefig(file_out)

    # show figure
    if show == True:
        plt.show()




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
        scale=0.8,
        width=0.005,
        color="indigo",
    )
    if key:
        ax.quiverkey(
            q, 0.94, 0.935, 0.1, r"$0.1 m/s$", labelpos="E", coordinates="figure"
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
                q = ax.quiver(set_up["X"][lat, lon], set_up["Y"][lat, lon], slope_u[lat, lon], slope_v[lat, lon], color='white', scale = 25, width=0.008)
            elif scipy.stats.linregress(u[:, lat, lon], time).pvalue < 0.05 or scipy.stats.linregress(v[:, lat, lon], time).pvalue < 0.05:
                q = ax.quiver(set_up["X"][lat, lon], set_up["Y"][lat, lon], slope_u[lat, lon], slope_v[lat, lon], color='grey', scale = 25, width=0.008)
            else:
                q = ax.quiver(set_up["X"][lat, lon], set_up["Y"][lat, lon],slope_u[lat, lon], slope_v[lat, lon], color='grey', scale = 25, width=0.003)
    if key:
        ax.quiverkey(q, 0.938, 0.495, 3, r'$3 \frac{m/s}{cent.}$', color='k', labelpos='E', coordinates='figure')
