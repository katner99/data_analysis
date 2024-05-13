from config_options import config_comparison
from plots import read_mask, pretty_labels, zoom_shelf
from directories_and_paths import output_path, grid_filepath
from mitgcm_python.grid import Grid
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from plots_2d import contour_func
import sys
import xarray as xr
import numpy as np


def comparison_plot(data, set_up, graph_params, graph_params_anom, experiment, file_out, save=True, show=False, linearity=False, zoom = None):
    fig, axs = plt.subplots(nrows=4, ncols=2, gridspec_kw={"hspace": 0.05, "wspace": 0.04}, figsize=(14, 20))
    axs = axs.flatten()

    main_figs = [0, 2, 4, 6]
    anom_figs = [3, 5, 7]
    axs[1].axis("off")

    pval = xr.open_dataset("CL_pval.nc")
    lens = pval.pvalue.values
    pval = xr.open_dataset("CW_pval.nc")
    wind = pval.pvalue.values
    pval = xr.open_dataset("CT_pval.nc")
    temp = pval.pvalue.values
    
    for count, i in enumerate(main_figs):
        if i == 6:
            hide_ticks_x = False
        else:
            hide_ticks_x = True
        main_plt = contour_func(axs[i], data[count], set_up, graph_params, hide_ticks_x, hide_ticks_y = False)
        zoom_shelf(axs[i], zoom)

        axs[i].text(
            0.5,
            0.85,
            experiment[count],
            FontSize = 16,
            horizontalalignment="center",
            transform=axs[i].transAxes,
            bbox=dict(facecolor="white", alpha=0.9),
        )
        pretty_labels(axs[i])

    hide_ticks_y = True
    anom_plt = contour_func(axs[3], data[1] - data[0], set_up, graph_params_anom, hide_ticks_x = True, hide_ticks_y = hide_ticks_y)
    axs[3].contourf(set_up["X"], set_up["Y"], lens, levels=[-np.inf, 0.05], colors='none', hatches=['....'], alpha=0)
    contour_func(axs[5], data[2] - data[0], set_up, graph_params_anom, hide_ticks_x = True, hide_ticks_y = hide_ticks_y)
    axs[5].contourf(set_up["X"], set_up["Y"], wind, levels=[-np.inf, 0.05], colors='none', hatches=['....'], alpha=0)
    contour_func(axs[7], data[3] - data[0], set_up, graph_params_anom, hide_ticks_x = False, hide_ticks_y = hide_ticks_y)
    axs[7].contourf(set_up["X"], set_up["Y"], temp, levels=[-np.inf, 0.05], colors='none', hatches=['....'], alpha=0)
                     
    for count, i in enumerate(anom_figs):
        zoom_shelf(axs[i], zoom)
        pretty_labels(axs[i])
        
        axs[i].text(
            0.5,
            0.85,
            experiment[count+1]+" - "+experiment[0],
            FontSize = 16,
            horizontalalignment="center",
            transform=axs[i].transAxes,
            bbox=dict(facecolor="white", alpha=0.9),
        )

    ticks=np.arange(graph_params["low_val"], graph_params["high_val"]+0.1, graph_params["interval"])
    cbar_ax = fig.add_axes([0.05, 0.04, 0.425, 0.02])  # Bottom color bar
    cbar = plt.colorbar(main_plt, cax=cbar_ax, orientation='horizontal')
    cbar.set_ticks(ticks)
    cbar.ax.xaxis.set_ticks_position('bottom')

    ticks=np.arange(graph_params_anom["low_val"], graph_params_anom["high_val"]+0.1, graph_params_anom["interval"])
    cbar_ax = fig.add_axes([0.525, 0.04, 0.425, 0.02])  # Bottom color bar
    cbar = plt.colorbar(anom_plt, cax=cbar_ax, orientation='horizontal')
    cbar.set_ticks(ticks)
    cbar.ax.xaxis.set_ticks_position('bottom')

    fig.suptitle(graph_params["title"], fontsize=20)

    plt.subplots_adjust(top=0.95)

    # save figure
    if save == True:
        fig.savefig(file_out, bbox_inches='tight')

    # show figure
    if show == True:
        plt.show()

def main():
    # set up the variables you need
    var = "THETA"
    save = True
    show = True
    period = "2070-2100"
    filepaths = [
        output_path + "average_" + ens + "_" + period + ".nc"
        for ens in ["CTRL", "LENS", "WIND", "TEMP"]
    ]
    file_out = f"mega_comparison{var}_{period}.png"

    for filepath in filepaths:
        try:
            open(filepath)
        except FileNotFoundError:
            sys.exit(f"Stopped - Could not find input file {filepath}")

    experiment = ["pre-industrial forcing", "high emissions forcing", "wind forcing", "thermodynamic forcing"]

    input_data = [
        xr.open_dataset(filepath, decode_times=False) for filepath in filepaths
    ]

    grid = Grid(grid_filepath)

    set_up = read_mask(input_data[0])
    [data, graph_params, graph_params_anom] = config_comparison(var, input_data, grid, period = period)

    comparison_plot(data,
        set_up,
        graph_params,
        graph_params_anom,
        experiment,
        file_out,
        save,
        show,
        True,
        zoom = "cont_shelf",
        )

if __name__ == "__main__":
    main()  # run the program
