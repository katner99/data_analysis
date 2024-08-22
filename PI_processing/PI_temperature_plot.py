from config_options import config_comparison
from .tools.plots import read_mask, pretty_labels, zoom_shelf
from directories_and_paths import output_path, grid_filepath
from mitgcm_python.grid import Grid
import matplotlib.pyplot as plt
from .tools.plots_2d import contour_func
import sys
import xarray as xr
import numpy as np


def comparison_plot(data, set_up, graph_params, graph_params_anom, experiment, file_out, save=True, show=False, linearity=True, zoom = None):
    fig, axs = plt.subplots(nrows=4, ncols=3, gridspec_kw={"hspace": 0.05, "wspace": 0.04}, figsize=(20, 20))
    axs = axs.flatten()

    main_figs = [0, 3, 6, 9]
    anom_figs = [4, 7, 8, 10, 11]
    axs[1].axis("off")
    #axs[2].axis("off")
    axs[5].axis("off")

    pval = xr.open_dataset("LENS_pvalue.nc")
    LC = pval.pvalue.values
    pval = xr.open_dataset("WIND_pvalue.nc")
    WC = pval.pvalue.values
    pval = xr.open_dataset("TEMP_pvalue.nc")
    TC = pval.pvalue.values
    pval = xr.open_dataset("WIND_LENS_pvalue.nc")
    WL = pval.pvalue.values
    pval = xr.open_dataset("TEMP_LENS_pvalue.nc")
    TL = pval.pvalue.values
    
    
    for count, i in enumerate(main_figs):
        if i == 9:
            hide_ticks_x = False
        else:
            hide_ticks_x = True
        main_plt = contour_func(axs[i], data[count], set_up, graph_params, hide_ticks_x, hide_ticks_y = False)
        zoom_shelf(axs[i], zoom)

        axs[i].text(
            0.5,
            0.9,
            experiment[count],
            FontSize = 19,
            horizontalalignment="center",
            transform=axs[i].transAxes,
            bbox=dict(facecolor="white", alpha=0.9),
        )
        pretty_labels(axs[i])

    hide_ticks_y = True
    anom_plt = contour_func(axs[4], data[1] - data[0], set_up, graph_params, hide_ticks_x = True, hide_ticks_y = hide_ticks_y)
    axs[4].contourf(set_up["X"], set_up["Y"], LC, levels=[-np.inf, 0.05], colors='none', hatches=['//'], alpha=0)
    contour_func(axs[7], data[2] - data[0], set_up, graph_params, hide_ticks_x = True, hide_ticks_y = hide_ticks_y)
    axs[7].contourf(set_up["X"], set_up["Y"], TC, levels=[-np.inf, 0.05], colors='none', hatches=['//'], alpha=0)
    contour_func(axs[10], data[3] - data[0], set_up, graph_params, hide_ticks_x = False, hide_ticks_y = hide_ticks_y)
    axs[10].contourf(set_up["X"], set_up["Y"], WC, levels=[-np.inf, 0.05], colors='none', hatches=['//'], alpha=0)
    contour_func(axs[8], data[1] - data[2], set_up, graph_params, hide_ticks_x = True, hide_ticks_y = hide_ticks_y)
    axs[8].contourf(set_up["X"], set_up["Y"], TL, levels=[-np.inf, 0.05], colors='none', hatches=['//'], alpha=0)
    contour_func(axs[11], data[1] - data[3], set_up, graph_params, hide_ticks_x = False, hide_ticks_y = hide_ticks_y)
    axs[11].contourf(set_up["X"], set_up["Y"], WL, levels=[-np.inf, 0.05], colors='none', hatches=['//'], alpha=0)

    if linearity:
        residual = data[0] + data[1] - data[2] - data[3]
        position = 2
        cs = contour_func(axs[position], residual, set_up, graph_params)

        pretty_labels(axs[position])
        zoom_shelf(axs[position], zoom)

        axs[position].text(
            0.5,
            0.9,
            "Non-linearity",
            FontSize = 18,
            horizontalalignment="center",
            transform=axs[position].transAxes,
            bbox=dict(facecolor="white", alpha=0.9),
        )

    for count, i in enumerate([4, 7, 10]):
        zoom_shelf(axs[i], zoom)
        pretty_labels(axs[i])
        
        axs[i].text(
            0.5,
            0.9,
            experiment[count+1]+" - "+experiment[0],
            FontSize = 18,
            horizontalalignment="center",
            transform=axs[i].transAxes,
            bbox=dict(facecolor="white", alpha=0.9),
        )

    for count, i in enumerate([8, 11]):
        zoom_shelf(axs[i], zoom)
        pretty_labels(axs[i])
        
        axs[i].text(
            0.5,
            0.9,
            experiment[1]+" - "+experiment[count+2],
            FontSize = 18,
            horizontalalignment="center",
            transform=axs[i].transAxes,
            bbox=dict(facecolor="white", alpha=0.9),
        )

    ticks=np.arange(graph_params["low_val"], graph_params["high_val"]+0.1, 0.5)
    cbar_ax = fig.add_axes([0.1, 0.04, 0.8, 0.02])  # Bottom color bar
    cbar = plt.colorbar(main_plt, cax=cbar_ax, orientation='horizontal')
    cbar.set_ticks(ticks)
    cbar.set_label('Potential Temperature (degC)', fontsize = 15)
    cbar.ax.xaxis.set_ticks_position('bottom')
    cbar.ax.tick_params(labelsize=16)

    fig.suptitle(graph_params["title"], fontsize=20, weight="bold")

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
    show = False
    period = "2070-2100"
    filepaths = [
        output_path + "average_" + ens + "_" + period + ".nc"
        for ens in ["CTRL", "LENS", "TEMP", "WIND"]
    ]
    file_out = f"mega_comparison{var}_{period}.png"

    for filepath in filepaths:
        try:
            open(filepath)
        except FileNotFoundError:
            sys.exit(f"Stopped - Could not find input file {filepath}")

    experiment = ["NONE", "ALL", "THERMO", "WIND"]

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
