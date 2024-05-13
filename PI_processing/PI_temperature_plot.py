import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from config_options import config_comparison
from plots_2d import contour_func
from plots import read_mask, pretty_labels, zoom_shelf
from directories_and_paths import output_path, grid_filepath
from mitgcm_python.grid import Grid

def comparison_plot(data, set_up, graph_params, graph_params_anom, experiment, file_out, save=True, show=False, zoom=None):
    fig, axs = plt.subplots(nrows=4, ncols=2, gridspec_kw={"hspace": 0.05, "wspace": 0.04}, figsize=(14, 20))
    axs = axs.flatten()

    pvals = [xr.open_dataset(f"{letter}.nc").pvalue.values for letter in ["CL", "CW", "CT"]]

    for count, i in enumerate([0, 2, 4, 6]):
        main_plt = contour_func(axs[i], data[count], set_up, graph_params, i != 6, hide_ticks_y=False)
        zoom_shelf(axs[i], zoom)
        axs[i].text(0.5, 0.85, experiment[count], FontSize=16, horizontalalignment="center", transform=axs[i].transAxes, bbox=dict(facecolor="white", alpha=0.9))
        pretty_labels(axs[i])

    for count, i in enumerate([3, 5, 7]):
        anom_plt = contour_func(axs[i], data[count] - data[0], set_up, graph_params_anom, hide_ticks_x=True, hide_ticks_y=True)
        zoom_shelf(axs[i], zoom)
        axs[i].text(0.5, 0.85, f"{experiment[count+1]} - {experiment[0]}", FontSize=16, horizontalalignment="center", transform=axs[i].transAxes, bbox=dict(facecolor="white", alpha=0.9))

        # Overlay p-values
        axs[i].contourf(set_up["X"], set_up["Y"], pvals[count], levels=[-np.inf, 0.05], colors='none', hatches=['....'], alpha=0)

    cbar_ax_params = [0.05, 0.04, 0.425, 0.02]
    for i, graph_param, cbar_ax_param in zip([main_plt, anom_plt], [graph_params, graph_params_anom], cbar_ax_params):
        ticks = np.arange(graph_param["low_val"], graph_param["high_val"] + 0.1, graph_param["interval"])
        cbar_ax = fig.add_axes(cbar_ax_param)  # Bottom color bar
        cbar = plt.colorbar(i, cax=cbar_ax, orientation='horizontal')
        cbar.set_ticks(ticks)
        cbar.ax.xaxis.set_ticks_position('bottom')

    fig.suptitle(graph_params["title"], fontsize=20)
    plt.subplots_adjust(top=0.95)

    if save:
        fig.savefig(file_out, bbox_inches='tight')
    if show:
        plt.show()

def main():
    var = "THETA"
    save = True
    show = True
    period = "2070-2100"
    filepaths = [output_path + f"average_{ens}_{period}.nc" for ens in ["CTRL", "LENS", "WIND", "TEMP"]]
    file_out = f"mega_comparison{var}_{period}.png"

    for filepath in filepaths:
        if not os.path.exists(filepath):
            sys.exit(f"Stopped - Could not find input file {filepath}")

    experiment = ["pre-industrial forcing", "high emissions forcing", "wind forcing", "thermodynamic forcing"]
    input_data = [xr.open_dataset(filepath, decode_times=False) for filepath in filepaths]
    grid = Grid(grid_filepath)
    set_up = read_mask(input_data[0])
    data, graph_params, graph_params_anom = config_comparison(var, input_data, grid, period=period)

    comparison_plot(data, set_up, graph_params, graph_params_anom, experiment, file_out, save, show, zoom="cont_shelf")

if __name__ == "__main__":
    main()
