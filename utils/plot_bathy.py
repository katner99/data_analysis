import matplotlib.pyplot as plt

import numpy as np
import xarray as xr

from plots_2d import plot_contour
from directories_and_paths import output_path
from plots import read_mask


def read_bathy(filepath):
    input_data = xr.open_dataset(filepath, decode_times=False)
    set_up = read_mask(input_data)
    graph_params = {
        "font_size": 12,
        "low_val": 0,
        "high_val": 1500,
        "interval": 500,
        "step": 15,
        "color_scheme": "bone_r",
    }
    return set_up, graph_params

def plot_bathy(
    graph_params,
    set_up,
):
    fig, axs = plt.subplots(
        figsize=(10, 8), sharex=True
    )
    plt.subplots_adjust(wspace=0.05, hspace=0.05)

    cs_bathy = plot_contour(
        axs,
        set_up["depth"],
        set_up,
        graph_params,
        "Bathymetry (m)",
        hide_ticks_x=False,
        hide_ticks_y=False,
    )

    ticks = np.arange(
        graph_params["low_val"],
        graph_params["high_val"] + 0.1,
        graph_params["interval"],
    )
    cbar = plt.colorbar(cs_bathy)
    cbar.set_ticks(ticks)

    # set various colors
    axs.spines['bottom'].set_color('white')
    axs.spines['top'].set_color('white') 
    axs.spines['right'].set_color('white')
    axs.spines['left'].set_color('white')
    axs.xaxis.label.set_color('white')
    axs.yaxis.label.set_color('white')
    axs.tick_params(colors='white', which='both')  # 'both' refers to minor and major axes

    fig.savefig("bathymetry_plot.png", transparent=True)

    plt.show()


def main():
    [set_up, graph_params] = read_bathy(
        output_path + "average_CTRL_1920-1950.nc"
    )

    plot_bathy(
        graph_params,
        set_up
    )


if __name__ == "__main__":
    main()