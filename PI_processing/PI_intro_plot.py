import matplotlib

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from plots_2d import contour_func, quiver_func
from directories_and_paths import output_path
from plots import read_mask, read_u_and_v


def main():
    filepath = output_path + "average_CTRL_1920-1950.nc"
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

    fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
    plt.tight_layout()
    axs = axs.flatten()
    contour_func(fig, axs[1], speed, set_up, graph_params)
    quiver_func(axs[1], u, v, set_up["lat"], set_up["lon"], chunk)
    axs[1].set_xlim([230, 260])
    axs[1].set_ylim([-75.5, -68])
    axs[1].text(
        0.5,
        0.9,
        "Bottom Velocity (m/s)",
        horizontalalignment="center",
        transform=axs[1].transAxes,
        bbox=dict(facecolor="white", alpha=0.9),
    )
    plt.show()


if __name__ == "__main__":
    main()
