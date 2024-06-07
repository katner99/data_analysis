from config_options import config_comparison
from plots import read_mask
from plots_2d import comparison
from directories_and_paths import output_path, grid_filepath
from mitgcm_python.grid import Grid

import sys
import xarray as xr


def main():
    # set up the variables you need
    var = "trend"
    save = True
    show = True

    if var in ["trend", "pvalue"]:
        var_name = "SALT"
        filepaths = [
            f"{output_path}{exp}_files_temp/{var_name}_trend.nc"
            for exp in ["CTRL", "LENS", "WIND", "TEMP"]
        ]
        file_out = f"mega_comparison_{var}_{var_name}.png"
    else:
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

    experiment = ["pre-industrial", "RCP 8.5", "wind f.", "thermo. f."]

    input_data = [
        xr.open_dataset(filepath, decode_times=False) for filepath in filepaths
    ]

    grid = Grid(grid_filepath)

    if var in ["trend", "pvalue"]:
        set_up_data = xr.open_dataset(f"{output_path}average_CTRL_1920-1950.nc", decode_times=False)
        [data, graph_params, graph_params_anom] = config_comparison(var, input_data, grid, var_name = var_name)

    


if __name__ == "__main__":
    main()  # run the program
