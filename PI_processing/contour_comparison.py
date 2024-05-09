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
    period = "2070-2100"

    if var in ["trend", "pvalue"]:
        val = "oceFWflx"
        filepaths = [
            f"{output_path}{exp}_files_temp/{val}_trend.nc"
            for exp in ["CTRL", "LENS", "WIND", "TEMP"]
        ]
    else:
        # load up the file paths for the monster
        filepaths = [
            output_path + "average_" + ens + "_" + period + ".nc"
            for ens in ["CTRL", "LENS", "WIND", "TEMP"]
        ]
    # check if the input files exist
    for filepath in filepaths:
        try:
            open(filepath)
        except FileNotFoundError:
            sys.exit(f"Stopped - Could not find input file {filepath}")

    experiment = ["pre-industrial", "RCP 8.5", "wind f.", "thermo. f."]

    # load up the input data
    input_data = [
        xr.open_dataset(filepath, decode_times=False) for filepath in filepaths
    ]

    # read in the general variables (these should be the same between the ensembles
    if var in ["trend", "pvalue"]:
        set_up_data = xr.open_dataset(f"{output_path}average_CTRL_1920-1950.nc", decode_times=False)
        set_up = read_mask(set_up_data)
    else:
        set_up = read_mask(input_data[0])

    grid = Grid(grid_filepath)

    [data, graph_params, graph_params_anom] = config_comparison(
        var, input_data, grid, period
    )

    file_out = "mega_comparison" + var + "_" + val + ".png"

    residual = data[0] + data[1] - data[2] - data[3]

    comparison(
        data,
        set_up,
        graph_params,
        graph_params_anom,
        experiment,
        graph_params["title"],
        file_out,
        save,
        show,
        True,
        residual,
        pvalue = graph_params["pvalue"]
    )


if __name__ == "__main__":
    main()  # run the program
