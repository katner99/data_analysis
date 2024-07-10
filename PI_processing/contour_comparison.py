from config_options import config_comparison
from plots import read_mask
from plots_2d import comparison
from directories_and_paths import output_path, grid_filepath
from mitgcm_python.grid import Grid

import sys
import xarray as xr

def read_data(var, var_name = None):
    """ read in analysed data """
    if var in ["trend", "pvalue"]:
        filename = f"{var_name}_spread.nc"
        file_out = f"mega_comparison_{var}_{var_name}.png"
    elif var == "mean":
        filename = f"{var_name}.nc"
        file_out = f"mega_comparison_{var}_{var_name}_mean.png"
    filepaths = [
        f"{output_path}{exp}_files_temp/{filename}"
        for exp in ["CTRL", "LENS", "WIND", "TEMP"]
    ]
    return filepaths, file_out

def main():
    """ plot function to create the ugly 4x4 beast of a plot """
    var = "mean"             # either set to the parameter you will be using or set to mean or trend (tells it where to read the values in)
    save = True              # save figure
    show = False             # print the figure on script completion

    # will follow the first option if going through analysed data, will otherwise get the climatology data for the corresponding var
    if var in ["trend", "pvalue", "mean"]:
        var_name = "oveFWflx"
        [filepaths, file_out] = read_data(var, var_name)
    else:
        period = "2070-2100"
        filepaths = [
            f"{output_path}average_{ens}_{period}.nc"
            for ens in ["CTRL", "LENS", "WIND", "TEMP"]
        ]
        file_out = f"mega_comparison{var}_{period}.png"

    for filepath in filepaths:
        try:
            open(filepath)
        except FileNotFoundError:
            sys.exit(f"Stopped - Could not find input file {filepath}")

    # experiment neames to add to the plot
    experiment = ["NONE", "ALL", "WIND", "THERMO"]

    input_data = [
        xr.open_dataset(filepath, decode_times=False) for filepath in filepaths
    ]

    grid = Grid(grid_filepath)

    # set up the data to be plotted
    if var in ["trend", "pvalue", "mean"]:
        set_up_data = xr.open_dataset(f"{output_path}average_CTRL_1920-1950.nc", decode_times=False)
        set_up = read_mask(set_up_data)
        [data, graph_params, graph_params_anom] = config_comparison(var, input_data, grid, var_name = var_name)
    else:
        set_up = read_mask(input_data[0])
        [data, graph_params, graph_params_anom] = config_comparison(var, input_data, grid, period = period)

    # plot the data
    comparison(
        data,
        set_up,
        graph_params,
        graph_params_anom,
        experiment,
        file_out,
        save,
        show,
        False,       # option to calculate the linearity in the plot
        "cont_shelf" # option to zoom into a certain area, choose between "cont_shelf" and "ice_shelf" or leave blank
    )

if __name__ == "__main__":
    main() 
