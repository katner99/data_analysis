from funcs import read_variable, find_nearest
from plots import read_mask
from plots_2d import comparison
from directories_and_paths import output_path, grid_filepath
from mitgcm_python.grid import Grid

import sys
import numpy as np
import xarray as xr
       
def main():
    """
    Main function that reads the command line arguments, checks them for validity,
    loads the input data and calls the appropriate plotting function.
    """
        
    # set up the variables you need
    var = "THETA"
    save = True
    show = True
    period = "2070-2100"
    
    # load up the file paths for the monster, needed 4
    filepaths = [output_path + "average_" + ens + "_" + period + ".nc" for ens in ["CTRL", "LENS", "WIND", "TEMP"]]
    # check if the input files exist
    for filepath in filepaths:
        try:
            open(filepath)
        except FileNotFoundError:
            sys.exit(f"Stopped - Could not find input file {filepath}")

    experiment = ["pre-industrial", "RCP 8.5", "wind f.", "thermo. f."]

    # load up the input data
    input_data = [xr.open_dataset(filepath, decode_times = False) for filepath in filepaths]

    # read in the general variables (these should be the same between the ensembles
    set_up = read_mask(input_data[0])

    grid = Grid(grid_filepath)

    # temperature
    if var == "THETA":
        depth_range = [find_nearest(input_data[0]["Z"].values, -200), find_nearest(input_data[0]["Z"].values, -700)]
        data = [read_variable(input, var, grid, depth_range) for input in input_data]
        color_scheme = "coolwarm"
        anom = 1.5
        min_val = -2
        max_val =  2.1
        title = f"Average temperature between 200 and 700m {period}"
        interval = 1
        interval_anom = 0.5

    # sea ice tickness
    elif var in ["SIheff", "oceFWflx", "SIfwmelt", "SIfwfrz", "EXFvwind", "oceQnet"]:
        data = [read_variable(input, var, grid)*3600*24*365/1000 for input in input_data]
        #color_scheme = "YlGnBu_r"
        #color_scheme = "seismic"
        color_scheme = "PRGn_r"
        anom = 1.5
        min_val = -5
        max_val = 5
        title = f"Freshwater fluxes m/yr {period}"

    # salinity
    elif var == "SALT":
        depth_range = [find_nearest(input_data[0]["Z"].values, -200), find_nearest(input_data[0]["Z"].values, -700)]
        data = [read_variable(input, var, grid, depth_range) for input in input_data]
        color_scheme = "PRGn_r"
        anom = 0.25
        min_val = int(np.min(data[0]))
        max_val =  int(np.max(data[0]))+1
        
    elif var == "SHIfwFlx":
        data = [-read_variable(input, var, grid)*3600*24*30*10**(-3) for input in input_data]
        color_scheme = "rainbow"
        anom = 25
        min_val = 0
        max_val = 76
        title = f"ice shelf basal melt rate ({period}) (m.w.e./yr)"

    # graph parameters 1.333112e-05 -1.0214189e-06
    graph_params = {
        "font_size": 12,
        "low_val": min_val,
        "high_val": max_val,
        "interval": interval,
        "color_scheme": color_scheme,
        "step": 15,
    }

    graph_params_anom = {
        "font_size": 12,
        "low_val": -anom,
        "high_val": anom,
        "interval": interval_anom,
        "color_scheme" : "PRGn_r",
        "step": 15,
    }

    
    file_out = "mega_comparison"+var+"_"+period+"_test.png"

    residual = data[0]+data[1]-data[2]-data[3]

    comparison(data, set_up, graph_params, graph_params_anom, experiment, title, file_out, save, show, True, residual, "cont_shelf")
          
if __name__ == '__main__':
    main() # run the program

