"""
Plotting script for generating comparison plots of oceanographic data.
"""

from .config_options import config_comparison
from .tools.plots import read_mask
from .tools.plots_2d import comparison
from .directories_and_paths import output_path

import sys
import xarray as xr

def read_data(var, var_name=None):
    """
    Reads in analyzed data based on the variable type.

    Parameters:
    - var (str): Type of data to read ("trend", "pvalue", "mean").
    - var_name (str): Name of the variable (default is None).

    Returns:
    - filepaths (list): List of file paths to the data files.
    - file_out (str): Output filename for the plot.
    """
    file_prefix = f"{output_path}{{}}_files_temp/{{}}"
    filename = f"{var_name}_spread.nc" if var in ["trend", "pvalue"] else f"{var_name}.nc"
    file_out = f"mega_comparison_{var}_{var_name}.png" if var in ["trend", "pvalue"] else f"mega_comparison_{var}_{var_name}_mean.png"
    
    experiments = ["CTRL", "LENS", "WIND", "TEMP"]
    filepaths = [file_prefix.format(exp, filename) for exp in experiments]
    
    return filepaths, file_out

def read_climatology_data(period):
    """
    Generates file paths and output filename for climatology data.

    Parameters:
    - period (str): Time period for climatology data.

    Returns:
    - filepaths (list): List of file paths to the data files.
    - file_out (str): Output filename for the plot.
    """
    experiments = ["CTRL", "LENS", "WIND", "TEMP"]
    filepaths = [f"{output_path}average_{exp}_{period}.nc" for exp in experiments]
    file_out = f"mega_comparison_climatology_{period}.png"
    
    return filepaths, file_out

def main():
    """
    Main function to create the 4x4 comparison plot.
    """
    var = "SIheff"  
    save = True     # Save figure to file
    show = False    # Display figure after script completion
    loc = "trend"   # Set to the parameter or "mean" or "trend" to determine data source

    if loc in ["trend", "mean"]:
        filepaths, file_out = read_data(loc, var)
    else:
        period = "2070-2100"
        filepaths, file_out = read_climatology_data(period)

    # Check for missing files
    for filepath in filepaths:
        if not open(filepath, 'r'):
            sys.exit(f"Stopped - Could not find input file {filepath}")

    # Experiment names for the plot
    experiments = ["NONE", "ALL", "WIND", "THERMO"]

    # Load datasets
    input_data = [xr.open_dataset(filepath, decode_times=False) for filepath in filepaths]
   
    # Set up data for plotting
    set_up = read_mask()

    # Configure plot parameters
    if loc in ["trend", "pvalue", "mean"]:
        data, graph_params, graph_params_anom = config_comparison(var, input_data, loc=loc)
    else:
        data, graph_params, graph_params_anom = config_comparison(var, input_data, period=period)

    # Generate the comparison plot
    comparison(
        data,
        set_up,
        graph_params,
        graph_params_anom,
        experiments,
        file_out,
        save,
        show,
        linearity=False,   # Option to calculate linearity in the plot
        zoom="cont_shelf"  # Option to zoom into a specific area: "cont_shelf", "ice_shelf", or leave blank
    )

if __name__ == "__main__":
    main()
