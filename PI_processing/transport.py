import sys
import os
import numpy as np
import xarray as xr
import pandas as pd
from funcs import find_nearest, find_years
from mitgcm_python.grid import Grid
from directories_and_paths import *
from config_options import *
from make_transport import append_transport


def main():
    '''main script to define or read in the experiment and ensemble member
    and calculate the timeseries transport'''
    
    start_date = None
    
    if len(sys.argv) == 4:
        start_date = str(sys.argv[3])
    
    if len(sys.argv) > 2:
        experiment = str(sys.argv[1])
        ens_member = str(sys.argv[2])       
    else:
        experiment = input("Enter the experiment (eg. WIND, CTRL, TEMP, LENS): ").strip()
        ens_member = input("Enter the ensemble member (1-9): ").strip()
       
    if experiment == "LENS":
        if int(ens_member) > 5:
            filepath = output_path + "PAS_LENS00" + str(ens_member) + "_noOBC/"
        else:
            filepath = lens_path + "PAS_LENS00" + str(ens_member) + "_noOBC/"
    else:
        filepath = output_path + experiment + "_ens0" + str(ens_member) + "_noOBC/"

    if not os.path.exists(filepath) or not os.path.isdir(filepath):
        sys.exit(f"Stopped - Could not find directory {filepath}")

    if start_date is None:
        prev_timeseries = [filename for filename in os.listdir(filepath) if filename.startswith("transport")]
        if not prev_timeseries:
            print("No previous transport timeseries detected")
        else:
            print(f"Previous transport timeseries detected: {prev_timeseries}")

        start_date = input("Enter the date code to start or add to the timeseries at (eg 199201): ").strip()
    
    if experiment == "LENS":
        n_years = 181
    else:
        filepath_years = filepath + "/output/"
        n_years = find_years(start_date, filepath_years)

    grid = Grid(grid_filepath)

    # set up output file
    start_year = int(start_date[:4])
    end_year = start_year + n_years

    if experiment == "LENS" and int(ens_member) < 6:
        output_file = f"{output_path}transport{str(end_year)}_experiment{ens_member}.nc"
    else:
        output_file = f"{filepath}transport{str(end_year)}.nc"

    if start_date == "192001":
        print(f"First year in the dataset. Creating new timeseries.")

        lat_range = lat_range_73
        lon_range = lon_range_73

        [transport, transport_south] = append_transport(n_years, start_year, filepath, grid, lat_range, lon_range)

        time = pd.date_range(start="1920-01-01", periods=len(transport), freq="M")

        # Create the xarray dataset
        dataset = xr.Dataset(
            {
                "transport": (["time"], transport),
                "transport_south": (["time"], transport_south),
            },
            coords={
                "time": time,
            },
        )
        
    else:
        # to make
        print("not made yet")
    dataset.to_netcdf(output_file)

if __name__ == "__main__":
    main()
