import sys
import os
import numpy as np
import pandas as pd
import xarray as xr
from dateutil.relativedelta import relativedelta
from funcs import find_nearest, append_years
from mitgcm_python.grid import Grid
from directories_and_paths import *

def find_years(start_year, filepath):
    year_directories = [name for name in os.listdir(filepath)]
    valid_directories = [year for year in year_directories if year >= start_year]
    num_years = len(valid_directories)
    print('Number of years detected:' + str(num_years) + ' equivalent to the following years: ')
    print(valid_directories)
    return valid_directories, num_years
    

def main():
    # ask for the start year
    start_date = input('Enter the date code to start or add to the timeseries at (eg 199201): ').strip()
    # Make sure it's a date
    valid_date = len(start_date)==6
    try:
        int(valid_date)
    except(ValueError):
        valid_date = False
    if not valid_date:
        print('Error: invalid date code ' + start_date)
        sys.exit()

    # ask for the ensemble member and experiment
    ens_member = input('Enter the ensemble member (1-3): ').strip()
    experiment = input('Enter the experiment (eg. WIND, CTRL, TEMP): ').strip()
    
    #set up filepath and check that it exists
    filepath = output_path + experiment + "_ens0" + str(ens_member)+ "_noOBC/output/"
    if not os.path.exists(filepath) or not os.path.isdir(filepath):
        sys.exit(f"Stopped - Could not find directory {filepath}")

    year_directories, n_years = find_years(start_date, filepath)

    # Set up a new xarray dataset with the required dimensions (replace this with your data)
    # set up grid
    grid = Grid(grid_filepath)
    grid_file = xr.open_dataset(grid_filepath)
    
    # set lat, lon, and depth range
    depth_range = [find_nearest(grid_file.Z.values, -200), find_nearest(grid_file.Z.values, -700)]
    lon_range = [find_nearest(grid_file.XC.values, 250), find_nearest(grid_file.XC.values, 260)]
    lat_range = [find_nearest(grid_file.YC.values, -76), find_nearest(grid_file.YC.values, -72)]
    
    # set up output file
    output_file = output_path + experiment + "_ens0" + str(ens_member)+ "_noOBC/" + 'timeseries.nc'
    filename = "output.nc"
    start_year = int(start_date[:4])
    end_year = start_year + n_years

    if start_date == "192001":
        print(f"First year in the dataset. Creating new timeseries.")
            
        [theta, salt, seaice] = append_years(n_years, start_year, filepath, filename, grid, lat_range, lon_range, depth_range)
            
        # Create the xarray dataset
        dataset = xr.Dataset(
            {
                'theta': (['time'], theta),
                'salt': (['time'], salt),
                'sea_ice': (['time'], seaice),
            },
            coords={
                'time': np.arange(12*n_years),
            }
        )
       
    else:
        # check for previous timeseries
        try:
            open(output_file)
            print(f"Previous dataset created. Appending new values to the previous timeseries.")
        except FileNotFoundError:
            sys.exit(f"Stopped - Could not find previous timeseries, go back to the 1920s {filepath}")
       
        dataset = xr.open_dataset(output_file)

        months_to_skip = (int(start_date[:4])-1920)*12
        
        # Add more variables (theta, salt, sea_ice) to the existing dataset
        [theta, salt, seaice] = append_years(n_years, start_year, filepath, filename, grid, lat_range, lon_range, depth_range)

        # Create the xarray dataset
        dataset_new = xr.Dataset(
            {
                'theta': (['time'], theta),
                'salt': (['time'], salt),
                'sea_ice': (['time'], seaice),
            },
            coords={
                'time': np.arange(months_to_skip, months_to_skip+12*n_years),
            }
        )

        # Concatenate the new data for each variable along the time dimension
        variables_to_concat = ['theta', 'salt', 'sea_ice']
        for variable in variables_to_concat:
            existing_data = dataset[variable]
            new_data = dataset_new[variable]
            combined_data = xr.concat([existing_data, new_data], dim='time')
            dataset[variable] = combined_data
        
        
            
    # Save the dataset to a NetCDF file
    print(int(start_date[:4])+n_years)
    dataset.to_netcdf(output_file) 
    with open(output_path + experiment + "_ens0" + str(ens_member)+ "_noOBC/" + 'calendar.txt', 'w') as f:
        f.write(str(end_year))
   

if __name__ == "__main__":
    main()
