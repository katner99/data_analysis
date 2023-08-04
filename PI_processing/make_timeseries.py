import sys
import os
import numpy as np
import xarray as xr
from funcs import find_nearest, append_years
from mitgcm_python.grid import Grid
from directories_and_paths import *

def find_years(start_year, filepath):
    year_directories = [name for name in os.listdir(filepath)]
    valid_directories = [year for year in year_directories if year >= start_year]
    num_years = len(valid_directories)
    print('Number of years detected:' + str(num_years) + ' equivalent to the following years: ')
    print(valid_directories)
    return num_years
    

def main():
    # ask for the ensemble member and experiment
    experiment = input('Enter the experiment (eg. WIND, CTRL, TEMP): ').strip()
    ens_member = input('Enter the ensemble member (1-3): ').strip()

    #set up filepath and check that it exists
    if experiment == "LENS":
        filepath = lens_path + "PAS_LENS00" + str(ens_member)+ "_noOBC/"
    else:
        filepath = output_path + experiment + "_ens0" + str(ens_member)+ "_noOBC/"

    if not os.path.exists(filepath) or not os.path.isdir(filepath):
        sys.exit(f"Stopped - Could not find directory {filepath}")

    # search for any previous timeseries:
    prev_timeseries = [filename for filename in os.listdir(filepath) if filename.startswith("timeseries")]

    if not prev_timeseries:
        print('No previous timeseries detected')
    else:
        print(f'Previous timeseries detected: {prev_timeseries}')
    
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

    filepath_years = filepath + "/output/"
    n_years = find_years(start_date, filepath_years)
    #n_years = 181

    # Set up a new xarray dataset with the required dimensions (replace this with your data)
    # set up grid
    grid = Grid(grid_filepath)
    grid_file = xr.open_dataset(grid_filepath)
    
    # set lat, lon, and depth range
    depth_range = [find_nearest(grid_file.Z.values, -200), find_nearest(grid_file.Z.values, -700)]
    lon_range = [find_nearest(grid_file.XC.values, 250), find_nearest(grid_file.XC.values, 260)]
    lat_range = [find_nearest(grid_file.YC.values, -76), find_nearest(grid_file.YC.values, -72)]

    # set up output file
    filename = "output.nc"
    start_year = int(start_date[:4])
    end_year = start_year + n_years
    if experiment == 'LENS':
        output_file = output_path + 'timeseries'+str(end_year)+'_experiment'+ens_member+'.nc'
    else:
        output_file = filepath + 'timeseries'+str(end_year)+'.nc'

    if start_date == "192001":
        print(f"First year in the dataset. Creating new timeseries.")
            
        [theta, salt, seaice] = append_years(n_years, start_year, filepath_years, filename, grid, lat_range, lon_range, depth_range)
            
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
        input_file = filepath + 'timeseries'+str(start_year)+'.nc'

        try:
            open(input_file)
            print(f"Previous dataset created. Appending new values to the previous timeseries.")
        except FileNotFoundError:
            sys.exit(f"Stopped - Could not find previous timeseries, go back to the 1920s {filepath}")
       
        dataset_old = xr.open_dataset(input_file)
        months_to_skip = (int(start_date[:4])-1920)*12

        # Add more variables (theta, salt, sea_ice) to the existing dataset
        [theta, salt, seaice] = append_years(n_years, start_year, filepath_years, filename, grid, lat_range, lon_range, depth_range)

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

        combined_data = []
        # Concatenate the new data for each variable along the time dimension
        variables_to_concat = ['theta', 'salt', 'sea_ice']
        for variable in variables_to_concat:
            existing_data = dataset_old[variable]
            new_data = dataset_new[variable]
            combined_data.append(xr.concat([existing_data, new_data], dim='time'))

        dataset = xr.Dataset(
            {
                'theta': (['time'], combined_data[0]),
                'salt': (['time'], combined_data[1]),
                'sea_ice': (['time'], combined_data[2]),
            },
            coords={
                'time': np.arange(len(combined_data[0])),
            }
        )
        # remove the previous file
        os.remove(input_file)
            
    # Save the dataset to a NetCDF file
    dataset.to_netcdf(output_file) 
   

if __name__ == "__main__":
    main()
