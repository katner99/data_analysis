import sys
import xarray as xr
import numpy as np
import pandas as pd
from scipy.stats import linregress
from directories_and_paths import output_path, lens_path
from multiprocessing import Pool, cpu_count

def concatenate_var_years(args):
    experiment, var, ens_member = args
    filename = "/MITgcm/output.nc"
    
    if experiment == "LENS":
        if int(ens_member) > 5:
            filepath = f"{output_path}PAS_LENS00{ens_member}_noOBC/"
        else:
            filepath = f"{lens_path}PAS_LENS00{ens_member}_noOBC/"
    else:
        filepath = f"{output_path}{experiment}_ens0{ens_member}_noOBC/"
    
    filepaths = [f"{filepath}output/{year}01{filename}" for year in range(1920, 2101)]
    
    sea_ice_data = []
    for file in filepaths:
        data = xr.open_dataset(file, decode_times=False)
        sea_ice_data.append(data[var])
    
    sea_ice_dataset = xr.concat(sea_ice_data, dim='time')
    time_coords = pd.date_range(start="1920-01-01", periods=len(sea_ice_dataset), freq="M")
    sea_ice_dataset['time'] = time_coords
    
    return sea_ice_dataset

def calculate_trend(data, time):
    slopes = np.empty(data.shape[1:])
    r_values = np.empty(data.shape[1:])
    p_values = np.empty(data.shape[1:])

    for y in range(data.shape[1]):
        for x in range(data.shape[2]):
            slope, intercept, r_value, p_value, std_err = linregress(time, data[:, y, x])
            slopes[y, x] = slope
            r_values[y, x] = r_value
            p_values[y, x] = p_value

    return slopes, r_values, p_values

def main():
    experiment = str(sys.argv[1])
    ensemble_members = list(range(1, 10))  
    var = "SHIfwFlx"

    args_list = [(experiment, var, str(ens_member)) for ens_member in ensemble_members]
    
    num_processes = min(22, cpu_count())
    
    with Pool(processes=num_processes) as pool:
        results = pool.map(concatenate_var_years, args_list)
    
    ensemble_mean = xr.concat(results, dim='ensemble_member').mean(dim='ensemble_member')
    
    ensemble_mean.to_netcdf(f"{output_path}{experiment}_files_temp/{var}.nc")

    filename = f"{output_path}{experiment}_files_temp/{var}.nc"
    dataset = xr.open_dataset(filename, decode_times=False)
    
    lat = dataset.YC.values
    lon = dataset.XC.values
    time = dataset.time.values
    data = dataset[var].values

    slopes, r_values, p_values = calculate_trend(data, time)

    # Create a new Dataset to store trend, p-value, and r-value
    trend_ds = xr.Dataset(
        {
            'trend': (['lat', 'lon'], slopes),
            'pvalue': (['lat', 'lon'], p_values),
            'rvalue': (['lat', 'lon'], r_values)
        },
        coords={'lat': lat, 'lon': lon}
    )

    # Save the dataset to a netCDF file
    trend_ds.to_netcdf(f"{output_path}{experiment}_files_temp/{var}_trend.nc")

if __name__ == "__main__":
    main()
