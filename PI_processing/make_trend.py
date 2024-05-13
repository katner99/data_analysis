import sys
import xarray as xr
import numpy as np
import pandas as pd
from scipy.stats import linregress
from directories_and_paths import output_path, lens_path
from multiprocessing import Pool, cpu_count

def concatenate_var_years(args):
    """ reads and concatentates years """
    experiment, var, option, ens_member = args
    filename = "/MITgcm/output.nc"
    
    if experiment == "LENS":
        if int(ens_member) > 5:
            filepath = f"{output_path}PAS_LENS00{ens_member}_noOBC/"
        else:
            filepath = f"{lens_path}PAS_LENS00{ens_member}_noOBC/"
    else:
        filepath = f"{output_path}{experiment}_ens0{ens_member}_noOBC/"
    
    filepaths = [f"{filepath}output/{year}01{filename}" for year in range(1920, 2101)]
    
    total_data = []
    for file in filepaths:
        data = xr.open_dataset(file, decode_times=False)
        if option == "total":
            total_data.append(data[var])
        elif option == "slice":
            from config_options import lat_slices, lon_slices
            sliced_data = data[var].sel(XC = lon_slices, method = "nearest")
            sliced_data = sliced_data.sel(YC = slice(lat_slices[0], lat_slices[1]))
            total_data.append(sliced_data)
    
    total_dataset = xr.concat(total_data, dim='time')
    time_coords = pd.date_range(start="1920-01-01", periods=len(total_dataset), freq="M")
    total_dataset['time'] = time_coords
    
    return total_dataset

def total_trend(filename, var, experiment):
    """calculates trend for simple 2d variable with dims time, lat, lon"""
    dataset = xr.open_dataset(filename, decode_times=False)
        
    lat = dataset.YC.values
    lon = dataset.XC.values
    time = dataset.time.values
    data = dataset[var].values

    slopes, r_values, p_values = calculate_trend(data, time)

    trend_ds = xr.Dataset(
        {
            'trend': (['lat', 'lon'], slopes),
            'pvalue': (['lat', 'lon'], p_values),
            'rvalue': (['lat', 'lon'], r_values)
        },
        coords={'lat': lat, 'lon': lon}
    )

    trend_ds.to_netcdf(f"{output_path}{experiment}_files_temp/{var}_trend.nc")

def slice_trend(filename, var, experiment):
    """ calculates trend over specific slices of the ocean, dims time, depth, lat, lon """
    dataset = xr.open_dataset(filename, decode_times=False)
        
    lat = dataset.YC.values
    lon = dataset.XC.values
    depth = dataset.Z.values
    time = dataset.time.values
    data = dataset[var].values

    def calculate_trend_slice(z):
        slopes = np.empty(len(lat))
        r_values = np.empty(len(lat))
        p_values = np.empty(len(lat))
        for y in range(len(lat)):
            for x in lon:
                slope, intercept, r_value, p_value, std_err = linregress(time, data[:,z, y, x])
                slopes[y] = slope
                r_values[y] = r_value
                p_values[y] = p_value
        
        return slopes, r_values, p_values

    num_processes = min(cpu_count(), 22)

    with Pool(processes=num_processes) as pool:
        results = pool.map(calculate_trend_slice, range(len(depth)))

    slopes = np.stack([result[0] for result in results])
    r_values = np.stack([result[1] for result in results])
    p_values = np.stack([result[2] for result in results])

    trend_ds = xr.Dataset(
        {
            'trend': (['depth', 'lat', 'lon'], slopes),
            'pvalue': (['depth', 'lat', 'lon'], p_values),
            'rvalue': (['depth', 'lat', 'lon'], r_values)
        },
        coords={'depth': depth, 'lat': lat, 'lon': lon}
    )

    trend_ds.to_netcdf(f"{output_path}{experiment}_files_temp/{var}_trend.nc")

def calculate_trend(data, time):
    """calculates trend for the slices (attempt at parallelisation)"""
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
    """ 
    script to read in the MITgcm output files and calculate grid trends with time
    INPUT: 
    experiment: options avail CTRL, LENS, WIND, TEMP
    
    OUTPUT:
    concatenated data
    trend data
    """
    # read in the experiment, set up ensemble members and the variable to run
    experiment = str(sys.argv[1])
    ensemble_members = list(range(1, 10))  
    var = "UVEL"
    option = "slice" # current options available: total, slice

    args_list = [(experiment, var, option, str(ens_member)) for ens_member in ensemble_members]
    
    num_processes = min(22, cpu_count())
    
    # I want to calculate trend with time so need all the years to be concatenated
    with Pool(processes=num_processes) as pool:
        results = pool.map(concatenate_var_years, args_list)
    
    ensemble_mean = xr.concat(results, dim='ensemble_member').mean(dim='ensemble_member')
    filename = f"{output_path}{experiment}_files_temp/{var}.nc"
    # save incase the next bit dies RIP
    ensemble_mean.to_netcdf(filename)

    filename = f"{output_path}{experiment}_files_temp/{var}.nc"

    # calculate trend
    if option == "total":
        total_trend(filename, var, experiment)
    else:
        slice_trend(filename, var, experiment)

if __name__ == "__main__":
    main()
