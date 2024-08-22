import sys
import xarray as xr
import numpy as np
import pandas as pd
from scipy.stats import linregress
from directories_and_paths import output_path, lens_path

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
            total_data.append(data[var][:,0,...])
        elif option == "slice":
            from config_options import lat_slices, lon_slices
            sliced_data = data[var].sel(XC = lon_slices, method = "nearest")
            sliced_data = sliced_data.sel(YC = slice(lat_slices[0], lat_slices[1]))
            total_data.append(sliced_data)
       
    total_dataset = xr.concat(total_data, dim='time')
    
    time_coords = pd.date_range(start="1920-01-01", periods=len(total_dataset), freq="M")
    total_dataset['time'] = time_coords

    return total_dataset

def calculate_trend_total(data, time):
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

def calculate_trend_slice(data, time):
    """calculates trend for the slices (attempt at parallelisation)"""
    slopes = np.empty(data.shape[1:])
    r_values = np.empty(data.shape[1:])
    p_values = np.empty(data.shape[1:])

    for z in range(data.shape[1]):
        for y in range(data.shape[2]):
            for x in range(data.shape[3]):
                slope, intercept, r_value, p_value, std_err = linregress(time, data[:, z, y, x])
                slopes[z, y, x] = slope
                r_values[z, y, x] = r_value
                p_values[z, y, x] = p_value

    return slopes, r_values, p_values

def calculate_trend_profile(data, time):
    """calculates trend for the slices (attempt at parallelisation)"""
    slopes = np.empty(data.shape[1:])
    r_values = np.empty(data.shape[1:])
    p_values = np.empty(data.shape[1:])

    for z in range(data.shape[1]):
        slope, intercept, r_value, p_value, std_err = linregress(time, data[:, z])
        slopes[z] = slope
        r_values[z] = r_value
        p_values[z] = p_value

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
    option = "total" # current options available: total, slice

    args_list = [(experiment, var, option, str(ens_member)) for ens_member in ensemble_members]
    
    #num_processes = min(22, cpu_count())
    
    # I want to calculate trend with time so need all the years to be concatenated
    results = map(concatenate_var_years, args_list)

    ensemble_spread = xr.concat(results, dim='ensemble_member').mean(dim='ensemble_member')

    filename = f"{output_path}{experiment}_files_temp/{var}.nc"
    # save incase the next bit dies RIP
    ensemble_spread.to_netcdf(filename)

def calc_trend():
    experiment = str(sys.argv[1]) 
    var = "DENSITY"
        
    if var == "DENSITY":
        from mitgcm_python.diagnostics import density
        dataset = xr.open_dataset(f"{output_path}{experiment}_files_temp/THETA_spread.nc", decode_times=False)
        temp = dataset.THETA.values
        dataset = xr.open_dataset(f"{output_path}{experiment}_files_temp/SALT_spread.nc", decode_times=False)
        salt = dataset.SALT.values
        depth = dataset.Z.values
        array_1_reshaped = depth.reshape(1, 50, 1, 1)
        press = np.broadcast_to(array_1_reshaped, salt.shape)
        data = density("MDJWF", salt, temp, press)
        
    else:
        filename = f"{output_path}{experiment}_files_temp/{var}_spread.nc"
        dataset = xr.open_dataset(filename, decode_times=False)
        data = dataset[var].values
        depth = dataset.Z.values

    time = (dataset.time.values)/365 # convert from days to year
    lat = dataset.YC.values
    lon = dataset.XC.values

    slopes = np.empty(data.shape[1:])
    r_values = np.empty(data.shape[1:])
    p_values = np.empty(data.shape[1:])

    for z in range(data.shape[1]):
        for y in range(data.shape[2]):
            for x in range(data.shape[3]):
                slope, intercept, r_value, p_value, std_err = linregress(time, data[:, z, y, x])
                slopes[z, y, x] = slope
                r_values[z, y, x] = r_value
                p_values[z, y, x] = p_value

    trend_ds = xr.Dataset(
        {
            'trend': (['z', 'lat', 'lon'], slopes),
            'pvalue': (['z', 'lat', 'lon'], p_values),
            'rvalue': (['z', 'lat', 'lon'], r_values)
        },
        coords={'z': depth, 'lat': lat, 'lon': lon}
    )

    filename = f"{output_path}{experiment}_files_temp/{var}_trend.nc"
    # save incase the next bit dies RIP
    trend_ds.to_netcdf(filename)

if __name__ == "__main__":
    calc_trend()
