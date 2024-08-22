import sys
import xarray as xr
import numpy as np
import pandas as pd
from scipy.stats import linregress
from tools.directories_and_paths import output_path, lens_path

def get_filepath(experiment, ens_member):
    """Get the appropriate file path based on the experiment and ensemble member."""
    base_path = output_path if experiment != "LENS" or int(ens_member) > 5 else lens_path
    return f"{base_path}PAS_LENS00{ens_member}_noOBC/" if experiment == "LENS" else f"{output_path}{experiment}_ens0{ens_member}_noOBC/"
    
def open_and_slice_data(file, var, option):
    """Open dataset and slice it based on the option."""
    data = xr.open_dataset(file, decode_times=False)
    var_data = data[var]
    
    if option == "total":
        if len(np.shape(var_data)) == 4:
            sliced_data = var_data[:, 0, ...]
        else:
            sliced_data = var_data
    elif option == "slice":
        from config_options import lat_slices, lon_slices
        sliced_data = var_data.sel(XC=lon_slices, method="nearest").sel(YC=slice(lat_slices[0], lat_slices[1]))
    else:
        raise ValueError(f"Invalid option: {option}")
    
    return sliced_data

def concatenate_var_years(experiment, var, option, ens_member):
    """Concatenate data over years for a given experiment and variable."""
    filepath = get_filepath(experiment, ens_member)
    file_template = f"{filepath}output/{{}}01/MITgcm/output.nc"
    
    total_data = []
    var_attrs = None  # To store the attributes of the variable
    
    for year in range(1920, 2101):
        data = open_and_slice_data(file_template.format(year), var, option)
        if var_attrs is None:
            var_attrs = data.attrs  # Store the attributes from the first file
        total_data.append(data)
    
    total_dataset = xr.concat(total_data, dim='time')
    total_dataset['time'] = pd.date_range(start="1920-01-01", periods=len(total_dataset), freq="M")
    
    # Reassign the original attributes to the concatenated data
    total_dataset.attrs = var_attrs
    
    return total_dataset

def create_ensemble_average(experiment, var, option, ensemble_members):
    """Create ensemble average for a variable across ensemble members."""
    args_list = [(experiment, var, option, str(ens_member)) for ens_member in ensemble_members]
    results = map(lambda args: concatenate_var_years(*args), args_list)
    
    ensemble_spread = xr.concat(results, dim='ensemble_member')
    ensemble_spread['ensemble_member'].attrs['long_name'] = "Ensemble member 01 to 09"
    filename = f"{output_path}{experiment}_files_temp/{var}.nc"
    ensemble_spread.to_netcdf(filename)

def calculate_trend(data, time):
    """Calculate trend, r-value, and p-value using linear regression."""
    
    if len(np.shape(data)) < 5:
        slopes, r_values, p_values = (np.empty([9, 384, 600]) for _ in range(3))
        for ens in range(data.shape[0]):
            for y in range(data.shape[2]):
                for x in range(data.shape[3]):
                    slope, _, r_value, p_value, _ = linregress(time, data[ens, :, y, x])
                    slopes[ens, y, x], r_values[ens, y, x], p_values[ens, y, x] = slope, r_value, p_value
    else:
        slopes, r_values, p_values = (np.empty([9, 50, 384, 600]) for _ in range(3))
        for z in range(data.shape[1]):
            for y in range(data.shape[2]):
                for x in range(data.shape[3]):
                    slope, _, r_value, p_value, _ = linregress(time, data[:, z, y, x])
                    slopes[z, y, x], r_values[z, y, x], p_values[z, y, x] = slope, r_value, p_value
    
    return slopes, r_values, p_values

def calc_trend(experiment, var="DENSITY"):
    """Calculate and save the trend for a given variable."""
    if var == "DENSITY":
        from mitgcm_python.diagnostics import density
        temp = xr.open_dataset(f"{output_path}{experiment}_files_temp/THETA_spread.nc", decode_times=False).THETA.values
        salt = xr.open_dataset(f"{output_path}{experiment}_files_temp/SALT_spread.nc", decode_times=False).SALT.values
        depth = xr.open_dataset(f"{output_path}{experiment}_files_temp/SALT_spread.nc", decode_times=False).Z.values
        press = np.broadcast_to(depth.reshape(1, 50, 1, 1), salt.shape)
        data = density("MDJWF", salt, temp, press)
    else: 
        filename = f"{output_path}{experiment}_files_temp/{var}.nc"
        dataset = xr.open_dataset(filename, decode_times=False)
        if var == "oceFWflx":
            data = dataset[var].values * 3600 * 24 * 365 / 1000 # convert from kg/m^2/s to m/yr
        if var =="UVEL":
            depth = dataset.Z.values

    time = dataset.time.values / 365  # Convert from days to years

    slopes, r_values, p_values = calculate_trend(data, time)
    
    if var in ["DENSITY", "UVEL"]:
        trend_ds = xr.Dataset(
            {
                'trend': (['z', 'lat', 'lon'], slopes),
                'pvalue': (['z', 'lat', 'lon'], p_values),
                'rvalue': (['z', 'lat', 'lon'], r_values)
            },
            coords={'z': depth, 'lat': dataset.YC.values, 'lon': dataset.XC.values}
        )
    else:
        trend_ds = xr.Dataset(
            {
                'trend': (['ens', 'lat', 'lon'], slopes),
                'pvalue': (['ens', 'lat', 'lon'], p_values),
                'rvalue': (['ens', 'lat', 'lon'], r_values)
            },
            coords={'ens': range(1, 10), 'lat': dataset.YC.values, 'lon': dataset.XC.values},
            attrs={
                'description': f'trend in time per grid point for {var}',
                'units': 'm/yr',
                },
        )
    
    trend_ds.to_netcdf(f"{output_path}{experiment}_files_temp/{var}_trend.nc")

if __name__ == "__main__":
    experiment = sys.argv[1]
    var = "oceFWflx"
    #create_ensemble_average(experiment, var, option="total", ensemble_members=range(1, 10))
    calc_trend(experiment, var)
