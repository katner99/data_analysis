import multiprocessing

from directories_and_paths import output_path, lens_path
from config_options import lat_range_U1, lon_range_U1, lat_range_U2, lon_range_U2

import xarray as xr

def read_data(ensemble, exp, year, lat_range, lon_range, loc):
    if exp == "LENS":
        if ensemble < 6:
            input_data = xr.open_dataset(f"{lens_path}PAS_LENS00{ensemble}_noOBC/output/{year}01/MITgcm/output.nc", decode_times=False)
        else:
            input_data = xr.open_dataset(f"{output_path}/PAS_{exp}00{ensemble}_noOBC/output/{year}01/MITgcm/output.nc", decode_times=False)
    else:
        input_data = xr.open_dataset(f"{output_path}/{exp}_ens0{ensemble}_noOBC/output/{year}01/MITgcm/output.nc", decode_times=False)
    
    salt = input_data.SALT.values[:, :, lat_range[0]:lat_range[1], lon_range]
    uvel = input_data.UVEL.values[:, :, lat_range[0]:lat_range[1], lon_range]
    time = input_data.time
    lat = input_data.YC.values[lat_range[0]:lat_range[1]]
    z = input_data.Z.values 
    
    data = xr.Dataset(
        {
            f"salt_{loc}": (["time", "z", "lat"], salt),
            f"uvel_{loc}": (["time", "z", "lat"], uvel),
        },
        coords={"time": time, "z": z, "lat": lat},
    )
    return data

def append_years_slice(ensemble, exp, start_year, n_years):
    from config_options import lat_range_U1, lat_range_U2, lat_range_U3, lat_range_U4, lon_range_U1, lon_range_U2, lon_range_U3, lon_range_U4
    slice = []
    for year in range(start_year, start_year+n_years):
        datasets = []
        data_U1 = read_data(ensemble, exp, year, lat_range_U1, lon_range_U1, "U1")
        datasets.append(data_U1)
        data_U2 = read_data(ensemble, exp, year, lat_range_U2, lon_range_U2, "U2")
        datasets.append(data_U2)
        data_U3 = read_data(ensemble, exp, year, lat_range_U3, lon_range_U3, "U3")
        datasets.append(data_U3)
        data_U4 = read_data(ensemble, exp, year, lat_range_U4, lon_range_U4, "U4")
        datasets.append(data_U4)
        data = xr.combine_by_coords([data_U1, data_U2, data_U3, data_U4])
        slice.append(data)   
    return xr.concat(slice, dim="time")


def main():
    experiments = ["WIND", "TEMP", "LENS"]
    ensemble = 9
    print(ensemble)

    for exp in experiments:
        data = append_years_slice(ensemble, exp, 1920, 181)
        output_file = f"{output_path}{exp}_ens0{ensemble}_noOBC/ensemble_mean_{exp}.nc"
        data.to_netcdf(output_file)

if __name__ == '__main__':
    main() 

