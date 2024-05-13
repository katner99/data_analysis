import sys
import xarray as xr
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from mitgcm_python.grid import Grid
from directories_and_paths import output_path, lens_path, grid_filepath
from multiprocessing import Pool, cpu_count
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

import multiprocessing as mp

def average_depth(filepath):
    filename = "MITgcm/output.nc"
    file = f"{filepath}/210001/{filename}"
    data = xr.open_dataset(file, decode_times=False)
    mask = data["hFacC"]
    temp = data["THETA"]
    temp = temp.where(mask == 1)
    depth_average = temp.sel(Z=slice(-200, -700)).mean(dim='Z')
    return depth_average.mean(dim='time')

def ttest_over_lat_lon(data_1, data_2):
    p_values = np.empty(np.shape(data_1)[1:])  # Initialize array to store p-values

    for lat_idx in range(np.shape(p_values)[0]):
        for lon_idx in range(np.shape(p_values)[1]):
            test = ttest_ind(data_1[:,lat_idx, lon_idx], data_2[:, lat_idx, lon_idx])
            p_values[lat_idx, lon_idx] = test.pvalue

    return p_values

def read_experiment(exp, ensemble):
    data = []
    for ens in ensemble:
        if exp == "LENS":
            if ens > 5:
                filepath = f"{output_path}PAS_LENS00{ens}_noOBC/output/"
            else:
                filepath = f"{lens_path}PAS_LENS00{ens}_noOBC/output/"
        else:
            filepath = f"{output_path}{exp}_ens0{ens}_noOBC/output/"
        data.append(average_depth(filepath))
    return data


def main():
    ensemble = list(range(1, 10))  
    data_1 = read_experiment("CTRL", ensemble)
    data_2 = read_experiment("LENS", ensemble)

    data = xr.concat(data_1, dim="ensemble")
    lat = data.YC.values
    lon = data.XC.values

    print(np.nanmax(data))

    p_value = ttest_over_lat_lon(np.array(data_1), np.array(data_2))

    print(np.nanmin(p_value))
    print(np.shape(p_value))
    trends = xr.Dataset(
        {
            'pvalue' : (['lat', 'lon'], p_value),
        },
        coords={'lat': lat, 'lon': lon}
    )

    

    trends.to_netcdf("CL_pval.nc")

if __name__ == "__main__":
    main()