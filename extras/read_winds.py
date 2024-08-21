from mitgcm_python.file_io import read_binary
from mitgcm_python.grid import Grid

from funcs import binary_coords, find_nearest, global_to_amundsen
from directories_and_paths import *
from config_options import *
 
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import interp2d

def identify_lat_and_lon():
    grid = Grid(grid_filepath)
    dimensions = ('t','y','x')
    
    path = '/data/oceans_input/raw_input_data/CESM/LENS/daily/'
    path += "PRECT/"
    data_binary = xr.open_dataset(path + "b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.04020101-04991231.nc", decode_times=False)
    lat = data_binary.lat.values
    lon = data_binary.lon.values

    #print(lon)

    data_relative = xr.open_dataset(f"{lens_noOBC_path}/192001/MITgcm/output.nc", decode_times=False)
    lat = data_relative.YC.values
    lon = data_relative.XC.values
    
    #print(lon[0], lon[-1])
    
    #print(grid.lon_1d[0]+360, grid.lon_1d[-1]+360)
    
    return data_relative
    
def cut_lat_and_lon():
    grid = Grid(grid_filepath)
    
    [data_binary, data_relative] = identify_lat_and_lon()
    
    lat = data_binary.lat.values
    lon = data_binary.lon.values
    
    lat_index = [find_nearest(lat, grid.lat_1d[0]), find_nearest(lat, grid.lat_1d[-1])]
    lon_index = [find_nearest(lon, grid.lon_1d[0]+360), find_nearest(lon, grid.lon_1d[-1]+360)]
    
    precip_binary = np.mean(data_binary.PRECT.values[:31, lat_index[0]:lat_index[1], lon_index[0]:lon_index[1]], axis = 0)
    precip_relative = data_relative.EXFpreci.values[0,:,:]
    
    fig, ax = plt.subplots(nrows=1, ncols=2,figsize=(12, 5))
    
    ax[0].contourf(precip_binary)
    ax[1].contourf(precip_relative)
    #plt.colorbar()
    fig.savefig("test_precip.png")

def compare_binary():
    
    grid = Grid(grid_filepath)
    
    [data_raw, data_relative] = identify_lat_and_lon()
    
    lat = data_raw.lat.values
    lon = data_raw.lon.values

    lat_len, lon_len, time = 192, 288, 365
    dimensions = ('t','y','x')

    lat_index = [find_nearest(lat, grid.lat_1d[0]), find_nearest(lat, grid.lat_1d[-1])]
    lon_index = [find_nearest(lon, grid.lon_1d[0]+360), find_nearest(lon, grid.lon_1d[-1]+360)]

    filepath = "/data/oceans_input/processed_input_data/CESM/LENS/"

    filename = f"{filepath}LENS_ens001_PRECT_1920"
    data_temporary = read_binary(filename, [lon_len, lat_len, time], dimensions)
    precip_binary = np.mean(data_temporary[:31,lat_index[0]:lat_index[1], lon_index[0]:lon_index[1]], axis = 0)
    #precip_binary = np.mean(data_temporary[:31,:,:], axis = 0)
    #precip_raw = np.mean(data_raw.PRECT.values[:31,:,:], axis = 0)
    precip_relative = data_relative.EXFpreci.values[0,:,:]

    fig, ax = plt.subplots(nrows=1, ncols=2,figsize=(12, 5))
    
    c = ax[0].contourf(precip_binary, levels=np.linspace(0, 3*10**(-8), 8))
    fig.colorbar(c, ax=ax[0])
    c = ax[1].contourf(precip_relative,levels=np.linspace(0, 3*10**(-8), 8))
    fig.colorbar(c, ax=ax[1])
    fig.savefig("test_precip.png")
    #plt.show()
  
def new_grid():
    grid = Grid(grid_filepath)
    var = "PRECT"
    [lat, lon] = binary_coords(var)
    
    dimensions = ('t','y','x')
    
    lat_index = [find_nearest(lat, grid.lat_1d[0]), find_nearest(lat, grid.lat_1d[-1])]
    lon_index = [find_nearest(lon, grid.lon_1d[0]+360), find_nearest(lon, grid.lon_1d[-1]+360)]

    filename = f"{binary_path}LENS_ens001_{var}_1920"

    data_global = read_binary(filename, [lon_bin, lat_bin, time_bin], dimensions)
    data_amundsen = data_global[:,lat_index[0]:lat_index[1], lon_index[0]:lon_index[1]]

    y_lo_res = lat[lat_index[0]:lat_index[1]]
    x_lo_res = lon[lon_index[0]:lon_index[1]]

    filename_example = f"{lens_O_path}192001/MITgcm/output.nc"
    high_res_data = xr.open_dataset(filename_example, decode_times = False)
    x_hi_res = high_res_data.XC.values
    y_hi_res = high_res_data.YC.values
    
    f = interp2d(x_lo_res, y_lo_res, data_amundsen, kind = 'return')
    data_amundsen_hi_res = f(x_hi_res, y_hi_res)
    
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))
    
    ax[0].contourf(np.mean(data_amundsen[:31,...], axis = 0), levels=np.linspace(0, np.max(data_amundsen), 10))
    ax[0].set_title("low res")
    
    ax[1].contourf(np.mean(data_amundsen_hi_res[:31,...], axis = 0), levels=np.linspace(0, np.max(data_amundsen), 10))
    ax[1].set_title("high res")
    
    c = ax[2].contourf(high_res_data.EXFpreci.values[0,:,:],levels=np.linspace(0, np.max(data_amundsen), 10))
    fig.colorbar(c, ax=ax[2])
    ax[2].set_title("simulation")

    #fig.savefig("test_precip_res_qui.png")

    #plt.show()

def calculate_monthly_average():
    grid = Grid(grid_filepath)
    dimensions = ('t','y','x')
    var = "UBOT"
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    
    filename = f"{binary_path}LENS/LENS_ens001_{var}_1920"

    data_global = read_binary(filename, [lon_bin, lat_bin, time_bin], dimensions)
    
    day_in_year = 0
    
    data = []

    for month in days_in_month:
        monthly_mean = np.mean(data_global[day_in_year:day_in_year + month,:,:], axis=(0)) 
        data_grid = global_to_amundsen(monthly_mean, var, grid)
        data.append(data_grid)
        day_in_year += month

    data = np.array(data)
    
    simulation = xr.open_dataset(f"{lens_O_path}192001/MITgcm/output.nc", decode_times = False)
    
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))
    
    ax[0].contourf(data[0,:,:], levels=np.linspace(np.min(data[0,:,:]), np.max(data[0,:,:]), 20))
    ax[0].set_title("averaged binary data")
    
    c = ax[1].contourf(simulation.EXFuwind.values[0,:,:],levels=np.linspace(np.min(data[0,:,:]), np.max(data[0,:,:]), 20))
    fig.colorbar(c, ax=ax[1])
    ax[1].set_title("simulation")

    fig.savefig("test_Vwind_res_month_av.png")

    plt.show()
    
def calculate_ensemble_mean(plot = False):
    grid = Grid(grid_filepath)
    dimensions = ('t','y','x')
    var = "UBOT"
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    
    data_ens = []
    for ens in range (1):
        data_year = []
        for year in range(1920, 1922):
            # read year by year and create data for each year
            filename = f"{binary_path}LENS/LENS_ens00{ens}_{var}_{year}"
            data_global = read_binary(filename, [lon_bin, lat_bin, time_bin], dimensions)
            day_in_year = 0
            for month in days_in_month:
                monthly_mean = np.mean(data_global[day_in_year:day_in_year + month,:,:], axis=(0)) 
                data_grid = global_to_amundsen(monthly_mean, var, grid)
                data_year.append(data_grid)
                day_in_year += month
        data_year = np.array(data_year)
        data_ens.append(data_year)

    data = np.array(np.mean(data_ens, axis = 0))
    
    simulation = xr.open_dataset(f"{output_path}average_LENS_1920-1950.nc", decode_times = False)
    
    if plot:
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))
    
        ax[0].contourf(data[0,:,:], levels=np.linspace(np.min(data[0,:,:]), np.max(data[0,:,:]), 20))
        ax[0].set_title("averaged binary data")
    
        c = ax[1].contourf(simulation.EXFvwind.values[0,:,:],levels=np.linspace(np.min(data[0,:,:]), np.max(data[0,:,:]), 20))
        fig.colorbar(c, ax=ax[1])
        ax[1].set_title("simulation")

        fig.savefig("test_uwind_ens_mean.png")

        plt.show()
    
    num_days = 24

    # Create a time coordinate using xarray with the appropriate frequency (assuming daily data)
    time_coord = pd.date_range(start="1920-01-01", periods=num_days, freq='M')

    
    # Create an xarray dataset
    ds = xr.Dataset(
        {
            "data": (["time", "y", "x"], data)
        },
        coords={
            "time": time_coord
        }
    )

    # Save the monthly average dataset to a NetCDF file
    ds.to_netcdf("vbot_lens.nc")
    
def main():
    grid = Grid(grid_filepath)
    dimensions = ('t','y','x')
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    exp = "PIctrl"
    data = []
    for var in ["VBOT", "UBOT"]:
        data_ens = []
        for ens in range (1, 10):
            data_year = []
            for year in range(1920, 2101):
                # read year by year and create data for each year
                filename = f"{binary_path}{exp}/{exp}_ens0{ens}_{var}_{year}"
                data_global = read_binary(filename, [lon_bin, lat_bin, time_bin], dimensions)
                day_in_year = 0
                for month in days_in_month:
                    monthly_mean = np.mean(data_global[day_in_year:day_in_year + month,:,:], axis=(0)) 
                    data_grid = global_to_amundsen(monthly_mean, var, grid)
                    data_year.append(data_grid)
                    day_in_year += month
            data_year = np.array(data_year)
            data_ens.append(data_year)

        data.append(np.array(np.mean(data_ens, axis = 0)))
    
    num_days = 12*181

    # Create a time coordinate using xarray with the appropriate frequency (assuming daily data)
    time_coord = pd.date_range(start="1920-01-01", periods=num_days, freq='M')
    
    # Create an xarray dataset
    ds = xr.Dataset(
        {
            "v_wind": (["time", "y", "x"], data[0]),
            "u_wind": (["time", "y", "x"], data[1])
        },
        coords={
            "time": time_coord
        }
    )

    # Save the monthly average dataset to a NetCDF file
    ds.to_netcdf(f"{output_path}pre_industrial_winds.nc")

if __name__ == '__main__':
    calculate_monthly_average()
