from mitgcm_python.file_io import read_binary
from mitgcm_python.grid import Grid
from plots import create_mask

from funcs import binary_coords, find_nearest, read_variable
from directories_and_paths import *
 
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import xarray as xr
import scipy.stats
from scipy.ndimage import zoom

def make_test_grid(plot = False):
    array_size = 100
    x, y, t = np.meshgrid(np.arange(array_size), np.arange(array_size), np.arange(12))
    distance = np.sqrt((x - array_size // 2)**2 + (y - array_size // 2)**2)
    decay_factor = 10e-3
    test_grid = 5 * np.exp(-decay_factor * distance**2)
    test_grid = np.clip(test_grid, 0, 5)
    
    for timestep in np.arange(1,12):
        test_grid[50:,50:,timestep]=test_grid[50:,50:,timestep-1]+0.5
        test_grid[:50,50:,timestep]=test_grid[:50,50:,timestep-1]-0.5
        
    if plot:
    
        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))
        
        ax[0].contourf(test_grid[:,:,0], levels=np.linspace(np.min(test_grid), np.max(test_grid), 20))
        ax[0].set_title("t = 0")
        
        ax[1].contourf(test_grid[:,:,4], levels=np.linspace(np.min(test_grid), np.max(test_grid), 20))
        ax[1].set_title("t = 5")
        
        c = ax[2].contourf(test_grid[:,:,9], levels=np.linspace(np.min(test_grid), np.max(test_grid), 20))
        ax[2].set_title("t = 10")
        fig.colorbar(c, ax=ax[2])
        
        fig.savefig("test_grid.png")
        #plt.show()
    
    return test_grid
    
def calc_trend():
    test_grid = make_test_grid(plot = True)
    
    # calculate the slope
    time = range(12)
    slope_test_grid = np.empty((100, 100))

    for lat in range(100):
        for lon in range(100):
            #print(scipy.stats.linregress(ubot[:, lat, lon], time).pvalue, scipy.stats.linregress(vbot[:, lat, lon], time).pvalue)
            if scipy.stats.linregress(test_grid[lat, lon, :], time).pvalue < 0.05:
                slope_test_grid[lat, lon] = scipy.stats.linregress(test_grid[lat, lon, :], time).slope
     
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 5))
    c = plt.contourf(slope_test_grid)
    fig.colorbar(c)
    
    fig.savefig("test_grid_trend.png")
    
    plt.show()

def old():
    exp = "WIND"
    if exp == "CTRL" or exp == "TEMP":
        domain = "pictrl"
    if exp == "WIND" or exp == "LENS":
        domain = "lens"

    # Load wind velocities
    ubot = xr.open_dataset(f"ubot_{domain}.nc").data.values[:50*12, :, :]
    vbot = xr.open_dataset(f"vbot_{domain}.nc").data.values[:50*12, :, :]

    lat_slope, lon_slope = binary_coords("VBOT")

    [LON, LAT] = np.meshgrid(lon_slope, lat_slope)

    # calculate the slope
    time = range(50*12)
    slope_u = np.empty((14, 48))
    slope_v = np.empty((14, 48))

    for lat in range(14):
        for lon in range(48):
            #print(scipy.stats.linregress(ubot[:, lat, lon], time).pvalue, scipy.stats.linregress(vbot[:, lat, lon], time).pvalue)
            if scipy.stats.linregress(ubot[:, lat, lon], time).pvalue < 0.05 or scipy.stats.linregress(vbot[:, lat, lon], time).pvalue < 0.05:
                slope_u[lat, lon] = scipy.stats.linregress(ubot[:, lat, lon], time).slope
                slope_v[lat, lon] = scipy.stats.linregress(vbot[:, lat, lon], time).slope
                #print("in")

    new_shape = (384, 600)
    print("end")

    # Create a grid of coordinates for the original array
    zoom_factors = (new_shape[0] / slope_u.shape[0], new_shape[1] / slope_u.shape[1])

    # Interpolate and resize the array to the new shape
    u = zoom(slope_u, zoom_factors, order=1)
    v = zoom(slope_v, zoom_factors, order=1)
    
    # download the top of the water column
    filepaths = output_path + exp + "_1920-1950_average.nc"
    input_data = xr.open_dataset(filepaths, decode_times=False)
    lat, lon, ice_mask_temp, depth = [input_data[param].values for param in ["YC", "XC", "maskC", "Depth"]]
    grid = Grid(grid_filepath)
    depth_range = [find_nearest(input_data["Z"].values, -200), find_nearest(input_data["Z"].values, -700)]
    ice_mask = ice_mask_temp[0,:,:]
    [land_mask, mask, colors] = create_mask(depth, ice_mask)
    data = read_variable(input_data, "SIfwfrz", grid)*3600*24*365

    [X, Y] = np.meshgrid(lon, lat)

    fig, ax = plt.subplots(figsize=(18,12))

    #print(np.shape(data))
    
    cs = ax.contourf(X, Y, data, cmap="YlGnBu_r", extend="both", levels=np.linspace(-0.5, 0, 15))
    ax.contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
    ax.contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
    q = ax.quiver(lon[0:-1:20], lat[0:-1:20], u[0:-1:20,0:-1:20], v[0:-1:20,0:-1:20], color = "black")
    qk = ax.quiverkey(q, 0.9, 0.9, 20, r'$20 \frac{m}{s}$', labelpos='E', coordinates='figure')
    ticks=np.arange(-0.5, 0.1, 0.1)
    cbar = fig.colorbar(cs)
    cbar.set_ticks(ticks)
    ax.set_title(exp + " wind and freezing", weight="bold")

    fig.savefig(exp+"_wind_trend_50.png")

    plt.show()

def main():
    domain = "high_emissions"
    data_1 = xr.open_dataset(f"{output_path}{domain}_winds.nc", decode_times = False)
    domain = "pre_industrial"
    data_2 = xr.open_dataset(f"{output_path}{domain}_winds.nc", decode_times = False)

    period_start = 2000
    start_index = (period_start - 1920) * 12

    uwinds = data_2.u_wind.values[start_index:,...]
    vwinds = data_2.v_wind.values[start_index:,...]
    time = data_1.time.values[start_index:,...]
    significance = np.empty(np.shape(uwinds[0,...]))
    slope_u = np.empty(np.shape(uwinds[0,...]))
    slope_v = np.empty(np.shape(uwinds[0,...]))

    # download the top of the water column
    filepaths = f"{output_path}average_LENS_1920-1950.nc"
    input_data = xr.open_dataset(filepaths, decode_times=False)
    lat, lon, ice_mask_temp, depth = [input_data[param].values for param in ["YC", "XC", "maskC", "Depth"]]
    grid = Grid(grid_filepath)
    depth_range = [find_nearest(input_data["Z"].values, -200), find_nearest(input_data["Z"].values, -700)]
    ice_mask = ice_mask_temp[0,:,:]
    [land_mask, mask, colors] = create_mask(depth, ice_mask)

    [X, Y] = np.meshgrid(lon, lat)

    fig, ax = plt.subplots(figsize=(18,12))

    ax.contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
    ax.contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
    
    for lat in range(0, uwinds.shape[1], 40):
        for lon in range(0, uwinds.shape[2], 40):
            slope_u[lat, lon] = scipy.stats.linregress(uwinds[:, lat, lon], time).slope / (100*12)
            slope_v[lat, lon] = scipy.stats.linregress(vwinds[:, lat, lon], time).slope / (100*12)
            if scipy.stats.linregress(uwinds[:, lat, lon], time).pvalue < 0.05 and scipy.stats.linregress(vwinds[:, lat, lon], time).pvalue < 0.05:
                q = ax.quiver(X[lat, lon], Y[lat, lon], slope_u[lat, lon], slope_v[lat, lon], color='black', width=0.005)
            elif scipy.stats.linregress(uwinds[:, lat, lon], time).pvalue < 0.05 or scipy.stats.linregress(vwinds[:, lat, lon], time).pvalue < 0.05:
                q = ax.quiver(X[lat, lon], Y[lat, lon], slope_u[lat, lon], slope_v[lat, lon], color='gray', width=0.005)
            else:
                q = ax.quiver(X[lat, lon], Y[lat, lon],slope_u[lat, lon], slope_v[lat, lon], color='lightgray', width=0.005)

    print(np.max(slope_u))
    ax.quiverkey(q, 0.9, 0.9, 3, r'$3 m/s/century$', labelpos='E', coordinates='figure')

    ax.set_title("wind trend in pre-industrial 2000-2100", weight="bold")

    fig.savefig("wind_trend_pre_industrial_2000-2100.png")

    plt.show()
    
            

if __name__ == '__main__':
    main() # run the program
