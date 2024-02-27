import os
import numpy as np
import xarray as xr
from funcs import find_nearest
from directories_and_paths import *
from config_options import *
from mitgcm_python.grid import Grid
from mitgcm_python.utils import add_time_dim

def basic_transport():
    sv = 10**(-6)
    vel = np.array([[1,1,1],[1,1,1]]) #transport should be 6 (?)
    dx = np.ones([3,1])
    dz = np.ones([1,2])
    dX, dZ = np.meshgrid(dx, dz)
    transport = np.sum(vel*dX*dZ)
    print(transport)
    
def grid_transport():
    sv = 10**(-6)   
    input_data = xr.open_dataset(grid_filepath, decode_times = False)
    lat_range = find_nearest(input_data.YC.values, -73)
    lon_range = [find_nearest(input_data.XC.values, 252.8),find_nearest(input_data.XC.values, 255)]
    grid = Grid(grid_filepath)
    dx = grid.dx_s[lat_range, lon_range[0]:lon_range[1]]
    dz = grid.dz
    print(np.shape(dx))
    dX, dZ = np.meshgrid(dx, dz)
    vel = np.ones(np.shape(dZ))
    print(np.shape(dX), dX)
    transport = sv*np.sum(vel*dX*dZ)
    print(transport)

def masked_transport():
    sv = 10**(-6)   
    input_data = xr.open_dataset(grid_filepath, decode_times = False)
    lat_range = find_nearest(input_data.YC.values, -73)
    lon_range = [find_nearest(input_data.XC.values, 252.8),find_nearest(input_data.XC.values, 255)]
    grid = Grid(grid_filepath)
    dx = grid.dx_s[lat_range, lon_range[0]:lon_range[1]]
    dz = grid.dz
    dX, dZ = np.meshgrid(dx, dz)
    vel = np.ones(np.shape(dZ))
    hfac = grid.hfac[:,lat_range, lon_range[0]:lon_range[1]]==0
    VEL = np.ma.masked_where(hfac, vel)
    #print(np.shape(VEL.mask), VEL)
    #print(np.shape(dX), dX)
    transport = sv*np.sum(VEL*dX*dZ)
    print(transport)
    
def year_transport():
    sv = 10**(-6)   
    input_data = xr.open_dataset(grid_filepath, decode_times = False)
    lat_range = find_nearest(input_data.YC.values, -73)
    lon_range = [find_nearest(input_data.XC.values, 252.8),find_nearest(input_data.XC.values, 255)]
    grid = Grid(grid_filepath)
    dx = grid.dx_s[lat_range, lon_range[0]:lon_range[1]]
    dz = grid.dz
    dX, dZ = np.meshgrid(dx, dz)
    dX = add_time_dim(dX, 12)
    dZ = add_time_dim(dZ, 12)
    vel = np.ones(np.shape(dZ))
    hfac = grid.hfac[:,lat_range, lon_range[0]:lon_range[1]]==0
    hfac = add_time_dim(hfac, 12)
    VEL = np.ma.masked_where(hfac, vel)
    transport = sv*np.sum(VEL*dX*dZ, axis=(-2,-1))
    print(np.shape(transport), transport)
    
def test_transport():
    transport = []
    sv = 10**(-6)  
    filepath = "/data/oceans_output/shelf/katner33/PIctrl_output/CTRL_ens01_noOBC/output/192001/MITgcm/output.nc"
    input_data = xr.open_dataset(filepath, decode_times = False)
    lat_range = find_nearest(input_data.YC.values, -73)
    lon_range = [find_nearest(input_data.XC.values, 252.8),find_nearest(input_data.XC.values, 255)]
    grid = Grid(grid_filepath)
    dx = grid.dx_s[lat_range, lon_range[0]:lon_range[1]]
    dz = grid.dz
    dX, dZ = np.meshgrid(dx, dz)
    dX = add_time_dim(dX, 12)
    dZ = add_time_dim(dZ, 12)
    vel = input_data.VVEL.values[...,lat_range, lon_range[0]:lon_range[1]]
    hfac = grid.hfac[:,lat_range, lon_range[0]:lon_range[1]]==0
    hfac = add_time_dim(hfac, 12)
    VEL = np.ma.masked_where(hfac, vel)
    VEL = np.ma.masked_where(VEL>0, VEL)
    transport = sv*np.sum(-VEL*dX*dZ, axis=(-2,-1))
    print(np.shape(transport), transport)

def append_transport(n_years, start_year, filepath, grid, lat_range, lon_range):
    transport = []
    transport_south = []
    
    for i in range(n_years):
        # read file of that year
        fileyear=str(start_year+i)
        print(fileyear)
        input_file = f"{filepath}output/{fileyear}01/MITgcm/output.nc"
        try:
            input_data = xr.open_dataset(input_file, decode_times=False)
        except FileNotFoundError:
            sys.exit(f"Stopped - Could not find directory {input_file}")
        transport_temp = calc_transport(input_data, grid, lat_range, lon_range)
        transport_south_temp = calc_transport(input_data, grid, lat_range, lon_range, option = "south")
        transport = np.ma.append(transport, transport_temp)
        transport_south = np.ma.append(transport_south, transport_south_temp)

    return transport, transport_south

def calc_transport(input_data, grid, lat_range, lon_range, option = "total"):
    dx = grid.dx_s[lat_range, lon_range[0]:lon_range[1]]
    dz = grid.dz
    dX, dZ = np.meshgrid(dx, dz)
    dX = add_time_dim(dX, 12)
    dZ = add_time_dim(dZ, 12)
    vel = input_data.VVEL.values[:12,:,lat_range, lon_range[0]:lon_range[1]]
    hfac = grid.hfac[:,lat_range, lon_range[0]:lon_range[1]]==0
    hfac = add_time_dim(hfac, 12)
    VEL = np.ma.masked_where(hfac, vel)
    if option == "south":
        VEL = np.ma.masked_where(VEL>0, VEL)
    return sv*np.sum(-VEL*dX*dZ, axis=(-2,-1))
    
if __name__ == "__main__":
    test_transport()
