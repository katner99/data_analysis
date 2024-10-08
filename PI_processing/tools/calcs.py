import sys
import numpy as np
import xarray as xr
from mitgcm_python.utils import add_time_dim
from tools.directories_and_paths import *
from config_options import *

def monthly_average(data):
    time = np.shape(data)[0]
    reshape_time = int(time/12)
    reshaped_array = data.reshape(reshape_time, 12, 384, 600)
    result_array = reshaped_array.mean(axis=1)
    return result_array

def calc_transport(input_data, grid, lat_range, lon_range, option="total"):
    """
    function to calculate the transport over a given longitudinal range
    """
    dx = grid.dx_s[lat_range, lon_range[0] : lon_range[1]]
    dz = grid.dz
    dX, dZ = np.meshgrid(dx, dz)
    dX = add_time_dim(dX, 12)
    dZ = add_time_dim(dZ, 12)
    vel = input_data.VVEL.values[:12, :, lat_range, lon_range[0] : lon_range[1]]
    mask = grid.hfac[:, lat_range, lon_range[0] : lon_range[1]] == 0
    mask = add_time_dim(mask, 12)
    VEL = np.ma.masked_where(mask, vel)
    if option == "south":
        VEL = np.ma.masked_where(VEL > 0, VEL)
    return SV * np.sum(-VEL * dX * dZ, axis=(-2, -1))


def append_transport(n_years, start_year, filepath, grid, lat_range, lon_range):
    """
    function to create a transport timeseries
    """
    transport = []
    transport_south = []

    for i in range(n_years):
        # read file of that year
        fileyear = str(start_year + i)
        print(fileyear)
        input_file = f"{filepath}output/{fileyear}01/MITgcm/output.nc"
        try:
            input_data = xr.open_dataset(input_file, decode_times=False)
        except FileNotFoundError:
            sys.exit(f"Stopped - Could not find directory {input_file}")
        transport_temp = calc_transport(input_data, grid, lat_range, lon_range)
        transport_south_temp = calc_transport(
            input_data, grid, lat_range, lon_range, option="south"
        )
        transport = np.ma.append(transport, transport_temp)
        transport_south = np.ma.append(transport_south, transport_south_temp)

    return transport, transport_south


def calc_melt(input_data, grid, lon_range, lat_range):
    """
    function to calculate the total melt over the area

    INPUT:
    dA (grid): grid cell area (m^2)
    SHIfwflx (iput_data): ice shelf freshwater flux (kg/m^2/s)

    OUTPUT:
    melt: melt over one year (Gt/yr)

    """
    dA = grid.dA
    dA_cut = add_time_dim(dA, 12)[:, lat_range[0] : lat_range[1], lon_range[0] : lon_range[1]]
    melt_flux = (
        input_data.SHIfwFlx.values[:12, :lat_range[1], lon_range[0] : lon_range[1]] * 3600 * 24 * 365
    )  # convert from seconds to years 
    return 10 ** (-12) * np.sum((-melt_flux) * dA_cut, axis=(-2, -1))


def append_melt(n_years, start_year, filepath, grid):
    """
    function to create a transport timeseries
    """
    melt = []

    for i in range(n_years):
        # read file of that year
        fileyear = str(start_year + i)
        print(fileyear)
        input_file = f"{filepath}output/{fileyear}01/MITgcm/output.nc"
        try:
            input_data = xr.open_dataset(input_file, decode_times=False)
        except FileNotFoundError:
            sys.exit(f"Stopped - Could not find directory {input_file}")

        [lat, lon] = [input_data[param].values for param in ["YC", "XC"]]
        lon_range = [find_nearest(lon, 250), find_nearest(lon, 260)]
        lat_range = [find_nearest(lat, -76), find_nearest(lat, -72)]

        melt_temp = calc_melt(input_data, grid, lon_range, lat_range)
        melt = np.ma.append(melt, melt_temp)
    print(np.shape(melt))
    return melt

def moving_average(a, n=3):
    """calculate a centred moving average"""
    if n < 2:
        raise ValueError(
            "Window size (n) must be at least 2 for a centered moving average."
        )

    data = np.empty_like(a, dtype=float)
    data.fill(np.nan)

    # Calculate the cumulative sum
    cumsum = np.cumsum(np.insert(a, 0, 0))

    # Calculate the centered moving average
    half_n = n // 2
    if n % 2 == 0:
        data[half_n - 1 : -half_n] = (cumsum[n:] - cumsum[:-n]) / n
    else:
        data[half_n:-half_n] = (cumsum[n:] - cumsum[:-n]) / n

    return data