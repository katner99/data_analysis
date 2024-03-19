import sys
import numpy as np
import xarray as xr
from mitgcm_python.utils import add_time_dim
from directories_and_paths import *
from config_options import *

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
    hfac = grid.hfac[:, lat_range, lon_range[0] : lon_range[1]] == 0
    hfac = add_time_dim(hfac, 12)
    VEL = np.ma.masked_where(hfac, vel)
    if option == "south":
        VEL = np.ma.masked_where(VEL > 0, VEL)
    return sv * np.sum(-VEL * dX * dZ, axis=(-2, -1))


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


def calc_melt(input_data, grid):
    """
    function to calculate the total melt over the area

    INPUT:
    dA (grid): grid cell area (m^2)
    SHIfwflx (iput_data): ice shelf freshwater flux (kg/m^2/s)

    OUTPUT:
    melt: melt over one year (Gt/yr)

    """
    dA = grid.dA
    dA = add_time_dim(dA, 12)
    melt_flux = (
        input_data.SHIfwFlx.values[:12, ...] * 3600 * 24 * 365
    )  # convert from seconds to years

    return 10 ** (-12) * np.sum(-melt_flux * dA)


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
        melt_temp = calc_melt(input_data, grid)
        melt = np.ma.append(melt, melt_temp)

    return melt
