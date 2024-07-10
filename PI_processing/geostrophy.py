import xarray as xr
import numpy as np
import math
import matplotlib.pyplot as plt
from funcs import chop_slice
from directories_and_paths import output_path
from config_options import lat_slices, lon_slices


def pressure_gradient(exp):
    """ function to calculate the pressure gradient given the experiment
    INPUT:
    exp: experiment name (expected: "LENS", "CTRL", "WIND", "TEMP")
    OUTPUT:
    pressure: calculated pressure at the cell centre
    z: depth at the cell centre
    """

    filepath = f"{output_path}{exp}_files_temp/DENSITY_trend.nc"
    dataset = xr.open_dataset(filepath)
    density = dataset.trend.values[:,:,0]
    
    g = 9.81                 # acceleration due to gravity (m/s^2)
    z = dataset.z.values     # depth at the cell centre (m)

    pressure = np.zeros(np.shape(density))

    # calculate the first point
    pressure[0,:] = g * density[0,:] * z[0] # equivalent to P(0) = g * rho * z

    # run through all the depths
    for depth in range(len(z)-1):
        # equivalent to P(n) = P(n-1) + g * rho(n) * [z(n) - z(n-1)]
        pressure[depth + 1,:] = pressure[depth,:] + (g * density[depth + 1,:] * (z[depth + 1] - z[depth])) 

    return pressure, z

def geostrophy(pressure, z):
    """ function to calculate the geostrophic component from the density gradients
    INPUT:
    pressure: hydrostatic pressure at each cell
    z: depth at cell centre
    
    OUTPUT:
    geostrophy_dataset: dataset containing the masked values for pressure and geostrophy
    """
    # original dataset for some reference values
    reference_dataset = xr.open_dataset(f"{output_path}/CTRL_ens01_noOBC/output/192001/MITgcm/output.nc", decode_times = False)
    
    omega = 7.2921 * (10*(-5))
    reference_density = reference_dataset.rhoRef.values
    vel_lat = reference_dataset.YG.sel(YG = slice(lat_slices[0], lat_slices[1])).values    # latitude at the cell interface
    pres_lat = reference_dataset.YC.sel(YC = slice(lat_slices[0], lat_slices[1])).values   # latitude at the cell centre
    distance = chop_slice(reference_dataset, "dyC")                                        # length of cell
    mask = chop_slice(reference_dataset, "hFacC", 2)                                       # mask to remove land

    f = 2 * omega * (np.sin(np.radians(vel_lat)))       # coriolis parameter calculated as f = 2*omega*sin(latitude)

    baroclinic_velocity = np.zeros(np.shape(pressure))

    # calculation of geostrophy progressing through latitudes
    
    for lat in range(len(vel_lat)):
        dP = pressure[:, lat + 1] - pressure[:, lat]
        dY = distance[lat]
        dPdY = np.divide(dP, dY)
        baroclinic_velocity[:, lat] = -dPdY / (np.multiply(reference_density,f[lat]))
    
    # mask our values
    baroclinic_velocity = np.ma.masked_where(mask == 0, baroclinic_velocity)
    pressure = np.ma.masked_where(mask == 0, pressure)

    geostrophy_dataset = xr.Dataset(
        {
            'pressure': (['Z', 'YC'], pressure),
            'geostrophy': (['Z', 'YG'], baroclinic_velocity),
        },
        coords={'Z': z, 'YC': pres_lat, 'YG': pres_lat}
    )

    return geostrophy_dataset
    

def main():
    exp = "NONE"
    [pressure, z] = pressure_gradient(exp)
    geostrophy_dataset = geostrophy(pressure, z)
    geostrophy_dataset.to_netcdf(f"{output_path}{exp}_files_temp/geostrophy.nc")

if __name__ == "__main__":
    main()