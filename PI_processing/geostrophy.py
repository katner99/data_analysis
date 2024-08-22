"""
Geostrophic Analysis Script: Calculates and visualizes the geostrophic component from oceanographic data.
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from .tools.funcs import chop_slice, read_data
from .directories_and_paths import output_path
from .config_options import LAT_SLICES

def pressure_gradient(exp):
    """
    Calculate the pressure gradient given the experiment.
    
    Parameters:
    - exp (str): Experiment name (e.g., "LENS", "CTRL", "WIND", "TEMP")
    
    Returns:
    - pressure (ndarray): Calculated pressure at the cell center.
    - z (ndarray): Depth at the cell center.
    """
    filepath = f"{output_path}{exp}_files_temp/DENSITY_trend.nc"
    dataset = xr.open_dataset(filepath)
    
    density = dataset.trend.values[:, :, 0]
    g = 9.81  # Acceleration due to gravity (m/s^2)
    z = dataset.z.values

    pressure = np.zeros_like(density)

    # Calculate pressure at each depth
    pressure[0, :] = g * density[0, :] * z[0]
    for depth in range(1, len(z)):
        pressure[depth, :] = pressure[depth - 1, :] + g * density[depth, :] * (z[depth] - z[depth - 1])

    return pressure, z

def geostrophy(pressure, z):
    """
    Calculate the geostrophic component from the pressure gradients.
    
    Parameters:
    - pressure (ndarray): Hydrostatic pressure at each cell.
    - z (ndarray): Depth at cell center.
    
    Returns:
    - geostrophy_dataset (xr.Dataset): Dataset containing the masked values for pressure and geostrophy.
    """
    reference_filepath = f"{output_path}/CTRL_ens01_noOBC/output/192001/MITgcm/output.nc"
    reference_dataset = xr.open_dataset(reference_filepath, decode_times=False)
    
    omega = 7.2921e-5  # Earth's rotation rate (rad/s)
    reference_density = reference_dataset.rhoRef.values
    vel_lat = reference_dataset.YG.sel(YG=slice(LAT_SLICES[0], LAT_SLICES[1])).values
    pres_lat = reference_dataset.YC.sel(YC=slice(LAT_SLICES[0], LAT_SLICES[1])).values
    distance = chop_slice(reference_dataset, "dyC")
    mask = chop_slice(reference_dataset, "hFacC", 2)

    f = 2 * omega * np.sin(np.radians(vel_lat))  # Coriolis parameter

    baroclinic_velocity = np.zeros_like(pressure)

    # Calculate geostrophic velocity
    for lat in range(len(vel_lat)):
        dP = pressure[:, lat + 1] - pressure[:, lat]
        dY = distance[lat]
        dPdY = dP / dY
        baroclinic_velocity[:, lat] = -dPdY / (reference_density * f[lat])

    # Apply mask to the results
    baroclinic_velocity = np.ma.masked_where(mask == 0, baroclinic_velocity)
    pressure = np.ma.masked_where(mask == 0, pressure)

    geostrophy_dataset = xr.Dataset(
        {
            'pressure': (['Z', 'YC'], pressure),
            'geostrophy': (['Z', 'YC'], baroclinic_velocity),
        },
        coords={'Z': z, 'YC': pres_lat}
    )

    return geostrophy_dataset

def shear(exp):
    """
    Calculate and save the surface velocity shear for the given experiment.
    
    Parameters:
    - exp (str): Experiment name (e.g., "LENS", "CTRL", "WIND", "TEMP")
    """
    data = xr.open_dataset(f"{output_path}{exp}_files_temp/UVEL_trend.nc")
    total_vel = data.trend.values[:, :, 0]
    z = data.z.values

    surface_vel = np.repeat(total_vel[0, :][np.newaxis, :], len(z), axis=0)

    dataset = xr.Dataset(
        {'surface_vel': (['z', 'lat'], surface_vel)},
        coords={'z': z, 'lat': data.lat.values}
    )

    dataset.to_netcdf(f"{output_path}{exp}_files_temp/surface_vel.nc")

def plot_experiments(axs, row, graph_params, lat, z, data):
    """
    Plot the data for different experiments in subplots.

    Parameters:
    - axs (array): Array of axes for subplots.
    - row (int): Row index for the subplot.
    - graph_params (dict): Parameters for the graph.
    - lat (ndarray): Latitude values.
    - z (ndarray): Depth values.
    - data (list): List of data arrays to plot.
    
    Returns:
    - cv (QuadContourSet): Contour set for color bar.
    """
    for col in range(3):
        print(row, col)
        cv = axs[row, col].contourf(
            lat, z, data[col], cmap=graph_params['color_scheme'], extend="both",
            levels=np.linspace(graph_params['low_val'], graph_params['high_val'], graph_params['step'])
        )
        axs[row, col].set_ylim([-1000, 0])
        axs[row, col].set_xlim([-74, -71.5])
        axs[row, col].tick_params(axis='x', rotation=45)

        if col > 0:
            axs[row, col].get_yaxis().set_visible(False)

        if row == 0:
            axs[row, col].set_title(graph_params['experiment'][col], weight="bold", fontsize=graph_params['font_size'])
            axs[row, col].get_xaxis().set_visible(False)
        else:
            axs[row, col].get_xaxis().set_visible(False)
    
    return cv

def mask_values(vel):
    """
    Apply mask to velocity data to remove land areas.
    
    Parameters:
    - vel (list): List of velocity arrays to mask.
    
    Returns:
    - masked_vel (list): List of masked velocity arrays.
    """
    reference_filepath = f"{output_path}/CTRL_ens01_noOBC/output/192001/MITgcm/output.nc"
    reference_dataset = xr.open_dataset(reference_filepath, decode_times=False)
    mask = chop_slice(reference_dataset, "hFacC", 2)
    masked_vel = [np.ma.masked_where(mask == 0, vel_exp) for vel_exp in vel]
    
    return masked_vel

def plot_geostrophy():
    """
    Generate and save geostrophy analysis plots.
    """
    options_to_plot = ["total_vel", "sheer_vel", "sheer_calc", "surface_vel", "total_calc"]
    
    # Create subplots with each variable on a new row
    fig, axs = plt.subplots(nrows=len(options_to_plot), ncols=3, 
                            gridspec_kw={"hspace": 0.05, "wspace": 0.05}, 
                            figsize=(12, 5 * len(options_to_plot)))

    for idx, option in enumerate(options_to_plot):
        if option == "total_vel":
            input_data, vel = read_data("UVEL", 0)
            z, lat = input_data[0].z.values, input_data[0]["lat"]
            title = "Modelled zonal velocity"

        elif option in ["surface_vel", "sheer_vel"]:
            filepaths = [f"{output_path}{exp}_files_temp/surface_vel.nc" for exp in ["LENS", "TEMP", "WIND"]]
            input_data = [xr.open_dataset(filepath, decode_times=False) for filepath in filepaths]
            surface_vel = [input.surface_vel.values * 100 for input in input_data]

            vel = mask_values(surface_vel)
            title = "Modelled surface velocity"
            if option == "sheer_vel":
                _, total_vel = read_data("UVEL", 0)
                sheer_vel = [total_vel[i] - surface_vel[i] for i in range(3)]
                vel = mask_values(sheer_vel)
                title = "Modelled shear velocity"

        elif option in ["sheer_calc", "total_calc"]:
            filepaths = [f"{output_path}{exp}_files_temp/geostrophy.nc" for exp in ["LENS", "TEMP", "WIND"]]
            input_data = [xr.open_dataset(filepath, decode_times=False) for filepath in filepaths]
            vel = [-input.geostrophy.values * 100 for input in input_data]
            title = "Geostrophic shear from \n density contours"

            if option == "total_calc":
                try:
                    surface_vel = mask_values(surface_vel)
                except NameError:
                    input_data = [xr.open_dataset(filepath, decode_times=False) for filepath in filepaths]
                    surface_vel = [input.surface_vel.values * 100 for input in input_data]
                    surface_vel = mask_values(surface_vel)
                    z, lat = input_data[0].z.values, input_data[0]["lat"]

                vel = [vel[i] + surface_vel[i] for i in range(3)]
                title = "Total velocity from \n density contours"

        # Graph parameters
        graph_params = {
            'color_scheme': "PiYG_r",
            'experiment': ["ALL", "THERMO", "WIND"],
            'low_val': -5e-05,
            'high_val': 5e-05,
            'step': 15,
            'font_size': 20,
        }

        # Plot the data
        cv = plot_experiments(axs, idx, graph_params, lat, z, vel)
        axs[idx, 0].set_ylabel(title, fontsize=graph_params['font_size'] - 4, weight="bold")

    # Adjust and save the plot
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
    fig.colorbar(cv, cax=cbar_ax)

    plt.suptitle("Geostrophy analysis (m s$^{-1}$ century$^{-1}$)", 
                 fontsize=graph_params['font_size'] + 2, weight="bold", y=0.925)

    fig.savefig("geostrophy.png")

def main():
    """
    Main function to calculate and save geostrophic data.
    """
    exp = "WIND"
    pressure, z = pressure_gradient(exp)
    geostrophy_dataset = geostrophy(pressure, z)
    geostrophy_dataset.to_netcdf(f"{output_path}{exp}_files_temp/geostrophy.nc")

if __name__ == "__main__":
    plot_geostrophy()
