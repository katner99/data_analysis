from funcs import read_variable, find_nearest, create_profile
from mitgcm_python.grid import Grid
from directories_and_paths import *

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import sys
import numpy as np
import xarray as xr

def main():
    ## check that enough arguments have been input by the user
    #if len(sys.argv) != 2:
    #    sys.exit("Stopped - Incorrect number of arguements. Use python PI_comparison.py <var>")
    
    # set up the variables you need
    #var = str(sys.argv[1])
    save = True
    show = True

    # load up the file paths for the monster, needed 4
    filepath_ctrl = "/data/oceans_output/shelf/katner33/PIctrl_output/PAS_wind09/output/209001/MITgcm/output.nc"
    #filepath_lens = "/data/oceans_output/shelf/katner33/PIctrl_output/LENS_ensemble_mean_2090.nc"
    filepath_wind = "/data/oceans_output/shelf/katner33/PIctrl_output/PAS_wind09test/output/209001/MITgcm/output.nc"
    #filepath_temp = "/data/oceans_output/shelf/katner33/PIctrl_output/THERM_ensemble_mean_2090.nc"
    #filepath_ctrl = "/data/oceans_output/shelf/katner33/PIctrl_output/climatology_runs/transient_average.nc"
    #filepath_lens = "/data/oceans_output/shelf/katner33/PIctrl_output/climatology_runs/presentday_climatology.nc"
    #filepath_wind = "/data/oceans_output/shelf/katner33/PIctrl_output/climatology_runs/historic_climatology.nc"
    #filepath_temp = "/data/oceans_output/shelf/katner33/PIctrl_output/climatology_runs/future_climatology.nc"

    # check if the input files exist
    for filepath in [filepath_ctrl, filepath_wind, ]:
        try:
            open(filepath)
        except FileNotFoundError:
            sys.exit(f"Stopped - Could not find input file {filepath}")

    # load up the input data    
    ctrl_input_data = xr.open_dataset(filepath_ctrl)
    #lens_input_data = xr.open_dataset(filepath_lens)
    #temp_input_data = xr.open_dataset(filepath_temp)
    wind_input_data = xr.open_dataset(filepath_wind)
    
    # read in the general variables (these should be the same between the ensembles
    # set lat, lon, and depth range
    grid_file = "/data/oceans_output/shelf/kaight/archer2_mitgcm/PAS_LENS001_O/output/192001/MITgcm/output.nc"
    grid = Grid(grid_file)
    lon_range = [find_nearest(ctrl_input_data.XC.values, 240), find_nearest(ctrl_input_data.XC.values, 260)]
    lat_range = [find_nearest(ctrl_input_data.YC.values, -80), find_nearest(ctrl_input_data.YC.values, -70)]
    
   
    # create the profile
    var = "THETA"
    theta_ctrl_profile = create_profile(ctrl_input_data, var, grid, lat_range, lon_range, timeseries =False)
    #theta_lens_profile = create_profile(lens_input_data, var, grid, lat_range, lon_range, timeseries =False)
    theta_wind_profile = create_profile(wind_input_data, var, grid, lat_range, lon_range, timeseries =False)
    #theta_temp_profile = create_profile(temp_input_data, var, grid, lat_range, lon_range, timeseries =False)
    
    # create the profile
    var = "SALT"
    ctrl_profile = create_profile(ctrl_input_data, var, grid, lat_range, lon_range,timeseries = False)
    #lens_profile = create_profile(lens_input_data, var, grid, lat_range, lon_range,timeseries = False)
    wind_profile = create_profile(wind_input_data, var, grid, lat_range, lon_range,timeseries = False)
    #temp_profile = create_profile(temp_input_data, var, grid, lat_range, lon_range,timeseries = False)
    
    z = ctrl_input_data.Z.values
    
    # create subplots with each variable on a new line
    fig, axs = plt.subplots(nrows=1, ncols=2, gridspec_kw={"hspace": 0.5, "wspace": 0.4})
      
    axs = axs.flatten()
    
    axs[0].plot(theta_ctrl_profile, z, color = 'seagreen')
    #axs[0].plot(theta_lens_profile, z, color = 'orchid')
    axs[0].plot(theta_wind_profile, z, '--', color = 'cornflowerblue')
    #axs[0].plot(theta_temp_profile, z, '--',color = 'orange')
    axs[0].set_title("THETA")
    axs[0].set_ylim([-800, 0])
    axs[0].set_ylabel("Depth")
    axs[0].set_xlabel("Theta (°C)")
    
    axs[1].plot(ctrl_profile, z, color = 'seagreen')
    #axs[1].plot(lens_profile, z, color = 'orchid')
    axs[1].plot(wind_profile, z, '--', color = 'cornflowerblue')
    #axs[1].plot(temp_profile, z, '--', color = 'orange')
    axs[1].set_title("Salinity")
    plt.ylim([-800, 0])
    plt.xlabel("Salinity")

    plt.subplots_adjust(right=0.7)
    plt.legend(["original","New Compiler"], bbox_to_anchor=(1.04, 1), loc="center left")

    fig.savefig("profile_compiled_"+var+".png")
    plt.show()
        

if __name__ == "__main__":
    main()
    
