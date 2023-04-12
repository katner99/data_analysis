import sys
import os.path
from typing import List

import numpy as np
import netCDF4 as nc
import xarray as xr

from funcs import find_nearest, read_variable_to_timeseries_cut, append_years
from mitgcm_python.grid import Grid

def main():
    ## check that enough arguments have been input by the user
    if len(sys.argv) != 2:
        sys.exit("Stopped - Incorrect number of arguements. Use python PI_append data.py <exp>")
    
    # set up the variables you need
    exp = sys.argv[1]
    start_year = 2070
    n_years = 29
    ensembles = range(1, 5) if exp == "lens" else range(1, 10)
        
    # initialise the final values
    theta_ensemble_mean = np.zeros((len(ensembles),(12*n_years)))
    salt_ensemble_mean = np.zeros((len(ensembles),(12*n_years)))
    melt_ensemble_mean = np.zeros((len(ensembles),(12*n_years)))
    
    # Set up the output data
    output_file = f"/data/oceans_output/shelf/katner33/PIctrl_output/ensemble_mean/{exp}_30_years.nc"
    with nc.Dataset(output_file, "w", format="NETCDF3_CLASSIC") as id_out:
        id_out.createDimension("t", None)
        nc_theta = id_out.createVariable("THETA", "f", ("t"))
        nc_salt = id_out.createVariable("SALT", "f", ("t"))
        nc_melt = id_out.createVariable("MELT", "f", ("t"))
    
        # initialise the timeseries
        theta_timeseries = []
        salt_timeseries = []
        melt_timeseries = []
        
        # save to out file
        id_out.variables['THETA'][:] = theta_timeseries
        id_out.variables['SALT'][:] = salt_timeseries  
        id_out.variables['MELT'][:] = melt_timeseries
        
        # set up grid
        grid_filepath = "/data/oceans_output/shelf/kaight/archer2_mitgcm/PAS_LENS001_O/output/192001/MITgcm/output.nc"
        grid = Grid(grid_filepath)
        grid_file = xr.open_dataset(grid_filepath)
        
        # set lat, lon, and depth range
        depth_range = [find_nearest(grid_file.Z.values, -200), find_nearest(grid_file.Z.values, -700)]
        lon_range = [find_nearest(grid_file.XC.values, 230), find_nearest(grid_file.XC.values, 260)]
        lat_range = [find_nearest(grid_file.YC.values, -75), find_nearest(grid_file.YC.values, -70)]
        

        # run through the ensembles
        for ens in ensembles:
            # load the path of the ensemble mamber
            if exp == "lens":
                ens_str = str(ens).zfill(3)
                filepath = "/data/oceans_output/shelf/kaight/archer2_mitgcm/PAS_LENS"+ens_str+"_noOBC/output/"
            else:
                ens_str = str(ens).zfill(2)
                filepath = "/data/oceans_output/shelf/katner33/PIctrl_output/PAS_"+exp+ens_str+"/output/"
            filename = "output.nc"
            
            # set up the input data
            input_file = filepath+str(start_year)+"01/MITgcm/"+filename
            input_data = xr.open_dataset(input_file)
            
            [theta_timeseries, salt_timeseries, melt_timeseries] = append_years(n_years, start_year, filepath, filename, grid, lat_range, lon_range, depth_range)
                
            # once you have run through all the years save the ensemble members and move on
            theta_ensemble_mean[ens - 1, :] = theta_timeseries.copy()
            salt_ensemble_mean[ens - 1, :] = salt_timeseries.copy()
            melt_ensemble_mean[ens - 1, :] = melt_timeseries.copy()

            theta_timeseries = []
            salt_timeseries = []
            melt_timeseries = []
                    
            print("ensemble member "+str(ens)+" complete")
        
        # save to the netcdf file
        id_out.variables['THETA'][:] = np.nanmean(theta_ensemble_mean, axis = 0)
        id_out.variables['SALT'][:] = np.nanmean(salt_ensemble_mean, axis = 0)
        id_out.variables['MELT'][:] = np.nanmean(melt_ensemble_mean, axis = 0)
    


if __name__ == "__main__":
    main()
    
