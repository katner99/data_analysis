import sys
import datetime
from dateutil.relativedelta import relativedelta
import numpy as np
import pandas as pd
import xarray as xr
from funcs import find_nearest, append_years
from mitgcm_python.grid import Grid

def main():
    ## check that enough arguments have been input by the user
    if len(sys.argv) != 3:
        sys.exit("Stopped - Incorrect number of arguements. Use python PI_append data.py <exp> <loc> <hovmoller>")
    
    # set up the variables you need
    exp = sys.argv[1]
    loc = sys.argv[2]
    start_year = 2070
    n_years = 31
    ensembles = range(1, 6) if exp == "lens" else range(1,10)
        
    # initialise the final values
    #theta_timeseries = np.zeros((len(ensembles),(12*n_years)))
    #salt_timeseries = np.zeros((len(ensembles),(12*n_years)))
    #seaice_timeseries = np.zeros((len(ensembles),(12*n_years)))
    theta_hovmoller = np.zeros((len(ensembles),(12*n_years), 50))
    
    # Set up the output data
    output_file = f"/data/oceans_output/shelf/katner33/PIctrl_output/ensemble_mean/{exp}_ensemble_mean_hov.nc"
    
    # initialise the single ensemble member variables
    hovmoller_member = np.zeros((12*n_years, 50))
    #theta_member = []
    #salt_member = []
    #seaice_member = []
        
    # set up grid
    grid_filepath = "/data/oceans_output/shelf/kaight/archer2_mitgcm/PAS_LENS001_O/output/192001/MITgcm/output.nc"
    grid = Grid(grid_filepath)
    grid_file = xr.open_dataset(grid_filepath)
    
    # set lat, lon, and depth range
    depth_range = [find_nearest(grid_file.Z.values, -200), find_nearest(grid_file.Z.values, -700)]
    print(depth_range)
    if loc == "PIG":
        lon_range = [find_nearest(grid_file.XC.values, 250), find_nearest(grid_file.XC.values, 260)]
    elif loc == "DG":
        lon_range = [find_nearest(grid_file.XC.values, 235), find_nearest(grid_file.XC.values, 250)]
    else:
        print("ERROR! The location you entered is not included choose Pine Island Bay (PIG) or Dotson-Getz (DG)")
        sys.exit()
    lat_range = [find_nearest(grid_file.YC.values, -76), find_nearest(grid_file.YC.values, -72)]
    
    # run through the ensembles
    for ens in ensembles:
        # load the path of the ensemble member
        ens_str = str(ens).zfill(3) if exp == "lens" else str(ens).zfill(2)
        filepath = "/data/oceans_output/shelf/kaight/archer2_mitgcm/PAS_LENS" + ens_str + "_noOBC/output/" if exp == "lens" else "/data/oceans_output/shelf/katner33/PIctrl_output/PAS_" + exp + ens_str + "/output/"
        filename = "output.nc"
            
        #[hovmoller_member, theta_member, salt_member, seaice_member] = append_years(n_years, start_year, filepath, filename, grid, lat_range, lon_range, depth_range)
        hovmoller_member = append_years(n_years, start_year, filepath, filename, grid, lat_range, lon_range, depth_range)
            
        # once you have run through all the years save the ensemble members and move on
        theta_hovmoller[ens - 1, :, :] = hovmoller_member[:, :]
        #theta_timeseries[ens - 1, :] = theta_member[:]
        #salt_timeseries[ens - 1, :] = salt_member[:]
        #seaice_timeseries[ens - 1, :] = seaice_member[:]
                    
        #print(np.shape(theta_timeseries))
        print("ensemble member "+str(ens)+" complete")

    # create time dimension
    start_date = datetime.date(start_year, 1, 1)
    end_date = datetime.date(2100, 12, 31)
    #time = [start_date + relativedelta(months=i) for i in range((end_date.year - start_date.year) * 12 + end_date.month - start_date.month + 1)]
    #times = pd.date_range("2070/01/01","2101/01/01",freq='M')
    #times2 = np.arange(len(theta_timeseries))
    #print(np.shape(times), np.shape(times2),len(theta_timeseries))
    
    ds = xr.Dataset(
        data_vars=dict(
            theta_hovmoller=(["time", "z"], np.nanmean(theta_hovmoller, axis = 0), {'units':'C'}),
            theta_std_hovmoller=(["time", "z"], np.nanstd(theta_hovmoller, axis = 0), {'units':'C'}),
            #theta_max=(["time"], np.nanmax(theta_timeseries, axis = 0), {'units':'C'}),
            #theta_min=(["time"], np.nanmin(theta_timeseries, axis = 0), {'units':'C'}),
            #theta_timeseries=(["time"], np.nanmean(theta_timeseries, axis = 0), {'units':'C'}),
            #salt_timeseries=(["time"], np.nanmean(salt_timeseries, axis = 0),{'units':''}),
            #seaice_max=(["time"], np.nanmax(seaice_timeseries, axis = 0), {'units':'C'}),
            #seaice_min=(["time"], np.nanmin(seaice_timeseries, axis = 0), {'units':'C'}),
            #seaice_timeseries=(["time"], np.nanmean(seaice_timeseries, axis = 0),{'units':'m'})
        ),
        coords=dict(
            depth=(["z"], grid_file.Z.values),
            time=np.arange(12*n_years),
            #reference_time=pd.Timestamp("2070-1-1")
        ),
        attrs=dict(description="Timeseries and hovmoller data."))
    
    #ds['time'] = ds['time'].apply(pd.to_datetime)

    ds.to_netcdf(output_file, 'w')

if __name__ == "__main__":
    main()
    
