import numpy as np
import sys
import xarray as xr

from mitgcm_python.utils import mask_land_ice, mask_3d, add_time_dim, z_to_xyz, apply_mask
from mitgcm_python.calculus import over_area
from mitgcm_python.file_io import read_binary

def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
    
def average_over_depth(grid, data, mask, depth_range):
    dz = z_to_xyz(grid.dz, grid)[depth_range[0]:depth_range[1],: ,:]
    dz = add_time_dim(dz, data.shape[0])
    hfac = grid.hfac[depth_range[0]:depth_range[1], :, :]
    hfac = add_time_dim(hfac, data.shape[0])
    return np.sum(data*dz*hfac*mask, axis=-3)/np.sum(dz*hfac*mask, axis=-3)

def read_variable(input_data, var, grid, depth_range = None):
    # Number of days in each month
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    if var == "THETA":
        data = mask_3d(input_data.THETA.values, grid, time_dependent=True)
        data_cut = data[:, depth_range[0]:depth_range[1], :, :]
        mask_cut = np.invert(data.mask).astype(float)[:, depth_range[0]:depth_range[1], :, :]
        return np.average(average_over_depth(grid, data_cut, mask_cut, depth_range), axis = 0, weights=days_in_month)

    elif var == "SALT":
        data = mask_3d(input_data.SALT.values, grid, time_dependent=True)
        return np.average(data[:,0,:,:], axis = 0, weights=days_in_month)
        
    else:
        hfac = grid.hfac[0,:,:]
        hfac = add_time_dim(hfac, input_data.time.values.shape[0])
        data = np.ma.masked_where(hfac == 0, input_data[var].values[..., :, :])
        return data[0,:,:]

def create_profile(input_data, var, grid, lat_range, lon_range, time = 12, timeseries = True):
    # Check if input_data contains 12 years of data
    if len(input_data.time.values) != 12:
        print('Error: input_data must contain {} years of data'.format(time))
        fillarray = np.full((12, 50), np.nan)
        return fillarray
    
    else:
        # Number of days in each month
        days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    
        # apply mask and cut over the area defined
        data = mask_3d(input_data[var].values, grid, time_dependent=True)
        data = apply_mask(data, (grid.get_hfac(gtype='t')==0) & (input_data.Depth.values >= 1500), time_dependent=True, depth_dependent=True)
        data_cut = data[:, :, lat_range[0]:lat_range[1], lon_range[0]:lon_range[1]]
        mask_cut = np.invert(data[:, :, :, :].mask).astype(float)[:, :, lat_range[0]:lat_range[1], lon_range[0]:lon_range[1]]
        dA = grid.dA
        dV = add_time_dim(dA, data.shape[1])
        dV = add_time_dim(dV, data.shape[0])    
        dV_cut = dV[:, :, lat_range[0]:lat_range[1],lon_range[0]:lon_range[1]]
    
        # calculate and return average, either a timeseries or average over the year
        if timeseries:
            return np.sum(data_cut*dV_cut*mask_cut, axis=(-2,-1))/np.sum(dV_cut*mask_cut, axis=(-2,-1))
        else:
            return np.average(np.sum(data_cut*dV_cut*mask_cut, axis=(-2,-1))/np.sum(dV_cut*mask_cut, axis=(-2,-1)), axis = 0, weights=days_in_month)
        
    
def interpolate_currents(V, zon_or_mer):
    V_interp = np.empty(np.shape(V))
    if zon_or_mer == "zonal":
        V_interp[...,:-1] = 0.5*(V[...,:-1] + V[...,1:]) 
        V_interp[...,-1] = V[...,-1]
    elif zon_or_mer == "meridional":
        V_interp[...,:-1,:] = 0.5*(V[...,:-1,:] + V[...,1:,:])
        V_interp[...,-1,:] = V[...,-1,:]
    return V_interp

# function to create timeseries over a given area slice (lat_range, lon_range)
def make_timeseries(var, input_data, grid, lat_range=None, lon_range=None, depth_range=None, time=12): 
    # Check if depth_range is set for THETA variable
    if var == 'THETA' and depth_range is None: 
        print('Error: depth_range must be set for THETA variable')
        sys.exit()
    
    # Check if input_data contains 12 years of data
    if len(input_data.time.values) != time:
        print('Error: input_data must contain {} years of data'.format(time))
        fillarray = np.full(12, np.nan)
        return fillarray

    # Prepare the data and grid for volume average
    if var == "THETA":
        data = mask_3d(input_data.THETA.values, grid, time_dependent=True)
        data_cut = data[:, depth_range[0]:depth_range[1], lat_range[0]:lat_range[1], lon_range[0]:lon_range[1]]
        dV = grid.dV
        dV = add_time_dim(dV, data.shape[0])
        dV_cut = dV[:, depth_range[0]:depth_range[1], lat_range[0]:lat_range[1], lon_range[0]:lon_range[1]]
        mask_cut = np.invert(data.mask).astype(float)[:, depth_range[0]:depth_range[1], lat_range[0]:lat_range[1], lon_range[0]:lon_range[1]]
        return np.sum(data_cut * dV_cut * mask_cut, axis=(-3,-2,-1)) / np.sum(dV_cut * mask_cut, axis=(-3,-2,-1))
    
    # Prepare the data and grid for 2D average

    elif var == "SALT":
        data = mask_3d(input_data[var].values, grid, time_dependent=True)
        data_cut = data[:, 0, lat_range[0]:lat_range[1], lon_range[0]:lon_range[1]]
        mask_cut = np.invert(data[:, 0, :, :].mask).astype(float)[:, lat_range[0]:lat_range[1], lon_range[0]:lon_range[1]]
        dA = grid.dA
        dA = add_time_dim(dA, data.shape[0])
        dA_cut = dA[:, lat_range[0]:lat_range[1],lon_range[0]:lon_range[1]]
        return np.sum(data_cut*dA_cut*mask_cut, axis=(-2,-1))/np.sum(dA_cut*mask_cut, axis=(-2,-1))
        
    else:
        hfac = grid.hfac[0,:,:]
        hfac = add_time_dim(hfac, input_data.time.values.shape[0])
        data = np.ma.masked_where(hfac == 0, input_data.SIheff.values)
        data_cut = data[:, lat_range[0]:lat_range[1], lon_range[0]:lon_range[1]]
        mask_cut = np.invert(data.mask).astype(float)[:, lat_range[0]:lat_range[1], lon_range[0]:lon_range[1]]
        dA = grid.dA
        dA = add_time_dim(dA, data.shape[0])
        dA_cut = dA[:, lat_range[0]:lat_range[1],lon_range[0]:lon_range[1]]       
        return np.nanmax(data_cut, axis=(-2, -1))
        #return np.sum(data_cut*dA_cut*mask_cut, axis=(-2,-1))/np.sum(dA_cut*mask_cut, axis=(-2,-1))

def mask(input_data, grid, time_dependent=True):
    hfac = grid.hfac[0,:,:]
    hfac = add_time_dim(hfac==0, 12)
    data = np.ma.masked_where(hfac, input_data)
    return data

def append_years(n_years, start_year, filepath, filename, grid, lat_range, lon_range, depth_range = None):
    # initialise the years
    #hovmoller_member = np.zeros((12 * n_years, 50))
    theta_member = []
    salt_member = []
    seaice_member = []

    print(start_year)

    #run through the years
    for i in range(n_years):
        # read file of that year
        fileyear=str(start_year+i)
        input_file = f"{filepath}{fileyear}01/MITgcm/{filename}"
        
        # check that the year exists
        try:
            input_data = xr.open_dataset(input_file)
        except FileNotFoundError:
            sys.exit(f"Stopped - Could not find directory {input_file}")
           
        # create timeseries
        theta = make_timeseries("THETA", input_data, grid, lat_range, lon_range, depth_range)
        salt = make_timeseries("SALT", input_data, grid, lat_range, lon_range)
        seaice = make_timeseries("SIheff", input_data, grid, lat_range, lon_range)
        #ovmoller = create_profile(input_data, "THETA", grid, lat_range, lon_range)
        
           
        #hovmoller_member[i * 12 : (i + 1) * 12,:] = hovmoller
        theta_member.extend(theta)
        salt_member.extend(salt)
        seaice_member.extend(seaice)

        # print progress
        percent = int(((i + 1) / n_years) * 100)
        print(f"{percent}% complete")

    return theta_member, salt_member, seaice_member
    #return hovmoller_member
