import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from netCDF4 import Dataset, MFDataset
import math

# These functions are used both for T point variables use in getData.py and in the velocity variables in 
# getData_vel.py
class forceVariableFile:
    def __init__(self, fileName, varId, nc_id, limits, tInterval, land, units,long_name):
        self.fileName = fileName      # name of file
        self.nc_id = nc_id            # set to -1 to know if the netcdf file is loading properly
        self.varId = varId            # variable name
        self.limits = limits          # upper limit so we can exclude nonsensical data
        self.tInterval = tInterval    # number of data points per day
        self.land = land              # does the variable account for land or not? if True there is data
                                      # over the land, if False the Interpolation will need to create land
                                      # values
        self.units = units            # variable units
        self.long_name = long_name    # brief explanation of the variable

def Interpolate(days_in_month, tInt, data_in, lons, lats, nc_lon, nc_lat, lim, land):
    
    # initialise the output variable 
    data_out = np.zeros((days_in_month, nc_lat.shape[0], nc_lat.shape[1]), dtype=float)
    for t in range(30):
        # create a daily average of the variable
        data_to_use = np.mean(data_in[t*tInt:(t*tInt)+tInt,:,:], axis = 0)
        
        # if the variable does not account for land we need to create dummie values so the coast can be
        # interpolated correctly
        
        if land == False:
            
            # set all unreasonable values to nan  
            data_to_use[data_to_use>lim] = np.nan
            
            # find location of all land values
            land_indices = []
            land_indices = np.isnan(data_to_use) 

            # set all land values to the longitudinal average     
            for iy, ix in np.ndindex(data_to_use.shape):
                if  np.isnan(data_to_use[iy,:]).all(): # if all longitudinal values are nan
                    data_to_use[iy,ix] = 270
                if  math.isnan(data_to_use[iy,ix]): 
                    data_to_use[iy,ix] = np.nanmean(data_to_use[iy,:])
            
            # smooth the longitudinal average values so the coastline is the least affected using a 
            # gaussian filter
            smoothed_data = np.zeros((144, 192), dtype=float)
            smoothed_data[:,:] = gaussian_filter(data_to_use[:,:], sigma = 2, mode = 'wrap')
                  
            # replace land points with the smoothed data
            data_to_use[ land_indices ] = smoothed_data[ land_indices ]
        
        # extend the grid longitudinally by 4 points each end to allow wrapping of data for fitting splines
        data_to_add_before = data_to_use[:,-4:]
        data_to_add_after = data_to_use[:,:4]
        data_in_extended = np.concatenate((data_to_add_before, data_to_use, data_to_add_after), axis = 1)

        # add nonsensical longtitudes to allow for the fit if needed
        lons_extended = np.append(lons[-4:]-360,np.append(lons,lons[:4]+360))

        # extend the grid in latitude and flip data as it moves over the pole
        data_to_add_before = np.flip(data_in_extended[:4,:],axis=0)
        data_to_add_after = np.flip(data_in_extended[-4:,:],axis=0)
        data_in_extended = np.concatenate((data_to_add_before,data_in_extended,data_to_add_after),axis = 0)

        #add nonsensical latitudes to allow for contraining fit if needed
        lats_extended = np.append(lats[-4:]-180,np.append(lats,lats[:4]+180))

        grid_lat,grid_lon = np.meshgrid(lats_extended,lons_extended)

        # then convert these to a list of points
        points = np.column_stack([grid_lat.ravel(), grid_lon.ravel()])

        # create a list of the values that match the points fron the input data
        vals = []
        for ix in range(0,len(lons_extended)):
            for iy in range(0,len(lats_extended)):
                vals.append(data_in_extended[iy,ix])

        #vals = [[vals.append(data_in_extended[iy,ix]) for ix in range(0,len(lons_extended))] for iy in range(0,len(lats_extended))]
        vals = np.array(vals) # and convert to array
        
        # interpolate the data into the ORCA grid   
        data_out[t,:,:] = griddata(points,vals,(nc_lat, nc_lon), method='cubic')
        
        # return the output data according to the number of days in month
        if (days_in_month == 28 and t == 27) or (days_in_month == 30 and t == 29):
            return data_out
    
    # if the month has 31 days we will need to create dummie data for the last day, we do this by copying 
    # the penultimate day
    data_out[30,:,:] = data_out[28,:,:]
    return data_out

# ------------------ 01000 NORTH POLE COHERANCE -------------------------------------------
def NorthPoleWindVectorCoherance(u_data_in, v_data_in, data_in_lon, data_in_lat):

    # Adjust north-pole row for cohereance
    
    for t in range(0,30): # each Met Office month is 30 days long, so run for each dat
        
        # take the top row at the north pole
        u_nth = u_data_in[t,-1,:]
        v_nth = v_data_in[t,-1,:]
        
        # respective co-ords
        lon_nth = np.array(data_in_lon)[:]
        lat_nth = np.array(data_in_lat)[-1:]
        u_nth_pol, v_nth_pol = ConvertRowOfVectorsToPolar(u_nth, v_nth, lon_nth)
        u_mean = np.mean(u_nth_pol)
        v_mean = np.mean(v_nth_pol)

        u_nth_corr, v_nth_corr = MakeARowFromAvgFromPolar(u_mean, v_mean, lon_nth)

        # these two rows are then used to over-write the last row (at the north pole)
        u_data_in[t,-1,:] = u_nth_corr
        v_data_in[t,-1,:] = v_nth_corr

    # returns the adjusted input data
    return u_data_in, v_data_in

def fx(glam, gphi):
    # 1. calculate the radius from the centre of the earth to get the correct curvature of the sphere
    # 2. using the current zonal components and the radius project onto the sphere and find the new polar
    #    components. Thus the following gives the angular coefficient of the zonal vestor
    res = 2 * np.cos( np.radians(glam) ) * np.tan( np.radians(180)/4. - np.radians(gphi)/2)
    return res
def fy(glam, gphi):
    # z = double( 2. * sin( rad*plam ) * tan( rpi/4. - rad*pphi/2. ))
    # does the same as the function above, but for the meridional component
    res = 2 * np.sin( np.radians(glam) ) * np.tan( np.radians(180)/4. - np.radians(gphi)/2)
    return res
    
def ConvertRowOfVectorsToPolar(vect_u, vect_v, lngs):
    new_u = vect_u * np.cos(np.radians(lngs)) - vect_v * np.sin(np.radians(lngs))
    new_v = vect_v * np.cos(np.radians(lngs)) + vect_u * np.sin(np.radians(lngs))
    return(new_u, new_v)

def MakeARowFromAvgFromPolar(u_val, v_val, lngs):
    new_u = u_val * np.cos(np.radians(lngs)) + v_val * np.sin(np.radians(lngs))
    new_v = v_val * np.cos(np.radians(lngs)) - u_val * np.sin(np.radians(lngs))
    return(new_u, new_v)