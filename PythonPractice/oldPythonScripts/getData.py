#!/usr/bin/env python
import sys
import os.path
import logging # print but better
import numpy as np
from netCDF4 import Dataset, MFDataset, date2num
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
import datetime as dt
from getDataFunctions import forceVariableFile, Interpolate
# -------------------------------- 00100 UTILITIES --------------------------------------
# Get command line arguments, arguments should contain the name of the program and the year 
# we are running
if len(sys.argv) != 2:
    sys.exit("Stopped - Incorrect arguements")

# extrapolate the year and create 12 months
year = sys.argv[1]
month = ["01","02","03","04","05","06","07","08","09","10","11","12"]

# set up logging to file 
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='forcing.'+str(year)+'.log',
                    filemode='w')

# define a Handler which writes INFO messages or higher
console = logging.StreamHandler()
console.setLevel(logging.INFO)

# set a format which is simpler for console use
formatter = logging.Formatter('%(levelname)-5s %(message)s')

# tell the handler to use this format
console.setFormatter(formatter)

# add the handler to the root logger
logging.getLogger('').addHandler(console)
log = logging.getLogger("Run")

# Get command line arguments, the following makes sure we have two arguments 
if len(sys.argv) != 2:
    log.error("Command line arguments not correct")
    log.error("Use: makeForcingData.py <year> ")
    sys.exit("Stopped")

# begin studying the year
log.info("Processing for month/year: "+str(year))

#----------------------------------- 00200 DOWNLOAD FILES ---------------------------------------------

# Load the file with the variable of interest for the first month, this will allow us to extrapolate 
# key information such as the variable name, its limits, the number of measurements made in 24 hours,
# the unit, and a general description of the variable. See getDataFunctions.py for more info on 
# how each variable class is organised

metOfficeFiles = []
# for the historic data use:
# path = "/gpfs/data/greenocean/software/resources/MetOffice/ozone/u-bc370_orgGRD/"
# and use the file name "bc370_<file_name>..."

# for the interactive ozone simulation from 2015 onwards:
path = "/gpfs/data/greenocean/software/resources/MetOffice/ozone/u-be682_orgGRD/"
metOfficeFiles.append(forceVariableFile(path+"be682_rad_y"+str(year)+"m","sw_down",-1,2000,8,True,"W m-2","Downward short-wave radiation"))
metOfficeFiles.append(forceVariableFile(path+"be682_bulk_y"+str(year)+"m","tair10m",-1,400,8,False,"K","Temperature"))
metOfficeFiles.append(forceVariableFile(path+"be682_rad_y"+str(year)+"m","lw_down",-1, 2000,8,True,"W m-2","Downward long-waver radiation"))
metOfficeFiles.append(forceVariableFile(path+"be682_bulk_y"+str(year)+"m","qair10m",-1,10,8,True,"kg kg-1","Specific Humidity"))
metOfficeFiles.append(forceVariableFile(path+"be682_bulk_y"+str(year)+"m","precip",-1, 10,8,True,"kg m-2 s-1","Precipitation Rate"))
metOfficeFiles.append(forceVariableFile(path+"be682_bulk_y"+str(year)+"m","snow",-1,10,8,True,"kg m-2 s-1","Snowfall Rate"))

#------------------------------------- 00300 OPEN FILES ---------------------------------------------
log.info("Running script to convert Met Office files to Forcing files for NEMO models")
log.info("Required files: ")

# run through the rad and bulk files and check that there is data available for each month
for v in range(2):
    for m in range(12):
        # checks that the variable appended exists
        if os.path.isfile(metOfficeFiles[v].fileName+month[m]+".nc"): 
            log.info(str(v)+": "+metOfficeFiles[v].fileName+month[m]+".nc"+" found")
        else:
            log.error(metOfficeFiles[v].fileName+month[m]+".nc"+" does not exist")

# import ORCA coordinates
ORCACoordinates=Dataset("mask.nc","r")
nc_lat = ORCACoordinates.variables['nav_lat'][:] 
nc_lon = ORCACoordinates.variables['nav_lon'][:]
log.info("Mesh opened.")

# import Met Office Coordinates
MetOfficeCoordinates=Dataset(metOfficeFiles[1].fileName+month[1]+".nc","r")
lons = MetOfficeCoordinates.variables['longitude'][:]
lats = MetOfficeCoordinates.variables['latitude'][:]

# change longitude to be 0 to 360, not -180 to 180 in both Met Office and ORCA coordinates
for i in range(len(lons)):
    if lons[i] < 0:
        lons[i] = 360 + lons[i]

for iy, ix in np.ndindex(nc_lon.shape):
    if nc_lon[iy,ix] < 0:
        nc_lon[iy,ix] = 360 + nc_lon[iy,ix]

#--------------------------------- 00400 CREATE OUTPUT FILES-------------------------------------

# create a file for all the variables and with the correct dimensions
for v in range(len(metOfficeFiles)):
    outputFile = "MetOffice_"+metOfficeFiles[v].varId+str(year)+".nc"
    nc_out_id = Dataset(outputFile, 'w', format='NETCDF3_CLASSIC')
    log.info("Creating: "+outputFile)
    
    # create dimensions
    out_xDim = nc_lat.shape[1]
    out_yDim = nc_lat.shape[0]

    # Set the available dimensions determined from meshmask file
    nc_out_id.createDimension("x", out_xDim) 
    nc_out_id.createDimension("y", out_yDim) 
    nc_out_id.createDimension("time_counter", None)
    
    # Set dimension attributes
    
    # Set the nav-lon and lat variables
    nc_v_nav_lon = nc_out_id.createVariable("nav_lon", nc_lon.dtype, ("y","x"))
    nc_v_nav_lat = nc_out_id.createVariable("nav_lat", nc_lat.dtype, ("y","x"))
    nc_out_id.variables['nav_lon'][:] = nc_lon
    nc_out_id.variables['nav_lat'][:] = nc_lat
    
    # the limits are duplicated from original Met Office forcing data
    nc_v_nav_lon.setncattr("units", "degrees_east")
    nc_v_nav_lon.setncattr("valid_min", np.array(-179.7507,'f'))
    nc_v_nav_lon.setncattr("valid_max", np.array(180.,'f'))
    nc_v_nav_lon.setncattr("long_name", "Longitude")
    nc_v_nav_lat.setncattr("units", "degrees_east")
    nc_v_nav_lat.setncattr("valid_min", np.array(-78.19058,'f'))
    nc_v_nav_lat.setncattr("valid_max", np.array(89.6139,'f'))
    nc_v_nav_lat.setncattr("long_name", "Latitude")
    
    # Set the time variable
    nc_v_time = nc_out_id.createVariable("time_counter", 'f', ('time_counter'))
    nc_v_time.setncattr("units", "seconds since 1800-01-01 00:00:00")
    nc_v_time.setncattr("calendar", "gregorian")
    nc_v_time.setncattr("title", "time")
    nc_v_time.setncattr("long_name", "   0-JAN-01 00:00:00")
    
    # create a daily time_counter using the year and creating points for 365 days:
    times_since_startOfYear = []
    start_date = dt.datetime(int(year),1,1,0,0)
    unit = "seconds since 1800-01-01"
    start_time = date2num(start_date, unit) # convert from date to seconds since 1800
    times = []
    times = np.arange(0, 365) # create an array of dates
    times = times * (24 * 60 * 60) # convert the array into seconds
    for t in range(0,365):
        times_since_startOfYear.append(start_time+times[t])
            
    # write out into variable
    nc_out_id.variables['time_counter'][:] = times_since_startOfYear
        
    # download variable from the Met Office File
    nc_v = nc_out_id.createVariable(metOfficeFiles[v].varId, float, ('time_counter', 'y', 'x'))
    nc_v.setncattr("long_name", metOfficeFiles[v].long_name)
    nc_v.setncattr("units", metOfficeFiles[v].units)
    nc_v.setncattr("short_name", metOfficeFiles[v].varId)
    nc_v.setncattr("axis", "TYX")
    nc_v.setncattr("interval_operation",  np.array(86400.,'f'))
    nc_v.setncattr("interval_write", np.array(86400.,'f'))
    nc_v.setncattr("associate", "time_counter nav_lat nav_lon")
    nc_v.setncattr("missing_value", np.array(1E20,'f'))

#--------------------------------- 00500 PROCESS DATA -----------------------------------------------
for v in range(len(metOfficeFiles)): 
    log.info("Parameter: "+metOfficeFiles[v].varId)
    outputFile = "MetOffice_"+metOfficeFiles[v].varId+str(year)+".nc"
    
    # append each month to create one file for each variable for the entire year
    nc_out_id = Dataset(outputFile, 'a', format='NETCDF3_CLASSIC')
   
    t=0 # the following time_counter is used to append the data correctly, it is reset for each variable       
    for m in range(12):
        # open the file containing the variable we want
        metOfficeFiles[v].nc_id = Dataset(metOfficeFiles[v].fileName+month[m]+".nc","r")
        
        # set the number of days in the month
        if m==0 or m==2 or m==4 or m==6 or m==7 or m==9 or m==11:
            days_in_month=31
        if m==1:
            days_in_month=28
        if m==3 or m==5 or m==8 or m==10:
            days_in_month=30
        
        # create the new variable dimensions
        new_var = np.zeros((days_in_month,out_yDim,out_xDim),dtype=float)
        data_in = np.zeros((30, len(lons), len(lats)),dtype=float)
            
        # check whether there is data available over land, import the data accordingly
        if metOfficeFiles[v].land == False:
            data_in = metOfficeFiles[v].nc_id.variables[metOfficeFiles[v].varId][:].data
        else:
            data_in = metOfficeFiles[v].nc_id.variables[metOfficeFiles[v].varId][:]
            
        # interpolate the new data
        new_var = Interpolate(days_in_month, metOfficeFiles[v].tInterval, data_in, lons, lats, nc_lon, nc_lat, metOfficeFiles[v].limits, metOfficeFiles[v].land)
        
        # the cubic interpolation creates negative values where there physically cannot be any, so correct as following:
        if metOfficeFiles[v].varId == "precip" or metOfficeFiles[v].varId == "snow" or metOfficeFiles[v].varId == "qair10m":
            new_var[new_var<0] = 0
        
        # save the variable to its new netCDF file
        nc_out_id.variables[metOfficeFiles[v].varId][t:t+days_in_month,:,:] = new_var[:]
        
        # increment the time counter to append the data at the right point
        t=t+days_in_month