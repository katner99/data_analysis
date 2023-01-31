###############################################################
# script to produce plots to compare data between two datasets
# created by Katherine Turner based on hovmoller.py script 
# 23 November 2022
###############################################################
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import animation
from PIL import Image,ImageFilter
import numpy as np
import sys
import datetime
import netCDF4 as nc
from mitgcm_python.utils import mask_land_ice
from mitgcm_python.calculus import over_area
from mitgcm_python.file_io import read_binary

################## AVAILABLE PLOTS #############################
# CONTOUR PLOTS:
# contourplots: plots a simple 2D plot
# animate_contour: animates contour plot
# compare_contour_plots: compares multiple plots for one variable
#
# TIMESERIES:
# make_timeseries_at_point: timesereis at a single point
# make timeseries over area: calculates the average over the area and pots as a timeseries
#
# OTHER PLOTS:
# quiver_plot: maps vector fields with or without overlay
#
# OTHER FUNCTIONS:
# find_nearest: finds the index with the value closest to the one wanted

################## CONTOUR PLOTS ################################

# plot contour plot:
# INPUT:
# lon, lat = longitude and latitude
# data = data to be input (requires two dimensions)
# title = graph title
# year = year of the contour plot
# exp = experiment name 
# var = variable looked at
# apply_land_mask, land_mask = do you need a land mask? expects a binary file, if not set assumes False and NaN
# apply_ice_mask, ice_mask = do you need an ice mask? expects a binary file, if not set assumes False and NaN
# show = shows the contour plot, if unset assumes False
# save = saves contour plot as a png, assumes True if unset

def contour_plots (lon, lat, data, title, year, exp, var, apply_land_mask=False, apply_ice_mask=False, land_mask=np.nan, ice_mask=np.nan, show=False, save=True):

    x = lon
    y = lat
    z = data

    # apply mask over land
    if apply_land_mask == True:
        data[land_mask == 0] = np.nan
    
    if apply_ice_mask == True:
        data[ice_mask < 1] = np.nan
    
    # create a 2D grid
    [X, Y] = np.meshgrid(lon, lat)

    # set up figure
    fig =  plt.figure(figsize=(15,10))
    cs = plt.contourf(X, Y, z, cmap = "ocean")
    plt.colorbar(cs)
    plt.title(title)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    
    # show figure
    if show == True:
        plt.show()
        
    # save figure
    if save == True:
        fig.savefig(exp+"_cont_"+year+"_"+var+".png")

# creates gif of contourplot:
# INPUT:
# lon, lat = longitude and latitude
# data = data to be input (requires two dimensions), expects monthly data
# year = year of the contour plot
# exp = experiment name 
# var = variable looked at
# cs = set color scheme
# apply_land_mask, land_mask = do you need a land mask? expects a binary file, if not set assumes False and NaN
# apply_ice_mask, ice_mask = do you need an ice mask? expects a binary file, if not set assumes False and NaN
# show = shows the contour plot, if unset assumes False
# save = saves contour plot as a png, assumes True if unset

def animate_contour (lon, lat, data, year, exp, var, cs, apply_land_mask=False, apply_ice_mask=False, land_mask=np.nan, ice_mask=np.nan, show=False, save=True):
    # prepare grid
    [X, Y] = np.meshgrid(lon, lat)
    
    # prepare values so all months have the same parameters
    low_val = np.nanmin(data)
    print(low_val)
    #low_val = -2
    #high_val = np.nanmax(data)
    #print(high_val)
    high_val = 0.002
    step = (high_val-low_val)/15
    print(step)
    
    # prepares the title
    labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    
    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure(figsize=(8,6))
    ax = plt.axes(xlim=(min(lon), max(lon)), ylim=(min(lat), max(lat)))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    
    # animation function
    def animate(i): 
        z = data[i,:,:]
        
        if apply_ice_mask == True:
            z[ice_mask[i,:,:] == 0] = np.nan
        # apply mask over land
        if apply_land_mask == True:
            z[land_mask == 0] = np.nan
              
        # create frame
        cont = plt.contourf(X, Y, z, np.arange(low_val, high_val, step), cmap=cs)
        #cont = plt.contourf(X, Y, z, cmap=cs)
        
        # write colorbar on the first or it will keep being copied to the figure
        
        plt.title(var+" "+labels[i]+" "+year)
        return cont  

    plt.colorbar(animate(0))
    
    # animate
    anim = animation.FuncAnimation(fig, animate, frames=12, interval = 200)

    if show == True:
        plt.show()
    
    if save == True:
        anim.save(exp+"_cont_"+var+"_"+year+".gif", fps = 2)

# givent two datasets compares the values of one variable
# INPUT:
# lon, lat = longitude and latitude
# data1, data2 = data to be input (requires two dimensions)
# month = month you are looking at
# apply_mask, mask = do you need a land mask? expects a binary file, if not set assumes False and NaN
# title1, title2 = titles ove the data you are looking at
# year = year of the contour plot
# var = variable looked at
# cm = set color scheme
# show = shows the contour plot, if unset assumes False
# save = saves contour plot as a png, assumes True if unset
def compare_contour_plots (lon, lat, data1, data2, month, apply_mask = False, mask = None, title1 = None, title2 = None, var = None, year = None, cm = "ocean", save = True, show = False):

    x = lon
    y = lat
    z1 = data1
    z2 = data2
    m = str(month)

    if apply_mask == True:
        z1[mask == 0] = np.nan
        z2[mask == 0] = np.nan
    
    # create a 2D grid
    [X, Y] = np.meshgrid(lon, lat)

    fig =  plt.figure(figsize=(17,5))

    plt.subplot(1,3,1) 
    cs = plt.contourf(X, Y, z1, cmap = cm)
    plt.colorbar(cs)
    plt.title(title1)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    plt.subplot(1,3,2) 
    cs = plt.contourf(X, Y, z2, cmap = cm)
    plt.colorbar(cs)
    plt.title(title2)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    plt.subplot(1,3,3) 
    cs = plt.contourf(X, Y, z1-z2, cmap = cm)
    plt.colorbar(cs)
    plt.title(title1+" - "+title2)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    # show figure
    if show == True:
        plt.show()
        
    # save figure
    if save == True:
        file_out = "PI_ctrl_comp_cont_"+year+"_"+m+"_"+var+".png"
        fig.savefig(file_out)
        

# givent two datasets compares the values of three variables
# INPUT:
# lon, lat = longitude and latitude
# data1, data2 = data to be input (requires two dimensions)
# month = month you are looking at
# apply_mask, mask = do you need a land mask? expects a binary file, if not set assumes False and NaN
# title1, title2 = titles ove the data you are looking at
# year = year of the contour plot
# var = variable looked at
# cm = set color scheme
# show = shows the contour plot, if unset assumes False
# save = saves contour plot as a png, assumes True if unset
def compare_contour_plots_TSM (lon, lat, temp1, temp2, salt1, salt2, melt1, melt2, month, ens1, ens2, apply_mask = False, mask = None, year = None, save = False, show = True):

    x = lon
    y = lat
    m = str(month)

    if apply_mask == True:
        temp1[mask == 0] = np.nan
        temp2[mask == 0] = np.nan
        salt1[mask == 0] = np.nan
        salt2[mask == 0] = np.nan
        melt1[mask == 0] = np.nan
        melt2[mask == 0] = np.nan
    
    print(np.shape(melt1))
    # create a 2D grid
    [X, Y] = np.meshgrid(lon, lat)

    fig =  plt.figure(figsize=(18,12))

    plt.subplot(3,3,1) 
    cs = plt.contourf(X, Y, temp1, levels=np.linspace(np.nanmin(temp1),np.nanmax(temp1),10), cmap = "coolwarm")
    plt.colorbar(cs)
    plt.title(ens1+" Theta")
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    plt.subplot(3,3,2) 
    cs = plt.contourf(X, Y, temp2, levels=np.linspace(np.nanmin(temp1),np.nanmax(temp1),10), cmap = "coolwarm")
    plt.colorbar(cs)
    plt.title(ens2+" Theta")
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    plt.subplot(3,3,3) 
    cs = plt.contourf(X, Y, temp1-temp2, cmap = "coolwarm")
    plt.colorbar(cs)
    plt.title("Anomaly")
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    plt.subplot(3,3,4) 
    cs = plt.contourf(X, Y, salt1, levels=np.linspace(np.nanmin(salt1),np.nanmax(salt1),10), cmap = "PRGn_r")
    plt.colorbar(cs)
    plt.title(ens1+" Salinity")
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    plt.subplot(3,3,5) 
    cs = plt.contourf(X, Y, salt2, levels=np.linspace(np.nanmin(salt1),np.nanmax(salt1),10), cmap = "PRGn_r")
    plt.colorbar(cs)
    plt.title(ens2+" Salinity")
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    plt.subplot(3,3,6) 
    cs = plt.contourf(X, Y, salt1-salt2, cmap = "PRGn_r")
    plt.colorbar(cs)
    plt.title("Anomaly")
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    plt.subplot(3,3,7) 
    cs = plt.contourf(X, Y, melt1, levels=np.linspace(0,0.002,10), cmap = "Blues_r")
    plt.colorbar(cs)
    plt.title(ens1+" Melt")
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    plt.subplot(3,3,8) 
    cs = plt.contourf(X, Y, melt2, levels=np.linspace(0,0.002,10), cmap = "Blues_r")
    plt.colorbar(cs)
    plt.title(ens2+" Melt")
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    plt.subplot(3,3,9) 
    cs = plt.contourf(X, Y, melt1-melt2, cmap = "Blues_r")
    plt.colorbar(cs)
    plt.title("Anomaly")
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
   
    fig.tight_layout(pad=2.0)

    # show figure
    if show == True:
        plt.show()
        
    # save figure
    if save == True:
        file_out = "PI_ctrl_comp_cont_"+year+"_"+m+"_TSM.png"
        fig.savefig(file_out)
######################### TIMESERIES #########################

# creates timeseries at point:
# INPUT:
# time = time variable, exoects monthly resolution
# data = data to be input (requires two dimensions), expects monthly data
# title = graph title
# year = year of the contour plot
# var = variable looked at
def make_timeseries_at_point (time, data, title, var, units, year):

    fig =  plt.figure(figsize=(12,10))

    x = time
    y = data

    labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

    plt.plot(x, y, color='green', linewidth = 2)
    plt.xticks(x, labels,rotation=45)
    plt.title(title)

    fig.savefig("temp"+year+".png")

# plot timeseries over the area:
# INPUT:
# time = time
# data = data to be input (requires two dimensions)
# title = graph title
# year = year of the contour plot
# var = variable looked at
# units = unit of the variable
# apply_mask,mask = do you need a land mask? expects a binary file, if not set assumes False and NaN
def make_timeseries_over_area (time, data, title, year, var, units, apply_mask=False, mask=np.nan):
    timeseries = []

    if apply_mask == True:
        for t in range(len(time)):
            temporary = data[t,:,:]
            temporary[mask == 0] = np.nan
            data[t,:,:] = temporary
    
    for t in range(data.shape[0]):
        timeseries.append(np.nanmean(data[t,:]))
    
    x = time
    y = np.array(timeseries) 

    fig =  plt.figure(figsize=(12,10))

    labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

    plt.plot(x, y, color='green', linewidth = 2)
    #plt.xticks(x, labels,rotation=45)
    plt.title(title)
    plt.xlabel('Days in the year')
    plt.ylabel(var+" "+units)
    
    fig.savefig("timeseries_"+var+"_"+year+".png")

# plots timeresies over multitude of years
# INPUT:
# filepath = initial path to file
# filename = netcdf file to read
# start_year = year to start from (expects string)
# n_years = how many years to run through
# var = what variable do you want to look at
# OUTPUT:
# y = array with the timeseries data
# ttimeseries_years = the years the timeseries runs through
# show = shows the contour plot, if unset assumes False
# save = saves contour plot as a png, assumes True if unset
def make_interannual_timeseries (filepath, filename, start_year, n_years, var, show = True, save = False):
    # set up the data
    input_data = filepath+str(start_year)+"01/MITgcm/"
    id = nc.Dataset(input_data+filename, 'r')
    time = id.variables["time"][:]
    lat = id.variables["YC"][:]
    lon = id.variables["XC"][:]
    mask = id.variables["maskC"][1,:,:]
    
    # run through the years
    timeseries = []
    timeseries_years = []
    for i in range(n_years):
        fileyear=str(start_year+i)
        timeseries_years.append(fileyear)
        input_data = filepath+fileyear+"01/MITgcm/"
        id = nc.Dataset(input_data+filename, 'r')
        data = id.variables[var][:,1,:,:]
        for t in range(len(time)):
            temporary = data[t,:,:]
            temporary[mask == 0] = np.nan
            data[t,:,:] = temporary
    
        for t in range(data.shape[0]):
            timeseries.append(np.nanmean(data[t,:]))
    
    #plot the data
    x = range(len(timeseries))
    y = np.array(timeseries) 

    fig =  plt.figure(figsize=(10,5))

    plt.plot(x, y, color='red', linewidth = 2)
    plt.xticks(np.arange(min(x),max(x), step=12), timeseries_years, rotation=45)
    plt.ylabel(var,fontsize=20)

    # show figure
    if show == True:
        plt.show()
        
    # save figure
    if save == True:
        file_out = "PI_ctrl_time_"+str(n_years)+"_years_"+var+".png"
        fig.savefig(file_out)

    return y, timeseries_years

######################## QUIVER PLOTS ##########################
#Plot vector fields on graph:
# INPUT:
# u, v = velocity vectors to be plotted
# lon, lat = associated longitude an latitude coordinates
# title = title of the plot
# exp = experiment name 
# var = variable looked at
# year = year of the contour plot
# overlay = do you want your quiver to be over another contour plot, set to False by default
# data = data to overlay
# x, y = coordinates associated with the data to overlay
# show = shows the contour plot, if unset assumes False
# save = saves contour plot as a png, assumes True if unset
def quiver_plot (u, v, lon, lat, title, exp, var, year, overlay = False, data = None, x = None, y = None, show = False, save = True):
    [LON,LAT] = np.meshgrid(lon, lat)

    fig, ax = plt.subplots()
    
    if overlay == True:
        [X,Y] = np.meshgrid(x, y)
        cp = plt.contourf(X, Y, data, cmap="coolwarm")
        cb = plt.colorbar(cp)
        quiv = plt.quiver(LON, LAT, u, v, color = "white")
        plt.title(title)
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
    else:
        plt.quiver(LON, LAT, u, v)
        plt.title(title)
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        
    if show == True:
        plt.show()
    
    if save == True:
        fig.save(exp+"_quiv_"+var+"_"+year+".pgn")

###################### OTHER FUNCTIONS ######################
# Fins the nearest value to a point
# INPUT:
# array = array you want to look through
# value = value you want to find
# OUTPUT:
# idx = returns the index with the closest value to the one chosen
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx




def compare_timeseries_at_point (time, data1, data2, data3, data4, title1, title2, year):

    fig =  plt.figure(figsize=(12,10))

    x = time
    y = data1
    z = data2
    i = data3
    j = data4
    labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

    plt.subplot(2,2,1)
    plt.plot(x, y, color='green', linewidth = 2)
    plt.plot(x, z, color='red', linewidth = 2)
    plt.xticks(x, labels,rotation=45)
    plt.title(title1)

    plt.subplot(2,2,3)
    plt.plot(x, y-z, color='black', linewidth = 2)
    plt.xticks(x, labels,rotation=45)
    plt.title("Difference between the two outputs")

    plt.subplot(2,2,2)
    plt.plot(x, i, color='blue', linewidth = 2)
    plt.plot(x, j, color='pink', linewidth = 2)
    plt.xticks(x, labels,rotation=45)
    plt.title(title1)

    plt.subplot(2,2,4)
    plt.plot(x, i-j, color='black', linewidth = 2)
    plt.xticks(x, labels,rotation=45)
    plt.title(title2)
   
    fig.tight_layout(pad=2.0)

    fig.savefig("temp"+year+".png")


if __name__ == "__main__":
    year = str(sys.argv[1])
    option = "binary"
    single_plot = False
    exp = "PI"
    if option == "NCfile":
        filepath = "/data/oceans_input/raw_input_data/CESM/LENS/daily/TREFHT/"
        filename = "b.e11.B1850C5CN.f09_g16.005.cam.h1.TREFHT.04020101-04991231.nc"
        var = "TREFHT"
        id = nc.Dataset(filepath+filename, 'r')
        temp = id.variables[var][0:365,:,:]
        lat = id.variables["lat"][:]
        lon = id.variables["lon"][:]
                                
    if option == "binary":
        grid_sizes = [192, 288]
        dimensions = ('t','y','x')
        filepath = "/data/oceans_input/processed_input_data/CESM/PIctrl/"
        filename = filepath+"PIctrl_ens03_PRECT_"+str(year)
        data = read_binary(filename, grid_sizes, dimensions)
        time = range(0, 365)
        print(np.shape(data))
   
    
    if option == "mask":
        id = nc.Dataset("/data/oceans_output/shelf/kaight/archer2_mitgcm/PAS_LENS001_O/output/189001/MITgcm/output.nc", 'r')
        mask = id.variables["maskC"][0,1:384:2,1:576:2]
        theta = [[np.nan for i in range(cols)] for j in range(rows)]
        stop_lat = find_nearest(lat, -62.389)
        start_lat = find_nearest(lat, -75.638)
        lat = lat[start_lat:stop_lat]
        stop_lon = find_nearest(lon, 279.94)
        start_lon = find_nearest(lon, 220.05)
        lon = lon[start_lon:stop_lon]
        temp = temp[0,start_lat:stop_lat,start_lon:stop_lon]
    
    #for t in time:
    lon = range(0, 288)
    lat = range(0, 192)
    data = np.transpose(data[0,:,:])
    #print(lat, lon)
    contourplots(lon,lat,data,"temp"+year, year, exp)
    #make_timeseries_at_point(time, thetaold, thetanew, saltold[:,380,590], saltnew[:,380,590], "Potential surface temperature at 62°S 100°W in "+year, "Salinity at 62°S 100°W in "+year, year)
    #theta_over_area=make_timeseries_over_area(time,theta,"temperature over the area "+year, year)
    


