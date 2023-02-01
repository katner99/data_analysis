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
# compare_contour_plots_TSM: compares multiple plots for temperature, salinity, and melt
#
# TIMESERIES:
# make_timeseries_at_point: timesereis at a single point
# make timeseries over area: calculates the average over the area and pots as a timeseries
# make_interannual_timeseries: timeseries over multiple years
# compare_timeseries_TSM: compare interannual timeseries for temperature, salinity, and melt
#
# OTHER PLOTS:
# quiver_plot: maps vector fields with or without overlay
#
# OTHER FUNCTIONS:
# find_nearest: finds the index with the value closest to the one wanted

################## CONTOUR PLOTS ################################

# CONTOUR PLOT:
# INPUT:
# lon, lat = longitude and latitude
# data     = data to be input (expects two dimensions)
# var      = variable name
# title    = graph title, if unset assumes None
# year     = year represented, if unset assumes None - NB needed to save
# cm       = colorscheme used, if unset assumes "ocean"
# low_val, high_val = max and min values for the colorbar, if unset assumes None
# apply_land_mask, land_mask = do you need a land mask? If unset assumes False and sets the mask to None
# apply_ice_mask, ice_mask   = do you need an ice mask? If unset assumes False and sets the mask to None
# show     = shows the contour plot, if unset assumes False
# save     = saves contour plot as a png, if unset assumes True

def contour_plots (lon, lat, data, var, title = None, year = None, cm = "ocean", low_val = None, high_val = None, apply_land_mask=False, apply_ice_mask=False, land_mask=None, ice_mask=None, show=False, save=True):

    x = lon
    y = lat
    z = data
    
    # prepare values so all months have the same parameters
    if low_val == None:
        low_val = np.nanmin(data)
    if high val == None:
        high_val = np.nanmax(data)
    step = (high_val-low_val)/15
    
    # apply mask over land
    if apply_land_mask == True:
        data[land_mask == 0] = np.nan
    
    # apply mask over ice
    if apply_ice_mask == True:
        data[ice_mask < 1] = np.nan
    
    # create a 2D grid
    [X, Y] = np.meshgrid(lon, lat)

    # set up figure
    fig =  plt.figure(figsize=(15,10))
    cs = plt.contourf(X, Y, z, np.arange(low_val, high_val, step), cmap = cm)
    plt.colorbar(cs)
    plt.title(title)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    
    # save figure
    if save == True:
        fig.savefig("PI_ctrl_cont_"+year+"_"+var+".png")

    # show figure
    if show == True:
        plt.show()
    
    return cs
    

# ANIMATE CONTOUR:
# INPUT:
# lon, lat = longitude and latitude
# data     = data to be input (expects monthly data)
# year     = represented
# var      = variable name
# cm       = set color scheme, if unset assumes "ocean"
# low_val, high_val = max and min values for the colorbar, if unset assumes None
# apply_land_mask, land_mask = do you need a land mask? If unset assumes False and sets the mask to None
# apply_ice_mask, ice_mask   = do you need an ice mask? If unset assumes False and sets the mask to None
# show = shows the contour plot, if unset assumes False
# save = saves contour plot as a png, if unset assumes True

def animate_contour (lon, lat, data, year, var, cm = "ocean", low_val = None, high_val = None, apply_land_mask=False, land_mask=None, show=False, save=True):
    # prepare grid
    [X, Y] = np.meshgrid(lon, lat)
    
    # prepare values so all months have the same parameters
    if low_val == None:
        low_val = np.nanmin(data)
    if high val == None:
        high_val = np.nanmax(data)
    step = (high_val-low_val)/15
    
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
              
        # create frame
        cont = contour_plots(lon, lat, data, var, cm = cm, low_val = low_val, high_val = high_val, apply_land_mask=apply_land_mask, land_mask=land_mask, show=False, save=False)
        plt.title(var+" "+labels[i]+" "+year)
        return cont  

    plt.colorbar(animate(0))
    
    # animate
    anim = animation.FuncAnimation(fig, animate, frames=12, interval = 200)
    
    if save == True:
        anim.save("PI_ctrl_cont_"+var+"_"+year+".gif", fps = 2)

    # show figure
    if show == True:
        plt.show()

# COMPARE CONTOUR PLOT:
# INPUT:
# lon, lat     = longitude and latitude
# data1, data2 = data to be input (expects two dimensions)
# month        = month represented
# apply_mask, mask = do you need a mask? If unset assumes False and sets the mask to None
# title1, title2   = graph titles, if unset assumes None
# year         = year represented, if unset assumes None
# var          = variable name
# cm           = set color scheme, if unset assumes None
# show         = shows the contour plot, if unset assumes False
# save         = saves contour plot as a png, if unset assumes True
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
    cs = plt.contourf(X, Y, z1, levels=np.linspace(np.nanmin(z1),np.nanmax(z1),10), cmap = cm)
    plt.colorbar(cs)
    plt.title(title1)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    plt.subplot(1,3,2) 
    cs = plt.contourf(X, Y, z2, levels=np.linspace(np.nanmin(z1),np.nanmax(z1),10), cmap = cm)
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
        

# COMPARE PLOTS FOR TEMPERATURE, SALINITY, AND MELT
# INPUT:
# lon, lat     = longitude and latitude
# temp1, temp2 = data for temperature
# salt1, salt2 = data for salinity
# melt1, melt2 = data for melt
# month        = represented
# apply_mask, mask = do you need a mask? If unset assumes False and sets the mask to None
# year         = year represented, if unset assumes None
# show         = shows the contour plot, if unset assumes False
# save         = saves contour plot as a png, if unset assumes True
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

    # create subplot with each variable on a new line
    # TEMPERATURE
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

    # SALINITY
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

    # MELT
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

    # save figure
    if save == True:
        file_out = "PI_ctrl_comp_cont_"+year+"_"+m+"_TSM.png"
        fig.savefig(file_out)

    # show figure
    if show == True:
        plt.show()
        
   
######################### TIMESERIES #########################

# TIMESERIES AT POINT
# INPUT:
# time  = time variable, exoects monthly resolution
# data  = data to be input (expects monthly data)
# title = graph title
# year  = year represented
# var   = variable name
# show         = shows the contour plot, if unset assumes False
# save         = saves contour plot as a png, if unset assumes True
def make_timeseries_at_point (time, data, title, var, units, year, save = True, show = False):

    fig =  plt.figure(figsize=(12,10))

    x = time
    y = data

    labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

    plt.plot(x, y, color='green', linewidth = 2)
    plt.xticks(x, labels,rotation=45)
    plt.title(title)

    # save figure
    if save == True:
        file_out = "PI_ctrl_time_"+year+"_"+var+".png"
        fig.savefig(file_out)

    # show figure
    if show == True:
        plt.show()
        
    return fig

# TIMESERIES OVER AREA:
# INPUT:
# time  = time
# data  = data to be input (requires two dimensions)
# title = graph title
# year  = year of the contour plot
# var   = variable looked at
# units = unit of the variable
# apply_mask, Mask = do you need a land mask? If unset assumes False and sets mask to None
# show  = shows the contour plot, if unset assumes False
# save  = saves contour plot as a png, if unset assumes True
def make_timeseries_over_area (time, data, title, year, var, units, apply_mask=False, mask=None, save = True, show = False):
    
    timeseries = []

    if apply_mask == True:
        for t in range(len(time)):
            temporary = data[t,:,:]
            temporary[mask == 0] = np.nan
            data[t,:,:] = temporary
    
    for t in range(data.shape[0]):
        timeseries.append(np.nanmean(data[t,:]))

    fig = make_timeseries_at_point(time, np.array(timeseries), title, var, units, year)

    # save figure
    if save == True:
        file_out = "PI_ctrl_time_"+year+"_"+var+".png"
        fig.savefig(file_out)

    # show figure
    if show == True:
        plt.show()

# INTERDECADAL TIMESERIES:
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
def make_interannual_timeseries (filepath, filename, start_year, n_years, var, plot = False, show = True, save = False):
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

        # read based on the variable read in
        if var == "THETA":
            data_d = id.variables[var][:,11:21,:,:]
            data = np.mean(data_d, axis=1)
        elif var == "SALT":
            data = id.variables[var][:,1,:,:]
            data[data == 0] = np.nan
        elif var == "SIfwmelt":
            data = id.variables[var][:,:,:]

        # apply mask
        for t in range(len(time)):
            temporary = data[t,:,:]
            temporary[mask == 0] = np.nan
            data[t,:,:] = temporary

        # create timeseries
        for t in range(data.shape[0]):
            timeseries.append(np.nanmean(data[t,:]))
            
    #plot the data
    x = range(len(timeseries))
    y = np.array(timeseries) 

    if plot == True:
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

# COMPARE TIMESERIES:
# INPUT:
# time         = time variable
# temp1, temp2 = data for temperature
# salt1, salt2 = data for salinity
# melt1, melt2 = data for melt
# ens1, ens2   = ensemble members
# show = shows the contour plot, if unset assumes False
# save = saves contour plot as a png, assumes True if unset
def compare_timeseries_TSM(time, temp1, temp2, salt1, salt2, melt1, melt2, ens1, ens2, show = False, save = True):

    x = range(len(time)*12)

    fig =  plt.figure(figsize=(20,15))

    labels = time[1:-1:2]

    # TEMPERATURE
    plt.subplot(2,3,1)
    plt.plot(x, temp1, color='green', linewidth = 2, label=ens1)
    plt.plot(x, temp2, color='red', linewidth = 2, label=ens2)
    plt.legend(loc="upper left")
    plt.xticks(np.arange(min(x),max(x), step=24), labels,rotation=45)
    plt.title("Theta")

    # SALINITY
    plt.subplot(2,3,2)
    plt.plot(x, salt1, color='green', linewidth = 2, label=ens1)
    plt.plot(x, salt2, color='red', linewidth = 2, label=ens2)
    plt.legend(loc="upper left")
    plt.xticks(np.arange(min(x),max(x), step=24), labels,rotation=45)
    plt.title("Salinity")

    # MELT
    plt.subplot(2,3,3)
    plt.plot(x, melt1, color='green', linewidth = 2, label=ens1)
    plt.plot(x, melt2, color='red', linewidth = 2, label=ens2)
    plt.legend(loc="upper left")
    plt.xticks(np.arange(min(x),max(x), step=24), labels,rotation=45)
    plt.title("Melt")

    # ANOMALIES
    plt.subplot(2,3,4)
    plt.plot(x, temp1-temp2, color='black', linewidth = 2)
    plt.xticks(np.arange(min(x),max(x), step=24), labels,rotation=45)
    plt.title("Temperature Anomaly")

    plt.subplot(2,3,5)
    plt.plot(x, salt1-salt2, color='black', linewidth = 2)
    plt.xticks(np.arange(min(x),max(x), step=24), labels,rotation=45)
    plt.title("Salinity Anomaly")

    plt.subplot(2,3,6)
    plt.plot(x, melt1-melt2, color='black', linewidth = 2)
    plt.xticks(np.arange(min(x),max(x), step=24), labels,rotation=45)
    plt.title("Melt Anomaly")
   

    # save figure
    if save == True:
        fig.savefig("PI_ctrl_comp_time.png")

    # show figure
    if show == True:
        plt.show()
        
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


if __name__ == "__main__":
    print("ERROR!! This is a file containing functions and cannot be run independently")



