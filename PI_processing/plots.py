###################################################
# script to produce controur plots to compare data
# between two datasets
# created by Katherine Turner based on hovmoller.py
# script 
# 23 November 2022
###################################################
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import animation
from PIL import Image,ImageFilter
import numpy as np
import sys
import datetime
import netCDF4 as nc
#from mitgcm_python.utils import mask_land_ice
#from mitgcm_python.calculus import over_area
#from mitgcm_python.file_io import read_binary

################## CONTOUR PLOTS ###############################

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

def contourplots (lon, lat, data, title, year, exp, var, apply_land_mask=False, apply_ice_mask=False, land_mask=np.nan, ice_mask=np.nan, show=False, save=True):

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
# apply_land_mask, land_mask = do you need a land mask? expects a binary file, if not set assumes False and NaN
# apply_ice_mask, ice_mask = do you need an ice mask? expects a binary file, if not set assumes False and NaN
# show = shows the contour plot, if unset assumes False
# save = saves contour plot as a png, assumes True if unset

def animate_contour (lon, lat, data, year, exp, var, apply_land_mask=False, apply_ice_mask=False, land_mask=np.nan, ice_mask=np.nan, show=False, save=True):
    # prepare grid
    [X, Y] = np.meshgrid(lon, lat)
    
    # prepare values so all months have the same parameters
    low_val = np.min(data)
    #low_val = -2
    high_val = np.max(data)+0.25
    #high_val = 1.5
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
        
        if apply_ice_mask == True:
            z[ice_mask[i,:,:] == 0] = np.nan
        # apply mask over land
        if apply_land_mask == True:
            z[land_mask == 0] = np.nan
              
        # create frame
        cont = plt.contourf(X, Y, z, np.arange(low_val, high_val,step), cmap="coolwarm")
        
        # write colorbar on the first or it will keep being copied to the figure
        
        plt.title(var+" "+labels[i]+" "+year)
        return cont  

    plt.colorbar(animate(0))
    
    # animate
    anim = animation.FuncAnimation(fig, animate, frames=12, interval=20)

    if show == True:
        plt.show()
    
    if save == True:
        anim.save(exp+"_cont_"+var+"_"+year+".gif")
        
    
    
def make_timeseries_at_point (time, data, title, var, units, year):

    fig =  plt.figure(figsize=(12,10))

    x = time
    y = data

    labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

    plt.plot(x, y, color='green', linewidth = 2)
    plt.xticks(x, labels,rotation=45)
    plt.title(title)

    fig.savefig("temp"+year+".png")


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

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def compare_contour_plots (lon, lat, data1, data2, title, var, year):

    x = lon
    y = lat
    z1 = data1
    z2 = data2

    # create a 2D grid
    [X, Y] = np.meshgrid(lon, lat)

    fig =  plt.figure(figsize=(8,5))

    # plot contour lines levels=np.linspace(-0.5,0.5,10)
    #plt.subplot(1,2,1) 
    cs = plt.contourf(X, Y, z1-z2, cmap = "BrBG")
    plt.colorbar(cs)
    #ax.set_ylim(-1100,0)
    plt.title(title1)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    #plt.subplot(1,2,2) 
    #cs = plt.contourf(X, Y, z2, cmap = "BrBG")
    #plt.colorbar(cs)
    #ax.set_ylim(-1100,0)
    #plt.title("1850")
    #plt.xlabel('Longitude')
    #plt.ylabel('Latitude')

    fig.savefig("compare_cont_"+var+"_"+year+".png")

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
    


