import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np
import sys
import datetime
import netCDF4 as nc
from plots import contourplots, make_timeseries_over_area
import numpy as np


if __name__ == "__main__":
    #year = str(sys.argv[1])
    option="thirty"
    filepath = "/data/oceans_output/shelf/katner33/PIctrl_output/PAS_ctrl/output/207001/MITgcm/"
    filename = "output.nc"
    id = nc.Dataset(filepath+filename, 'r')
    time = id.variables["time"][:]
    lat = id.variables["YC"][:]
    lon = id.variables["XC"][:]
    mask = id.variables["maskC"][1,:,:]
    labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    
    if option=="animate":
        filepath = "/data/oceans_output/shelf/katner33/PIctrl_output/PAS_ctrl/output/209901/MITgcm/"
        filename = "output.nc"
        var = "SALT"
        id = nc.Dataset(filepath+filename, 'r')
        temp = id.variables[var][:,1,:,:]
        [X, Y] = np.meshgrid(lon, lat)
    
        # First set up the figure, the axis, and the plot element we want to animate
        fig = plt.figure(figsize=(10,8))
        ax = plt.axes(xlim=(min(lon), max(lon)), ylim=(min(lat), max(lat)))
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
    
        # animation function
        def animate(i): 
            z = temp[i,:,:]
            z[mask == 0] = np.nan
            cont = plt.contourf(X, Y, z, 25)
            plt.title(labels[i]+" 2099")
            return cont  

        anim = animation.FuncAnimation(fig, animate, frames=12, interval=20)

        plt.show()
    
    year=2070
    years=[]
    timeseries = []
    if option=="thirty":
        for i in range(0,29):
            fileyear=str(year+i)
            years.append(fileyear)
            filepath = "/data/oceans_output/shelf/katner33/PIctrl_output/PAS_ctrl/output/"+fileyear+"01/MITgcm/"
            filename = "output.nc"
            var = "THETA"
            id = nc.Dataset(filepath+filename, 'r')
            data = id.variables[var][:,1,:,:]
            for t in range(len(time)):
                temporary = data[t,:,:]
                temporary[mask == 0] = np.nan
                data[t,:,:] = temporary
    
            for t in range(data.shape[0]):
                timeseries.append(np.nanmean(data[t,:]))

        x = range(len(timeseries))
        y = np.array(timeseries) 

        fig =  plt.figure(figsize=(10,5))

        #labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

        plt.plot(x, y, color='red', linewidth = 2)
        plt.xticks(np.arange(min(x),max(x), step=12), years, rotation=45)
        plt.ylabel("Temperature (°C)",fontsize=20)
        plt.show()
        
        #plt.xticks(x, labels,rotation=45)
        #plt.title(title)
        #plt.xlabel('Days in the year')
        #plt.ylabel(var+" "+units)
    
        #fig.savefig("timeseries_"+var+"_"+year+".png")
    
        #print(temp[1:10])
        
