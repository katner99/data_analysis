import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np
import sys
import datetime
import netCDF4 as nc


if __name__ == "__main__":
    year = str(sys.argv[1])
    filename = "output.nc"
    var = "SALT"
    id = nc.Dataset(filename, 'r')
    temp = id.variables[var][:,1,:,:]
    time = id.variables["time"][:]
    lat = id.variables["YC"][:]
    lon = id.variables["XC"][:]
    mask = id.variables["maskC"][1,:,:]
    labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    
    [X, Y] = np.meshgrid(lon, lat)
    
    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure(figsize=(10,8))
    ax = plt.axes(xlim=(min(lon), max(lon)), ylim=(min(lat), max(lat)))
    plt.xlabel('Longitude', fontsize=20)
    plt.ylabel('Latitude', fontsize=20)
    
    # animation function
    def animate(i): 
    	#cb.cla()
        z = temp[i,:,:]
        z[mask == 0] = np.nan
        cont = plt.contourf(X, Y, z, 25, cmap="plasma")
        #cb = plt.colorbar(cont, ax=ax)
        plt.title("Salinity at level 1 "+labels[i]+" 2099", fontsize=20)
        return cont  

    anim = animation.FuncAnimation(fig, animate, frames=12, interval=10)

    #plt.show()
    anim.save("salt_2099.gif", fps=1)
    #contourplots(lon, lat, temp, "Pre-industrial potential temperature for 2070", year, "PIctrl",True,mask)

    #make_timeseries_over_area(time, temp, "Timeseries of theta in 2070", year,"Theta","(°C)", True, mask)

    #frames = []
    #for i in range(1, len(time)):
    #    frams.append(contourplots(lon, lat, temp[i,:,:], "Pre-industrial potential temperature for 2070", year, "PIctrl",True,mask,False))

    #gif.save(frames, "theta_2070.gif", duration=15,unit="s",between="startend")
