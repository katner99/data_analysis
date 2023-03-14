''' script to animate a set of contour plots '''
import datetime
import sys

import matplotlib
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
from PIL import Image, ImageFilter
from matplotlib import animation
from mpl_toolkits.axes_grid1 import AxesGrid
from plots import animate_contour, create_mask

def make_animation(fig, ax, params, plot_date, plot_name):
    """Create an animated contour plot of a 2D dataset.

    Args:
        fig (matplotlib.figure.Figure): The figure object to use for the plot.
        ax (matplotlib.axes.Axes): The axes object to use for the plot.
        params (dict): A dictionary containing the following parameters:
            - lon_plot (array): 1D array of longitudes.
            - lat_plot (array): 1D array of latitudes.
            - time_plot (array): 1D array of datetime objects.
            - data_plot (array): 3D array of data to animate.
            - depth_plot (array): 2D array of depth.
            - ocean_plto (bool): boolean for wet grid cells.
            - color_scheme_plot (str): Matplotlib colormap for the data.
            - speedy (int, optional): Number of frames per second in the animation.
        plot_date (Date): series of dates represented in the figure
        plot_name (str): what to save the animation as

    Returns:
        None

    Raises:
        ValueError: If the dimensions of the input arrays do not match.
    """
    # Check input array dimensions
    if len(params['lon_plot']) != params['data_plot'].shape[2]:
        raise ValueError("Length of lon_plot does not match data_plot dimension")
    if len(params['lat_plot']) != params['data_plot'].shape[1]:
        raise ValueError("Length of lat_plot does not match data_plot dimension")
    if len(params['time_plot']) != params['data_plot'].shape[0]:
        raise ValueError("Length of time_plot does not match data_plot dimension")

    # Create meshgrid
    X, Y = np.meshgrid(params['lon_plot'], params['lat_plot'])

    # Set up the plot
    ax.set_xlim(min(params['lon_plot']), max(params['lon_plot']))
    ax.set_ylim(min(params['lat_plot']), max(params['lat_plot']))
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

    # set colorbar limits
    low_lim = np.nanmin(params['data_plot'])
    high_lim = np.nanmax(params['data_plot'])
    
    # calculate the mask
    [land_mask, shelf_mask, colors] = create_mask(params['depth_plot'], params['ocean_plot'])
    
    # Animation function
    def animate(i):
        z = params['data_plot'][i, :, :]


        cont = ax.contourf(
            X, Y, z, 
            levels=np.linspace(low_lim, high_lim, 15),
            extend="both", 
            cmap=params['color_scheme_plot'])
            
        ax.contourf(
            X, Y, land_mask, 
            cmap=matplotlib.colors.ListedColormap(params['colors_plot'])) # mask out the land
        ax.contour(X, Y, shelf_mask, 
        2, 
        cmap = "Greys",
        linestyles='dashed')
        
        plt.title(plot_date[i].strftime("%m-%Y"))
            
        return cont  
    
    # plot the colorbar once
    plt.colorbar(animate(0), ticks = [np.linspace(0, 20, 10)])
    
    
    # animate
    anim = animation.FuncAnimation(fig, animate, frames=len(params['time_plot'])
        
    # save and display
    anim.save(plot_name, fps = speedy)
    plt.show()
    
def main():
    
    """Main function to run the program"""
    
    # check that enough arguments have been input by the user
    if len(sys.argv) != 4:
        sys.exit("Stopped - Incorrect number of arguements. Use python PI_animation.py <year or slice> <var> <exp>")
    
    # set the year, variable name, and experiment name    
    year = str(sys.argv[1])
    var = str(sys.argv[2])
    exp = str(sys.argv[3])
    
    # run over a number of years (a slice)
    if year == "slice":
        # set the filepath and colour scheme
        filepath = "/data/oceans_output/shelf/katner33/PIctrl_output/ensemble_mean/WIND_ensemble_mean_slice.nc"
    
    else:
         # set up the filepath
        filepath = "/data/oceans_output/shelf/katner33/PIctrl_output/PAS_ctrl08/output/"+year+"01/MITgcm/output.nc"
        
    # load up the different variables
    id = nc.Dataset(filepath, 'r')
    
    # set the colour scheme and load up the variables
    if var == "THETA":
        data = id.variables[var][:,11:21,:,:]
        color_scheme = "coolwarm"
    if var == "SALT":
        data = id.variables[var][:,1,:,:]
        data[data == 0] = np.nan
        print(np.nanmax(data))
        color_scheme = "PRGn_r"
    if var == "SHIfwFlx" or var == "SIheff":
        data = id.variables[var][:,:,:]
        color_scheme = "YlGnBu"
    
    # load up the parameters
    time = id.variables["time"][:]
    lat = id.variables["YC"][:]
    lon = id.variables["XC"][:]
    depth = id.variables["Depth"][:, :]
    ocean = id.variables["maskC"][1, :, :]
    
    # set mask
    #

    params = {
        lon_plot : lon
        lat_plot : lat
        time_plot : time
        data_plot : data
        depth_plot : depth
        ocean_plot : ocean
        color_scheme_plot : color_scheme
        speeedy : 5
    }
    
    # create list of dates shown in function
    start = datetime.datetime.strptime("01-2070", "%m-%Y")
    end = datetime.datetime.strptime("12-2101", "%m-%Y")
    date_generated = [start + datetime.timedelta(days=x) for x in range(0, (end - start).days, 30)]
    
    # write out the file name
    name = exp_plot +"_cont_"+var+"_.gif"
    
    # pass through the figure and axes you want
    fig, ax = plt.subplot(figsize=(8,6))
    
    # run the animation
    make_animation(fig, ax, params, date_generated, name)

    #animate_contour(lon, lat, data, year, exp, var, cs, True, False, land_mask, ice_mask, False, True)
        
    
if __name__ == '__main__':
    main() # run the program
