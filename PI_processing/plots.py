import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib import animation

from PIL import Image,ImageFilter
import numpy as np
import sys
import datetime
import time
import netCDF4 as nc
import xarray as xr
from mitgcm_python.utils import mask_land_ice
from mitgcm_python.calculus import over_area
from mitgcm_python.file_io import read_binary
from funcs import interpolate_currents, moving_average

def plot_timeseries_comparison(ax, var, data, experiments, ensemble_members, plot_info):
    """
    This function generates a plot comparing different experiments based on the provided data and plot information.

    Parameters:
    - var: Variable information (not used directly in the function)
    - data: Dictionary containing the data for different experiments
    - experiments: List of experiment names to be compared
    - ensemble_members: List of ensemble members
    - plot_info: Dictionary containing plot configuration parameters including:
      - ylabel: Label for the y-axis
      - xlabel: Label for the x-axis
      - time: Time information
      - file_out: Output file name for saving the plot
      - smooth: Smoothing factor for the data (optional, default is 0)
      - linearity: Flag indicating whether to plot linearity (optional, default is False)
      - warming: Flag indicating whether to adjust for warming (optional, default is False)
      - percentage: Flag indicating whether to calculate percentage (optional, default is False)
      - shade_range: Flag indicating whether to shade the range (optional, default is True)
      - title: Title for the plot (optional)
      - experiment_full: List containing full names of experiments (optional)

    Returns:
    - ax: Axis object representing the main plot
    - all_means: List containing the means of all experiments' data
    """
    colors = ["dodgerblue", "purple", "forestgreen", "orangered"]
    experiment_full = plot_info["experiment_full"]
    all_means = []

    if plot_info.get("warming", False):
        ctrl_mean = np.nanmean([data["CTRL"][ens][ts_idx][:] for ens in ensemble_members for ts_idx in range(len(data["CTRL"][ens]))], axis=0)
        warming_diff = -np.nanmean(ctrl_mean) / 10
    
    for exp_idx, exp in enumerate(experiments):
        experiment_mean = np.nanmean([data[exp][ens][ts_idx][:] for ens in ensemble_members for ts_idx in range(len(data[exp][ens]))], axis=0) / 10
        experiment_max = np.nanmax([data[exp][ens][ts_idx][:] for ens in ensemble_members for ts_idx in range(len(data[exp][ens]))], axis=0) / 10
        experiment_min = np.nanmin([data[exp][ens][ts_idx][:] for ens in ensemble_members for ts_idx in range(len(data[exp][ens]))], axis=0) / 10
        
        if plot_info.get("smooth", 0) > 0:
            smoothed_mean = moving_average(experiment_mean, plot_info["smooth"] * 12)
            if plot_info.get("shade_range", True):
                smoothed_max = moving_average(experiment_max, plot_info["smooth"] * 12)
                smoothed_min = moving_average(experiment_min, plot_info["smooth"] * 12)
        else:
            smoothed_mean, smoothed_max, smoothed_min = experiment_mean, experiment_max, experiment_min
        
        if plot_info.get("warming", False):
            smoothed_mean += warming_diff
            smoothed_max += warming_diff
            smoothed_min += warming_diff
        
        ax.plot(smoothed_mean, color=colors[exp_idx], label=experiment_full[exp_idx], alpha=1)
        
        if plot_info.get("shade_range", True):
            ax.fill_between(range(len(smoothed_min)), smoothed_min, smoothed_max, color=colors[exp_idx], alpha=0.3)
        
        all_means.append(experiment_mean[:])

    if plot_info.get("linearity", False) or plot_info.get("percentage", False):
        linearity = moving_average(all_means[3] + all_means[2] - all_means[0] - all_means[1], 12 * plot_info.get("smooth", 0))
        ax.plot(linearity, color="black", label="linearity")
        ax.axhline(0, color='grey')

    ax.set_ylabel(plot_info['ylabel'], fontsize=14)
    ax.set_xticks(np.arange(6 * plot_info.get("smooth", 0), plot_info['time'] - (6 * plot_info.get("smooth", 0)), 120))
    ax.set_xticklabels(plot_info['xlabel'].astype(str), rotation=45)
    ax.set_xlim(plot_info['x_lim'])
    ax.set_title(plot_info.get("title"))

    if plot_info.get("warming", False):
        ax.set_ylim([-1.6, 1.6])
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend()
    ax.grid(alpha=0.8)

    return ax, all_means


def create_mask (depth, ice_mask):
    land_mask = np.zeros(np.shape(depth))
    land_mask[depth == 0] = 1

    # apply mask over the ice shelf (determiend by the ice-mask) and the continental shelf (roughly where the depth is less than 1500m)
    mask = np.zeros(np.shape(depth))
    mask[depth < 1500] = 1
    mask[ice_mask == 0] = 2
        
    # set the colors to block the continent (set to grey)
    colors = [
        (1.0, 1.0, 1.0, 0),
        (0.7, 0.7, 0.7, 1),
        (0.6, 0.6, 0.6, 1)]
    return land_mask, mask, colors

# CONTOUR PLOT: plots a simple contour plot

def contour_plots (x, y, data, year, var, depth, exp, color_scheme = "ocean", apply_mask = False, mask = None, lonlatplot = True, font_size = 12, low_val = None, high_val = None, xlabel = None, ylabel = False, title = None, save_as = None, show=False, save=True):
    
    # prepare values so all months have the same parameters
    if low_val == None:
        low_val = np.nanmin(data)
    if high_val == None:
        high_val = np.nanmax(data)    

    # set up title
    if title == None:
        title = var+" "+year
    
    # apply mask over land, i.e. where ocean depth is zero
    if mask is not None:
        land_mask = np.zeros(np.shape(depth))
        land_mask[depth == 0] = 1

        # apply mask over the ice shelf (determiend by the ice-mask) and the continental shelf (roughly where the depth is less than 1500m)
        mask = np.zeros(np.shape(depth))
        mask[depth < 1500] = 1
        mask[ice_mask == 0] = 2
    
        # set the colors to block the continent (set to grey)
        colors = [
            (1.0, 1.0, 1.0, 0),
            (0.7, 0.7, 0.7, 1),
            (0.6, 0.6, 0.6, 1)]
    
    # create a 2D grid
    [X, Y] = np.meshgrid(x, y)

    # set up figure
    fig = plt.figure(figsize=(15,10))
    cont = plt.contourf(X, Y, z, levels=np.linspace(low_val, high_val,15), extend = "both", cmap = color_scheme)
    if apply_mask == True:
        plt.contourf(X, Y, land_mask, cmap = matplotlib.colors.ListedColormap(colors))
        plt.contour(X, Y, mask, 2, cmap = "Greys",linestyles='dashed')
    plt.colorbar(cont)
    plt.title(title, fontsize = font_size)
    plt.xlabel(xlabel, fontsize = font_size)
    plt.ylabel(yalbel, fontsize = font_size)
    
    # save figure
    if save == True:
        if save_as is not None:
            fig.savefig(save_as)
        else:
            fig.savefig(exp+"_cont_"+year+"_"+var+".png")

    # show figure
    if show == True:
        plt.show()
    
    return cont
    

# ANIMATE CONTOUR: animate contour plot
def animate_contour (x, y, data, year, var, depth, exp, color_scheme = "ocean", mask = None, lonlatplot = True, font_size = 12, low_val = None, high_val = None, xlabel = None, ylabel = False, save_as = None, show=False, save=True):
    # prepare grid
    [X, Y] = np.meshgrid(x, y)
    
    # prepare values so all months have the same parameters
    if low_val == None:
        low_val = np.nanmin(data)
    if high_val == None:
        high_val = np.nanmax(data)
    step = (high_val-low_val)/15
    
    # prepares the title
    labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    
    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure(figsize=(8,6))
    ax = plt.axes(xlim=(min(x), max(y)), ylim=(min(y), max(y)))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    # animation function
    def animate(i): 
        z = data[i,:,:]  
              
        cont = contour_plots(x = x, y = y, data = z, year = year, var = var, depth = depth, exp = exp, color_scheme = color_scheme, mask = mask, lonlatplot = lonlatplot, font_size = font_size, low_val = low_val, high_val = high_val, xlabel = xlabel, ylabel = ylabel, title = None, save_as = None, show = False, save = False)

        plt.title(var+" "+labels[i]+" "+year)
        return cont  

    plt.colorbar(animate(0))
    
    # animate
    anim = animation.FuncAnimation(fig, animate, frames=12, interval = 200)
    
    if save == True:
        if save_as is not None:
            anim.save(save_as)
        else:
            anim.save(exp +"_cont_"+year+"_"+var+".gif", fps = 2)

    # show figure
    if show == True:
        plt.show()
        

# COMPARE PLOTS 3X3: creates subplots showign two datasets and their anomaly 
def compare_contour_plots_LATLON (exp, lon, lat, temp1, temp2, salt1, salt2, stress1, stress2, month, ens1, ens2, depth, ice_mask, font_size = 10, year = None, save = False, show = True):

    x = lon
    y = lat
    m = str(month)

    # apply mask over land, i.e. where ocean depth is zero
    land_mask = np.zeros(np.shape(depth))
    land_mask[depth == 0] = 1

    # apply mask over the ice shelf (determiend by the ice-mask) and the continental shelf (roughly where the depth is less than 1500m)
    mask = np.zeros(np.shape(depth))
    mask[depth < 1500] = 1
    mask[ice_mask == 0] = 2
    
    # set the colors to block the continent (set to grey)
    colors = [
        (1.0, 1.0, 1.0, 0),
        (0.7, 0.7, 0.7, 1),
        (0.6, 0.6, 0.6, 1)]
    
    [X, Y] = np.meshgrid(lon, lat)

    fig =  plt.figure(figsize=(18,12))

    # create subplot with each variable on a new line
    # TEMPERATURE
    plt.subplot(3,3,1) 
    cs = plt.contourf(X, Y, temp1, levels=np.linspace(np.nanmin(temp1),np.nanmax(temp2),10), cmap = "coolwarm")
    plt.contourf(X, Y, land_mask, cmap = matplotlib.colors.ListedColormap(colors))
    plt.contour(X, Y, mask, 2, cmap = "Greys",linestyles='dashed')
    plt.colorbar(cs)
    plt.title(ens1+" Theta", fontsize = font_size)
    plt.xlabel('Longitude', fontsize = font_size)
    plt.ylabel('Latitude', fontsize = font_size)

    plt.subplot(3,3,2) 
    cs = plt.contourf(X, Y, temp2, levels=np.linspace(np.nanmin(temp1),np.nanmax(temp2),10), cmap = "coolwarm")
    plt.contourf(X, Y, land_mask, cmap = matplotlib.colors.ListedColormap(colors))
    plt.contour(X, Y, mask, 2, cmap = "Greys",linestyles='dashed')
    plt.colorbar(cs)
    plt.title(ens2+" Theta", fontsize = font_size)
    plt.xlabel('Longitude', fontsize = font_size)
    plt.ylabel('Latitude', fontsize = font_size)

    plt.subplot(3,3,3) 
    cs = plt.contourf(X, Y, temp2-temp1, levels=np.linspace(-np.nanmax(np.abs(temp1-temp2)), np.nanmax(np.abs(temp1-temp2)), 10), cmap = "coolwarm")
    plt.contourf(X, Y, land_mask, cmap = matplotlib.colors.ListedColormap(colors))
    plt.contour(X, Y, mask, 2, cmap = "Greys",linestyles='dashed')
    plt.colorbar(cs)
    plt.title("Anomaly", fontsize = font_size)
    plt.xlabel('Longitude', fontsize = font_size)
    plt.ylabel('Latitude', fontsize = font_size)

    # SALINITY
    plt.subplot(3,3,4) 
    cs = plt.contourf(X, Y, salt1, levels=np.linspace(np.nanmin(salt2),np.nanmax(salt1),10), cmap = "PRGn_r")
    plt.contourf(X, Y, land_mask, cmap = matplotlib.colors.ListedColormap(colors))
    plt.contour(X, Y, mask, 2, cmap = "Greys",linestyles='dashed')
    plt.colorbar(cs)
    plt.title(ens1+" Salinity", fontsize = font_size)
    plt.xlabel('Longitude', fontsize = font_size)
    plt.ylabel('Latitude', fontsize = font_size)

    plt.subplot(3,3,5) 
    cs = plt.contourf(X, Y, salt2, levels=np.linspace(np.nanmin(salt2),np.nanmax(salt1),10), cmap = "PRGn_r")
    plt.contourf(X, Y, land_mask, cmap = matplotlib.colors.ListedColormap(colors))
    plt.contour(X, Y, mask, 2, cmap = "Greys",linestyles='dashed')
    plt.colorbar(cs)
    plt.title(ens2+" Salinity", fontsize = font_size)
    plt.xlabel('Longitude', fontsize = font_size)
    plt.ylabel('Latitude', fontsize = font_size)

    plt.subplot(3,3,6) 
    cs = plt.contourf(X, Y, salt2-salt1, levels=np.linspace(-np.nanmax(np.abs(salt1-salt2)), np.nanmax(np.abs(salt1-salt2)),10), cmap = "PRGn_r")
    plt.contourf(X, Y, land_mask, cmap = matplotlib.colors.ListedColormap(colors))
    plt.contour(X, Y, mask, 2, cmap = "Greys",linestyles='dashed')
    plt.colorbar(cs)
    plt.title("Anomaly", fontsize = font_size)
    plt.xlabel('Longitude', fontsize = font_size)
    plt.ylabel('Latitude', fontsize = font_size)

    # STRESS
    plt.subplot(3,3,7) 
    cs = plt.contourf(X, Y, stress1, levels=np.linspace(np.nanmin(stress2),np.nanmax(stress2),10), cmap = "BrBG")
    plt.contourf(X, Y, land_mask, cmap = matplotlib.colors.ListedColormap(colors))
    plt.contour(X, Y, mask, 2, cmap = "Greys",linestyles='dashed')
    plt.colorbar(cs)
    plt.title(ens1+" Stress", fontsize = font_size)
    plt.xlabel('Longitude', fontsize = font_size)
    plt.ylabel('Latitude', fontsize = font_size)

    plt.subplot(3,3,8) 
    cs = plt.contourf(X, Y, stress2, levels=np.linspace(np.nanmin(stress2),np.nanmax(stress2),10), cmap = "BrBG")
    plt.contourf(X, Y, land_mask, cmap = matplotlib.colors.ListedColormap(colors))
    plt.contour(X, Y, mask, 2, cmap = "Greys",linestyles='dashed')
    plt.colorbar(cs)
    plt.title(ens2+" Stress", fontsize = font_size)
    plt.xlabel('Longitude', fontsize = font_size)
    plt.ylabel('Latitude', fontsize = font_size)

    plt.subplot(3,3,9) 
    cs = plt.contourf(X, Y, stress2-stress1, levels=np.linspace(-np.nanmax(stress2 - stress1),np.nanmax(stress2 - stress1),10),cmap = "BrBG")
    plt.contourf(X, Y, land_mask, cmap = matplotlib.colors.ListedColormap(colors))
    plt.contour(X, Y, mask, 2, cmap = "Greys",linestyles='dashed')
    plt.colorbar(cs)
    plt.title("Anomaly", fontsize = font_size)
    plt.xlabel('Longitude', fontsize = font_size)
    plt.ylabel('Latitude', fontsize = font_size)
   
    fig.tight_layout(pad=2.5)

    # save figure
    if save == True:
        file_out = exp+"_comp_cont_"+year+"_TSM.png"
        fig.savefig(file_out)

    # show figure
    if show == True:
        plt.show()

######################## SLICES ##############################

# COMPARE_SLICE: compare slices of two variables (A and B) from two datasets (1 and 2) with the length (lon or lat) and the depth (z)
def compare_slice (length, z, varA_1, varA_2, varB_1, varB_2, varA_min, varA_max, varB_min, varB_max, varA_name, varB_name, loc_name, len_name, exp1_name, exp2_name, font_size = 12, len_limit = -800, save = 1, show = 1, save_as = None):

    # create grid from the length chosen and the depth
    [X, Y] = np.meshgrid(length, z)

    print(np.linspace(-np.nanmax(abs(varA_2 - varA_1)), np.nanmax(abs(varA_2 - varA_1))))

    # create plot
    fig = plt.figure(figsize=(15,10))
    
    # First line looks at var A and the anomaly
    plt.subplot(2,3,1)
    cs = plt.contourf(X, Y, varA_1, levels=np.linspace(varA_min, varA_max, 10), extend = "both", cmap = "RdBu_r")
    plt.ylim(len_limit, 0)
    plt.title(varA_name+" at "+loc_name+" "+exp1_name, fontsize = font_size)
    plt.xticks(rotation=45)
    plt.xlabel(len_name, fontsize = font_size)
    plt.ylabel('Depth', fontsize = font_size)
    plt.colorbar(cs)
    
    plt.subplot(2,3,2)
    cs = plt.contourf(X, Y, varA_2, levels=np.linspace(varA_min, varA_max, 10), extend = "both", cmap = "RdBu_r")
    plt.ylim(len_limit, 0)
    plt.title(varA_name+" at "+loc_name+" "+exp2_name, fontsize = font_size)
    plt.xticks(rotation=45)
    plt.xlabel(len_name, fontsize = font_size)
    plt.ylabel('Depth', fontsize = font_size)
    plt.colorbar(cs)
    
    #  levels=np.linspace(-np.nanmax(abs(varA_2 - varA_1)), np.nanmax(abs(varA_2 - varA_1)), 10)
    plt.subplot(2,3,3)
    cs = plt.contourf(X, Y, varA_2 - varA_1, cmap = "RdBu_r")
    plt.ylim(len_limit, 0)
    plt.title(varA_name+" at "+loc_name+" "+exp2_name+" - "+exp1_name, fontsize = font_size)
    plt.xticks(rotation=45)
    plt.xlabel(len_name, fontsize = font_size)
    plt.ylabel('Depth', fontsize = font_size)
    plt.colorbar(cs)
    
    # Second line looks at var B and the anomaly
    plt.subplot(2,3,4)
    cs = plt.contourf(X, Y, varB_1, levels=np.linspace(varB_min, varB_max, 10), extend = "both", cmap = "PRGn_r")
    plt.ylim(len_limit, 0)
    plt.title(varB_name+" at "+loc_name+" "+exp1_name, fontsize = font_size)
    plt.xticks(rotation=45)
    plt.xlabel(len_name, fontsize = font_size)
    plt.ylabel('Depth', fontsize = font_size)
    plt.colorbar(cs)
    
    plt.subplot(2,3,5)
    cs = plt.contourf(X, Y, varB_2, levels=np.linspace(varB_min, varB_max, 10), extend = "both", cmap = "PRGn_r")
    plt.ylim(len_limit, 0)
    plt.title(varB_name+" at "+loc_name+" "+exp2_name, fontsize = font_size)
    plt.xticks(rotation=45)
    plt.xlabel(len_name, fontsize = font_size)
    plt.ylabel('Depth', fontsize = font_size)
    plt.colorbar(cs)
    
    plt.subplot(2,3,6)
    cs = plt.contourf(X, Y, varB_2 - varB_1, cmap = "PRGn_r")
    plt.ylim(len_limit, 0)
    plt.title(varB_name+" at "+loc_name+" "+exp2_name+" - "+exp1_name, fontsize = font_size)
    plt.xticks(rotation=45)
    plt.xlabel(len_name, fontsize = font_size)
    plt.ylabel('Depth', fontsize = font_size)
    plt.colorbar(cs)
    
    fig.tight_layout(pad=2.5)
    
    # save figure
    if save == 1:
        if save_as is not None:
            fig.savefig(save_as)
        else:
            file_out = "THERMO_vs_LENS__slice.png"
            fig.savefig(file_out)

    # show figure
    if show == 1:
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
    if save == 1:
        file_out = "PI_ctrl_time_"+year+"_"+var+".png"
        fig.savefig(file_out)

    # show figure
    if show == 1:
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
def make_interannual_timeseries (filepath, filename, var, plot = False, show = True, save = False):
    # set up the data
    input_data = xr.open_dataset(filepath+str(start_year)+"01/MITgcm/")
    
    
    # run through the years
    timeseries = []
    timeseries_years = []
    for year in range(2070, 2100):
        fileyear=str(start_year+i)
        timeseries_years.append(fileyear)
        input_data = xr.open_dataset(filepath+str(start_year)+"01/MITgcm/")
        # read based on the variable read in
        if var == "THETA":
            data_d = id.variables[var][:,11:21,:,:]
            data = np.mean(data_d, axis=1)
        elif var == "SALT":
            data = id.variables[var][:,1,:,:]
            data[data == 0] = np.nan
        elif var == "SHIfwFlx":
            data = id.variables[var][:,:,:]

        # apply mask
        for t in range(len(time)):
            temporary = data[t,:,:]
            temporary[(mask>0) & (mask<1500)] = np.nan
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
def compare_timeseries_TSM(time, temp1, temp2, salt1, salt2, melt1, melt2, ens1, ens2, show = True, save = True):

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
        fig.savefig("PIctrlVSLENS_1920_comp_time.png")

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
def quiver_plot (u, v, lon, lat, exp, var, depth, ice_mask, year, interp = False, overlay = False, func = None, show = False, save = True):
    [LON,LAT] = np.meshgrid(lon, lat)
    [LON_short, LAT_short] = (lon[0:-1:20], lat[0:-1:20])

    if interp == True:
        u_interp = interpolate_currents(u, "zonal")
        v_interp = interpolate_currents(v, "meridional")
        [U, V] = np.meshgrid(u_interp, v_interp)
    else:
        [U,V] = np.meshgrid(u, v)

    # apply mask over land, i.e. where ocean depth is zero
    land_mask = np.zeros(np.shape(depth))
    land_mask[depth == 0] = 1

    # apply mask over the ice shelf (determiend by the ice-mask) and the continental shelf (roughly where the depth is less than 1500m)
    mask = np.zeros(np.shape(depth))
    mask[depth < 1500] = 1
    mask[ice_mask == 0] = 2
    
    # set the colors to block the continent (set to grey)
    colors = [
        (1.0, 1.0, 1.0, 0),
        (0.7, 0.7, 0.7, 1),
        (0.6, 0.6, 0.6, 1)]
    lim = np.max(np.abs(func))

    fig, ax = plt.subplots()
    
    if overlay == True:
        cp = plt.contourf(LON, LAT, func, levels = np.linspace(-lim, lim,15), cmap="seismic")
        cb = plt.colorbar(cp)
        plt.contourf(LON, LAT, land_mask, cmap = matplotlib.colors.ListedColormap(colors))
        plt.contour(LON, LAT, mask, 2, cmap = "Greys",linestyles='dashed')
        quiv = plt.quiver(LON_short, LAT_short, u[0:-1:20,0:-1:20], v[0:-1:20,0:-1:20], color = "white")
        plt.title(exp+" "+var)
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
    else:
        plt.contourf(LON, LAT, land_mask,cmap = matplotlib.colors.ListedColormap(colors))
        plt.contour(LON, LAT, mask, 2, cmap = "Greys",linestyles='dashed')
        plt.quiver(LON_short, LAT_short, U[0:-1:20,0:-1:20], V[0:-1:20,0:-1:20])
        plt.title(exp+" "+var)
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
    
    if show == True:
        plt.show()
    
    if save == True:
        fig.savefig(exp+"_quiv_"+var+"_"+year+".png")



def compare_quiver_plots (U1, V1, U2, V2, func1, func2, exp1, exp2, lon, lat, var, depth, ice_mask, year, interp = False, show = False, save = True):
    start_time = time.time()
    [LON,LAT] = np.meshgrid(lon, lat, copy=False)
    [LON_short, LAT_short] = np.meshgrid(lon[0:-1:20], lat[0:-1:20], copy=False)
    [u1,v1] = np.meshgrid(U1, V1, copy = False)
    [u2,v2] = np.meshgrid(U2, V2, copy=False)
    print("--- %s seconds ---" % (time.time() - start_time))
    # apply mask over land, i.e. where ocean depth is zero
    land_mask = np.zeros(np.shape(depth))
    land_mask[depth == 0] = 1

    # apply mask over the ice shelf (determiend by the ice-mask) and the continental shelf (roughly where the depth is less than 1500m)
    mask = np.zeros(np.shape(depth))
    mask[depth < 1500] = 1
    mask[ice_mask == 0] = 2
    
    # set the colors to block the continent (set to grey)
    colors = [
        (1.0, 1.0, 1.0, 0),
        (0.7, 0.7, 0.7, 1),
        (0.6, 0.6, 0.6, 1)]
    lim = np.max(np.abs(func1))
    anlim = np.max(np.abs(func1-func2))
    print("--- %s seconds ---" % (time.time() - start_time))
    
    fig =  plt.figure(figsize=(17,5))

    sc = 10

    plt.subplot(1,3,1) 
    #cp = plt.contourf(LON, LAT, func1, levels = np.linspace(-lim, lim,15), cmap="seismic")
    #cb = plt.colorbar(cp)
    #plt.contourf(LON, LAT, land_mask, cmap = matplotlib.colors.ListedColormap(colors))
    #plt.contour(LON, LAT, mask, 2, cmap = "Greys",linestyles='dashed')
    quiv = plt.quiver(LON_short, LAT_short, u1[0:-1:20,0:-1:20], v1[0:-1:20,0:-1:20], color = "k")
    plt.title(exp1+" "+var)
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")

    plt.subplot(1,3,2) 
    #cp = plt.contourf(LON, LAT, func2, levels = np.linspace(-lim, lim,15), cmap="seismic")
    #cb = plt.colorbar(cp)
    #plt.contourf(LON, LAT, land_mask, cmap = matplotlib.colors.ListedColormap(colors))
    #plt.contour(LON, LAT, mask, 2, cmap = "Greys",linestyles='dashed')
    quiv = plt.quiver(LON_short, LAT_short, u2[0:-1:20,0:-1:20], v2[0:-1:20,0:-1:20])
    plt.title(exp2+" "+var)
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")

    plt.subplot(1,3,3) 
    #cp = plt.contourf(LON, LAT, func1-func2, levels = np.linspace(-anlim, anlim,15), cmap="seismic")
    #cb = plt.colorbar(cp)
    #plt.contourf(LON, LAT, land_mask, cmap = matplotlib.colors.ListedColormap(colors))
    #plt.contour(LON, LAT, mask, 2, cmap = "Greys",linestyles='dashed')
    quiv = plt.quiver(LON_short, LAT_short, u1[0:-1:20,0:-1:20]-u2[0:-1:20,0:-1:20], v1[0:-1:20,0:-1:20]-v2[0:-1:20,0:-1:20])
    plt.title(exp1+" - "+exp2+" "+var)
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")

    # show figure
    if show == True:
        plt.show()
        
    # save figure
    if save == True:
        file_out = "PI_ctrl_comp_quiv_"+year+"_"+var+".png"
        fig.savefig(file_out)

if __name__ == "__main__":
    print("ERROR!! This is a file containing functions and cannot be run independently")



