import sys

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import netCDF4 as nc
import numpy as np
import xarray as xr

def main():
    # check that enough arguments have been input by the user
    if len(sys.argv) != 2:
        sys.exit("Stopped - Incorrect number of arguements. Use python PI_slice.py <option>")
    
    # check the option cosen is valid
    lemenu = ["compare_slice", "animate_slice", "LENS_WND_CTRL"]

    option = str(sys.argv[1])

    if option not in lemenu:
        sys.exit("Stopped - Invalid option. Please choose from my amazing menu selection of <compare_slice> <animate_slice> <LENS_WND_CTRL>")

    if option == "compare_slice":
        # load up the file paths for the monster, needed 4
        filepath_ctrl = "/data/oceans_output/shelf/katner33/PIctrl_output/average_historic.nc"
        filepath_lens = "/data/oceans_output/shelf/katner33/PIctrl_output/average_current.nc"
        filepath_wind = "/data/oceans_output/shelf/katner33/PIctrl_output/average_future_2.nc"

        # check if the input files exist
        for filepath in [filepath_ctrl, filepath_lens, filepath_wind]:
            try:
                open(filepath)
            except FileNotFoundError:
                sys.exit(f"Stopped - Could not find input file {filepath}")

        # load up the input data    
        ctrl_input_data = xr.open_dataset(filepath_ctrl)
        lens_input_data = xr.open_dataset(filepath_lens)
        wind_input_data = xr.open_dataset(filepath_wind)
        #temp_input_data = xr.open_dataset(filepath_temp)
        
        # read in the general variables (these should be the same between the ensembles
        z = ctrl_input_data.Z.values
        lon = ctrl_input_data.XC.values[329:350]
        ice_mask = ctrl_input_data.maskC.values[:,98,329:350]
        color_scheme = "PRGn"
       
        # load the variables and calculate the average over the year
        days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        
        ctrl_vvel = np.average(ctrl_input_data.VVEL.values[:,:,98,329:350], axis = 0, weights = days_in_month)
        lens_vvel = np.average(lens_input_data.VVEL.values[:,:,98,329:350], axis = 0, weights = days_in_month)
        wind_vvel = np.average(wind_input_data.VVEL.values[:,:,98,329:350], axis = 0, weights = days_in_month)

        ctrl_salt = np.average(ctrl_input_data.SALT.values[:,:,98,329:350], axis = 0, weights = days_in_month)
        lens_salt = np.average(lens_input_data.SALT.values[:,:,98,329:350], axis = 0, weights = days_in_month)
        wind_salt = np.average(wind_input_data.SALT.values[:,:,98,329:350], axis = 0, weights = days_in_month)
        #temp_vvel = np.average(temp_input_data.VVEL.values[:,:,98,329:350], axis = 0, weights = days_in_month)
        
        ctrl_vvel[ice_mask == 0] = np.nan
        lens_vvel[ice_mask == 0] = np.nan
        wind_vvel[ice_mask == 0] = np.nan
        #temp_vvel[ice_mask == 0] = np.nan
        
        low_val = np.nanmin(ctrl_vvel)
        high_val = -low_val
        low_salt = 33
        high_salt = np.nanmax(ctrl_salt)
        print(low_salt, high_salt)
        step = 15

        font_size = 16
        
        # create subplots with each variable on a new line
        fig, axs = plt.subplots(nrows=1, ncols=3, gridspec_kw={"hspace": 0.5, "wspace": 0.4}, figsize=(15,5))
        
        axs = axs.flatten()

        # PI_ctrl
        cs = axs[0].contourf(lon, z, ctrl_vvel, cmap=color_scheme, extend="both", levels=np.linspace(low_val, high_val, step))
        salt = axs[0].contour(lon, z, ctrl_salt, levels=np.linspace(low_salt, high_salt, 10), colors = "k")
        axs[0].clabel(salt, fontsize=8, inline=1, fmt = '%0.1f')
        axs[0].set_title("Pre-industrial", weight="bold", fontsize = font_size)
        axs[0].set_ylabel("Depth (m)", fontsize = font_size)
        axs[0].set_xlabel("Longitude", fontsize = font_size)
        axs[0].set_ylim([-800, 0])
        
        # LENS
        cs = axs[1].contourf(lon, z, lens_vvel, cmap=color_scheme, extend="both", levels=np.linspace(low_val, high_val, step))
        salt = axs[1].contour(lon, z, lens_salt, levels=np.linspace(low_salt, high_salt, 10), colors = "k")
        axs[1].clabel(salt, fontsize=8, inline=1, fmt = '%0.1f')
        axs[1].set_title("Current Winds", weight="bold", fontsize = font_size)
        axs[1].set_xlabel("Longitude", fontsize = font_size)
        axs[1].set_ylim([-800, 0])
        
        # WIND
        cs = axs[2].contourf(lon, z, wind_vvel, cmap=color_scheme, extend="both", levels=np.linspace(low_val, high_val, step))
        salt = axs[2].contour(lon, z, wind_salt, levels=np.linspace(low_salt, high_salt, 10), colors = "k")
        axs[2].clabel(salt, fontsize=8, inline=1, fmt = '%0.1f')
        axs[2].set_title("Future Winds", weight="bold", fontsize = font_size)
        axs[2].set_xlabel("Longitude", fontsize = font_size)
        axs[2].set_ylim([-800, 0])
        
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        fig.colorbar(cs, cax=cbar_ax, ticks = np.arange(-0.1, 0.12, 0.05))
        
        plt.suptitle("Meridional velocity across the Pine Island Bay Trough at 73°S", fontsize = font_size, weight = "bold")
        
        fig.savefig("poster_slice.png")
        plt.show()
        
        
    if option == "animate_slice":

        # load the data
        filepath = "/users/katner33/ensemble_mean/"
        filename = year+"_ensemble_mean.nc"
        id = nc.Dataset(filepath+filename, 'r')
        theta = id_ctrl.variables["THETA"][:,:,98,329:350]

        z = id.variables["Z"][:]
        lon = id.variables["XC"][329:350]
        mask = id.variables["maskC"][:,98,329:350]
        
        # animate_contour (lon, z, theta, year, var = "theta", depth, exp = "PI_ctrl", color_scheme = "jet", mask = mask, lonlatplot = False, xlabel = "Longitude", ylabel = "Depth", show=True, save=True)
    
    # compare between the different simulations
    if option == "LENS_WND_CTRL":
        # set the filepath
        filepath_ctrl = "/data/oceans_output/shelf/katner33/PIctrl_output/average_historic.nc"
        filepath_lens = "/data/oceans_output/shelf/katner33/PIctrl_output/average_current.nc"
        filepath_wind = "/data/oceans_output/shelf/katner33/PIctrl_output/average_future_2.nc"
        

        # load up the files
        input_ctrl = xr.open_dataset(filepath_ctrl)
        input_lens = xr.open_dataset(filepath_lens)
        input_wind = xr.open_dataset(filepath_wind)
        
        
        # read in the general variables (these should be the same between the ensembles
        z = input_ctrl.Z.values
        lon = input_ctrl.XC.values[329:350]
        mask = input_ctrl.maskC.values[:,98,329:350]
        
        color_scheme = "coolwarm"

        # load up the datasets
        [data_ctrl, salt_ctrl] = read_xarray(var, input_ctrl, mask)
        [data_lens, salt_lens] = read_xarray(var, input_lens, mask)
        [data_wind, salt_wind] = read_xarray(var, input_wind, mask)

        # set up graph
        [X, Y] = np.meshgrid(lon, z)
        
        font_size = 10
        low_val = np.nanmin(data_ctrl)
        high_val = np.nanmax(data_ctrl)
        tick = np.arange(-2, 1.6, 0.5)

        low_val_anom = -1.5
        high_val_anom = 1.5
        ticks_anom = np.arange(low_val_anom, high_val_anom+0.0001, 0.05)
        #ticks_anom = None

        save = True
        show = True
        
        fig =  plt.figure(figsize=(18,12))
        
        salt_range = np.arange(34, 35, 0.1)
        salt_anom_range = np.arange(-0.5, 0.5, 0.1)
        
        #fig.subplots_adjust(top=1)

        # create subplot with each variable on a new line
        # PI_ctrl
        plt.subplot(4,4,1) 
        cs = plt.contourf(X, Y, data_ctrl, levels=np.linspace(low_val, high_val, 10), extend = "both", cmap = color_scheme)
        salt = plt.contour(X, Y, salt_ctrl, levels=Contourrange, colors = "k")
        plt.clabel(salt, fontsize=8, inline=1,fmt = '%0.1f',ticks=Contourrange)
        plt.colorbar(cs, ticks = tick)
        #plt.colorbar(cs)
        plt.ylim(-800, 0)
        plt.title("PI_ctrl", fontsize = font_size, weight='bold')
        #plt.xlabel('Longitude', fontsize = font_size)
        plt.ylabel('Depth', fontsize = font_size)
        
        # LENS - CTRL
        plt.subplot(4,4,5) 
        cs = plt.contourf(X, Y, data_lens - data_ctrl, levels=np.linspace(low_val_anom, high_val_anom, 10), cmap = color_scheme)
        salt= plt.contour(X, Y, salt_lens - salt_ctrl, levels = Anomrange, colors = "k")
        plt.clabel(salt, fontsize=8, inline=1,fmt = '%0.01f',ticks=Anomrange)
        plt.colorbar(cs, ticks=np.arange(-2, 2.3, 0.5))
        #plt.colorbar(cs)
        plt.ylim(-800, 0)
        plt.title("LENS - PI_ctrl", fontsize = font_size)
        #plt.xlabel('Longitude', fontsize = font_size)
        plt.ylabel('Depth', fontsize = font_size)

        # LENS
        plt.subplot(4, 4, 6) 
        cs = plt.contourf(X, Y, data_lens, levels=np.linspace(low_val, high_val, 10), extend = "both", cmap = color_scheme)
        salt = plt.contour(X, Y, salt_lens, levels=Contourrange, colors = "k")
        plt.clabel(salt, fontsize=8, inline=1,fmt = '%0.1f',ticks=Contourrange)
        plt.colorbar(cs, ticks = tick)
        plt.ylim(-800, 0)
        #plt.colorbar(cs)
        plt.title("LENS", fontsize = font_size, weight='bold')
        #plt.xlabel('Longitude', fontsize = font_size)
        #plt.ylabel('Depth', fontsize = font_size)
        
        # WIND - CTRL
        plt.subplot(4, 4, 9) 
        cs = plt.contourf(X, Y, data_wind - data_ctrl, levels=np.linspace(low_val_anom, high_val_anom, 10),cmap = color_scheme)
        #cs = plt.contourf(X, Y, data_wind - data_ctrl, levels=np.linspace(-2.5, 2.5,10),extend = "both", cmap = color_scheme2)
        #plt.contourf(X, Y, land_mask, cmap = matplotlib.colors.ListedColormap(colors))
        salt = plt.contour(X, Y, salt_wind - salt_ctrl, colors = "k")
        #plt.colorbar(cs)
        plt.clabel(salt, fontsize=8, inline=1,fmt = '%0.01f',ticks=Anomrange)
        plt.colorbar(cs, ticks=np.arange(-1, 1.3, 0.5))
        plt.ylim(-800, 0)
        plt.title("WIND - PI_ctrl", fontsize = font_size)
        #plt.xlabel('Longitude', fontsize = font_size)
        plt.ylabel('Depth', fontsize = font_size)

        # WIND - LENS
        plt.subplot(4, 4, 10) 
        cs = plt.contourf(X, Y, data_wind - data_lens, levels=np.linspace(low_val_anom, high_val_anom, 10),cmap = color_scheme)
        #cs = plt.contourf(X, Y, data_wind - data_lens, levels=np.linspace(-2.5, 2.5, 10),extend = "both", cmap = color_scheme2)
        #plt.contourf(X, Y, land_mask, cmap = matplotlib.colors.ListedColormap(colors))
        salt = plt.contour(X, Y, salt_wind - salt_lens, colors = "k")
        plt.colorbar(cs, ticks = np.arange(-2, 2.3, 0.5))
        plt.clabel(salt, fontsize=8, inline=1,fmt = '%0.01f',ticks=Anomrange)
        #plt.colorbar(cs)
        plt.ylim(-800, 0)
        plt.title("WIND - LENS", fontsize = font_size)
        #plt.xlabel('Longitude', fontsize = font_size)
        #plt.ylabel('Latitude', fontsize = font_size)
        
        # WIND
        plt.subplot(4, 4, 11) 
        cs = plt.contourf(X, Y, data_wind, levels=np.linspace(low_val, high_val, 10), extend = "both", cmap = color_scheme)
        #plt.contourf(X, Y, land_mask, cmap = matplotlib.colors.ListedColormap(colors))
        salt=plt.contour(X, Y, salt_wind, colors = "k")
        plt.colorbar(cs, ticks = tick)
        plt.clabel(salt, fontsize=8, inline=1,fmt = '%0.1f',ticks=Contourrange)
        plt.ylim(-800, 0)
        #plt.colorbar(cs)
        plt.title("WIND", fontsize = font_size, weight='bold')
        #plt.xlabel('Longitude', fontsize = font_size)
        #plt.ylabel('Depth', fontsize = font_size)

         # THERMO - CTRL
        plt.subplot(4, 4, 13) 
        cs = plt.contourf(X, Y, data_temp - data_ctrl, levels=np.linspace(low_val_anom, high_val_anom, 10), cmap = color_scheme)
        #cs = plt.contourf(X, Y, data_temp - data_ctrl, levels=np.linspace(-2.5,2.5, 10),extend = "both", cmap = color_scheme2)
        #plt.contourf(X, Y, land_mask, cmap = Matplotlib.colors.ListedColormap(colors))
        salt = plt.contour(X, Y, salt_temp - salt_ctrl,colors = "k")
        plt.colorbar(cs, ticks = np.arange(-1, 1.3, 0.5))
        plt.clabel(salt, fontsize=8, inline=1,fmt = '%0.01f',ticks=Anomrange)
        plt.ylim(-800, 0)
        #plt.colorbar(cs)
        plt.title("THERMO - PI_ctrl", fontsize = font_size)
        plt.xlabel('Longitude', fontsize = font_size)
        plt.ylabel('Depth', fontsize = font_size)
        
        # TEMP - LENS
        plt.subplot(4, 4, 14) 
        cs = plt.contourf(X, Y, data_temp - data_lens,levels=np.linspace(low_val_anom, high_val_anom, 10), cmap = color_scheme)
        #cs = plt.contourf(X, Y, data_temp - data_lens,levels=np.linspace(-2.5,2.5, 10), extend = "both", cmap = color_scheme2)
        #plt.contourf(X, Y, land_mask, cmap = matplotlib.colors.ListedColormap(colors))
        salt = plt.contour(X, Y, salt_temp - salt_lens, colors = "k")
        plt.colorbar(cs, ticks=np.arange(-1, 1.3, 0.5))
        plt.clabel(salt, fontsize=8, inline=1,fmt = '%0.01f',ticks=Anomrange)
        plt.ylim(-800, 0)
        #plt.colorbar(cs)
        plt.title("THERMO - LENS", fontsize = font_size)
        plt.xlabel('Longitude', fontsize = font_size)
        #plt.ylabel('Depth', fontsize = font_size)

        # TEMP - WIND
        plt.subplot(4, 4, 15) 
        cs = plt.contourf(X, Y, data_temp - data_wind,levels=np.linspace(low_val_anom, high_val_anom, 10),cmap = color_scheme)
        #cs = plt.contourf(X, Y, data_temp - data_wind,levels=np.linspace(-2.5,2.5, 10), extend = "both",cmap = color_scheme2)
        #plt.contourf(X, Y, land_mask, cmap = matplotlib.colors.ListedColormap(colors))
        salt = plt.contour(X, Y, salt_temp - salt_wind, colors = "k")
        plt.colorbar(cs, ticks = np.arange(-1.5, 1.8, 0.5))
        plt.clabel(salt, fontsize=8, inline=1,fmt = '%0.01f',ticks=Anomrange)
        #plt.colorbar(cs)
        plt.ylim(-800, 0)
        plt.title("THERMO - WIND", fontsize = font_size)
        plt.xlabel('Longitude', fontsize = font_size)
        #plt.ylabel('Depth', fontsize = font_size)
        
        # THERMO
        plt.subplot(4, 4, 16) 
        cs = plt.contourf(X, Y, data_temp, levels=np.linspace(low_val, high_val, 10), extend = "both", cmap = color_scheme)
        #plt.contourf(X, Y, land_mask, cmap = matplotlib.colors.ListedColormap(colors))
        salt = plt.contour(X, Y, salt_temp, colors = "k")
        #plt.colorbar(cs)
        plt.ylim(-800, 0)
        plt.colorbar(cs, ticks = tick)
        plt.clabel(salt,fontsize=8, inline=1,fmt = '%0.1f',ticks=Contourrange)
        plt.title("THERMO", fontsize = font_size, weight='bold')
        plt.xlabel('Longitude', fontsize = font_size)
        #plt.ylabel('Depth', fontsize = font_size)
        
        #fig.suptitle(var+" ensemble mean 2090")
       
        fig.tight_layout(pad=2.5)

        # save figure
        if save == True:
            file_out ="comp_slice_2090_"+var+".png"
            fig.savefig(file_out)

        # show figure
        if show == True:
            plt.show()
if __name__ == '__main__':
    main() # run the program
