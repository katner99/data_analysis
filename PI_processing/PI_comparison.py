from plots import compare_contour_plots_LATLON, create_mask
from funcs import read_variable, find_nearest
from mitgcm_python.grid import Grid

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import sys
import numpy as np
import xarray as xr


def main():
    """
    Main function that reads the command line arguments, checks them for validity,
    loads the input data and calls the appropriate plotting function.
    """
    ## check that enough arguments have been input by the user
    if len(sys.argv) != 4:
        sys.exit("Stopped - Incorrect number of arguements. Use python PI_comparison.py <year> <option> <var>")
    
    # check the option chosen is valid
    lemenu = ["oneVSone", "the_monster"]
    option = str(sys.argv[2])
    if option not in lemenu:
        sys.exit("Stopped - Invalid option. Please choose from my amazing menu selection of <oneVSone>, <the_monster>")
    
    # set up the variables you need
    year = sys.argv[1]
    var = str(sys.argv[3])
    save = True
    show = True
    
    if option == "the_monster":
        # load up the file paths for the monster, needed 4
        filepath_ctrl = "/data/oceans_output/shelf/katner33/PIctrl_output/CTRL_ensemble_mean_2090.nc"
        filepath_lens = "/data/oceans_output/shelf/katner33/PIctrl_output/LENS_ensemble_mean_2090.nc"
        filepath_wind = "/data/oceans_output/shelf/katner33/PIctrl_output/WIND_ensemble_mean_2090.nc"
        filepath_temp = "/data/oceans_output/shelf/katner33/PIctrl_output/THERM_ensemble_mean_2090.nc"

        # check if the input files exist
        for filepath in [filepath_ctrl, filepath_lens, filepath_wind, filepath_temp]:
            try:
                open(filepath)
            except FileNotFoundError:
                sys.exit(f"Stopped - Could not find input file {filepath}")

        # load up the input data    
        ctrl_input_data = xr.open_dataset(filepath_ctrl)
        lens_input_data = xr.open_dataset(filepath_lens)
        wind_input_data = xr.open_dataset(filepath_wind)
        temp_input_data = xr.open_dataset(filepath_temp)
        
        # read in the general variables (these should be the same between the ensembles
        lat = ctrl_input_data.YC.values
        lon = ctrl_input_data.XC.values
        ice_mask = ctrl_input_data.maskC.values[0,:,:]
        depth = ctrl_input_data.Depth.values
        grid_filepath = "/data/oceans_output/shelf/kaight/archer2_mitgcm/PAS_LENS001_O/output/192001/MITgcm/output.nc"
        grid = Grid(grid_filepath)
   
        # temperature
        if var == "THETA":
            depth_range = [find_nearest(ctrl_input_data.Z.values, -200), find_nearest(ctrl_input_data.Z.values, -700)]
            data_ctrl = read_variable(ctrl_input_data, var, grid, depth_range)
            data_lens = read_variable(lens_input_data, var, grid, depth_range)
            data_wind = read_variable(wind_input_data, var, grid, depth_range)
            data_temp = read_variable(temp_input_data, var, grid, depth_range)
            color_scheme = "coolwarm"

        # sea ice tickness
        elif var in ["SIheff", "oceFWflx", "SIfwmelt", "SIfwfrz"]:
            data_ctrl = read_variable(ctrl_input_data, var, grid)
            data_lens = read_variable(lens_input_data, var, grid)
            data_wind = read_variable(wind_input_data, var, grid)
            data_temp = read_variable(temp_input_data, var, grid)
            #color_scheme = "YlGnBu_r"
            color_scheme = "BrBG"

        # salinity
        elif var == "SALT":
            data_ctrl = read_variable(ctrl_input_data, var, grid)
            data_lens = read_variable(lens_input_data, var, grid)
            data_wind = read_variable(wind_input_data, var, grid)
            data_temp = read_variable(temp_input_data, var, grid)
            color_scheme = "PRGn_r"

        # set mask
        [land_mask, mask, colors] = create_mask(depth, ice_mask)

        # set up the grid
        [X, Y] = np.meshgrid(lon, lat)
        
        
        # graph parameters
        graph_params = {
            "figsize": (18, 12),
            "font_size": 12,
            "low_val": int(np.nanmin(data_lens)),
            "high_val": int(np.nanmax(data_ctrl)),
            "step": 15,
            "low_val_anom": -1,
            "high_val_anom": 1,
            "ticks_anom": np.arange(-1, 1.5, 0.5)
        }

        # create subplots with each variable on a new line
        fig, axs = plt.subplots(nrows=4, ncols=4, gridspec_kw={"hspace": 0.5, "wspace": 0.4}, figsize=graph_params["figsize"])
        
        axs = axs.flatten()
        
        axs[1].axis("off")
        axs[2].axis("off")
        axs[3].axis("off")
        axs[6].axis("off")
        axs[7].axis("off")
        axs[11].axis("off")

        # PI_ctrl
        cs = axs[0].contourf(X, Y, data_ctrl, cmap=color_scheme, extend="both", levels=np.linspace(graph_params["low_val"], graph_params["high_val"], graph_params["step"]))
        axs[0].contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
        axs[0].contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
        fig.colorbar(cs, ax=axs[0], ticks=np.arange(int(graph_params["low_val"]), int(graph_params["high_val"]), 0.5))
        axs[0].set_title("PI_ctrl", fontsize=graph_params["font_size"], weight="bold")
        axs[0].set_ylabel("Latitude", fontsize=graph_params["font_size"])
        
        # LENS - CTRL
        cs = axs[4].contourf(X, Y, data_lens - data_ctrl, cmap=color_scheme, extend="both", levels=np.linspace(graph_params["low_val_anom"], graph_params["high_val_anom"], graph_params["step"]))
        axs[4].contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
        axs[4].contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
        fig.colorbar(cs, ax=axs[4], ticks=graph_params["ticks_anom"])
        axs[4].set_title("LENS  - PI_ctrl", fontsize=graph_params["font_size"])
        axs[4].set_ylabel("Latitude", fontsize=graph_params["font_size"])
        
        # LENS
        cs = axs[5].contourf(X, Y, data_lens, cmap=color_scheme, extend="both", levels=np.linspace(graph_params["low_val"], graph_params["high_val"], graph_params["step"]))
        axs[5].contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
        axs[5].contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
        fig.colorbar(cs, ax=axs[5], ticks=np.arange(int(graph_params["low_val"]), int(graph_params["high_val"]), 0.5))
        axs[5].set_title("LENS", fontsize=graph_params["font_size"], weight="bold")
        
        
        # WIND - CTRL
        cs = axs[8].contourf(X, Y, data_wind - data_ctrl, cmap=color_scheme, extend="both", levels=np.linspace(graph_params["low_val_anom"], graph_params["high_val_anom"], graph_params["step"]))
        axs[8].contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
        axs[8].contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
        fig.colorbar(cs, ax=axs[8], ticks=graph_params["ticks_anom"])
        axs[8].set_title("WIND  - PI_ctrl", fontsize=graph_params["font_size"])
        axs[8].set_ylabel("Latitude", fontsize=graph_params["font_size"])

        # WIND - LENS
        cs = axs[9].contourf(X, Y, data_wind - data_lens, cmap=color_scheme, extend="both", levels=np.linspace(graph_params["low_val_anom"], graph_params["high_val_anom"], graph_params["step"]))
        axs[9].contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
        axs[9].contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
        fig.colorbar(cs, ax=axs[9], ticks=graph_params["ticks_anom"])
        axs[9].set_title("WIND  - LENS", fontsize=graph_params["font_size"])
        
        # WIND
        cs = axs[10].contourf(X, Y, data_wind, cmap=color_scheme, extend="both", levels=np.linspace(graph_params["low_val"], graph_params["high_val"], graph_params["step"]))
        axs[10].contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
        axs[10].contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
        fig.colorbar(cs, ax=axs[10], ticks=np.arange(int(graph_params["low_val"]), int(graph_params["high_val"]), 0.5))
        axs[10].set_title("WIND", fontsize=graph_params["font_size"], weight="bold")
        
        # THERMO - CTRL
        cs = axs[12].contourf(X, Y, data_temp - data_ctrl, cmap=color_scheme, extend="both", levels=np.linspace(graph_params["low_val_anom"], graph_params["high_val_anom"], graph_params["step"]))
        axs[12].contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
        axs[12].contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
        fig.colorbar(cs, ax=axs[12], ticks=graph_params["ticks_anom"])
        axs[12].set_title("THERMO  - CTRL", fontsize=graph_params["font_size"])
        axs[12].set_xlabel("Longitude", fontsize=graph_params["font_size"])
        axs[12].set_ylabel("Latitude", fontsize=graph_params["font_size"])
        
        # THERMO - LENS
        cs = axs[13].contourf(X, Y, data_temp - data_lens, cmap=color_scheme, extend="both", levels=np.linspace(graph_params["low_val_anom"], graph_params["high_val_anom"], graph_params["step"]))
        axs[13].contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
        axs[13].contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
        fig.colorbar(cs, ax=axs[13], ticks=graph_params["ticks_anom"])
        axs[13].set_title("THERMO  - LENS", fontsize=graph_params["font_size"])
        axs[13].set_xlabel("Longitude", fontsize=graph_params["font_size"])

        # THERMO - WIND
        cs = axs[14].contourf(X, Y, data_temp - data_wind, cmap=color_scheme, extend="both", levels=np.linspace(graph_params["low_val_anom"], graph_params["high_val_anom"], graph_params["step"]))
        axs[14].contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
        axs[14].contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
        fig.colorbar(cs, ax=axs[14], ticks=graph_params["ticks_anom"])
        axs[14].set_title("THERMO  - WIND", fontsize=graph_params["font_size"])
        axs[14].set_xlabel("Longitude", fontsize=graph_params["font_size"])
        
        # THERMO
        cs = axs[15].contourf(X, Y, data_temp, cmap=color_scheme, extend="both", levels=np.linspace(graph_params["low_val"], graph_params["high_val"], graph_params["step"]))
        axs[15].contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
        axs[15].contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
        fig.colorbar(cs, ax=axs[15], ticks=np.arange(int(graph_params["low_val"]), int(graph_params["high_val"]), 0.5))
        axs[15].set_title("THERMO", fontsize=graph_params["font_size"], weight="bold")
        axs[15].set_xlabel("Longitude", fontsize=graph_params["font_size"])
        
        fig.suptitle(var, fontsize=16)

        # save figure
        if save == True:
            file_out ="comp_cont_2090_"+var+".png"
            fig.savefig(file_out)

        # show figure
        if show == True:
            plt.show()
            
if __name__ == '__main__':
    main() # run the program


