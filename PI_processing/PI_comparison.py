from funcs import read_variable, find_nearest
from plots import create_mask
from directories_and_paths import output_path, grid_filepath
from mitgcm_python.grid import Grid
from mitgcm_python.plot_utils.labels import lat_label, lon_label

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import sys
import numpy as np
import xarray as xr

def comparison(graph_params, data, X, Y, color_scheme, land_mask, colors, mask, experiment, title, file_out, save=True, show=False, linearity=False, residual = None, ice_shelf = False):
    total = len(data)
    
    fig, axs = plt.subplots(nrows=total, ncols=total, gridspec_kw={"hspace": 0.05, "wspace": 0.04}, figsize=graph_params["figsize"])
    axs = axs.flatten()

    if linearity:
        position = total-1
        cs = axs[position].contourf(X, Y, residual, cmap= "PRGn_r", extend="both", levels=np.linspace(-1, 1, 25))
        axs[position].contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
        axs[position].contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
        if ice_shelf:
            axs[position].set_ylim([-75.6, -70])
            axs[position].set_xlim([235, 265])
        else:
            axs[position].set_xlim([230, 265])
            axs[position].set_ylim([-75.5, -68])
        fig.colorbar(cs, ax=axs[position], ticks=np.arange(-1, 1.1, 1))
        axs[position].set_title("Residual", fontsize=graph_params["font_size"], weight="bold")
        axs[position].get_yaxis().set_visible(True)
        axs[position].get_xaxis().set_visible(True)

        lon_ticks = axs[position].get_xticks() - 360
        lon_labels = []
        for x in lon_ticks:
            lon_labels.append(lon_label(x,2))
        axs[position].set_xticklabels(lon_labels)
        axs[position].tick_params(axis='x', labelrotation=45)
        
        lat_ticks = axs[position].get_yticks()
        lat_labels = []
        for y in lat_ticks:
            lat_labels.append(lat_label(y,2))
        axs[position].set_yticklabels(lat_labels)
    
    for i in range(total):

        # MAIN GRAPH
        diagonal = (total+1)*i

        cs_diag = axs[diagonal].contourf(X, Y, data[i], cmap=color_scheme, extend="both", levels=np.linspace(graph_params["low_val"], graph_params["high_val"], graph_params["step"]))
        axs[diagonal].contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
        axs[diagonal].contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
        if ice_shelf:
            axs[diagonal].set_ylim([-75.6, -70])
            axs[diagonal].set_xlim([235, 265])
        else: 
            axs[diagonal].set_xlim([230, 265])
            axs[diagonal].set_ylim([-75.5, -68])
        #ticks=np.arange(graph_params["low_val"], graph_params["high_val"]+0.1, 2.5)
        #cbar = fig.colorbar(cs, ax=axs[diagonal])
        #cbar.set_ticks(ticks)
        axs[diagonal].set_title(experiment[i], fontsize=graph_params["font_size"], weight="bold")
        if diagonal > 0:
            axs[diagonal].get_yaxis().set_visible(False)
        else:
            lat_ticks = axs[diagonal].get_yticks()
            lat_labels = []
            for y in lat_ticks:
                lat_labels.append(lat_label(y,2))
            axs[diagonal].set_yticklabels(lat_labels)

        if diagonal != ((total*total)-1):
            axs[diagonal].get_xaxis().set_visible(False)
        else:
            lon_ticks = axs[diagonal].get_xticks() - 360
            lon_labels = []
            for x in lon_ticks:
                lon_labels.append(lon_label(x,2))
            axs[diagonal].set_xticklabels(lon_labels)
            axs[diagonal].tick_params(axis='x', labelrotation=45)

        # ANOMALY
        for j in range(i+1, total):
            anomaly = (total*j)+i

            cs_anom = axs[anomaly].contourf(X, Y, data[j] - data[i], cmap= "PRGn_r", extend="both", levels=np.linspace(graph_params["low_val_anom"], graph_params["high_val_anom"], graph_params["step"]))
            axs[anomaly].contourf(X, Y, land_mask, cmap=matplotlib.colors.ListedColormap(colors))
            axs[anomaly].contour(X, Y, mask, 2, cmap="Greys", linestyles="dashed")
            #ticks=graph_params["ticks_anom"]
            if ice_shelf:
                axs[anomaly].set_ylim([-75.6, -70])
                axs[anomaly].set_xlim([235, 265])
            else:
                axs[anomaly].set_xlim([230, 265])
                axs[anomaly].set_ylim([-75.5, -68])
            #cbar = fig.colorbar(cs, ax=axs[anomaly])
            #cbar.set_ticks(ticks)
            #axs[anomaly].set_title(experiment[j]+" - "+experiment[i], fontsize=graph_params["font_size"])
            axs[anomaly].text(
                0.5,
                0.85,
                experiment[j]+" - "+experiment[i],
                horizontalalignment="center",
                transform=axs[anomaly].transAxes,
                bbox=dict(facecolor="white", alpha=0.9),
            )
            if anomaly % total != 0:
                axs[anomaly].get_yaxis().set_visible(False)
            else:
                lat_ticks = axs[anomaly].get_yticks()
                lat_labels = []
                for y in lat_ticks:
                    lat_labels.append(lat_label(y,2))
                axs[anomaly].set_yticklabels(lat_labels)

            if anomaly < (total * (total - 1)):
                axs[anomaly].get_xaxis().set_visible(False)
            else:
                lon_ticks = axs[anomaly].get_xticks() - 360
                lon_labels = []
                for x in lon_ticks:
                    lon_labels.append(lon_label(x,2))
                axs[anomaly].set_xticklabels(lon_labels)
                axs[anomaly].tick_params(axis='x', labelrotation=45)
        # discarded graphs
        for k in range(i):
            axs[i + (total*k)].axis("off")
    
    ticks=np.arange(graph_params["low_val"], graph_params["high_val"]+0.1, 1)
    cbar_ax = fig.add_axes([0.05, 0.525, 0.02, 0.4])
    cbar = plt.colorbar(cs_diag, cax=cbar_ax, orientation='vertical')
    cbar.set_ticks(ticks)
    cbar.ax.yaxis.set_ticks_position('left')

    ticks=graph_params["ticks_anom"]
    cbar_ax = fig.add_axes([0.05, 0.1, 0.02, 0.4])
    cbar = plt.colorbar(cs_anom, cax=cbar_ax, orientation='vertical')
    cbar.set_ticks(ticks)
    cbar.ax.yaxis.set_ticks_position('left')

    fig.suptitle(title, fontsize=16)

    # save figure
    if save == True:
        fig.savefig(file_out)

    # show figure
    if show == True:
        plt.show()
       

def main():
    """
    Main function that reads the command line arguments, checks them for validity,
    loads the input data and calls the appropriate plotting function.
    """
        
    # set up the variables you need
    var = "THETA"
    save = True
    show = True
    period = "2070-2100"
    
    # load up the file paths for the monster, needed 4
    filepaths = [output_path + "average_" + ens + "_" + period + ".nc" for ens in ["CTRL", "LENS", "WIND", "TEMP"]]
    # check if the input files exist
    for filepath in filepaths:
        try:
            open(filepath)
        except FileNotFoundError:
            sys.exit(f"Stopped - Could not find input file {filepath}")

    experiment = ["pre-industrial", "RCP 8.5", "wind f.", "thermo. f."]

    # load up the input data
    input_data = [xr.open_dataset(filepath, decode_times = False) for filepath in filepaths]

    # read in the general variables (these should be the same between the ensembles
    [lat, lon, ice_mask_temp, depth] = [input_data[0][param].values for param in ["YC", "XC", "maskC", "Depth"]]
    ice_mask = ice_mask_temp[0,:,:]

    grid = Grid(grid_filepath)

    # temperature
    if var == "THETA":
        depth_range = [find_nearest(input_data[0]["Z"].values, -200), find_nearest(input_data[0]["Z"].values, -700)]
        data = [read_variable(input, var, grid, depth_range) for input in input_data]
        color_scheme = "coolwarm"
        anom = 1.5
        min_val = -2
        max_val =  2.1
        title = f"Average temperature between 200 and 700m {period}"

    # sea ice tickness
    elif var in ["SIheff", "oceFWflx", "SIfwmelt", "SIfwfrz", "EXFvwind", "oceQnet"]:
        data = [read_variable(input, var, grid)*3600*24*365/1000 for input in input_data]
        #data =3600*24*365
        #color_scheme = "YlGnBu_r"
        color_scheme = "PRGn_r"
        anom = 1.5
        min_val = -5
        max_val = 5
        #color_scheme = "seismic"
        title = f"Freshwater fluxes m/yr {period}"

    # salinity
    elif var == "SALT":
        depth_range = [find_nearest(input_data[0]["Z"].values, -200), find_nearest(input_data[0]["Z"].values, -700)]
        data = [read_variable(input, var, grid, depth_range) for input in input_data]
        color_scheme = "PRGn_r"
        anom = 0.25
        min_val = int(np.min(data[0]))
        max_val =  int(np.max(data[0]))+1
        
    elif var == "SHIfwFlx":
        data = [-read_variable(input, var, grid)*3600*24*30*10**(-3) for input in input_data]
        color_scheme = "rainbow"
        anom = 25
        min_val = 0
        max_val = 76
        title = f"ice shelf basal melt rate ({period}) (m.w.e./yr)"

    # set mask
    [land_mask, mask, colors] = create_mask(depth, ice_mask)

    # set up the grid
    [X, Y] = np.meshgrid(lon, lat)

    # graph parameters 1.333112e-05 -1.0214189e-06
    graph_params = {
        "figsize": (18, 12),
        "font_size": 12,
        "low_val": min_val,
        "high_val": max_val,
        "step": 15,
        "low_val_anom": -anom,
        "high_val_anom": anom,
        "ticks_anom": np.arange(-anom, anom +0.2, 0.5)
    }
    
    file_out = "mega_comparison"+var+"_"+period+"_test.png"

    print(np.shape(data[1]))

    residual = data[0]+data[1]-data[2]-data[3]

    print(np.max(residual), np.min(residual))

    comparison(graph_params, data, X, Y, color_scheme, land_mask, colors, mask, experiment, title, file_out, save, show, linearity = True, residual = residual, ice_shelf = False)
             
if __name__ == '__main__':
    main() # run the program

