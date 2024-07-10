import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mitgcm_python.plot_utils.labels import lat_label
import numpy as np
from directories_and_paths import output_path
from funcs import read_data
import xarray as xr

def paper():
    UC = 0 # undercurrent section
    [input_data, vel] = read_data("UVEL", UC)
    [input_data, salt] = read_data("DENSITY", UC)
      
    z = input_data[0].z.values
    lat = input_data[0]["lat"]
    color_scheme = "PiYG_r"
    experiment = ["NONE", "ALL", "WIND", "THERMO"]
    low_val = -e-13
    high_val = -low_val
    print(low_sal, high_sal)
    step = 15
    font_size = 20
        
    # create subplots with each variable on a new line
    fig, axs = plt.subplots(nrows=1, ncols=3, gridspec_kw={"hspace": 0.08, "wspace": 0.05}, figsize=(20,8)) 
    axs = axs.flatten()

    fmt = ticker.ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((0, 0))
    for i in range(1, 4):
        cv = axs[i-1].contourf(lat, z, salt[i]-salt[0], cmap=color_scheme, extend="both", levels=np.linspace(low_sal, high_sal, step))
        cs = axs[i-1].contour(lat, z, vel[i]-vel[0], colors='black')
        
        axs[i-1].set_title(experiment[i], weight="bold", fontsize = font_size)
        axs[i-1].set_xlabel("Latitude", fontsize = font_size)
        axs[i-1].set_ylim([-1000, 0])
        axs[i-1].set_xlim([-74, -71.5])
        axs[i-1].tick_params(axis='x', labelsize=16, rotation=45)
        axs[i-1].tick_params(axis='y', labelsize=16)
        
        if i > 0:
            axs[i-1].get_yaxis().set_visible(False)
            axs[i-1].set_ylabel("Depth (m)", fontsize = font_size)
        
        axs[i-1].locator_params(axis='x', nbins=6)
        lat_ticks = axs[i-1].get_xticks()
        lat_labels = []
        for x in lat_ticks[:-1]:
            lat_labels.append(lat_label(x, 2))
        axs[i-1].set_xticklabels(lat_labels, size = 12)

        axs[i-1].set_aspect("auto", adjustable="box")
        plt.clabel(cs, cs.levels, inline = True,  fmt=lambda x: f'{x:.0e}', fontsize = 16) 

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(cv, cax=cbar_ax, ticks = np.arange(-2e-04, 2.1e-04, 1e-04))
    tick_labels = [-2, -1, 0, 1, 2]
    cbar.ax.set_yticklabels(tick_labels)
    cbar.ax.text(1.4, 1, r'$\times 10^{-4}$', fontsize = font_size, transform=cbar.ax.transAxes)
    cbar.set_label("kg m$^{-3}$ century$^{-1}$", fontsize = font_size)
   
    plt.show()
        
def main():
    filepaths = [
        f"{output_path}{exp}_files_temp/geostrophy.nc"
        for exp in ["LENS", "WIND", "TEMP"]
    ]
    input_data = [xr.open_dataset(filepath, decode_times = False) for filepath in filepaths]
    vel = [input.geostrophy.values for input in input_data]
    z = input_data[0].Z.values
    lat = input_data[0]["YC"]

    color_scheme = "PiYG_r"
    experiment = ["ALL", "WIND", "THERMO"]
    low_val = -1e-13
    high_val = -low_val
    print(low_val)
    step = 15
    font_size = 20
        
    # create subplots with each variable on a new line
    fig, axs = plt.subplots(nrows=1, ncols=3, gridspec_kw={"hspace": 0.08, "wspace": 0.05}, figsize=(20,8)) 
    axs = axs.flatten()

    fmt = ticker.ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((0, 0))
    for i in range(3):
        cv = axs[i].contourf(lat, z, vel[i], cmap=color_scheme, extend="both", levels=np.linspace(low_val, high_val, step))
        axs[i].set_title(experiment[i], weight="bold", fontsize = font_size)
        axs[i].set_xlabel("Latitude", fontsize = font_size)
        axs[i].set_ylim([-1000, 0])
        axs[i].set_xlim([-74, -71.5])
        axs[i].tick_params(axis='x', labelsize=16, rotation=45)
        axs[i].tick_params(axis='y', labelsize=16)
        
        if i > 1:
            axs[i].get_yaxis().set_visible(False)
            axs[i].set_ylabel("Depth (m)", fontsize = font_size)
        
        axs[i].locator_params(axis='x', nbins=6)
        lat_ticks = axs[i].get_xticks()
        lat_labels = []
        for x in lat_ticks[:-1]:
            lat_labels.append(lat_label(x, 2))
        axs[i].set_xticklabels(lat_labels, size = 12)

        axs[i].set_aspect("auto", adjustable="box")
        
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(cv, cax=cbar_ax)
    cbar.set_label("m s$^{-1}$ century$^{-1}$", fontsize = font_size)
    plt.suptitle("geostrophic anomaly from the trends in density anomaly (m s$^{-1}$ century$^{-1}$)", fontsize = font_size+2, weight = "bold")
    
    fig.savefig("geostrophy.png")
    plt.show()
        
if __name__ == '__main__':
    main() # run the program

