import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import numpy as np
from funcs import read_data

def main():
    UC = 0 # undercurrent section
    [input_data, vel] = read_data("UVEL", UC)
    [input_data, salt] = read_data("SALT", UC)
      
    z = input_data[0].z.values
    lat = input_data[0]["lat"]
    color_scheme = "PiYG_r"
    experiment = ["pre-industrial", "high-emissions forcing", "wind forcing", "thermodynamic forcing"]
    low_val = -np.max(vel[1])
    high_val = -low_val
    low_sal = -np.max(salt[1])
    high_sal = -low_val
    print(low_sal, high_sal)
    step = 15
    font_size = 20
        
    # create subplots with each variable on a new line
    fig, axs = plt.subplots(nrows=3, ncols=1, gridspec_kw={"hspace": 0.08, "wspace": 0.05}, figsize=(8,20)) 
    axs = axs.flatten()

    fmt = ticker.ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((0, 0))
    for i in range(1, 4):
        cv = axs[i-1].contourf(lat, z, vel[i]-vel[0], cmap=color_scheme, extend="both", levels=np.linspace(low_val, high_val, step))
        cs = axs[i-1].contour(lat, z, salt[i]-salt[0], colors='black')
        
        axs[i-1].set_title(experiment[i], weight="bold", fontsize = font_size)
        axs[i-1].set_xlabel("Latitude", fontsize = font_size)
        axs[i-1].set_ylim([-1000, 0])
        axs[i-1].tick_params(axis='x', labelsize=16)
        axs[i-1].tick_params(axis='y', labelsize=16)
        if i < 3:
            axs[i-1].get_xaxis().set_visible(False)
        axs[i-1].set_aspect("auto", adjustable="box")
        plt.clabel(cs, cs.levels, inline = True,  fmt=lambda x: f'{x:.0e}', fontsize = 16) 
    
        axs[i-1].set_ylabel("Depth (m)", fontsize = font_size)  

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.82, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(cv, cax=cbar_ax, ticks = np.arange(-5e-05, 5.1e-05, 2.5e-05))
    cbar.set_label("m s$^{-1}$ century$^{-1}$", fontsize = font_size)

    
    plt.suptitle("Zonal velocity and salinity trends per century", fontsize = font_size+2, weight = "bold")
        
    fig.savefig(f"uvel_trend_undercurrent_{UC}.png", bbox_inches='tight', transparent=True)
    plt.show()
        
    
if __name__ == '__main__':
    main() # run the program
