from config_options import config_comparison, SLICE_RANGES
from tools.plots import read_mask, pretty_labels, read_var_fluxes, read_var_profile
from tools.plots_2d import contour_func
from tools.funcs import create_profile, read_data
from tools.directories_and_paths import grid_filepath
from mitgcm_python.grid import Grid
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mitgcm_python.plot_utils.labels import lat_label
import numpy as np


#def plot_fluxes(axs, fig, data, data_SH, exp_names, set_up, graph_params, graph_params_SH):
    

def plot_contour(ax, data, set_up, graph_params, title, hide_ticks_x=True, hide_ticks_y=True):
    cs = contour_func(ax, data, set_up, graph_params, hide_ticks_x, hide_ticks_y)
    ax.set_xlim([230, 265])
    ax.set_ylim([-75.5, -67])
    ax.set_ylabel(title, weight = 'bold')
    pretty_labels(ax)

    return cs

def plot_profiles(axs, input_data):
    grid = Grid(grid_filepath)
    lon_range = SLICE_RANGES["lon_range_cont"]
    lat_range = SLICE_RANGES["lat_range_cont"]
    
    profiles_SALT = [create_profile(data, "SALT", grid, lat_range, lon_range, timeseries=False) for data in input_data]
    profiles_TEMP = [create_profile(data, "THETA", grid, lat_range, lon_range, timeseries=False) for data in input_data]
    z = input_data[0].Z.values
      
    temp_handles = []  # to store handles for temperature legend
    temp_labels = []   # to store labels for temperature legend
    salt_handles = []  # to store handles for salinity legend
    salt_labels = []   # to store labels for salinity legend
  
    for i, position in enumerate([1, 3, 5]):
        axs[position].get_xaxis().set_visible(True)

        [ax1, ax2, a ,b ,c, d] = plot_profile(axs[position], profiles_TEMP[i + 1], profiles_TEMP[0], profiles_SALT[i + 1], profiles_SALT[0], z)

        if i > 0:
            ax2.get_xaxis().set_ticklabels([])
        if i < 2:
            ax1.get_xaxis().set_ticklabels([])  

    axs[1].set_title("End of century salinity and temperature \n profiles over the continental shelf")
    axs[5].legend(['Forced Temperature', 'Control Temperature'], loc='lower center', bbox_to_anchor=(0.55, -0.3), ncol=2)
    ax2.legend(['Forced Salinity', 'Control Salinity'], loc='lower center', bbox_to_anchor=(0.55, -0.45), ncol=2)
    
        
def plot_profile(ax, TEMP_sim, TEMP_ctr, SALT_sim, SALT_ctr, z):
    ax2 = ax.twiny()
    a = ax.plot(TEMP_sim, z, color="orchid", label = 'Forced Temperature')
    b = ax.plot(TEMP_ctr, z, color="orchid", linestyle="--", label = 'Control Temperature')
    ax.set_xlim([-2, 2])
    ax.set_xticks(np.arange(-2, 2, 1))
    ax.set_ylim([-600, 0])
    ax.tick_params(axis='x', colors='orchid', labelsize = 16)

    c = ax2.plot(SALT_sim, z, color="skyblue", label = 'Forced Salinity')
    d = ax2.plot(SALT_ctr, z, color="skyblue", linestyle="--", label = 'Control Salinity')
    ax2.set_xlim([33, 35])
    ax2.set_xticks(np.arange(33.5, 35, 0.5))
    ax2.tick_params(axis='x', colors='skyblue', labelsize = 16)

    ax.grid(alpha=0.8)

    return ax, ax2

def plot_slice(axs, fig, col):
    UC = 0  # undercurrent section
    [input_data, vel] = read_data("UVEL", UC)
    [input_data, salt] = read_data("DENSITY", UC)
      
    z = input_data[0].z.values
    lat = input_data[0]["lat"]
    color_scheme = "PiYG_r"
    low_val = np.min(vel[1])
    high_val = -low_val
    low_sal = np.min(salt[1])
    high_sal = -low_sal
    print(low_sal, high_sal)
    step = 15

    fmt = ticker.ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((0, 0))
    axs[0][col].set_title("Trends in zonal velocity (m s$^{-1}$ century$^{-1}$) \nand density (kg m$^{-3}$ century$^{-1}$)")
        
    levels = np.arange(-4, 5, 1)  # Adjusted to include the 4 level
    cut = 70
    for i in range(3):
        cv = axs[i][col].contourf(lat[:-cut], z, vel[i][:,:-cut], cmap=color_scheme, extend="both", levels=np.linspace(low_val, high_val, step))
        cs = axs[i][col].contour(lat[:-cut], z, salt[i][:,:-cut], colors='black', levels=np.linspace(low_sal, high_sal, step))
        
        axs[i][col].set_ylim([-1000, 0])
        #axs[position].set_xlim([-74, -71.5])
        axs[i][col].tick_params(axis='x', rotation=45)
        axs[i][col].tick_params(axis='y')
        axs[i][col].set_ylabel("Depth (m)")
        axs[i][col].yaxis.set_label_position("right")
        axs[i][col].yaxis.tick_right()
        
        if i == 2:
            axs[i][col].locator_params(axis='x', nbins=6)
            lat_ticks = axs[i][col].get_xticks()
            lat_labels = []
            for x in lat_ticks[:-1]:
                lat_labels.append(lat_label(x, 2))
            axs[i][col].set_xticklabels(lat_labels, size=12)
        else:
            axs[i][col].get_xaxis().set_visible(False)
        
        axs[i][col].set_aspect("auto", adjustable="box")
        axs[i][col].clabel(cs, cs.levels, inline=True, fmt = lambda x: f'{x:.0e}', fontsize=10)
    
    ticks = np.arange(-0.02, 0.03, 0.01)
    cbar_ax_SI = fig.add_axes([0.55, 0.075, 0.4, 0.02])
    cbar_SI = plt.colorbar(cv, cax=cbar_ax_SI, orientation='horizontal')
    cbar_SI.set_ticks(ticks) 
    cbar_SI.set_label("Velocity trend")

    #cbar_SI.ax.tick_params(rotation=45)
    



def just_profiles():
    input_data = read_var_profile()

    grid = Grid(grid_filepath)
    lon_range = SLICE_RANGES["lon_range_cont"]
    lat_range = SLICE_RANGES["lat_range_cont"]
    
    experiment = ["ALL", "THERMO", "WIND"]

    profiles_SALT = [create_profile(data, "SALT", grid, lat_range, lon_range, timeseries=False) for data in input_data]
    profiles_TEMP = [create_profile(data, "THETA", grid, lat_range, lon_range, timeseries=False) for data in input_data]
    z = input_data[0].Z.values

    fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(18, 10))
  
    for i in range(3):
        axs[i].get_xaxis().set_visible(True)
        axs[i].set_title(experiment[i], fontsize = 16)
        [ax1, ax2] = plot_profile(axs[i], profiles_TEMP[i + 1], profiles_TEMP[0], profiles_SALT[i + 1], profiles_SALT[0], z)

        if i > 0:
            ax2.get_yaxis().set_ticklabels([])
        else:
            axs[0].set_ylabel("Depth (m)", fontsize = 16)
            axs[0].legend(['Forced Temperature', 'NONE Temperature'], fontsize = 16, loc = 'lower left')
            ax2.legend(['Forced Salinity', 'NONE Salinity'], fontsize = 16, loc = 'upper right', bbox_to_anchor=(0.66, 0.22))



    plt.suptitle("End of century salinity and temperature profiles over the continental shelf", fontsize = 18, weight = 'bold')
    
    fig.savefig("profiles.png", transparent=False)
    
def plot_fluxes(axs, fig, col, data_SI, data_SH, exp_names, grid, set_up, graph_params_SI, graph_params_SH, pvalue = "False", hide_y=False):
    dA = grid.dA

    for i in range(len(exp_names)):
        print(np.sum(data_SH[i]*dA/(10**9))) # print the total melt rate 
        
        hide_x = False if i == len(exp_names)-1 else True # don't plot x axis unless you are at the last plot to save space

        cs_SI = plot_contour(axs[i][col], data_SI[i], set_up, graph_params_SI, exp_names[i], hide_ticks_x=hide_x, hide_ticks_y=hide_y)
        cs_SH = plot_contour(axs[i][col], data_SH[i], set_up, graph_params_SH, exp_names[i], hide_ticks_x=hide_x, hide_ticks_y=hide_y)
        
        if pvalue:
            axs[i][col].contourf(set_up["X"], set_up["Y"], graph_params_SI["pvalue"][i], levels=[-np.inf, 0.05], colors='none', hatches=['///'], alpha=0)
            axs[i][col].contourf(set_up["X"], set_up["Y"], graph_params_SH["pvalue"][i], levels=[-np.inf, 0.05], colors='none', hatches=['///'], alpha=0)
    
    axs[0][col].set_title("Trends in freshwater fluxes \n from sea-ice and ice shelf \n melting (m yr$^{-1}$ century$^{-1}$)")

    ticks_SI = np.arange(-1.5, 1.6, 0.5)
    cbar_ax_SI = fig.add_axes([0.1, 0.075, 0.4, 0.02])
    cbar_SI = plt.colorbar(cs_SI, cax=cbar_ax_SI, orientation='horizontal')
    cbar_SI.set_ticks(ticks_SI)
    cbar_SI.ax.xaxis.set_ticks_position('bottom')  # Move the ticks to the bottom
    cbar_SI.ax.tick_params(axis='x', direction='inout', length=0)
    cbar_SI.ax.xaxis.set_tick_params(pad=-12) 
    cbar_SI.set_label("Sea Ice FW flux trend")

    ticks_SH = np.arange(-5, 5.1, 2.5)
    cbar_ax_SH = fig.add_axes([0.1, 0.035, 0.4, 0.02])
    cbar_SH = plt.colorbar(cs_SH, cax=cbar_ax_SH, orientation='horizontal')
    cbar_SH.set_ticks(ticks_SH) 
    tick_labels = [f'{tick:.1f}' for tick in ticks_SH]  # Format tick labels as desired
    cbar_SH.set_ticklabels(tick_labels)
    cbar_SH.ax.xaxis.set_ticks_position('bottom')  # Move the ticks to the bottom
    cbar_SH.ax.tick_params(axis='x', direction='inout', length=0)
    cbar_SH.ax.xaxis.set_tick_params(pad=-12) 
    cbar_SH.set_label("Ice shelf FW flux trend")

def plot_all():
    """
    Function creates Figure 4 from Turner et al. (2024)
    """

    file_out = f"mega_comparison_trend_freshwater.png"
    exp_names = ["ALL", "THERMO", "WIND"]        # experiment names

    ### PLOT LEFT COLUMN - THE TOTAL FRESHWATER FLUXES ###   
    input_data_SI = read_var_fluxes("oceFWflx")   # total freshwater flux
    input_data_SH = read_var_fluxes("SHIfwFlx")   # ice shelf melt fluxes
   
    grid = Grid(grid_filepath)                    # MITgcm grid info
    set_up = read_mask()                          # read X, Y and land and ice masks

    [data_SI, graph_params_SI, _] = config_comparison("oceFWflx", input_data_SI, loc="trend", calc_pval="ctrl")
    [data_SH, graph_params_SH, _] = config_comparison("SHIfwFlx", input_data_SH, loc="trend", calc_pval="ctrl")

    fig, axs = plt.subplots(3, 2, figsize=(10, 12))
    plt.subplots_adjust(wspace=0.05, hspace=0.05)  

    grid = Grid(grid_filepath)

    plot_fluxes(axs, fig, 0, data_SI, data_SH, exp_names, grid, set_up, graph_params_SI, graph_params_SH, pvalue = True)
    axs[0][0].vlines(x=238, ymin=-74, ymax=-71.5, color='black', linewidth = 3) # if plotting the velocity profiles, indicate where on the map
    plot_slice(axs, fig, 1)
    plt.subplots_adjust(bottom = 0.15)
    #plt.text(0, 0.01, "Ice shelf trends")
    #plt.text(0, -0.01, "Sea ice trends")
    #plot_profiles(axs, input_data_ST)
    fig.savefig(file_out, transparent=False)
    #plt.show()

def main():
    """
    Function that creates plot 4 and S3 from the paper Turner et al. (2024)
    plot 4 is created using  the plot_all function
    plot S3 is created using the just_profiles function
    """
    plot_all()
    #just_profiles()

if __name__ == "__main__":
    main()