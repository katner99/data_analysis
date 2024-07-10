from config_options import config_comparison, slice_ranges
from plots import read_mask, pretty_labels, read_var_fluxes, read_var_profile
from plots_2d import contour_func
from funcs import create_profile, read_data
from directories_and_paths import output_path, grid_filepath
from mitgcm_python.grid import Grid
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

def plot_fluxes(axs, fig, data, data_SH, exp_names, var, set_up, graph_params, graph_params_SH):
    grid = Grid(grid_filepath)
    dA = grid.dA

    for i, place in enumerate([0, 2, 4]):
        print(np.sum(data_SH[i]*dA/(10**9)))
        hide_x = False if place == 4 else True

        cs_SI = plot_contour(axs[place], data[i], set_up, graph_params, exp_names[i], hide_ticks_x=hide_x, hide_ticks_y=False)
        cs_SH = plot_contour(axs[place], data_SH[i], set_up, graph_params_SH, exp_names[i], hide_ticks_x=hide_x, hide_ticks_y=False)
        
        axs[place].contourf(set_up["X"], set_up["Y"], graph_params["pvalue"][i], levels=[-np.inf, 0.05], colors='none', hatches=['///'], alpha=0)
        axs[place].contourf(set_up["X"], set_up["Y"], graph_params_SH["pvalue"][i], levels=[-np.inf, 0.05], colors='none', hatches=['///'], alpha=0)
    axs[0].set_title("Trends in freshwater fluxes \n from sea-ice and ice shelf \n melting (m$^{2}$ yr$^{-1}$ century$^{-1}$)")

    ticks_SI = np.arange(-1.5, 1.6, 0.5)
    cbar_ax_SI = fig.add_axes([0.1, 0.01, 0.4, 0.02])
    cbar_SI = plt.colorbar(cs_SI, cax=cbar_ax_SI, orientation='horizontal')
    cbar_SI.set_ticks(ticks_SI)
    cbar_SI.ax.xaxis.set_ticks_position('bottom')  # Move the ticks to the bottom
    cbar_SI.ax.tick_params(axis='x', direction='inout', length=0)
    cbar_SI.ax.xaxis.set_tick_params(pad=-12) 


    ticks_SH = np.arange(-5, 5.1, 2.5)
    cbar_ax_SH = fig.add_axes([0.1, 0.035, 0.4, 0.02])
    cbar_SH = plt.colorbar(cs_SH, cax=cbar_ax_SH, orientation='horizontal')
    cbar_SH.set_ticks(ticks_SH) 
    tick_labels = [f'{tick:.1f}' for tick in ticks_SH]  # Format tick labels as desired
    cbar_SH.set_ticklabels(tick_labels)
    cbar_SH.ax.xaxis.set_ticks_position('bottom')  # Move the ticks to the bottom
    cbar_SH.ax.tick_params(axis='x', direction='inout', length=0)
    cbar_SH.ax.xaxis.set_tick_params(pad=-12) 

def plot_contour(ax, data, set_up, graph_params, title, hide_ticks_x=True, hide_ticks_y=True):
    cs = contour_func(ax, data, set_up, graph_params, hide_ticks_x, hide_ticks_y)
    ax.set_xlim([230, 265])
    ax.set_ylim([-75.5, -67])
    ax.set_ylabel(title, weight = 'bold')
    pretty_labels(ax)

    return cs

def plot_profiles(axs, input_data):
    grid = Grid(grid_filepath)
    lon_range = slice_ranges["lon_range_cont"]
    lat_range = slice_ranges["lat_range_cont"]
    
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
    ax.tick_params(axis='x', colors='orchid')

    c = ax2.plot(SALT_sim, z, color="skyblue", label = 'Forced Salinity')
    d = ax2.plot(SALT_ctr, z, color="skyblue", linestyle="--", label = 'Control Salinity')
    ax2.set_xlim([33, 35])
    ax2.set_xticks(np.arange(33.5, 35, 0.5))
    ax2.tick_params(axis='x', colors='skyblue')

    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.spines['right'].set_position(('outward', 0))
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('black')
    ax.set_ylabel("Depth (m)")
    ax.grid(alpha=0.8)

    return ax, ax2, a, b, c, d

def plot_slice(axs, fig):
    import matplotlib.ticker as ticker
    from plots import lat_label
    UC = 0 # undercurrent section
    [input_data, vel] = read_data("UVEL", UC)
    [input_data, salt] = read_data("DENSITY", UC)
      
    z = input_data[0].z.values
    lat = input_data[0]["lat"]
    color_scheme = "PiYG_r"
    experiment = ["NONE", "ALL", "WIND", "THERMO"]
    low_val = np.min(vel[1])
    high_val = -low_val
    low_sal = -2.5e-4
    high_sal = -low_sal
    print(low_sal, high_sal)
    step = 15
    font_size = 20

    fmt = ticker.ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((0, 0))
    axs[1].set_title("Zonal velocity (m s$^{-1}$ century$^{-1}$) and \ndensity (kg m$^{-3}$ century$^{-1}$) trends \nper century")
        
    for idx, position in enumerate([1, 3, 5]):
        cv = axs[position].contourf(lat, z, vel[idx+1]-vel[0], cmap=color_scheme, extend="both", levels=np.linspace(low_val, high_val, step))
        cs = axs[position].contour(lat, z, salt[idx+1]-salt[0], colors='black')
        axs[position].set_ylim([-1000, 0])
        axs[position].set_xlim([-74, -71.5])
        axs[position].tick_params(axis='x', rotation=45)
        axs[position].tick_params(axis='y')
        axs[position].set_ylabel("Depth (m)")
        axs[position].yaxis.set_label_position("right")
        axs[position].yaxis.tick_right()
        
        if idx == 2:
            axs[position].locator_params(axis='x', nbins=6)
            lat_ticks = axs[position].get_xticks()
            lat_labels = []
            for x in lat_ticks[:-1]:
                lat_labels.append(lat_label(x, 2))
            axs[position].set_xticklabels(lat_labels, size = 12)
        else:
            axs[position].get_xaxis().set_visible(False)

        axs[position].set_aspect("auto", adjustable="box")
        plt.clabel(cs, cs.levels, inline = True,  fmt=lambda x: f'{x:.0e}', fontsize = 16) 

    cbar_ax_SI = fig.add_axes([0.55, 0.035, 0.4, 0.02])
    cbar_SI = plt.colorbar(cv, cax=cbar_ax_SI, orientation='horizontal')
    #cbar_SI.set_ticks(ticks_SI)
    #cbar_SI.ax.xaxis.set_ticks_position('bottom')  # Move the ticks to the bottom
    #cbar_SI.ax.tick_params(axis='x', direction='inout', length=0)
    #cbar_SI.ax.xaxis.set_tick_params(pad=-12)
    #cbar = fig.colorbar(cv, cax=cbar_ax, ticks = np.arange(-2e-04, 2.1e-04, 1e-04))
    #tick_labels = [-2, -1, 0, 1, 2]
    #cbar.ax.set_yticklabels(tick_labels)
    #cbar.ax.text(1.4, 1, r'$\times 10^{-4}$', fontsize = font_size, transform=cbar.ax.transAxes)
    #cbar.set_label("kg m$^{-3}$ century$^{-1}$", fontsize = font_size)

    
    #plt.suptitle("Zonal velocity (m s$^{-1}$ century$^{-1}$) and density (kg m$^{-3}$ century$^{-1}$) trends per century", fontsize = font_size+2, weight = "bold")
        
    #fig.savefig(f"uvel_trend_undercurrent_{UC}.png", bbox_inches='tight', transparent=False)
    #plt.show()

def main():
    var = "oceFWflx"
    file_out = f"mega_comparison_trend_{var}.png"
    experiment = ["ALL", "WIND", "THERMO"]

    input_data_SI = read_var_fluxes(var)          # sea_ice freshwater flux
    input_data_SH = read_var_fluxes("SHIfwFlx")   # ice shelf melt fluxes
    input_data_ST = read_var_profile()            # salinity and temperature profiles

    grid = Grid(grid_filepath)

    set_up_data = xr.open_dataset(f"{output_path}average_CTRL_1920-1950.nc", decode_times=False)
    set_up = read_mask(set_up_data)
    [data_SI, graph_params_SI, graph_params_anom_SI] = config_comparison("trend", input_data_SI, grid, var_name=var)
    [data_SH, graph_params_SH, graph_params_anom_SH] = config_comparison("trend", input_data_SH, grid, var_name="SHIfwFlx")

    fig, axs = plt.subplots(3, 2, figsize=(10, 12))
    plt.subplots_adjust(wspace=0.05, hspace=0.05)
    axs = axs.flatten()

    plot_fluxes(axs, fig, data_SI, data_SH, experiment, var, set_up, graph_params_SI, graph_params_SH)
    plot_slice(axs, fig)
    #plot_profiles(axs, input_data_ST)
    fig.savefig(file_out, transparent=False)
    #plt.show()

if __name__ == "__main__":
    main()